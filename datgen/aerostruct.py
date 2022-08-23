# Importing classes requried for openmdao and mphys
import openmdao.api as om
from mphys.multipoint import Multipoint
from  mphys.scenario_aerostructural import ScenarioAeroStructural

# Importing builders for required tools
from adflow.mphys import ADflowBuilder
from tacs.mphys import TacsBuilder
from funtofem.mphys import MeldBuilder
from pygeo.mphys import OM_DVGEOCOMP

# Importing other python packages
import numpy as np
from pygeo import geo_utils
from baseclasses import AeroProblem
from mpi4py import MPI
from pyDOE2 import lhs

comm = MPI.COMM_WORLD

class DefaultOptions():
    """
        Class creates a default option list (for solvers and other settings)
        which is later edited/appended with user provided options.
    """

    def __init__(self):

        self.printAllOptions = True
        self.dataFileName = "training_data"

        # Parameters used in aero problem instantiation
        self.flowParameters = {
            "mach" : 0.8,
            "altitude" : 10000,
            "areaRef" : 45.5,
            "chordRef" : 3.25 
        }

        # Setting up the design variable dict
        self.designVariables = {}

        # Sampling Method
        self.sampleType = "lhs"
        self.numberOfSamples = 2

        # Solver Options
        self.aeroSolver = "adflow"
        self.structSolver = "tacs"

        if self.aeroSolver == "adflow":
            self.aeroSolverOptions = {
                "printAllOptions": False,
                "printIntro": False
            }

        self.structSolverOptions = {}

        self.ffdFile = "./ffd.xyz"

class Top(Multipoint):
    """
        Class contains the openmdao model for aero-struct case.
    """

    def setOptions(self, aero_options, struct_options):
        """
            Method is used to set aero and struct solver options provided by user.
        """
        self.aero_options = aero_options
        self.struct_options = struct_options

    def setup(self):
        """
            This method is run by when user calls setup method of problem class. It builds
            the openmdao model, ready for analysis or optimization.
        """

        # Initialize the aerodynamics solver using builder class object
        aero_builder = ADflowBuilder(self.aero_options, scenario="aerostructural")
        aero_builder.initialize(self.comm)

        # Adding an aero mesh component which will output the initial surface 
        # mesh coordinates.
        self.add_subsystem("mesh_aero", aero_builder.get_mesh_coordinate_subsystem())

        # Initialize the structural solver using builder class object
        struct_builder = TacsBuilder(self.struct_options)
        struct_builder.initialize(self.comm)

        # Adding an struct mesh component which will output the initial surface 
        # mesh coordinates.
        self.add_subsystem("mesh_struct", struct_builder.get_mesh_coordinate_subsystem())

        # Initialize MELD using builder class object
        # Note: isym is an important parameter. Value depends on the orientation of span
        # of wing.
        # Note: Why to set check partials as true?
        ldxfer_builder = MeldBuilder(aero_builder, struct_builder, isym=2)
        ldxfer_builder.initialize(self.comm)

        # Adding an IVC for design variables / parameters of the wing
        # It promotes all the variables present in it.
        dvs = self.add_subsystem("dvs", om.IndepVarComp(), promotes=["*"])

        # Adding a subsystem for pygeo. It will take design variables and initial surface aero mesh
        # and structural mesh to output deformed surface aero mesh and struct mesh.
        ############################################################## TO DO: Create the options list properly.
        self.add_subsystem("geometry", OM_DVGEOCOMP(ffd_file="ffd.xyz"))

        # Setting solvers for coupled primal and adjoint solution
        # Nonlinear Block Gauss-Seidel is for primal solution
        # Linear Block Gauss-Seidel is for adjoint solution
        nonlinear_solver = om.NonlinearBlockGS(maxiter=25, iprint=2, use_aitken=True, rtol=1e-7, atol=1e-8)
        linear_solver = om.LinearBlockGS(maxiter=25, iprint=2, use_aitken=True, rtol=1e-6, atol=1e-8)

        # Adding the coupling sub-group to the model which aerostructal scenario
        ############################################################## TO DO: Create the option list properly.
        # Need to refactor such that user-provided name can be used.
        aero_struct_scenario = ScenarioAeroStructural(aero_builder=aero_builder, 
                                        struct_builder=struct_builder, ldxfer_builder=ldxfer_builder)
        self.mphys_add_scenario("scenario", aero_struct_scenario, nonlinear_solver, linear_solver)

        # Manually connecting vars in the geo component to cruise
        for discipline in ["aero", "struct"]:
            self.connect("geometry.x_%s0" % discipline, "scenario.x_%s0" % discipline)

        # Add the structural thickness DVs
        ndv_struct = struct_builder.get_ndv()
        dvs.add_output("dv_struct", np.array(ndv_struct * [0.01]))
        self.connect("dv_struct", "scenario.dv_struct")

        # Connecting output of mesh subsystems to the input of geometry subsystem
        # Data being transfered is initial aerodynamics surface mesh and structural mesh coordinates.
        self.connect("mesh_aero.x_aero0", "geometry.x_aero_in")
        self.connect("mesh_struct.x_struct0", "geometry.x_struct_in")

    def configure(self):

        # Call this to configure the coupling solver
        ########################################################## Do I need to do this?
        # super().configure()

        # Add the objective function to the scenario
        ########################################################## Do I need to do this?
        # self.scenario.aero_post.mphys_add_funcs()

        # Get the surface coordinates from the mesh component
        points = self.mesh_aero.mphys_get_surface_mesh()

        # Adding the pointset (contains the mesh coordinates - aero only) for both aero and struct
        self.geometry.nom_add_discipline_coords("aero", points)
        self.geometry.nom_add_discipline_coords("struct")

        ########### Twist Variable ###########

        # Create reference axis for the twist variable
        nRefAxPts = self.geometry.nom_addRefAxis(name="wingAxis", xFraction=0.25, alignIndex="k")

        # Set up global design variables. No change in root twist
        def twist(val, geo):
            for i in range(1, nRefAxPts):
                geo.rot_z["wingAxis"].coef[i] = -val[i - 1]

        # Add the twist design variable as the dvs component's output
        ############################################################# TO DO: need to change the value of twist: Provided by sample.
        self.dvs.add_output("twist", val=np.array([0] * (nRefAxPts - 1)))

        # Add twist variable in the pygeo's model
        self.geometry.nom_addGeoDVGlobal(dvName="twist", value=np.array([0] * (nRefAxPts - 1)), func=twist)

        # Manually connect twist from dvs to geometry
        self.connect("twist", "geometry.twist")

        ########### Shape Variable ###########

        # Adding shape as a variable in pygeo's model
        pts = self.geometry.DVGeo.getLocalIndex(0)
        indexList = pts[:, :, :].flatten()
        PS = geo_utils.PointSelect("list", indexList)
        nShapes = self.geometry.nom_addGeoDVLocal(dvName="shape", pointSelect=PS)

        # Adding shape as an output of dvs component
        ############################################################# TO DO: val must come from sample.
        self.dvs.add_output("shape", val=np.array([0] * nShapes))

        # Manually connect shape from dvs to geometry
        self.connect("shape", "geometry.shape")

        ########### AOA Variable ###########

        aoa0 = 5.0 # Initial value of aoa

        # Adding aoa as an output of dvs component
        ############################################################# TO DO: val must come from sample.
        self.dvs.add_output("aoa0", val=aoa0, units="deg")

        ########### Mach Variable ###########

        mach0 = 0.8 # Initial value of mach

        # Adding aoa as an output of dvs component
        ############################################################# TO DO: val must come from sample.
        self.dvs.add_output("mach0", val=mach0)

        ########### Setting up the aeroproblem for Adflow ###########

        aero_problem = AeroProblem(
            name="scenario",
            mach=0.5,
            altitude=10000,
            alpha=aoa0,
            areaRef=45.5,
            chordRef=3.25,
            evalFuncs=["lift", "drag", "cl", "cd"],
        )

        # Adding alpha as an input to the scenario group (I think it promotes the variables)
        aero_problem.addDV("alpha", value=aoa0, name="aoa", units="deg")

        # Adding mach as an input to the scenario group (I think it promotes the variables)
        aero_problem.addDV("mach", value=mach0, name="mach")

        # Adds the alpha as design variables
        self.scenario.coupling.aero.mphys_set_ap(aero_problem)
        self.scenario.aero_post.mphys_set_ap(aero_problem)

        # Manually connect aoa from dvs to the scenario group
        self.connect("aoa0", ["scenario.coupling.aero.aoa", "scenario.aero_post.aoa"])

        # Manually connect mach from dvs to the scenario group
        self.connect("mach0", ["scenario.coupling.aero.mach", "scenario.aero_post.mach"])

class AeroStruct():

    def __init__(self, options):

        self.aero_options = options["aeroSolverOptions"]
        self.struct_options = options["structSolverOptions"]
        
        # If 'options' is None, notify the user
        if options is None:
            self._error("The 'options' argument not provided.")
        
        # If 'options' is not a dictionary, notify the user
        if not isinstance(options, dict):
            self._error("The 'options' argument provided is not a dictionary.")

        self.options = {}

        # Setting up default options
        self._getDefaultOptions()

        # Updating/Appending the default option list with user provided options
        self._setOptions(options)

        self._bounds()

        self._lhs()

    def _getDefaultOptions(self):
        """
            Setting up the initial values of options which are common across all functions.
        """
        
        defaultOptions = DefaultOptions()

        for key in vars(defaultOptions):
            value = getattr(defaultOptions, key)
            self.options[key] = value

    def _setOptions(self, options):
        """
            Method for assigning user provided options. This method should be called only after checks.
        """

        dvlist = ["aoa", "mach", "twist", "shape"]

        if "designVariables" not in options or type(options["designVariables"]) == {}:
            self._error("Design variable option is not provided.")
        
        if type(options["designVariables"]) == dict:
            if not set(options["designVariables"].keys()).issubset(set(dvlist)):
                self._error("One or more design variable(s) is not recognized.")
        else:
            self._error("Design variable option is not a dictionary.")
        

        for key in options.keys():
            # Checking whether the provided option is valid
            if key in self.options.keys():
                if isinstance(options[key], dict): 
                    self.options[key].update(options[key]) 
                else:
                    self.options[key] = options[key]
            else:
                self._error(key + " is not a valid option. Please remove/edit it.")

    def _bounds(self):

        lowerBound = []
        upperBound = []

        dummy = []

        for key in self.options["designVariables"]:

            if key == "aoa":
                lowerBound.append(self.options["designVariables"][key]["lowerBound"])
                upperBound.append(self.options["designVariables"][key]["upperBound"])

                dummy.extend([key])

            if key == "mach":
                lowerBound.append(self.options["designVariables"][key]["lowerBound"])
                upperBound.append(self.options["designVariables"][key]["upperBound"])

                dummy.extend([key])

            if key == "twist":
                try:
                    sections = self.options["designVariables"][key]["noOfSections"]
                except:
                    self._error(key + " option provided in twist variable is not recognized")
                lowerBound.extend([self.options["designVariables"][key]["lowerBound"]] * sections)
                upperBound.extend([self.options["designVariables"][key]["upperBound"]] * sections)

                dummy.extend([key] * sections)

            if key == "shape":
                ffdPoints = self.options["designVariables"][key]["noOfFFDPoints"]
                lowerBound.extend([self.options["designVariables"][key]["lowerBound"]] * ffdPoints)
                upperBound.extend([self.options["designVariables"][key]["upperBound"]] * ffdPoints)

                dummy.extend([key] * ffdPoints) 

        self.lowerBound = np.array(lowerBound)
        self.upperBound = np.array(upperBound)
        self.dummy = dummy
        
    def _error(self, message):
        """
            Method for printing a user mistake in formatted manner.
        """

        msg = "\n+" + "-" * 78 + "+" + "\n" + "| Datgen Error: "
        i = 19
        for word in message.split():
            if len(word) + i + 1 > 78:  # Finish line and start new one
                msg += " " * (78 - i) + "|\n| " + word + " "
                i = 1 + len(word) + 1
            else:
                msg += word + " "
                i += len(word) + 1
        msg += " " * (78 - i) + "|\n" + "+" + "-" * 78 + "+" + "\n"

        if comm.rank == 0:
            print(msg, flush=True)

        exit()
        
    def _lhs(self):
        """
            Method for creating a lhs sample.
        """

        dim = len(self.lowerBound)

        samples = lhs(dim, samples=self.options["numberOfSamples"], criterion='cm', iterations=50)

        self.samples = self.lowerBound + (self.upperBound - self.lowerBound) * samples
        
    def _getSample(self, sampleNo, variable):
        """
            Method to get a specific sample
        """

        dummy = np.array(self.dummy)
        samples = self.samples[sampleNo,:]

        return samples[(dummy == variable)]        

    def generateSamples(self):

        prob = om.Problem()
        prob.model = Top()

        prob.model.setOptions(self.options["aeroSolverOptions"], self.options["structSolverOptions"])

        prob.setup()

        for i in range(self.options["numberOfSamples"]):

            # prob['aoa0'] = self._getSample(i, "aoa")
            prob['mach0'] = self._getSample(i, "mach")
            # prob['shape'] = self._getSample(i, "shape")
            # prob['twist'] = self._getSample(i, "twist")

            prob.run_model()

            if prob.model.comm.rank == 0:
           	    print("cl = ", prob["scenario.aero_post.cl"])
                # print("cd = ", prob["scenario.aero_post.cd"])
                # print("km =  ", prob["scenario.struct_post.km"])

            # om.n2(prob, show_browser=False)
