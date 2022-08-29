# Importing classes requried for openmdao and mphys
import openmdao.api as om
from mphys.multipoint import Multipoint
from  mphys.scenario_aerostructural import ScenarioAeroStructural

# Importing builders for required tools
from adflow.mphys import ADflowBuilder
from tacs.mphys import TacsBuilder
from funtofem.mphys import MeldBuilder
from tacs import elements, constitutive, functions

# Importing other python packages
from baseclasses import AeroProblem
from mpi4py import MPI
import numpy as np
import pickle

comm = MPI.COMM_WORLD

# Material properties
rho = 2780.0  # density, kg/m^3
E = 73.1e9  # elastic modulus, Pa
nu = 0.33  # poisson's ratio
kcorr = 5.0 / 6.0  # shear correction factor
ys = 324.0e6  # yield stress, Pa

# Shell thickness
t = 0.003  # m
tMin = 0.002  # m
tMax = 0.05  # m

def element_callback(dvNum, compID, compDescript, elemDescripts, specialDVs, **kwargs):
    # Setup (isotropic) property and constitutive objects
    prop = constitutive.MaterialProperties(rho=rho, E=E, nu=nu, ys=ys)
    # Set one thickness dv for every component
    con = constitutive.IsoShellConstitutive(prop, t=t, tNum=dvNum, tlb=tMin, tub=tMax)

    # For each element type in this component,
    # pass back the appropriate tacs element object
    transform = None
    elem = elements.Quad4Shell(transform, con)

    return elem

def problem_setup(scenario_name, fea_assembler, problem):
    """
    Helper function to add fixed forces and eval functions
    to structural problems used in tacs builder
    """
    # Add TACS Functions
    # Only include mass from elements that belong to pytacs components (i.e. skip concentrated masses)
    problem.addFunction("mass", functions.StructuralMass)
    problem.addFunction("ks_vmfailure", functions.KSFailure, safetyFactor=1.0, ksWeight=50.0)

    # Add gravity load
    g = np.array([0.0, 0.0, -9.81])  # m/s^2
    problem.addInertialLoad(g)

structOptions = {
    "element_callback": element_callback,
    "problem_setup": problem_setup,
    "mesh_file": "wingbox.bdf",
}

class Top(Multipoint):
    """
        Class containing aerostruct openmdao model
    """

    def setOptions(self):
        """
            Method used to set various parameters need by openmdao model.
        """

        filehandler = open("input.pickle", 'rb') 
        input = pickle.load(filehandler)
        filehandler.close()

        self.aero_options = input["aeroSolverOptions"]
        self.struct_options = structOptions

        self.aero_options["gridFile"] = "../../" + self.aero_options["gridFile"]
        self.struct_options["mesh_file"] = "../../" + self.struct_options["mesh_file"]

        self.sample = input["sample"]

    def setup(self):
        """
            This method is executed when user calls setup method of problem class. It builds
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
        ldxfer_builder = MeldBuilder(aero_builder, struct_builder, isym=2)
        ldxfer_builder.initialize(self.comm)

        # Adding an IVC for design variables / parameters of the wing
        # It promotes all the variables present in it.
        dvs = self.add_subsystem("dvs", om.IndepVarComp(), promotes=["*"])

        # Adding design variables in the design variable IVC
        if "aoa" in self.sample.keys():
            self.dvs.add_output("aoa", val=self.sample["aoa"], units="deg")
        if "mach" in self.sample.keys():
            self.dvs.add_output("mach", val=self.sample["mach"])

        # Setting solvers for coupled primal and adjoint solution
        # Nonlinear Block Gauss-Seidel is for primal solution
        # Linear Block Gauss-Seidel is for adjoint solution
        nonlinear_solver = om.NonlinearBlockGS(maxiter=25, iprint=2, use_aitken=True, rtol=1e-7, atol=1e-8)
        linear_solver = om.LinearBlockGS(maxiter=25, iprint=2, use_aitken=True, rtol=1e-6, atol=1e-8)

        aerostruct_scenario = ScenarioAeroStructural(aero_builder=aero_builder, 
                                                     struct_builder=struct_builder, 
                                                     ldxfer_builder=ldxfer_builder)

        # Adding the group called scenario (contains the entire coupling)
        self.mphys_add_scenario("scenario", aerostruct_scenario, nonlinear_solver, linear_solver)

        # Connecting aero and struct mesh component's output to the inputs in teh above added scenario group
        for discipline in ["aero", "struct"]:
            self.mphys_connect_scenario_coordinate_source("mesh_%s" % discipline, "scenario", discipline)

        # Add the structural thickness DVs
        ndv_struct = struct_builder.get_ndv()
        dvs.add_output("dv_struct", np.array(ndv_struct * [0.01]))

        # Connecting output of dvs to all the inputs in scenario group
        self.connect("dv_struct", "scenario.dv_struct") 

    def configure(self):
        """
            This method is executed after setup method of problem class is executed.
        """

        aerostruct_problem = AeroProblem(
            name="ap",
            mach=0.8,
            altitude=10000,
            alpha=2.0,
            areaRef=45.5,
            chordRef=3.25,
            evalFuncs=["cl", "cd"],
        )

        if "aoa" in self.sample.keys():
            # You need to add name while adding dv. It is used for output/input.
            aerostruct_problem.addDV("alpha", name="aoa", units="deg")
        if "mach" in self.sample.keys():
            aerostruct_problem.addDV("mach", name="mach")

        self.scenario.coupling.aero.mphys_set_ap(aerostruct_problem)
        self.scenario.aero_post.mphys_set_ap(aerostruct_problem)

        if "aoa" in self.sample.keys():
            self.connect("aoa", ["scenario.coupling.aero.aoa", "scenario.aero_post.aoa"])

        if "mach" in self.sample.keys():
            self.connect("mach", ["scenario.coupling.aero.mach", "scenario.aero_post.mach"])

prob = om.Problem()
prob.model = Top()

prob.model.setOptions()

prob.setup()

prob.run_model()

om.n2(prob, outfile="n2.html", show_browser=False)

if prob.model.comm.rank == 0:
    print("Scenario 0")
    print("cl = ", prob["scenario.aero_post.cl"])
    print("cd = ", prob["scenario.aero_post.cd"])
    print("failure = ", prob["scenario.struct_post.eval_funcs.ks_vmfailure"])

if prob.model.comm.rank == 0:
    output = {"cl" : prob["scenario.aero_post.cl"],
              "cd" : prob["scenario.aero_post.cd"],
              "failure" : prob["scenario.struct_post.eval_funcs.ks_vmfailure"]}
    
    filehandler = open("output.pickle", "xb")
    pickle.dump(output, filehandler)
    filehandler.close()
