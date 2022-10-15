# Importing classes requried for openmdao and mphys
import openmdao.api as om
from mphys.multipoint import Multipoint
from  mphys.scenario_aerostructural import ScenarioAeroStructural

# Importing builders for required tools
from adflow.mphys import ADflowBuilder
from tacs.mphys import TacsBuilder
from funtofem.mphys import MeldBuilder

# Importing structural solver setup file
import tacsSetup

# Importing other python packages
from baseclasses import AeroProblem
from mpi4py import MPI
import numpy as np
import pickle

comm = MPI.COMM_WORLD

class Top(Multipoint):
    """
        Class containing aerostruct openmdao model
    """

    def setOptions(self):
        """
            Method used to set various parameters for a single analysis.
        """

        # Reading input file for the analysis
        filehandler = open("input.pickle", 'rb') 
        input = pickle.load(filehandler)
        filehandler.close()

        self.aero_options = input["aeroSolverOptions"]
        self.aero_options["gridFile"] = "grid.cgns"

        self.sample = input["sample"]

        self.parameters = input["parameters"]

        self.objectives = input["objectives"]

        self.struct_options = {
            "element_callback": tacsSetup.element_callback,
            "problem_setup": tacsSetup.problem_setup,
            "mesh_file": "mesh.bdf",
        }

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
        if "altitude" in self.sample.keys():
            self.dvs.add_output("altitude", val=self.sample["altitude"], units="m")

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

        # Connecting aero and struct mesh component's output to the inputs in the above added scenario group
        for discipline in ["aero", "struct"]:
            self.mphys_connect_scenario_coordinate_source("mesh_%s" % discipline, "scenario", discipline)

        # Add the structural thickness - they are not DVs yet
        ndv_struct = struct_builder.get_ndv()
        dvs.add_output("dv_struct", np.array(ndv_struct * [0.01]))

        # Connecting output of dvs to all the inputs in scenario group
        self.connect("dv_struct", "scenario.dv_struct") 

    def configure(self):
        """
            This method is executed after setup method of problem class is executed.
        """

        # Checking if the three entities are design variable or not
        # If yes, then initialize their value, it will be over-ridden by the ivc.
        # If no, then assign the parameter value obtained from user.
        if "altitude" not in self.sample.keys():
            altitude = self.parameters["altitude"]
        else:
            altitude = 10000 # dummy value for initialization, in m
        
        if "mach" not in self.sample.keys():
            mach = self.parameters["mach"]
        else:
            mach = 0.8 # dummy value for initialization
        
        if "aoa" not in self.sample.keys():
            alpha = self.parameters["aoa"]
        else:
            alpha = 2.0 # dummy value for initialization, in deg

        # Assign parameters directly here, since these cannot be design variables.
        areaRef = self.parameters["areaRef"]
        chordRef = self.parameters["chordRef"]

        # Assigning objectives obtained from user
        evalFuncs = self.objectives

        aerostruct_problem = AeroProblem(
            name="asp",
            mach=mach,
            altitude=altitude,
            alpha=alpha,
            areaRef=areaRef,
            chordRef=chordRef,
            evalFuncs=evalFuncs
        )

        # You need to add name while adding dv. It is the same name used for output/input.
        if "aoa" in self.sample.keys():
            aerostruct_problem.addDV("alpha", name="aoa", units="deg")
        if "mach" in self.sample.keys():
            aerostruct_problem.addDV("mach", name="mach")
        if "altitude" in self.sample.keys():
            aerostruct_problem.addDV("altitude", name="altitude", units="m")

        self.scenario.coupling.aero.mphys_set_ap(aerostruct_problem)
        self.scenario.aero_post.mphys_set_ap(aerostruct_problem)

        if "aoa" in self.sample.keys():
            self.connect("aoa", ["scenario.coupling.aero.aoa", "scenario.aero_post.aoa"])

        if "mach" in self.sample.keys():
            self.connect("mach", ["scenario.coupling.aero.mach", "scenario.aero_post.mach"])

        if "altitude" in self.sample.keys():
            self.connect("altitude", ["scenario.coupling.aero.altitude", "scenario.aero_post.altitude"])

prob = om.Problem()
prob.model = Top()

prob.model.setOptions()

prob.setup()

prob.run_model()

if prob.model.comm.rank == 0:

    output = {}

    for value in prob.model.objectives:
        if "cl" == value:
            print("cl = ", prob["scenario.aero_post.cl"])
            output["cl"] = prob["scenario.aero_post.cl"]
        if "cd" == value:
            print("cd = ", prob["scenario.aero_post.cd"])
            output["cd"] = prob["scenario.aero_post.cd"]
        if "lift" == value:
            print("lift = ", prob["scenario.aero_post.lift"])
            output["lift"] = prob["scenario.aero_post.lift"]
        if "drag" == value:
            print("drag = ", prob["scenario.aero_post.drag"])
            output["drag"] = prob["scenario.aero_post.drag"]
        if prob["scenario.struct_post.eval_funcs.ks_vmfailure"]:
            print("failure = ", prob["scenario.struct_post.eval_funcs.ks_vmfailure"])
            output["failure"] = prob["scenario.struct_post.eval_funcs.ks_vmfailure"]
        if prob["scenario.struct_post.mass_funcs.mass"]:
            print("mass = ", prob["scenario.struct_post.mass_funcs.mass"])
            output["mass"] = prob["scenario.struct_post.mass_funcs.mass"]

    # print(prob.model.list_outputs(val=False))
    # print(type(prob.model.list_outputs(val=False)))
    # print("ks_vmfailure" in prob.model.list_outputs(val=False))
    # print("mass" in prob.model.list_outputs(val=False))
    
    filehandler = open("output.pickle", "xb")
    pickle.dump(output, filehandler)
    filehandler.close()
