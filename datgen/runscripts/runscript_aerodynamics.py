# Importing classes requried for openmdao and mphys
import openmdao.api as om
from mphys.multipoint import Multipoint
from  mphys.scenario_aerodynamic import ScenarioAerodynamic

# Importing builders for required tools
from adflow.mphys import ADflowBuilder

# Importing other python packages
from baseclasses import AeroProblem
from mpi4py import MPI
import pickle

comm = MPI.COMM_WORLD

class Top(Multipoint):
    """
        Class containing aerodynamics openmdao model
    """

    def setOptions(self):
        """
            Method used to set various parameters need by openmdao model.
        """

        filehandler = open("input.pickle", 'rb') 
        input = pickle.load(filehandler)
        filehandler.close()

        self.aero_options = input["aeroSolverOptions"]

        self.sample = input["sample"]

        self.parameters = input["parameters"]

        self.objectives = input["objectives"]

    def setup(self):

        adflow_builder = ADflowBuilder(self.aero_options, scenario="aerodynamic")
        adflow_builder.initialize(self.comm)

        # Adding mesh component which is IVC and outputs surface mesh coordinates
        self.add_subsystem("mesh", adflow_builder.get_mesh_coordinate_subsystem())

        # IVC to keep the top level DVs
        self.add_subsystem("dvs", om.IndepVarComp(), promotes=["*"])

        if "aoa" in self.sample.keys():
            self.dvs.add_output("aoa", val=self.sample["aoa"], units="deg")
        if "mach" in self.sample.keys():
            self.dvs.add_output("mach", val=self.sample["mach"])
        if "altitude" in self.sample.keys():
            self.dvs.add_output("altitude", val=self.sample["altitude"], units="m")

        # Adding an aerodynamics group which does analysis, pre and post postprocessing
        self.mphys_add_scenario("scenario", ScenarioAerodynamic(aero_builder=adflow_builder))

        # Connecting output of mesh component to input of scenario group
        self.connect("mesh.x_aero0", "scenario.x_aero")

    def configure(self):

        # Checking if the three entities are design variable or not
        # If yes, then initialize their value, it will be over-riden by the ivc.
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

        aero_problem = AeroProblem(
            name="ap",
            mach=mach,
            altitude=altitude,
            alpha=alpha,
            areaRef=areaRef,
            chordRef=chordRef,
            evalFuncs=evalFuncs
        )

        if "aoa" in self.sample.keys():
            # You need to add name while adding dv. It is the same name used for output/input.
            aero_problem.addDV("alpha", name="aoa", units="deg")
        if "mach" in self.sample.keys():
            aero_problem.addDV("mach", name="mach")
        if "altitude" in self.sample.keys():
            aero_problem.addDV("altitude", name="altitude", units="m")

        # set the aero problem in the coupling and post coupling groups
        self.scenario.coupling.mphys_set_ap(aero_problem)
        self.scenario.aero_post.mphys_set_ap(aero_problem)

        if "aoa" in self.sample.keys():
            self.connect("aoa", ["scenario.coupling.aoa", "scenario.aero_post.aoa"])

        if "mach" in self.sample.keys():
            self.connect("mach", ["scenario.coupling.mach", "scenario.aero_post.mach"])

        if "altitude" in self.sample.keys():
            self.connect("altitude", ["scenario.coupling.altitude", "scenario.aero_post.altitude"])

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
    
    filehandler = open("output.pickle", "xb")
    pickle.dump(output, filehandler)
