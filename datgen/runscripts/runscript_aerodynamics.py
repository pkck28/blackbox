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
            Method is used to set aero solver options provided by user.
        """

        filehandler = open("input.pickle", 'rb') 
        input = pickle.load(filehandler)

        self.aero_options = input["aeroSolverOptions"]

        self.aero_options["gridFile"] = "../../" + self.aero_options["gridFile"]

        self.sample = input["sample"]

    def setup(self):

        self.designVariables = ["mach"]

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

        # Adding an aerodynamics group which does analysis, pre and post postprocessing
        self.mphys_add_scenario("scenario", ScenarioAerodynamic(aero_builder=adflow_builder))

        # Connecting output of mesh component to input of scenario group
        self.connect("mesh.x_aero0", "scenario.x_aero")

    def configure(self):

        aero_problem = AeroProblem(
            name="ap", mach=0.8, altitude=10000, alpha=2.0, areaRef=45.5, chordRef=3.25, evalFuncs=["cl", "cd"]
        )

        if "aoa" in self.sample.keys():
            # You need to add name while adding dv. It is what is used for output/input
            aero_problem.addDV("alpha", name="aoa", units="deg")
        if "mach" in self.sample.keys():
            aero_problem.addDV("mach", name="mach")

        # set the aero problem in the coupling and post coupling groups
        self.scenario.coupling.mphys_set_ap(aero_problem)
        self.scenario.aero_post.mphys_set_ap(aero_problem)

        if "aoa" in self.sample.keys():
            self.connect("aoa", ["scenario.coupling.aoa", "scenario.aero_post.aoa"])

        if "mach" in self.sample.keys():
            self.connect("mach", ["scenario.coupling.mach", "scenario.aero_post.mach"])

prob = om.Problem()
prob.model = Top()

prob.model.setOptions()

prob.setup()

prob.run_model()

if prob.model.comm.rank == 0:
    print("Scenario 0")
    print("cl =", prob["scenario.aero_post.cl"])
    print("cd =", prob["scenario.aero_post.cd"])

if prob.model.comm.rank == 0:
    output = {"cl":prob["scenario.aero_post.cl"], "cd":prob["scenario.aero_post.cd"]}
    
    filehandler = open("output.pickle", "xb")
    pickle.dump(output, filehandler)
