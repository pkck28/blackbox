# Importing classes requried for openmdao and mphys
import openmdao.api as om
from mphys.multipoint import Multipoint
from  mphys.scenario_aerodynamic import ScenarioAerodynamic

# Importing builders for required tools
from adflow.mphys import ADflowBuilder
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
        Class creates a default option list (for solver and other settings)
        which is later edited/appended with user provided options.
    """

    def __init__(self):

        # Setting up the design variable dict
        self.designVariables = {}

        # Sampling options
        self.samplingMethod = "lhs"
        self.numberOfSamples = 2

        # Solver Options
        self.aeroSolver = "adflow"

        if self.aeroSolver == "adflow":
            self.aeroSolverOptions = {
                "printAllOptions": False,
                "printIntro": False
            }

class Top(Multipoint):
    """
        Class containing aerodynamics openmdao model
    """

    def setOptions(self, aero_options):
        """
            Method is used to set aero solver options provided by user.
        """
        self.aero_options = aero_options

    def setup(self):

        adflow_builder = ADflowBuilder(self.aero_options, scenario="aerodynamic")
        adflow_builder.initialize(self.comm)

        # Adding mesh component which is IVC and outputs surface mesh coordinates
        self.add_subsystem("mesh", adflow_builder.get_mesh_coordinate_subsystem())

        # IVC to keep the top level DVs
        self.add_subsystem("dvs", om.IndepVarComp(), promotes=["*"])

        # Adding an aerodynamics group which does analysis, pre and post postprocessing
        self.mphys_add_scenario("cruise", ScenarioAerodynamic(aero_builder=adflow_builder))

        # Connecting output of mesh component to input of cruise group
        self.connect("mesh.x_aero0", "cruise.x_aero")


    def configure(self):
        aoa = 1.5
        ap0 = AeroProblem(
            name="ap0", mach=0.8, altitude=10000, alpha=aoa, areaRef=45.5, chordRef=3.25, evalFuncs=["cl", "cd"]
        )
        ap0.addDV("alpha", value=aoa, name="aoa", units="deg")

        # set the aero problem in the coupling and post coupling groups
        self.cruise.coupling.mphys_set_ap(ap0)
        self.cruise.aero_post.mphys_set_ap(ap0)

        # add dvs to ivc and connect
        self.dvs.add_output("aoa", val=aoa, units="deg")

        # call the promote inputs to propagate aoa dvs
        # TODO does not work now
        # self.cruise._mphys_promote_inputs()
        # so connect manually
        self.connect("aoa", ["cruise.coupling.aoa", "cruise.aero_post.aoa"])

class Aerodynamics():
    
    def __init__(self, options):
        
        # If 'options' is None, notify the user
        if options is None:
            self._error("The 'options' argument not provided.")
        
        # If 'options' is not a dictionary, notify the user
        if not isinstance(options, dict):
            self._error("The 'options' argument provided is not a dictionary.")

        # Creating an empty options dictionary
        self.options = {}

        # Setting up default options
        self._getDefaultOptions()

        # Updating/Appending the default option list with user provided options
        self._setOptions(options)

        # Setting up OpenMDAO model
        self._setupModel()

    def _getDefaultOptions(self):
        """
            Setting up the initial values of options.
        """
        
        defaultOptions = DefaultOptions()

        for key in vars(defaultOptions):
            value = getattr(defaultOptions, key)
            self.options[key] = value

    def _setOptions(self, options):
        """
            Method for assigning user provided options.
        """

        # List of allowed design variables
        dvlist = ["aoa", "mach"]

        # Design variables checking and assignment
        if "designVariables" not in options or type(options["designVariables"]) == {}:
            self._error("Design variable option is not provided.")
        
        if type(options["designVariables"]) == dict:
            if not set(options["designVariables"].keys()).issubset(set(dvlist)):
                self._error("One or more design variable(s) is not recognized.")
        else:
            self._error("Design variable option is not a dictionary.")
        
        # Checking whether the provided option is valid
        for key in options.keys():
            if key in self.options.keys():
                # If the value is dictionary, update the default dictionary.
                # Otherwise, assign values.
                if isinstance(options[key], dict): 
                    self.options[key].update(options[key]) 
                else:
                    self.options[key] = options[key]
            else:
                self._error(key + " is not a valid option. Please remove/edit it.")

    def _setupModel(self):
        """
            Method to setup the openmdao model.
        """

        self.prob = om.Problem()
        self.prob.model = Top()

        self.prob.model.setOptions(self.options["aeroSolverOptions"])

        self.prob.setup()

        om.n2(self.prob, show_browser=False)

    def _lhs(self):
        """
            Method for creating a lhs sample.
        """

        # lower and upper bound are created as numpy arrays and dummy is a normal list here.
        lowerBound = np.array([])
        upperBound = np.array([])
        dummy = np.array([])
        self.samples = {}

        for key in self.options["designVariables"]:

            if key == "aoa" or "mach":
                lowerBound = np.append(lowerBound, self.options["designVariables"][key]["lowerBound"])
                upperBound = np.append(upperBound, self.options["designVariables"][key]["upperBound"])
            else:
                self._error("Unrecognized design variable")

            self.samples[key] = np.array([])
            dummy = np.append(dummy, key)

        dim = len(lowerBound)

        samples = lhs(dim, samples=self.options["numberOfSamples"], criterion='cm', iterations=50)

        samples = lowerBound + (upperBound - lowerBound) * samples

        for sampleNo in range(self.options["numberOfSamples"]):
            sample = samples[sampleNo,:]

            for key in self.options["designVariables"]:
                self.samples[key] = np.append(self.samples[key], sample[(dummy == key)])

    def generateSamples(self):
        """
            Method to generate samples.
        """

        if self.options["samplingMethod"] == "lhs":
            self._lhs()

        for sampleNo in range(self.options["numberOfSamples"]):
            self.prob['aoa'] = self.samples["aoa"][sampleNo]

            self.prob.run_model()

            if self.prob.model.comm == 0:
                print("cl = ", self.prob["cruise.aero_post.cl"])
                print("cd = ", self.prob["cruise.aero_post.cl"])

    def _error(self, message):
        """
            Method for printing a user mistake in nice manner.
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
