# Imports
import os, sys, shutil, pickle, time, psutil
import numpy as np
from scipy.io import savemat
from scipy import integrate
from pyDOE2 import lhs
from mpi4py import MPI
from baseclasses import AeroProblem
from pygeo import DVGeometry, geo_utils
from prefoil import Airfoil
from prefoil.utils import readCoordFile

# Trying to import matplotlib
try:
    import matplotlib.pyplot as plt
except ImportError:
    msg_matplotlib = "Matplotlib is not installed"
else:
    msg_matplotlib = None

comm = MPI.COMM_WORLD

class DefaultOptions():
    """
        Class creates a default option which are later 
        edited/appended with user provided options.
    """

    def __init__(self):
        
        # Aero solver Options
        self.aeroSolver = "adflow"

        # Other options
        self.directory = "output"
        self.noOfProcessors = 4
        self.writeSliceFile = False

class WingFFD():
    """
        This class provides methods for generating samples for a simple wing
        using FFD parameterization.
    """

    def __init__(self, options):

        # Partial checking of options argument provided by the user.
        if not isinstance(options, dict):
            self._error("The 'options' argument provided is not a dictionary.")
        elif options == {}:
            self._error("The 'options' argument provided is an empty dictionary.")

        # Creating an empty options dictionary
        self.options = {}

        # Setting up default options
        self._getDefaultOptions()

        # Setting up the required options list
        requiredOptions = ["solverOptions", "volumeMesh", "ffdFile", "aeroProblem"]

        # Validating user provided options
        self._checkOptions(options, requiredOptions)

        # Updating/Appending the default option list with user provided options
        self._setOptions(options)

        # Overiding/set some solver options
        self.options["solverOptions"]["printAllOptions"] = False
        self.options["solverOptions"]["printIntro"] = False
        self.options["solverOptions"]["outputDirectory"] = "."
        self.options["solverOptions"]["numberSolutions"] = False
        self.options["solverOptions"]["printTiming"] = False

        # Getting abs path for the storage directory
        self.options["directory"] = os.path.abspath(self.options["directory"])

        # Setting up the folder for saving the results
        directory = self.options["directory"]

        # Creating directory for storing the results
        if not os.path.isdir(directory):
            os.system("mkdir {}".format(directory))
        else:
            os.system("rm -r {}".format(directory))
            os.system("mkdir {}".format(directory))

        # Creating DVGeometry object
        self.DVGeo = DVGeometry(self.options["ffdFile"])

        # Some initializations which will be used later
        self.DV = []
        self.genSamples = 0

    # ----------------------------------------------------------------------------
    #                       Design Variable related methods
    # ----------------------------------------------------------------------------

    def addDV(self, name: str, lowerBound: list, upperBound: list) -> None:
        """
            Method for adding a DV for CST parameterization.
        """

        # Checking
        self._checkDV(name, lowerBound, upperBound)

        if name == "shape":

            # Getting the number of ffd points
            nffd = self.DVGeo.getLocalIndex(0).flatten().shape[0]

            if len(self.DV) == 0:
                self.upperBound = np.array([upperBound]*nffd)
                self.lowerBound = np.array([lowerBound]*nffd)
                self.locator = np.array(["{}".format(name)]*nffd)
            else:
                self.upperBound = np.append(self.upperBound, np.array([upperBound]*nffd))
                self.lowerBound = np.append(self.lowerBound, np.array([lowerBound]*nffd))
                self.locator = np.append(self.locator, np.array(["{}".format(name)]*nffd))

            # Adding ffd as dv
            self.DVGeo.addLocalDV("shape", lower=lowerBound, upper=upperBound, axis="y", scale=1.0)

        else:
            locator = np.array(["{}".format(name)])

            if len(self.DV) == 0:
                self.upperBound = np.array([upperBound])
                self.lowerBound = np.array([lowerBound])
                self.locator = np.array([locator])
            else:
                self.upperBound = np.append(self.upperBound, upperBound)
                self.lowerBound = np.append(self.lowerBound, lowerBound)
                self.locator = np.append(self.locator, locator)    

        # Adding the DV to the list
        self.DV.append(name)

    # ----------------------------------------------------------------------------
    #                       Methods related to validation
    # ----------------------------------------------------------------------------

    def _checkOptions(self, options: dict, requiredOptions: list) -> None:
        """
            This method validates user provided options.
        """

        defaultOptions = list(self.options.keys())
        allowedUserOptions = defaultOptions
        allowedUserOptions.extend(requiredOptions)
        userProvidedOptions = list(options.keys())

        # Checking if user provided option contains only allowed attributes
        if not set(userProvidedOptions).issubset(allowedUserOptions):
            self._error("Option dictionary contains unrecognized attribute(s): {}"\
                        .format(set(userProvidedOptions) - set(allowedUserOptions)))

        # Checking if user has mentioned all the requried attributes
        if not set(requiredOptions).issubset(userProvidedOptions):
            self._error("Option dictionary doesn't contain following attribute(s): {}"\
                        .format(set(requiredOptions) - set(userProvidedOptions)))

        ############ Validating aeroProblem
        if not isinstance(options["aeroProblem"], AeroProblem):
            self._error("\"aeroProblem\" attribute is not an aeroproblem.")

        ############ Validating solverOptions
        if not isinstance(options["solverOptions"], dict):
            self._error("\"solverOptions\" attribute is not a dictionary.")

        if "gridFile" in options["solverOptions"].keys():
            self._warning("\"gridFile\" attribute in solver options is not required.")

        ############ Validating ffdFile
        if not os.path.exists(os.path.abspath(options["ffdFile"])):
            self._error("Provided FFD file doesn't exists.")
        else:
            options["ffdFile"] = os.path.abspath(options["ffdFile"])

        ############ Validating volumeMesh
        if not os.path.exists(os.path.abspath(options["volumeMesh"])):
            self._error("Provided Volume mesh file doesn't exists.")
        else:
            options["volumeMesh"] = os.path.abspath(options["volumeMesh"])

        ############ Validating noOfProcessors
        if "noOfProcessors" in userProvidedOptions:
            if not isinstance(options["noOfProcessors"], int):
                self._error("\"noOfProcessors\" attribute is not an integer.")

            if psutil.cpu_count(False) < options["noOfProcessors"] + 1:
                self._error("\"noOfProcessors\" requested is more than available processors.")

        ############ Validating writeSliceFile
        if "writeSliceFile" in userProvidedOptions:
            if not isinstance(options["writeSliceFile"], bool):
                self._error("\"writeSliceFile\" attribute is not a boolean value.")

        ############ Validating directory attribute
        if "directory" in userProvidedOptions:
            if not isinstance(options["directory"], str):
                self._error("\"directory\" attribute is not string.")

    def _checkDV(self, name: str, lb: float or np.ndarray, ub: float or np.ndarray) -> None:
        """
            Method for validating DV user wants to add.
        """

        # List of possible DVs
        possibleDVs = ["shape", "alpha", "mach", "altitude"]

        # Validating name of the DV
        if not isinstance(name, str):
            self._error("Name argument is not a string.")

        # Checking if the DV is allowed
        if name not in possibleDVs:
            self._error("{} argument is not a valid DV.".format(name))

        # Checking if the DV is already added
        if name in self.DV:
            self._error("{} already exists.".format(name))

        # Checking if these variables are initialized through aero problem or not
        if name in ["mach", "altitude"]:
            if name not in self.options["aeroProblem"].inputs.keys():
                self._error("You need to initialize \"{}\" in the aero problem to set it as design variable.".format(name))

        if not isinstance(lb, float):
            self._error("Lower Bound argument is not a float.")

        if not isinstance(ub, float):
            self._error("Upper Bound argument is not a float.")

        if lb >= ub:
            self._error("Lower bound is greater than or equal to upper bound.")

    # ----------------------------------------------------------------------------
    #                               Other methods
    # ----------------------------------------------------------------------------

    def _getDefaultOptions(self) -> None:
        """
            Setting up the initial values of options.
        """

        defaultOptions = DefaultOptions()

        for key in vars(defaultOptions):
            value = getattr(defaultOptions, key)
            self.options[key] = value

    def _setOptions(self, options: dict) -> None:
        """
            Method for assigning user provided options.
        """

        for key in options.keys():
            # update strategy is different for a dictionary
            # and other type of variables
            if isinstance(options[key], dict):
                # if the option is already present in default dictionary, 
                # then append the user provided key-value pairs to the 
                # default dictionary. For example: solverOptions.
                if key in self.options.keys():
                    self.options[key].update(options[key]) 
                else:
                    self.options[key] = options[key]
            else:
                self.options[key] = options[key]

    def _lhs(self, numSamples) -> np.ndarray:
        """
            Method to generate the lhs samples.
        """

        # Number of dimensions
        dim = len(self.lowerBound)

        # Generating normalized lhs samples
        samples = lhs(dim, samples=numSamples, criterion='cm', iterations=1000)

        # Scaling the samples
        x = self.lowerBound + (self.upperBound - self.lowerBound) * samples

        return x

    def _warning(self, message: str) -> None:
        """
            Method for printing warnings in nice manner.
        """

        ############ To Do: Redundant - error and warning method can be combined.

        # Initial message - total len is 80 characters
        msg = "\n+" + "-" * 78 + "+" + "\n" + "| Blackbox Warning: "

        # Initial number of characters
        i = 16

        for word in message.split():
            if len(word) + i + 1 > 76:  # Finish line and start new one
                msg += " " * (76 - i) + " |\n| " + word + " " # Adding space and word in new line
                i = len(word) + 1 # Setting i value for new line
            else:
                msg += word + " " # Adding the word with a space
                i += len(word) + 1 # Increase the number of characters
        msg += " " * (76 - i) + " |\n" + "+" + "-" * 78 + "+" + "\n" # Adding last line
 
        print(msg, flush=True)

        exit()

    def _error(self, message: str) -> None:
        """
            Method for printing errors in nice manner.
        """

        # Initial message - total len is 80 characters
        msg = "\n+" + "-" * 78 + "+" + "\n" + "| Blackbox Error: "

        # Initial number of characters
        i = 16

        for word in message.split():
            if len(word) + i + 1 > 76:  # Finish line and start new one
                msg += " " * (76 - i) + " |\n| " + word + " " # Adding space and word in new line
                i = len(word) + 1 # Setting i value for new line
            else:
                msg += word + " " # Adding the word with a space
                i += len(word) + 1 # Increase the number of characters
        msg += " " * (76 - i) + " |\n" + "+" + "-" * 78 + "+" + "\n" # Adding last line
 
        print(msg, flush=True)

        exit()
