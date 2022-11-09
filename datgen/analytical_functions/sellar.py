# Importing python packages
import os, sys, shutil, pickle
import numpy as np
from scipy.io import savemat
from pyDOE2 import lhs, fullfact
from ..base import BaseClass

class DefaultOptions():
    """
        Class creates a default option list which is later 
        edited/appended with user provided options.
    """

    def __init__(self):
        self.directory = "output"

class Sellar(BaseClass):
    """
        Class contains essential methods for generating data
        of Sellar MDA problem (refer the main sellar paper):

            f(x1, x2, x3, y1, y2) = x1^2 + x3 + y1 + exp(-y2)
            y1(x1, x2, x3, y2) = x1^2 + x2 + x3 - 0.2y2
            y2(x1, x3, y1) = sqrt(y1) + x1 + x3
            g1(y1) = 1 - (y1/3.16) <= 0
            g2(y2) = (y2/24.0) - 1 <= 0

        -10 <= x1 <= 10, 0 <= x2 <= 10, 0 <= x3 <= 10

        There are two values possible for type option: "single" and "multi".
        Options argument is not needed for "single" type analysis, while it
        is needed for "multi" type analysis. Provide an appropriate options
        dict which contains necessary attributes. There are only three possible
        attributes:

        "directory" : Folder name where the data.mat file will be saved (string, optional).
        "numberOfSamples" : number of samples to be generated (integer).
        "samplingMethod" : name of the sampling method ("lhs" or "fullfactorial") (string).

        Note: For "single", options argument is not need.
    """

    def __init__(self, type="multi", options=None):

        if type == "multi":
            # If 'options' is None, notify the user
            if options is not None:
                if not isinstance(options, dict):
                    self._error("The 'options' argument provided is not a dictionary.")
                elif options == {}:
                    self._error("The 'options' argument provided is an empty dictionary.")
            else:
                self._error("Options argument not provided.")

            self._setupMultiAnalysis(options)

        elif type == "single":
            # If 'options' is not None, then raise an error
            if options is not None:
                self._error("Options argument is not needed when type is \"single\".")
            self._setupSingleAnalysis()

        else:
            self._error("Value of type argument not recognized.")


    # ----------------------------------------------------------------------------
    #                   All the methods for multi analysis
    # ----------------------------------------------------------------------------

    def _setupMultiAnalysis(self, options):
        """
            Method to perform initialization when the type is "multi".
        """

        # Creating an empty options dictionary and assign value of type
        self.options = {}
        self.options["type"] = "multi"

        # Setting up default options
        self._getDefaultOptions(DefaultOptions())

        # Validating user provided options
        self._checkOptionsForMultiAnalysis(options)

        self.options["lowerBound"] = [-10, 0, 0]
        self.options["upperBound"] = [10, 10, 10]
        self.d1_counts = 0
        self.d2_counts = 0

        # Updating/Appending the default option list with user provided options
        self._setOptions(options)

        # Setting up the folders for saving the results
        self._setDirectory()

    def _checkOptionsForMultiAnalysis(self, options):
        """
            This method validates user provided options for type = "multi".
        """

        # Creating list of various different options
        defaultOptions = list(self.options.keys())
        requiredOptions = ["numberOfSamples", "samplingMethod"]
        allowedUserOptions = defaultOptions
        allowedUserOptions.extend(requiredOptions)

        userProvidedOptions = list(options.keys())

        # Checking if user provided option contains only allowed attributes
        if not set(userProvidedOptions).issubset(allowedUserOptions):
            self._error("Option dictionary contains unrecognized attribute(s).")

        # Checking if user has mentioned all the requried attributes
        if not set(requiredOptions).issubset(userProvidedOptions):
            self._error("Option dictionary doesn't contain all the requried options. \
                        {} attribute(s) is/are missing.".format(set(requiredOptions) - set(userProvidedOptions)))

        # Validating number of samples attribute
        if type(options["numberOfSamples"]) is not int:
            self._error("\"numberOfSamples\" attribute is not an integer.")
        
        # Setting minimum limit on number of samples
        if options["numberOfSamples"] < 2:
            self._error("Number of samples need to be at-least 2.")

        # Validating sampling method
        if options["samplingMethod"] not in ["lhs", "fullfactorial"]:
            self._error("\"samplingMethod\" attribute is not correct. \"lhs\" and \"fullfactorial\" are only allowed.")

        # Validating directory attribute
        if "directory" in userProvidedOptions:
            if type(options["directory"]) is not str:
                self._error("\"directory\" attribute is not string.")

    def generateSamples(self):
        """
            Method to generate samples and save the data for further use.
        """

        if self.options["type"] != "multi":
            self._error("You cannot call generateSamples() method when type is not \"multi\".")

        # Generating x based on user provided method
        if self.options["samplingMethod"] == "lhs":
            self._lhs()
        elif self.options["samplingMethod"] == "fullfactorial":
            self._fullfactorial()
        else:
            self._error("Sampling method is not recognized.")

        for index, sample in enumerate(self.x):
            if index == 0:
                self.y = self._function(sample).reshape(1,-1)
            else:
                result = self._function(sample).reshape(1,-1)
                self.y = np.concatenate((self.y, result))

        data = {"x" : self.x, "y" : self.y }

        print("Total Discipline 1 counts: " + str(self.d1_counts))
        print("Total Discipline 2 counts: " + str(self.d2_counts))

        # Saving data file in the specified folder
        os.chdir(self.options["directory"])
        savemat("data.mat", data)
        os.chdir("../")

    def _lhs(self):
        """
            Method for creating a lhs sample.
            Stores a 2D numpy array of size (samples vs  dimensions).
            Each row represents a new sample and each column corresponds to
            a particular design variable.
        """

        lowerBound = np.array(self.options["lowerBound"])
        upperBound = np.array(self.options["upperBound"])

        dim = len(self.options["lowerBound"])

        samples = lhs(dim, samples=self.options["numberOfSamples"], criterion='cm', iterations=100)

        self.x = lowerBound + (upperBound - lowerBound) * samples

    def _fullfactorial(self):
        """
            Method to create fullfactorial samples
        """

        lowerBound = np.array(self.options["lowerBound"])
        upperBound = np.array(self.options["upperBound"])

        dim = len(self.options["lowerBound"])

        samplesInEachDimension = round(np.exp( np.log(self.options["numberOfSamples"]) / dim ))

        print("{} full-factorial samples are generated".format(samplesInEachDimension**dim))

        samples = fullfact([samplesInEachDimension]*dim)

        self.x = lowerBound + samples * (upperBound - lowerBound) / (samplesInEachDimension- 1)

    # ----------------------------------------------------------------------------
    #                      All the methods for single analysis.
    # ----------------------------------------------------------------------------

    def _setupSingleAnalysis(self):
        """
            Method to setup object for single analysis
        """
        
        self.options = {}
        self.options["type"]  = "single"

        # There is no other special initialization needed for 
        # Forrester function when the type is "single".

    def getObjectives(self, x):
        """
            Method to calculate y, g1, g2 based on input x.
            Input x should be a numpy array with three integer entries. 
            Output will be a numpy array of y, g1, g2 in order.
        """

        self.d1_counts = 0
        self.d2_counts = 0

        if self.options["type"] != "single":
            self._error("You cannot call getObjectives() method when type is not \"single\".")

        # Validating x provided by the user
        if type(x) != np.ndarray:
            self._error("Provided x is not a numpy array.")
        elif len(x) != 3:
            self._error("Provided x doesn't contain three entries.")

        for index, number in enumerate(x):
            if type(number) != np.float64:
                self._error("{} entry in x is not a float.".format(index+1))

        return self._function(x), self.d1_counts, self.d2_counts

    # ----------------------------------------------------------------------------
    #          Other required methods, irrespective of type of analysis.
    # ----------------------------------------------------------------------------

    def _function(self, x):
        """
            Sasena function. Note: Output of function should always be
            of size num_samples X 4 (f and three constraints).
        """

        pkgdir = sys.modules["datgen"].__path__[0]
        filepath = os.path.join(pkgdir, "runscripts/runscript_sellar.py")
        shutil.copy(filepath, ".")

        input = {}
        input["x"] = x

        filehandler = open("input.pickle", "xb")
        pickle.dump(input, filehandler)
        filehandler.close()

        os.system("python runscript_sellar.py")

        output = {}

        filehandler = open("output.pickle", "rb")
        output = pickle.load(filehandler)
        filehandler.close()

        os.system("rm -r input.pickle output.pickle runscript_sellar.py reports")

        self.d1_counts += output["d1_counts"]
        self.d2_counts += output["d2_counts"]

        return output["y"]
    
    def _setDirectory(self):
        """
            Method for setting up directory
        """

        directory = self.options["directory"]

        if not os.path.isdir(directory):
            os.system("mkdir {}".format(directory))
        else:
            os.system("rm -r {}".format(directory))
            os.system("mkdir {}".format(directory))
