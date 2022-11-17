# Importing python packages
import os, sys, shutil, pickle
import numpy as np
from scipy.io import savemat
from ..base import BaseClass

class Sellar(BaseClass):
    """
        Class contains essential methods for generating data
        of Sellar MDA problem (refer the main sellar paper):

            f(x1, x2, x3, y1, y2) = x1^2 + x3 + y1 + exp(-y2)
            y1(x1, x2, x3, y2) = x1^2 + x2 + x3 - 0.2y2
            y2(x1, x3, y1) = sqrt(y1) + x1 + x3
            g1(y1) = 1 - (y1/8.0) <= 0
            g2(y2) = (y2/10.0) - 1 <= 0

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
        """
            Initialization.
        """

        self._initialization(type, options)

    # ----------------------------------------------------------------------------
    #                   All the methods for multi analysis
    # ----------------------------------------------------------------------------

    def _setupMultiAnalysis(self, options):
        """
            Method to perform initialization when the type is "multi".
        """

        # Creating an empty options dictionary and assign value of type
        self.options = {}
        
        # Setting up default options
        self._getDefaultOptions()

        # Creating list of required options
        requiredOptions = ["numberOfSamples", "samplingMethod"]

        # Validating user provided options
        self._checkOptionsForMultiAnalysis(options, requiredOptions)

        self.options["type"] = "multi"
        self.options["lowerBound"] = [-10, 0, 0]
        self.options["upperBound"] = [10, 10, 10]
        self.d1_counts = 0
        self.d2_counts = 0

        # Updating/Appending the default option list with user provided options
        self._setOptions(options)

        # Setting up the folders for saving the results
        self._setDirectory()

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

        print("{} {} samples generated".format(self.options["numberOfSamples"], self.options["samplingMethod"]))
        print("Total Discipline 1 counts: " + str(self.d1_counts))
        print("Total Discipline 2 counts: " + str(self.d2_counts))

        # Saving data file in the specified folder
        os.chdir(self.options["directory"])
        savemat("data.mat", data)
        os.chdir("../")

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
            Sellar function. Note: Output of function should always be
            of size num_samples X 4 (f and three constraints).
        """

        pkgdir = sys.modules["blackbox"].__path__[0]
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
