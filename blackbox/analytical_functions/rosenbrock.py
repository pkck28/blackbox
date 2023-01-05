# Importing python packages
import numpy as np
from ..base import BaseClass

class Rosenbrock(BaseClass):
    """
        Class contains essential methods for generating data
        of Rosenbrock Function:

            y(x1,x2) = (1 - x1)^2 + 100*(x2-x1^2)^2

        There are two values possible for type: "single" and "multi" (default). For
        "multi", following is the list of possible attributes:

        "directory" : Folder name where the data.mat file will be saved (string, optional).
        "numberOfSamples" : number of samples to be generated (integer).
        "lowerBound" : lower bound (list with two integer entries).
        "upperBound" : upper bound (list with two integer entries).
        "samplingMethod" : name of the sampling method ("lhs" or "fullfactorial") (string).

        For "single", options argument is not need.
    """

    def __init__(self, type="multi", options=None):
        """
            Initialization.
        """

        self._initialization(type, options)

    # ----------------------------------------------------------------------------
    #               All the methods for multi analysis
    # ----------------------------------------------------------------------------

    def _setupMultiAnalysis(self, options):
        """
            Method to perform initialization when the type is "multi".
        """

        # Creating an empty options dictionary and assign value of type
        self.options = {}
        self.options["type"] = "multi"

        # Setting up default options
        self._getDefaultOptions()

        # Creating list of required options
        requiredOptions = ["numberOfSamples", "samplingMethod"]

        # Validating user provided options
        self._checkOptionsForMultiAnalysis(options, requiredOptions)

        self.options["lowerBound"] = [-2, -2]
        self.options["upperBound"] = [2, 2]

        # Updating/Appending the default option list with user provided options
        self._setOptions(options)

        # Setting up the folders for saving the results
        self._setDirectory()

    # ----------------------------------------------------------------------------
    #          All the methods for single analysis
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
            Method to generate a single output y, based on input x.
            Input x should be a numpy array of size (2,).
            Output y will be a float.
        """

        if self.options["type"] != "single":
            self._error("You cannot call getObjectives() method when type is not \"single\".")

        # Validating x provided by the user
        if type(x) != np.ndarray:
            self._error("Provided x is not a numpy array.")
        elif len(x) != 2:
            self._error("Provided x doesn't contain two entries.")

        return self._function(x)

    # ----------------------------------------------------------------------------
    #          Other required methods, irrespective of type of analysis.
    # ----------------------------------------------------------------------------

    def _function(self, x):
        """
            Rosenbrock function. Note: output of function should always be
            of size num_samples X 1, so reshape is used for x1 and x2.
        """

        a = 1
        b = 100

        if self.options["type"] == "multi":
            x1 = x[:, 0].reshape(-1,1)
            x2 = x[:, 1].reshape(-1,1)

        elif self.options["type"] == "single":
            x1 = x[0]
            x2 = x[1]

        y = (a-x1)**2 + b*(x2-x1**2)**2

        return y
