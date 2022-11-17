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

        # Validating user provided options
        self._checkOptionsForMultiAnalysis(options)

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
        requiredOptions = ["numberOfSamples", "lowerBound", "upperBound", "samplingMethod"]
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
            self._error("Number of samples need to least 2.")

        # Validating bounds provided by the user
        self._verifyBounds(options)

        # Validating sampling method
        if options["samplingMethod"] not in ["lhs", "fullfactorial"]:
            self._error("\"samplingMethod\" attribute is not correct. \"lhs\" and \"fullfactorial\" are only allowed.")

        # Validating directory attribute
        if "directory" in userProvidedOptions:
            if type(options["directory"]) is not str:
                self._error("\"directory\" attribute is not string.")

    def _verifyBounds(self, options):
        """
            Method for checking bounds provided by user.
        """

        if not type(options["lowerBound"]) == list:
            self._error("Lower bound option is not a list")

        if not len(options["lowerBound"]) == 2:
            self._error("Two entries in lower bounds list are need")

        if not type(options["upperBound"]) == list:
            self._error("Upper bound option is not a list")

        if not len(options["upperBound"]) == 2:
            self._error("Two entries in upper bounds list are need")

        for index, lb in enumerate(options["lowerBound"]):
            if lb >= options["upperBound"][index]:
                self._error("Lower bound for variable {} is greater than or equal upper bound.".format(index+1))

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
