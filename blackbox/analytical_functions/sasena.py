# Importing python packages
import numpy as np
from ..base import BaseClass

class Sasena(BaseClass):
    """
        Class contains essential methods for generating data
        of sasena Function:

            y(x1,x2) = -(x1 - 1)^2 - (x2 - 0.5)^2 
            g1(x1,x2) = (x1 - 3)^2 + (x2 + 2)^2 exp(-x2^7) - 12 <= 0
            g2(x1,x2) = 10x1 + x2 - 7 <= 0
            g3(x1,x2) = (x1 - 0.5)^2 + (x2 - 0.5)^2 - 0.2 <= 0

        0 <= x1, x2 <= 1

        There are two values possible for type: "single" and "multi" (default). For
        "multi", following is the list of possible attributes:

        "directory" : Folder name where the data.mat file will be saved (string, optional).
        "numberOfSamples" : number of samples to be generated (integer).
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

        # Setting up default options
        self._getDefaultOptions()

        # Creating list of required options
        requiredOptions = ["numberOfSamples", "samplingMethod"]

        # Validating user provided options
        self._checkOptionsForMultiAnalysis(options, requiredOptions)

        self.options["type"] = "multi"
        self.options["lowerBound"] = [0, 0]
        self.options["upperBound"] = [1, 1]

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
            Method to calculate y, g1, g2, and g3 based on input x. Input 
            x should be a numpy array of size (2,). Output will be a numpy 
            array of size (4,) containing y, g1, g2, and g3 in order.
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
            Sasena function. Note: Output of function should always be
            of size num_samples X 4 (f and three constraints).
        """

        if self.options["type"] == "multi":
            x1 = x[:, 0].reshape(-1,1)
            x2 = x[:, 1].reshape(-1,1)

            y = -(x1 - 1)**2 - (x2 - 0.5)**2 
            g1 = (x1 - 3)**2 + (x2 + 2)**2 * np.exp(-x2**7) - 12
            g2 = 10*x1 + x2 - 7
            g3 = (x1 - 0.5)**2 + (x2 - 0.5)**2 - 0.2

            y = np.concatenate((y, g1, g2, g3), axis = 1)

        elif self.options["type"] == "single":
            x1 = x[0]
            x2 = x[1]
            y = [0, 0, 0, 0]

            y[0] = -(x1 - 1)**2 - (x2 - 0.5)**2 
            y[1] = (x1 - 3)**2 + (x2 + 2)**2 * np.exp(-x2**7) - 12
            y[2] = 10*x1 + x2 - 7
            y[3] = (x1 - 0.5)**2 + (x2 - 0.5)**2 - 0.2

        return y
