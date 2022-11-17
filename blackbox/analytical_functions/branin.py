# Importing python packages
import numpy as np
from ..base import BaseClass

class Branin(BaseClass):
    """
        Class contains essential methods for generating data
        of branin Function:

            y(x1,x2) = a(x2 - bx1^2 + cx1 - r)^2 + s(1 - t)cos(x1) + s

            g1(x1,x2) = 1 - (x1*x2)/(x1'*x2')
            
        a = 1, b = 5.1/(4pi^2), c = 5/pi, r = 6, s = 10, t = 1/(8pi) 

        -5 <= x1 <= 10, 0 <= x2 <= 15, f(x') = 0.397887 with 
        x' = (-pi,12.275), x' = (pi,2.275), x' = (3pi,2.475)
    
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
        # Input is object of DefaultOptions class
        self._getDefaultOptions()

        # Creating list of required options
        requiredOptions = ["numberOfSamples", "samplingMethod"]

        # Validating user provided options
        self._checkOptionsForMultiAnalysis(options, requiredOptions)

        self.options["type"] = "multi"
        self.options["lowerBound"] = [-5, 0]
        self.options["upperBound"] = [10, 15]

        # Updating/Appending the default option list with user provided options
        self._setOptions(options)

        # Setting up the folders for saving the results
        self._setDirectory()

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
            Method to calculate y and g1 based on input x. Input x 
            should be a numpy array of size (2,). Output will be a 
            numpy array of size (2,) containing y and g1 in order.
        """

        if self.options["type"] != "single":
            self._error("You cannot call getObjectives() method when type is not \"single\".")

        # Validating of x provided by the user
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
            Branin function. Note: Output of function should always be
            of size num_samples X 2 (f and constraints).
        """
        
        a = 1
        b = 5.1/(4*np.pi**2)
        c = 5/np.pi
        r = 6
        s = 10
        t = 1/(8*np.pi)

        if self.options["type"] == "multi":
            x1 = x[:, 0].reshape(-1,1)
            x2 = x[:, 1].reshape(-1,1)

            g = 1 - (x1*x2)/(3*np.pi*2.475)
            y = a*(x2 - b*x1**2 + c*x1 - r)**2 + s*(1 - t)*np.cos(x1) + s

            result = np.concatenate((y, g), axis = 1)

        elif self.options["type"] == "single":
            x1 = x[0]
            x2 = x[1]

            g = 1 - (x1*x2)/(3*np.pi*2.475)
            y = a*(x2 - b*x1**2 + c*x1 - r)**2 + s*(1 - t)*np.cos(x1) + s

            result = np.array([y, g])

        return result
