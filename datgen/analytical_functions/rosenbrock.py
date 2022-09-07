# Importing python packages
import os
import shutil
import sys
import numpy as np
from pyDOE2 import lhs, fullfact
import pickle
from scipy.io import savemat
import math

class DefaultOptions():
    """
        Class creates a default option list which is later 
        edited/appended with user provided options.
    """

    def __init__(self):
        self.samplingMethod = "lhs"
        self.numberOfSamples = 2
        self.directory = "output"
        self.lowerBound = [-2, -2]
        self.upperBound = [2, 2]
        self.parameters = {
            "a" : 1,
            "b" : 100
        }

class Rosenbrock():
    """
        Class contains essential methods for generating data
        of Rosenbrock Function:

            y(x1,x2) = (a - x1)^2 + b*(x2-x1^2)^2

        User can generate data without providing any options.
        Data will be generated based on default values shown below:

            Lower Bound: [-2, -2]
            Upper Bound: [2, 2]
            Number of Samples: 2
            Directory: output
            Sampling Method: LHS
            Parameters: a = 1, b = 100
        
        Provide an appropriate options dict to change these options.
    """

    def __init__(self, type="single", options=None):

        if type == "batch":
            self._setupBatchAnalysis(options)

        elif type == "single":
            self._setupSingleAnalysis(options)

        else:
            self._error("Value of type argument not recognized.")

    # ----------------------------------------------------------------------------
    #               All the methods for batch analysis
    # ----------------------------------------------------------------------------

    def _setupBatchAnalysis(self, options):

            # If 'options' is None, notify the user
            if options is not None:
                if not isinstance(options, dict):
                    self._error("The 'options' argument provided is not a dictionary.")
                elif options == {}:
                    self._error("The 'options' argument provided is an empty dictionary.")

            # Creating an empty options dictionary
            self.options = {}

            # Setting up default options
            self._getDefaultOptions()

            # Updating/Appending the default option list with user provided options
            self._setOptions(options)

            # Setting up the folders for saving the results
            self._setDirectory()

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
            Method for checking and assigning user provided options.
        """

        if "lowerBound" in options.keys() or "upperBound" in options.keys():
            self._verifyBounds(options)
        
        # Checking whether the other provided options are valid
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

    def _verifyBounds(self, options):
        """
            Method for checking bounds provided by user.
        """

        if "lowerBound" not in options: 
            self._error("Lower bound option is not provided")

        if not type(options["lowerBound"]) == list:
            self._error("Lower bound option is not a list")

        if not len(options["lowerBound"]) == 2:
            self._error("Two entries in lower bounds list are need")

        if "upperBound" not in options: 
            self._error("Upper bound option is not provided")

        if not type(options["upperBound"]) == list:
            self._error("Upper bound option is not a list")

        if not len(options["upperBound"]) == 2:
            self._error("Two entries in upper bounds list are need")

        for index, lb in enumerate(options["lowerBound"]):
            if lb >= options["upperBound"][index]:
                self._error("Lower bound for variable {} is greater than upper bound.".format(index+1))

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

    def generateSamples(self):
        """
            Method to generate samples and save the data for further use.
        """

        if self.options["samplingMethod"] == "lhs":
            self._lhs()
        elif self.options["samplingMethod"] == "fullfactorial":
            self._fullfactorial()
        else:
            self._error("Sampling method is not recognized.")

        self.y = self._function(self.samples)

        data = {"x" : self.samples, "y" : self.y }

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

        samples = lhs(dim, samples=self.options["numberOfSamples"], criterion='cm', iterations=50)

        self.samples = lowerBound + (upperBound - lowerBound) * samples

    def _fullfactorial(self):
        """
            Method to create fullfactorial samples
        """

        lowerBound = np.array(self.options["lowerBound"])
        upperBound = np.array(self.options["upperBound"])

        dim = len(self.options["lowerBound"])

        samplesInEachDimension = round(math.exp( math.log(self.options["numberOfSamples"]) / dim ))

        print("{} full-factorial samples are generated".format(samplesInEachDimension**dim))

        samples = fullfact([samplesInEachDimension]*dim)

        self.samples = lowerBound + samples * (upperBound - lowerBound) / (samplesInEachDimension- 1)

    # ----------------------------------------------------------------------------
    #          All the methods for single analysis
    # ----------------------------------------------------------------------------

    def _setupSingleAnalysis(self, options):
        pass

    # ----------------------------------------------------------------------------
    #          Other required methods, irrespective of type of analysis.
    # ----------------------------------------------------------------------------

    def _function(self, x):
        """
            Rosenbrock function. Note: output of function should always be
            of size num_samples X num_features (?), so reshape is used for x1 and x2.
        """

        x1 = x[:, 0].reshape(-1,1)
        x2 = x[:, 1].reshape(-1,1)

        y = (self.options["parameters"]["a"]-x1)**2 + self.options["parameters"]["b"]*(x2-x1**2)**2

        return y

    def _error(self, message):
        """
            Method for printing errors in nice manner.
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
 
        print(msg, flush=True)

        exit()
