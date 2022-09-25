# Importing python packages
import os
import numpy as np
from pyDOE2 import lhs, fullfact
from scipy.io import savemat
import math

class DefaultOptions():
    """
        Class creates a default option list which is later 
        edited/appended with user provided options.
    """

    def __init__(self):
        self.directory = "output"
        self.parameters = {
            "a" : 1,
            "b" : 100
        }


class Rosenbrock():
    """
        Class contains essential methods for generating data
        of Rosenbrock Function:

            y(x1,x2) = (a - x1)^2 + b*(x2-x1^2)^2

        There are two values possible for type: "single" and "multi". For
        "multi", following is the list of possible attributes:

        "directory" : Folder name where the data.mat file will be saved (string, optional).
        "parameters" : Dictionary containing only 'a' and 'b' as attributes (dict, optional).
        "numberOfSamples" : number of samples to be generated (integer).
        "lowerBound" : lower bound (list with two integer entries).
        "upperBound" : upper bound (list with two integer entries).
        "samplingMethod" : name of the sampling method ("lhs" or "fullfactorial") (string).

        For "single", following is the list of possbile attributes:

        "parameters" : Dictionary containing only 'a' and 'b' as attributes (dict, optional).
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
            self._setupSingleAnalysis(options)

        else:
            self._error("Value of type argument not recognized.")

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

        # Validating parameter attribute
        if "parameters" in userProvidedOptions:
            if type(options["parameters"]) is not dict:
                self._error("\"parameters\" attribute is not a dictionary.")

            if options["parameters"] == {}:
                self._error("\"parameters\" attribute is an empty dictionary.")

            for key in options["parameters"].keys():
                if key not in ["a", "b"]:
                    self._error("\"parameters\" dictionary cannot contain attribute other than \"a\" and \"b\".")

                if type(options["parameters"][key]) != int:
                    self._error("\"{}\" attribute in \"parameters\" dictionary is not an integer.".format(key))

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
            if type(lb) is not int:
                self._error("Lower bound for variable {} is not an integer.".format(index+1))

            if type(options["upperBound"][index]) is not int:
                self._error("Upper bound for variable {} is not an integer.".format(index+1))

            if lb >= options["upperBound"][index]:
                self._error("Lower bound for variable {} is greater than or equal upper bound.".format(index+1))

    def _setOptions(self, options):
        """
            Method for checking and assigning user provided options.
        """
        
        # Checking whether the other provided options are valid
        for key in options.keys():
            # If the value is dictionary, update the default dictionary.
            # Otherwise, assign values.
            if isinstance(options[key], dict): 
                self.options[key].update(options[key]) 
            else:
                self.options[key] = options[key]

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
        
        if self.options["type"] != "multi":
            self._error("You cannot call generateSamples() method when type is not \"multi\".")

        # Generating x based on user provided method
        if self.options["samplingMethod"] == "lhs":
            self._lhs()
        elif self.options["samplingMethod"] == "fullfactorial":
            self._fullfactorial()
        else:
            self._error("Sampling method is not recognized.")

        self.y = self._function(self.samples)

        data = {"x" : self.samples, "y" : self.y }

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
        """
            Method to setup object for single analysis
        """

        self.options = {}
        self.options["type"]  = "single"

        # Setting up default options
        self._getDefaultOptions()

        # Removing directory element from the dictionary
        self.options.pop("directory")

        if options != None:
            if not isinstance(options, dict):
                self._error("The 'options' argument provided is not a dictionary.")
            elif options == {}:
                self._error("The 'options' argument provided is an empty dictionary.")

            # Validating user provided options
            self._checkOptionsForSingleAnalysis(options)

            # Updating/Appending the default option list with user provided options
            self._setOptions(options)

    def _checkOptionsForSingleAnalysis(self, options):
        """
            This method validates user provided options for type = "single".
        """

        # Creating list of various different options
        defaultOptions = list(self.options.keys())
        userProvidedOptions = list(options.keys())

        if not set(userProvidedOptions).issubset(set(defaultOptions)):
            self._error("Only \"parameters\" attribute is allowed in the dictionary.")

        if type(options["parameters"]) is not dict:
                self._error("\"parameters\" attribute is not a dictionary.")

        if options["parameters"] == {}:
            self._error("\"parameters\" attribute is an empty dictionary.")

        for key in options["parameters"].keys():
            if key not in ["a", "b"]:
                self._error("\"parameters\" dictionary cannot contain attribute other than \"a\" and \"b\".")

            if type(options["parameters"][key]) != int:
                self._error("\"{}\" attribute is not an integer.".format(key))

    def getObjectives(self, x):
        """
            Method to generate a single output y, based on input x.
            Input x should be a list with two integer entries.
            Output y will be an integer.
        """

        if self.options["type"] != "single":
            self._error("You cannot call getObjectives() method when type is not \"single\".")

        # Validating x provided by the user
        if type(x) != list:
            self._error("Provided x is not a list.")
        elif len(x) != 2:
            self._error("Provided x doesn't contain two entries.")

        for index, number in enumerate(x):
            if type(number) != int:
                self._error("{} entry in x is not an integer.".format(index+1))

        return self._function(x)

    # ----------------------------------------------------------------------------
    #          Other required methods, irrespective of type of analysis.
    # ----------------------------------------------------------------------------

    def _getDefaultOptions(self):
        """
            Setting up the initial values of options.
        """

        defaultOptions = DefaultOptions()

        for key in vars(defaultOptions):
            value = getattr(defaultOptions, key)
            self.options[key] = value

    def _function(self, x):
        """
            Rosenbrock function. Note: output of function should always be
            of size num_samples X num_features (?), so reshape is used for x1 and x2.
        """

        a = self.options["parameters"]["a"]
        b = self.options["parameters"]["b"]

        if self.options["type"] == "multi":
            x1 = x[:, 0].reshape(-1,1)
            x2 = x[:, 1].reshape(-1,1)

        elif self.options["type"] == "single":
            x1 = x[0]
            x2 = x[0]

        y = (a-x1)**2 + b*(x2-x1**2)**2

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
