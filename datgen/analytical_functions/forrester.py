# Importing python packages
import os
import numpy as np
from scipy.io import savemat

class DefaultOptions():
    """
        Class creates a default option list which is later 
        edited/appended with user provided options.
    """

    def __init__(self):
        self.directory = "output"


class Forrester():
    """
        Class contains essential methods for generating data
        of Forrester Function:

            y(x) = (6x - 2)**2 * sin(12x - 4)

        Data will be generated based on default values shown below.
        Provide an appropriate options dict to change these options.
        Note: There is no sampling method option for forrester function.
    """

    def __init__(self, type="multi", options=None):

        # If 'options' is None, notify the user
        if options is not None:
            if not isinstance(options, dict):
                self._error("The 'options' argument provided is not a dictionary.")
            elif options == {}:
                self._error("The 'options' argument provided is an empty dictionary.")
        else:
            self._error("Options argument not provided.")

        # Initializing based on the type
        if type == "multi":
            self._setupMultiAnalysis(options)
        elif type == "single":
            self._setupSingleAnalysis(options)
        else:
            self._error("Value of type argument is not recognized. Only \"multi\" and \"single\" are allowed.")

    # ----------------------------------------------------------------------------
    #               All the methods for multi analysis
    # ----------------------------------------------------------------------------

    def _setupMultiAnalysis(self, options):

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

        # Validating user provided options
        self._checkOptionsForMultiAnalysis(options)

        # Updating/Appending the default option list with user provided options
        self._setOptions(options)

        # Setting up the folder for saving the results
        self._setDirectory()

    def _checkOptionsForMultiAnalysis(self, options):
        """
            This method validates user provided options for type = "multi".
        """

        # Creating list of various different options
        defaultOptions = list(self.options.keys())
        requiredOptions = ["numberOfSamples", "lowerBound", "upperBound"]
        allowedUserOptions = defaultOptions
        allowedUserOptions.extend(requiredOptions)

        userProvidedOptions = list(options.keys())

        # Checking if the user provided option contains only allowed attributes
        if not set(userProvidedOptions).issubset(allowedUserOptions):
            self._error("Option dictionary contains unrecognized attribute(s).")

        # Checking if user has mentioned all the requried attributes
        if not set(requiredOptions).issubset(userProvidedOptions):
            self._error("Option dictionary doesn't contain all the requried options. \
                        {} attribute(s) is/are missing.".format(set(requiredOptions) - set(userProvidedOptions)))

        # Checking type of required options
        for attribute in requiredOptions:
            if type(options[attribute]) is not int:
                self._error("\"{}\" attribute is not an integer".format(attribute))

        # Setting minimum limit on number of samples
        if options["numberOfSamples"] < 2:
            self._error("Number of samples need to least 2.")

        # Checking correctness of bound
        if options["lowerBound"] >= options["upperBound"]:
            self._error("Lower bound is greater than upper bound.")

    def generateSamples(self):
        """
            Method to generate samples and save the data for further use.
        """

        self.samples = np.linspace(self.options["lowerBound"], self.options["upperBound"], self.options["numberOfSamples"])
        self.samples = np.reshape(self.samples, (-1,1))

        self.y = self._function(self.samples)

        data = {"x" : self.samples, "y" : self.y }

        os.chdir(self.options["directory"])
        savemat("data.mat", data)
        os.chdir("../")

    # ----------------------------------------------------------------------------
    #          All the methods for single analysis
    # ----------------------------------------------------------------------------

    def _setupSingleAnalysis(self, options):

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
        requiredOptions = ["numberOfSamples"]
        allowedUserOptions = defaultOptions
        allowedUserOptions.extend(requiredOptions)

        userProvidedOptions = list(options.keys())

        # Checking if the user provided option contains only allowed attributes
        if not set(userProvidedOptions).issubset(allowedUserOptions):
            self._error("Option dictionary contains unrecognized attribute(s).")

        # Checking if user has mentioned all the requried attributes
        if not set(requiredOptions).issubset(userProvidedOptions):
            self._error("Option dictionary doesn't contain all the requried options. \
                        {} attribute(s) is/are missing.".format(set(requiredOptions) - set(userProvidedOptions)))

        # Checking type of required options
        for attribute in requiredOptions:
            if type(options[attribute]) is not int:
                self._error("\"{}\" attribute is not an integer".format(attribute))

        # Setting minimum limit on number of samples
        if options["numberOfSamples"] < 2:
            self._error("Number of samples need to least 2.")

    def getObjectives(self, x):

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

        if not type(options["lowerBound"]) == int:
            self._error("Lower bound option is not an integer")

        # if not len(options["lowerBound"]) == 2:
        #     self._error("Two entries in lower bounds list are need")

        if "upperBound" not in options: 
            self._error("Upper bound option is not provided")

        if not type(options["upperBound"]) == int:
            self._error("Upper bound option is not an integer")

        # if not len(options["upperBound"]) == 2:
        #     self._error("Two entries in upper bounds list are need")

        if options["lowerBound"] >= options["upperBound"]:
            self._error("Lower bound is greater than upper bound.")

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

    def _function(self, x):
        """
            Forrester function. Note: Output of function should always be
            of size num_samples X 1. Here, input x is already of size num_samples X 1,
            so no need to use reshape.
        """

        y = (6*x - 2)**2 * np.sin(12*x - 4)

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
