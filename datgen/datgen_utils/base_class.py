import os
from collections import OrderedDict

from pyDOE2 import lhs, fullfact

class DefaultOptions():
    """
        This class is used to set initial options which are common across all functions.
        To set initial value of options which are specific to a function, use the respective module.
    """

    def __init__(self):

        self.folderName = "output"
        self.case = "aerostruct"


class BaseClass():
    """
        Base class for all the function. All the function classes should inherit from this class.
    """
    
    # ----------------------------------------------------------------------------
    #               Initialization for the class - Constructor function
    # ----------------------------------------------------------------------------
    def __init__(self, options):

        # If 'options' is None, notify the user
        if options is None:
            self._error("The 'options' argument not provided.")
        
        # If 'options' is not a dictionary, notify the user
        if not isinstance(options, dict):
            self._error("The 'options' argument provided is not a dictionary.")

        self.options = {}

        # Setting up initial value of options
        self._getDefaultOptions()

        # Setting up common options provided by user
        self._setOptions(options)

    def _getDefaultOptions(self):
        """
            Setting up the initial values of options which are common across all functions.
        """
        
        defaultOptions = DefaultOptions()

        for key in vars(defaultOptions):
            value = getattr(defaultOptions, key)
            self.options[key] = value

    def _setOptions(self, options):
        """
            Method for assigning user provided options. This methos should be called only after checks.
        """

        for key in options.keys():
            if key in self.options.keys():
                self.options[key] = options[key]
            else:
                self._error(key + " is not a recognized option. Please remove it.")

    def _lhs(self):
        """
            Method for creating a lhs sample.
        """
        lhs(self.num_dim, samples=self.sample_size, criterion='cm', iterations=50)

        pass

    def _fullfactorial(self):
        """
            Method for creating a fullfactorial sample.
        """

        pass

    def _error(self, message):
        """
            Method for printing a user mistake in formatted manner.
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
