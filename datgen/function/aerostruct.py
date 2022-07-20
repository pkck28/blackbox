import os

from ..datgen_utils.base_class import BaseClass
from ..datgen_utils.aerostruct_utils import *

class DefaultOptions():

    def __init__():
        pass

class AeroStruct(BaseClass):

    # ----------------------------------------------------------------------------
    #               Initialization for the class - Constructor function
    # ----------------------------------------------------------------------------
    def __init__(self, options=None):

        super().__init__()

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

    # def _getDefaultOptions(self):
    #     print("Hello")