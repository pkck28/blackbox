import os, sys, shutil, math
import numpy as np
from pyDOE2 import lhs, fullfact
from scipy.io import savemat

class DefaultOptions():
    """
        Class creates a default option list which is later 
        edited/appended with user provided options. This class
        can be used as a parent class for application specific
        default class. No need to define instance variables
        already defined here.
    """

    def __init__(self):
        self.directory = "output"

class BaseClass():
    """
        Base class is for providing methods which are 
        frequently used by many other class (a way to reduce redundancy). 
        Inherit this base class and use the methods provided. Methods can also
        be over-ridden in the child class or new methods can also be created.
    """

    def _error(self, message):
        """
            Method for printing errors in nice manner.
        """

        msg = "\n+" + "-" * 78 + "+" + "\n" + "| Datgen Error: "
        i = 19
        for word in message.split():
            if len(word) + i + 1 > 78:  # Finish line and start new one
                msg += " " * (82 - i) + "|\n| " + word + " "
                i = 1 + len(word) + 1
            else:
                msg += word + " "
                i += len(word) + 1
        msg += " " * (78 - i) + "|\n" + "+" + "-" * 78 + "+" + "\n"
 
        print(msg, flush=True)

        exit()

    def _initialization(self, type, options):
        """
            Method for performing initialization.
        """
        
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
            # If 'options' is not None, then raise an error
            if options is not None:
                self._error("Options argument is not needed when type is \"single\".")
            self._setupSingleAnalysis()

        else:
            self._error("Value of type argument not recognized.")

    def _getDefaultOptions(self, defaultOptions=DefaultOptions()):
        """
            Setting up the initial values of options.
        """

        for key in vars(defaultOptions):
            value = getattr(defaultOptions, key)
            self.options[key] = value

    def _setOptions(self, options):
        """
            Method for assigning user provided options.
        """

        for key in options.keys():
            # update strategy is different for a dictionary
            # and other type of variables
            if isinstance(options[key], dict):
                # if the option is already present in default dictionary, 
                # then add the user provided key-value pairs to the default dictionary.
                if key in self.options.keys():
                    self.options[key].update(options[key]) 
                else:
                    self.options[key] = options[key]
            else:
                self.options[key] = options[key]

    def _setDirectory(self):
        """
            Method for setting up directory.
        """

        directory = self.options["directory"]

        if not os.path.isdir(directory):
            os.system("mkdir {}".format(directory))
        else:
            os.system("rm -r {}".format(directory))
            os.system("mkdir {}".format(directory))

        # if self.options["type"] == "multi":
        #     FFD = False
        #     for key in self.options["varyingParameters"].keys():
        #         if key == "twist" or key == "shape":
        #             FFD = True
        #             break

        #     for sampleNo in range(self.options["numberOfSamples"]):
        #         os.system("mkdir {}/{}".format(directory,sampleNo))
        #         pkgdir = sys.modules["datgen"].__path__[0]
        #         filepath = os.path.join(pkgdir, "runscripts/runscript_aerodynamics.py")
        #         shutil.copy(filepath, "{}/{}".format(directory,sampleNo))
        #         os.system("cp -r {} {}/{}/grid.cgns".format(self.options["aeroSolverOptions"]["gridFile"],directory,sampleNo))

        #         # copying ffd file
        #         if FFD:
        #             pkgdir = sys.modules["datgen"].__path__[0]
        #             filepath = os.path.join(pkgdir, "runscripts/deform_mesh.py")
        #             shutil.copy(filepath, "{}/{}".format(directory,sampleNo))
        #             os.system("cp -r {} {}/{}/ffd.xyz".format(self.options["ffdFile"],directory,sampleNo))

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

        samples = lhs(dim, samples=self.options["numberOfSamples"], criterion='cm', iterations=100)

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

    # ----------------------------------------------------------------------------
    #       All the abstract methods which child class needs to implement.
    # ----------------------------------------------------------------------------

    def _setupMultiAnalysis(self):
        """
            Method for performing multi analysis initialization.
            Developer needs to write an implementation for this method.
        """

        raise NotImplementedError("Please implement setupMultiAnalysis method")

    def _setupSingleAnalysis(self):
        """
            Method for performing single analysis initialization.
            Developer needs to write an implementation for this method.
        """

        raise NotImplementedError("Please implement setupSingleAnalysis method")

    def getObjectives(self):
        """
            Method for performing single analysis.
            Developer needs to write an implementation for this method.
        """

        raise NotImplementedError("Please implement getObjectives method")
    