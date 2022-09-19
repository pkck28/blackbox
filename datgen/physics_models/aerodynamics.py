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
        Class creates a default option list (for solver and other settings)
        which is later edited/appended with user provided options.
    """

    def __init__(self):

        # Setting up the design variable, parameter, and objectives
        self.designVariables = {}
        self.parameters = {}
        self.objectives = ["cl", "cd", "lift", "drag"]

        # Sampling options
        self.samplingMethod = "lhs"
        self.numberOfSamples = 2

        # Solver Options
        self.aeroSolver = "adflow"
        self.aeroSolverOptions = {}
        if self.aeroSolver == "adflow":
            self.aeroSolverOptions = {
                "printAllOptions": False,
                "printIntro": False,
                "outputDirectory": "."
            }

        # Other options
        self.directory = "output"
        self.noOfProcessors = 4

class Aerodynamics():
    """
        Class contains essential methods for generating data in two different
        ways. One is for generating data in bulk and other is for running single analysis.
    """
    
    def __init__(self, type="multi", options=None):

        if type == "multi":
            self._setupMultiAnalysis(options)
        elif type == "single":
            self._setupSingleAnalysis(options)
        else:
            self._error("Value of type argument not recognized.")

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
            else:
                self._error("Options argument not provided.")

            # Creating an empty options dictionary
            self.options = {}
            self.options["type"] = "multi"

            # Setting up default options
            self._getDefaultOptions()

            # Updating/Appending the default option list with user provided options
            self._setOptions(options)

            # Setting up the folders for saving the results
            self._setDirectory()

    def generateSamples(self):
        """
            Method to generate multiple samples and save the data for further use.
        """

        if self.options["type"] != "multi":
            self._error("You can run generateSamples() method only when type is multi.")

        if self.options["samplingMethod"] == "lhs":
            self._lhs()
        elif self.options["samplingMethod"] == "fullfactorial":
            self._fullfactorial()
        else:
            self._error("Unrecognized sampling technique provided. Only two techniques are available: \"lhs\" and \"fullfactorial\".")

        # self._createInputFile()

        # y = {}

        # for value in self.options["objectives"]:
        #     y[value] = np.array([])

        # for sampleNo in range(self.options["numberOfSamples"]):
        #     os.chdir("{}/{}".format(self.options["directory"],sampleNo))
        #     print("Running analysis {} of {}".format(sampleNo + 1, self.options["numberOfSamples"]))
        #     os.system("mpirun -n {} python runscript_aerodynamics.py >> log.txt".format(self.options["noOfProcessors"]))
        #     os.system("rm -r input.pickle runscript_aerodynamics.py reports")

        #     filehandler = open("output.pickle", 'rb')
        #     output = pickle.load(filehandler)

        #     for value in self.options["objectives"]:
        #         y[value] = np.append(y[value], output[value])

        #     os.system("rm -r output.pickle {}".format(self.options["aeroSolverOptions"]["gridFile"]))
        #     os.chdir("../..")

        # for index, value in enumerate(self.options["objectives"]):
        #     if index == 0:
        #         Y = y[value].reshape(-1,1)
        #     else:
        #         Y = np.concatenate((Y, y[value].reshape(-1,1)), axis=1)

        # data = {'x' : self.x, 'y' : Y}

        # os.chdir("{}".format(self.options["directory"]))
        # savemat("data.mat", data)
        # os.chdir("../")

    def _lhs(self):
        """
            Method for creating a lhs sample.
            Stores a 2D numpy array of size (samples vs  dimensions).
            Each row represents a new sample and each column corresponds to
            a particular design variable.
        """

        # lower and upper bound are created as numpy arrays and dummy is a normal list here.
        lowerBound = np.array([])
        upperBound = np.array([])
        dummy = np.array([])
        self.samples = {}

        for key in self.options["designVariables"]:
            if key == "aoa" or "mach":
                lowerBound = np.append(lowerBound, self.options["designVariables"][key]["lowerBound"])
                upperBound = np.append(upperBound, self.options["designVariables"][key]["upperBound"])
            else:
                self._error("Unrecognized design variable")

            self.samples[key] = np.array([])
            dummy = np.append(dummy, key)

        dim = len(lowerBound)

        samples = lhs(dim, samples=self.options["numberOfSamples"], criterion='cm', iterations=50)

        samples = lowerBound + samples * (upperBound - lowerBound)

        self.x = samples

        for sampleNo in range(self.options["numberOfSamples"]):
            sample = samples[sampleNo,:]

            for key in self.options["designVariables"]:
                self.samples[key] = np.append(self.samples[key], sample[(dummy == key)])

    def _fullfactorial(self):
        """
            Method for creating a full-factorial sample.
        """

        # lower and upper bound are created as numpy arrays and dummy is a normal list here.
        lowerBound = np.array([])
        upperBound = np.array([])
        dummy = np.array([])
        self.samples = {}

        for key in self.options["designVariables"]:
            if key == "aoa" or "mach":
                lowerBound = np.append(lowerBound, self.options["designVariables"][key]["lowerBound"])
                upperBound = np.append(upperBound, self.options["designVariables"][key]["upperBound"])
            else:
                self._error("Unrecognized design variable")

            self.samples[key] = np.array([])
            dummy = np.append(dummy, key)

        dim = len(lowerBound)

        samplesInEachDimension = round(math.exp( math.log(self.options["numberOfSamples"]) / dim ))

        print("{} full-factorial samples are generated".format(samplesInEachDimension**dim))

        samples = fullfact([samplesInEachDimension]*dim)

        samples = lowerBound + samples * (upperBound - lowerBound) / (samplesInEachDimension- 1)

        self.x = samples

        for sampleNo in range(self.options["numberOfSamples"]):
            sample = samples[sampleNo,:]

            for key in self.options["designVariables"]:
                self.samples[key] = np.append(self.samples[key], sample[(dummy == key)])

        print(self.samples)


    # ----------------------------------------------------------------------------
    #               All the methods for single analysis
    # ----------------------------------------------------------------------------

    def _setupSingleAnalysis(self, options):

        # If 'options' is None, notify the user
            if options is not None:
                if not isinstance(options, dict):
                    self._error("The 'options' argument provided is not a dictionary.")
                elif options == {}:
                    self._error("The 'options' argument provided is an empty dictionary.")
            else:
                self._error("Options argument not provided.")

            # Forcing the name of directory irrespective of user input
            # options["directory"] = "additional_samples"

            # Creating an empty options dictionary
            self.options = {}
            self.options["type"] = "single"
            self.sampleNo = 0

            if "numberOfSamples" in options:
                self._error("Number of samples option should not be used or single analysis.")

            # Setting up default options
            self._getDefaultOptions()

            # Updating/Appending the default option list with user provided options
            self._setOptions(options)

            # Setting up the folder for saving the result
            directory = self.options["directory"]

            if not os.path.isdir(directory):
                os.system("mkdir {}".format(directory))
            else:
                os.system("rm -r {}".format(directory))
                os.system("mkdir {}".format(directory))

    def getObjectives(self, user_sample):
        """
            Method for running a single sample analysis.
        """

        # Checking whether the type is consistent
        if self.options["type"] != "single":
            self._error("You can run getObjectives() method only when type is single.")

        # Checking if the sample provided by user is consistent
        noOfDV = len(self.options["designVariables"].keys())

        if type(user_sample) == list:
            self._error("Sample provided by user is not a list")

        if len(user_sample) != noOfDV:
            self._error("No of values provided by user is not matching the number of design variables.")

        # Creating input dictionary for this analysis

        
        # input = {}
        # sample = {}
        # parameters = self.options["parameters"]
        # objectives = self.options["objectives"]

        

        # for key in self.options["designVariables"]:
        #     sample[key] = 

        # input = {
        #     "aeroSolverOptions" : self.options["aeroSolverOptions"],
        #     "sample" : sample,
        #     "parameters" : parameters,
        #     "objectives" : objectives
        # }
        
        # directory = self.options["directory"]
        # self.sampleNo = self.sampleNo + 1

        # os.system("mkdir {}/{}".format(directory, self.sampleNo))
        # pkgdir = sys.modules["datgen"].__path__[0]
        # filepath = os.path.join(pkgdir, "runscripts/runscript_aerodynamics.py")
        # shutil.copy(filepath, "{}/{}".format(directory, self.sampleNo))
        # os.system("cp -r {} {}/{}".format(self.options["aeroSolverOptions"]["gridFile"], directory, self.sampleNo))

        

        # os.chdir("{}/{}".format(directory, self.sampleNo))
        # filehandler = open("input.pickle", "xb")
        # pickle.dump(input, filehandler)
        # filehandler.close()
        # os.chdir("../..")

    def _createSample(self):
        """
            Method for creating a sample for single analysis.
        """

        dummy = np.array([])
        self.samples = {}

        for key in self.options["designVariables"]:

            if key == "aoa" or "mach":
                lowerBound = np.append(lowerBound, self.options["designVariables"][key]["lowerBound"])
                upperBound = np.append(upperBound, self.options["designVariables"][key]["upperBound"])
            else:
                self._error("Unrecognized design variable")

            self.samples[key] = np.array([])
            dummy = np.append(dummy, key)

        dim = len(lowerBound)

        samples = lhs(dim, samples=self.options["numberOfSamples"], criterion='cm', iterations=50)

        samples = lowerBound + samples * (upperBound - lowerBound)

        self.x = samples

        for sampleNo in range(self.options["numberOfSamples"]):
            sample = samples[sampleNo,:]

            for key in self.options["designVariables"]:
                self.samples[key] = np.append(self.samples[key], sample[(dummy == key)])

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
            Method for assigning user provided options.
            This method should be called only after setting
            up a default option dictionary.
        """

        # List of allowed design variables and parameters
        dvlist = ["aoa", "mach", "altitude"]
        parameters = ["aoa", "mach", "altitude", "areaRef", "chordRef"]

        ########################### Design variable checking
        if "designVariables" not in options or options["designVariables"] == {}:
            self._error("Design variable option is not provided.")
        
        if type(options["designVariables"]) == dict:
            if not set(options["designVariables"].keys()).issubset(set(dvlist)):
                self._error("One or more design variable(s) is not recognized.")
        else:
            self._error("Design variable option is not a dictionary.")

        # Removing design variables from parameters
        parameters = list(set(parameters) - set(options["designVariables"].keys()))
        
        ########################### Parameter checking
        if "parameters" not in options or options["parameters"] == {}:
            self._error("Parameter option is not provided.")

        if type(options["parameters"]) == dict:
            if not set(options["parameters"].keys()) == set(parameters):
                self._error("One or more parameter(s) are not provided or are part of design variable. Only these parameters are required: {}".format(parameters))
        else:
            self._error("Parameter option is not a dictionary.")

        ########################### Objectives checking
        if "objectives" not in options or options["objectives"] == []:
            self._error("Objectives option is not provided.")

        if type(options["objectives"]) == list:
            if not set(options["objectives"]).issubset(set(self.options["objectives"])):
                self._error("One or more objective(s) are not valid. Only these objectives are supported: {}".format(self.options["objectives"]))
        else:
            self._error("Objective option is not a list.")

        ########################### Other Checks

        if "outputDirectory" in options["aeroSolverOptions"]:
            self._error("You cannot input directory name in aerodynamics solver options dictionary. \
                         Mention it outside the solver options, using \"directory\" option.")

        ########################### Checking whether the other provided options are valid
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

        for sampleNo in range(self.options["numberOfSamples"]):
            os.system("mkdir {}/{}".format(directory,sampleNo))
            pkgdir = sys.modules["datgen"].__path__[0]
            filepath = os.path.join(pkgdir, "runscripts/runscript_aerodynamics.py")
            shutil.copy(filepath, "{}/{}".format(directory,sampleNo))
            os.system("cp -r {} {}/{}".format(self.options["aeroSolverOptions"]["gridFile"],directory,sampleNo))

    def _createInputFile(self):
        """
            Method to create an input file for analysis
        """

        directory = self.options["directory"]

        if self.options["type"] == "multi":
            for sampleNo in range(self.options["numberOfSamples"]):
                os.chdir("{}/{}".format(directory,sampleNo))

                input = {}
                sample = {}
                parameters = self.options["parameters"]
                objectives = self.options["objectives"]

                for key in self.samples:
                    sample[key] = self.samples[key][sampleNo]

                input = {
                    "aeroSolverOptions" : self.options["aeroSolverOptions"],
                    "sample" : sample,
                    "parameters" : parameters,
                    "objectives" : objectives
                }

                filehandler = open("input.pickle", "xb")
                pickle.dump(input, filehandler)
                filehandler.close()
                os.chdir("../..")

        elif self.options["type"] == "single":
            os.chdir("{}/{}".format(directory,sampleNo))

            input = {}
            sample = {}
            parameters = self.options["parameters"]
            objectives = self.options["objectives"]

            for key in self.samples:
                sample[key] = self.samples[key][sampleNo]

            input = {
                "aeroSolverOptions" : self.options["aeroSolverOptions"],
                "sample" : sample,
                "parameters" : parameters,
                "objectives" : objectives
            }

            filehandler = open("input.pickle", "xb")
            pickle.dump(input, filehandler)
            filehandler.close()

            os.chdir("../..")
            

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
