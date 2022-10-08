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

        # Aero solver Options
        self.aeroSolver = "adflow"

        # Other options
        self.directory = "output"
        self.noOfProcessors = 4

class Aerodynamics():
    """
        Class contains essential methods for generating aerodynamic data.
        There are two values possible for type: "single" and "multi" (default). 
        For "multi", following is the list of possible attributes:
        
        Optional: (.,.) shows datatype and defualt value respectively.
            "directory" : Folder name where the data.mat file will be saved (string, "output").
            "noOfProcessors" : Number of processors to use (integer, 4).
            "aeroSolver" : Name of the aerodynamics solver (string, "adflow").
            "structSolver" : Name of the structural solver (string, "tacs").
        Compulsory: (.) shows datatype.
            "numberOfSamples" : number of samples to be generated (integer).
            "fixedParameters" : List of all the valid fixed parameters (list of strings).
            "varyingParameters" : List of all the valid varying parameters (list of strings).
            "lowerBound" : List of lower bound values for varying parameters (list of integers).
            "upperBound" : List of upper bound values for varying parameters (list of integers).
            "samplingMethod" : name of the sampling method ("lhs" or "fullfactorial") (string).
            "objectives" : List of desired objectives in y (list of string).
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

        # Creating an empty options dictionary
        self.options = {}
        self.options["type"] = "multi"

        # Setting up default options
        self._getDefaultOptions()

        # Validating user provided options
        self._checkOptionsForMultiAnalysis(options)

        # Setting some default solver options
        self.aeroSolverOptions = {
                "printAllOptions": False,
                "printIntro": False,
                "outputDirectory": "."
            }

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
        requiredOptions = ["fixedParameters", "varyingParameters", "objectives", \
                            "numberOfSamples", "samplingMethod", "aeroSolverOptions"]
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

        ############ Checking parameters
        self._checkParameters(options)

        ############ Checking bounds
        self._checkBounds(options)

        ############ Checking objectives
        self._checkObjectives(options)

        ############ Checking numberOfSamples
        if type(options["numberOfSamples"]) is not int:
            self._error("\"numberOfSamples\" attribute is not an integer.")

        # Setting minimum limit on number of samples
        if options["numberOfSamples"] < 2:
            self._error("Number of samples need to least 2.")

        ############ Checking samplingMethod
        if options["samplingMethod"] not in ["lhs", "fullfactorial"]:
            self._error("\"samplingMethod\" attribute is not correct. \"lhs\" and \"fullfactorial\" are only allowed.")

        ############ Checking aeroSolverOptions
        if type(options["aeroSolverOptions"]) != dict:
            self._error("\"aeroSolverOptions\" attribute is not a dictionary.")

        ############ Checking directory
        if "directory" in userProvidedOptions:
            if type(options["directory"]) is not str:
                self._error("\"directory\" attribute is not string.")

        ############ Checking noOfProcessors
        if "noOfProcessors" in userProvidedOptions:
            if type(options["directory"]) is not int:
                self._error("\"noOfProcessors\" attribute is not integer.")

    def generateSamples(self):
        """
            Method to generate multiple samples and save the data for further use.
        """

        if self.options["type"] != "multi":
            self._error("You can call generateSamples() method only when type is multi.")

        if self.options["samplingMethod"] == "lhs":
            self._lhs()
        elif self.options["samplingMethod"] == "fullfactorial":
            self._fullfactorial()

        self._createInputFile()

        y = {}

        for value in self.options["objectives"]:
            y[value] = np.array([])

        for sampleNo in range(self.options["numberOfSamples"]):
            os.chdir("{}/{}".format(self.options["directory"],sampleNo))
            print("Running analysis {} of {}".format(sampleNo + 1, self.options["numberOfSamples"]))
            os.system("mpirun -n {} python runscript_aerodynamics.py >> log.txt".format(self.options["noOfProcessors"]))
            os.system("rm -r input.pickle runscript_aerodynamics.py reports")

            filehandler = open("output.pickle", 'rb')
            output = pickle.load(filehandler)

            for value in self.options["objectives"]:
                y[value] = np.append(y[value], output[value])

            os.system("rm -r output.pickle {}".format(self.options["aeroSolverOptions"]["gridFile"]))
            os.chdir("../..")

        for index, value in enumerate(self.options["objectives"]):
            if index == 0:
                Y = y[value].reshape(-1,1)
            else:
                Y = np.concatenate((Y, y[value].reshape(-1,1)), axis=1)

        data = {'x' : self.x, 'y' : Y}

        os.chdir("{}".format(self.options["directory"]))
        savemat("data.mat", data)
        os.chdir("../")

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

        for key in self.options["varyingParameter"]:
            lowerBound = np.append(lowerBound, self.options["varyingParameter"][key]["lowerBound"])
            upperBound = np.append(upperBound, self.options["varyingParameter"][key]["upperBound"])

            self.samples[key] = np.array([])
            dummy = np.append(dummy, key)

        dim = len(lowerBound)

        samples = lhs(dim, samples=self.options["numberOfSamples"], criterion='cm', iterations=50)

        samples = lowerBound + samples * (upperBound - lowerBound)

        self.x = samples

        for sampleNo in range(self.options["numberOfSamples"]):
            sample = samples[sampleNo,:]

            for key in self.options["varyingParameter"]:
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

        for key in self.options["varyingParameter"]:
            lowerBound = np.append(lowerBound, self.options["varyingParameter"][key]["lowerBound"])
            upperBound = np.append(upperBound, self.options["varyingParameter"][key]["upperBound"])

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

            for key in self.options["varyingParameter"]:
                self.samples[key] = np.append(self.samples[key], sample[(dummy == key)])

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

    def _checkParameters(self, options):
        """
            Method to check whether user provided correct parameters.
        """

        # Checking datatype of provided parameter option
        if not type(options["varyingParameters"]) == dict:
            self._error("\"varyingParameters\" option is not a dictionary.")
        if not type(options["fixedParameters"]) == dict:
            self._error("\"fixedParameters\" option is not a dictionary.")    

        # Defining list of parameters
        parameters = ["areaRef", "chordRef"]
        varyingParameters = ["aoa", "mach", "altitude"]
        parameters.extend(varyingParameters)

        userVaryingParameters = options["varyingParameters"].keys()
        userFixedParameters = options["fixedParameters"].keys()

        if not set(userVaryingParameters).issubset(set(varyingParameters)):
            self._error("\"varyingParameters\" dictionary contains attribute(s) which are not allowed vary. \
                Only \"aoa\", \"mach\", \"altitude\" are allowed.")

        requiredFixedParameters = list(set(parameters) - set(userVaryingParameters))

        if not set(requiredFixedParameters) == set(userFixedParameters):
            self._error("Fixed Parameter dictionary doesn't contain all the required attribute(s).\
                {} attribute(s) is/are missing.".format(set(requiredFixedParameters) - set(userFixedParameters)))
        
        for key in options["varyingParameters"]:
            if type(options["varyingParameters"][key]) != dict:
                self._error("Value of " + key + " in \"varyingParameters\" dictionary is not a dictionary.")
            
            if set(["lowerBound", "upperBoubnd"]) == set(options["varyingParameters"][key]):
                self._error(key + " dictionary can only have \"lowerBound\" and \"upperBound\" attributes.")

            if type(options["varyingParameters"][key]["lowerBound"]) != float:
                print(type(options["varyingParameters"][key]["lowerBound"]))
                self._error("Value of \"lowerBound\" in " + key + " dictionary is not a float.")
                
            if type(options["varyingParameters"][key]["upperBound"]) != float:
                self._error("Value of \"upperBound\" in " + key + " dictionary is not a float.")

            if not options["varyingParameters"][key]["upperBound"] > options["varyingParameters"][key]["lowerBound"]:
                self._error("Value of upper bound in " + key + " dictionary is smaller or equal to lower bound.")

        for key in options["fixedParameters"]:
            if type(options["fixedParameters"][key]) != int:
                self._error("Value of " + key + " in \"fixedParameters\" dictionary is not a dictionary.")

    # def _checkBounds(self, options):
    #     """
    #         Method for checking bounds provided by user.
    #     """

    #     dim = len(options["varyingParameters"].keys())

    #     if not type(options["lowerBound"]) == list:
    #         self._error("Lower bound option is not a list")

    #     if not len(options["lowerBound"]) == dim:
    #         self._error("Lower bound list size is not correct.")

    #     if not type(options["upperBound"]) == list:
    #         self._error("Upper bound option is not a list")

    #     if not len(options["upperBound"]) == dim:
    #         self._error("Two entries in upper bounds list are need")

    #     for index, lb in enumerate(options["lowerBound"]):
    #         if type(lb) is not int:
    #             self._error("Lower bound at index location {} is not an integer.".format(index+1))

    #         if type(options["upperBound"][index]) is not int:
    #             self._error("Upper bound at index location {} is not an integer.".format(index+1))

    #         if lb >= options["upperBound"][index]:
    #             self._error("Lower bound at index location {} is greater than or equal upper bound.".format(index+1))

    def _checkObjectives(self, options):
        """
            Checking the objectives provided by the user
        """

        allowedObjectives = ["cl", "cd", "lift", "drag"]

        if type(options["objectives"]) == list:
            if not set(options["objectives"]).issubset(allowedObjectives):
                self._error("One or more objective(s) are not valid. Only these \
                    objectives are supported: {}".format(allowedObjectives))
        else:
            self._error("\"objectives\" option is not a list.")

    def _setOptions(self, options):
        """
            Method for assigning user provided options.
        """

        for key in options.keys():
            # If the value is dictionary, update the default dictionary.
            # Otherwise, assign values.
            if isinstance(options[key], dict): 
                self.options[key].update(options[key]) 
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
                parameters = self.options["fixedParameters"]
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
