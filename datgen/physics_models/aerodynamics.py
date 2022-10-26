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
        self.ffdFile = None

class Aerodynamics():
    """
        Class contains essential methods for generating aerodynamic data.
        There are two values possible for type: "single" and "multi" (default). 
        For "multi", following is the list of possible attributes:
        
        Optional: (.,.) shows datatype and defualt value respectively.
            "directory" : Folder name where the data.mat file will be saved (string, "output").
            "noOfProcessors" : Number of processors to use (integer, 4).
            "aeroSolver" : Name of the aerodynamics solver (string, "adflow").
        Compulsory: (.) shows datatype.
            "numberOfSamples" : number of samples to be generated (integer).
            "fixedParameters" : Dictionary of all the valid fixed parameters (dict).
            "varyingParameters" : Dictionary of all the valid varying parameters (dict).
            "samplingMethod" : name of the sampling method ("lhs" or "fullfactorial") (string).
            "objectives" : List of desired objectives in y (list of string).
            "aeroSolverOptions" : Dictionary of containing various ADflow options (dict).

        For "single", following is the list of possible attributes:
        
        Optional: (.,.) shows datatype and defualt value respectively.
            "directory" : Folder name where the analysis output will be saved (string, "additional_samples").
            "noOfProcessors" : Number of processors to use (integer, 4).
            "aeroSolver" : Name of the aerodynamics solver (string, "adflow").
        Compulsory: (.) shows datatype.
            "fixedParameters" : Dictionary of all the valid fixed parameters (dict).
            "varyingParameters" : Dictionary of all the valid varying parameters (dict).
            "objectives" : List of desired objectives in y (list of string).
            "aeroSolverOptions" : Dictionary of containing various ADflow options (dict).
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
        self.options["aeroSolverOptions"] = {
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

        # Checking for not allowed aero solver options
        commonElements = set(options["aeroSolverOptions"].keys()).intersection(set(["printAllOptions", "printIntro", "outputDirectory"]))
        if len(commonElements) != 0:
            self._error("Please remove {} attribute from the \"aeroSolverOptions\"".format(commonElements))

        ############ Checking directory
        if "directory" in userProvidedOptions:
            if type(options["directory"]) is not str:
                self._error("\"directory\" attribute is not string.")

        ############ Checking noOfProcessors
        if "noOfProcessors" in userProvidedOptions:
            if type(options["noOfProcessors"]) is not int:
                self._error("\"noOfProcessors\" attribute is not an integer.")

        ############ Checking FFD file exists - This is done in checkParamters() method since this
        ############ attribute is required only when varying parameter contain twist or shape

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

        FFD = False

        for key in self.options["varyingParameters"]:
            if key == "twist" or key == "shape":
                FFD = True

        # Running analysis
        for sampleNo in range(self.options["numberOfSamples"]):
            os.chdir("{}/{}".format(self.options["directory"],sampleNo))
            print("Running analysis {} of {}".format(sampleNo + 1, self.options["numberOfSamples"]))

            # Deform the mesh if DV has twist or shape
            if FFD:
                os.system("python deform_mesh.py >> mesh_deformation_log.txt")
                os.system("rm deform_mesh.py")
            
            # Run the analysis
            os.system("mpirun -n {} --use-hwthread-cpus python runscript_aerodynamics.py >> analysis_log.txt".format(self.options["noOfProcessors"]))
            os.system("mv ap_000_surf.plt aero_analysis_result.plt")
            os.system("rm -r input.pickle runscript_aerodynamics.py reports grid.cgns")

            filehandler = open("output.pickle", 'rb')
            output = pickle.load(filehandler)
            filehandler.close()

            os.system("rm -r output.pickle")
            os.chdir("../..")

            for value in self.options["objectives"]:
                if sampleNo == 0:
                    y[value] = np.array([])
                y[value] = np.append(y[value], output[value])

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

        # Generating lower and upper bound
        for key in self.options["varyingParameters"]:

            if key == "twist" or key == "shape":
                number = self.options["varyingParameters"][key]["numberOfVariables"]
                lowerBound = np.append(lowerBound, [self.options["varyingParameters"][key]["lowerBound"]] * number)
                upperBound = np.append(upperBound, [self.options["varyingParameters"][key]["upperBound"]] * number)
                dummy = np.append(dummy, [key] * number)
            else:
                lowerBound = np.append(lowerBound, self.options["varyingParameters"][key]["lowerBound"])
                upperBound = np.append(upperBound, self.options["varyingParameters"][key]["upperBound"])
                dummy = np.append(dummy, key)

            self.samples[key] = np.array([])

        # Dimension of the variable
        dim = len(lowerBound)

        # Generating samples
        samples = lhs(dim, samples=self.options["numberOfSamples"], criterion='cm', iterations=50)

        # Scaling the samples
        samples = lowerBound + samples * (upperBound - lowerBound)

        # Saving the X in the data
        self.x = samples

        # Extracting the generated samples for each variable and storing them together
        for sampleNo in range(self.options["numberOfSamples"]):
            sample = samples[sampleNo,:]

            for key in self.options["varyingParameters"]:
                self.samples[key] = np.append(self.samples[key], [sample[(dummy == key)]])

        # Reshaping the array so that they can be divided using sample number
        for key in self.samples.keys():
            if key == "twist" or key == "shape":
                number = self.options["varyingParameters"][key]["numberOfVariables"]
                self.samples[key] = self.samples[key].reshape((-1,number))

    def _fullfactorial(self):
        """
            Method for creating a full-factorial sample.
            Stores a 2D numpy array of size (samples vs  dimensions).
            Each row represents a new sample and each column corresponds to
            a particular design variable.
        """

        for key in self.options["varyingParameters"]:
            if key == "twist" or key == "shape":
                self._error("Full-factorial sampling when twist or shape is a DV is not recommended. Please use LHS.")

        # lower and upper bound are created as numpy arrays and dummy is a normal list here.
        lowerBound = np.array([])
        upperBound = np.array([])
        dummy = np.array([])
        self.samples = {}

        # Generating lower and upper bound
        for key in self.options["varyingParameters"]:
            if key == "twist" or key == "shape":
                number = self.options["varyingParameters"][key]["numberOfVariables"]
                lowerBound = np.append(lowerBound, [self.options["varyingParameters"][key]["lowerBound"]] * number)
                upperBound = np.append(upperBound, [self.options["varyingParameters"][key]["upperBound"]] * number)
                dummy = np.append(dummy, [key] * number)
            else:
                lowerBound = np.append(lowerBound, self.options["varyingParameters"][key]["lowerBound"])
                upperBound = np.append(upperBound, self.options["varyingParameters"][key]["upperBound"])
                dummy = np.append(dummy, key)

            self.samples[key] = np.array([])

        # Dimension of the variable
        dim = len(lowerBound)

        # Calculating number of samples in each dimension
        samplesInEachDimension = round(math.exp( math.log(self.options["numberOfSamples"]) / dim ))

        print("{} full-factorial samples are generated".format(samplesInEachDimension**dim))

        # Generating samples
        samples = fullfact([samplesInEachDimension]*dim)

        # Scaling the samples
        samples = lowerBound + samples * (upperBound - lowerBound) / (samplesInEachDimension- 1)

        # Saving the X in the data
        self.x = samples

        # Reshaping the array so that they can be divided using sample number
        for key in self.samples.keys():
            if key == "twist" or key == "shape":
                number = self.options["varyingParameters"][key]["numberOfVariables"]
                self.samples[key] = self.samples[key].reshape((-1,number))

    # ----------------------------------------------------------------------------
    #               All the methods for single analysis
    # ----------------------------------------------------------------------------

    def _setupSingleAnalysis(self, options):

        # Creating an empty options dictionary
        self.options = {}
        self.options["type"] = "single"
        
        # Setting up default options
        self._getDefaultOptions()

        # Changing default directory value for single analysis
        self.options["directory"] = "additional_samples"

        # Sample Number under current instantiation
        self.sampleNo = 0

        # Validating user provided options
        self._checkOptionsForSingleAnalysis(options)

        # Setting some default solver options
        self.options["aeroSolverOptions"] = {
                "printAllOptions": False,
                "printIntro": False,
                "outputDirectory": "."
            }

        # Updating/Appending the default option list with user provided options
        self._setOptions(options)

        # Setting up the folder for saving the result
        self._setDirectory()

    def _checkOptionsForSingleAnalysis(self, options):
        """
            This method validates user provided options for type = "single".
        """

        # Creating list of various different options
        defaultOptions = list(self.options.keys())
        requiredOptions = ["fixedParameters", "varyingParameters", "objectives", "aeroSolverOptions"]
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

        ############ Checking objectives
        self._checkObjectives(options)

        ############ Checking aeroSolverOptions
        if type(options["aeroSolverOptions"]) != dict:
            self._error("\"aeroSolverOptions\" attribute is not a dictionary.")

        # Checking for common elements between fixed and varying parameters
        commonElements = set(options["aeroSolverOptions"].keys()).intersection(set(["printAllOptions", "printIntro", "outputDirectory"]))
        if len(commonElements) != 0:
            self._error("Please remove {} attribute from the \"aeroSolverOptions\"".format(commonElements))

        ############ Checking directory
        if "directory" in userProvidedOptions:
            if type(options["directory"]) is not str:
                self._error("\"directory\" attribute is not string.")

        ############ Checking noOfProcessors
        if "noOfProcessors" in userProvidedOptions:
            if type(options["noOfProcessors"]) is not int:
                self._error("\"noOfProcessors\" attribute is not an integer.")

    def getObjectives(self, user_sample):
        """
            Method for running a single sample analysis.
        """
        
        directory = self.options["directory"]

        # Checking whether the type is consistent
        if self.options["type"] != "single":
            self._error("You can run getObjectives() method only when type is single.")

        # if type(user_sample) != list:
        #     self._error("Sample provided by user is not a list")

        if len(user_sample) != len(self.options["varyingParameters"].keys()):
            self._error("No of values provided by user is not matching the number of design variables.")

        # Creating new directory for the analysis
        os.system("mkdir {}/{}".format(directory, self.sampleNo))

        # Setting up sample variable for single analysis
        self._createSample(user_sample)
        
        # # Creating input file for single analysis
        # self._createInputFile()

        # # Pasting essential files in the directory for running directory
        # pkgdir = sys.modules["datgen"].__path__[0]
        # filepath = os.path.join(pkgdir, "runscripts/runscript_aerodynamics.py")
        # shutil.copy(filepath, "{}/{}".format(directory, self.sampleNo))
        # os.system("cp -r {} {}/{}/grid.cgns".format(self.options["aeroSolverOptions"]["gridFile"], directory, self.sampleNo))

        # # Changing directory and running the analysis
        # os.chdir("{}/{}".format(self.options["directory"], self.sampleNo))
        # print("Running analysis {}".format(self.sampleNo))
        # os.system("mpirun -n {} --use-hwthread-cpus python runscript_aerodynamics.py >> log.txt".format(self.options["noOfProcessors"]))
        # os.system("mv ap_000_surf.plt aero_analysis_result.plt")

        # # Cleaning up the analysis directory
        # os.system("rm -r input.pickle runscript_aerodynamics.py reports grid.cgns")

        # filehandler = open("output.pickle", 'rb')
        # output = pickle.load(filehandler)
        # filehandler.close()

        # os.system("rm -r output.pickle")
        # os.chdir("../..")

        # y = []

        # for value in output:
        #     y.append(output[value][0])

        # self.sampleNo = self.sampleNo + 1

        # return y

    def _createSample(self, user_sample):
        """
            Method for creating a samples dictionary for single analysis.
        """

        self.samples = {}
        dummy = np.array([])

        for key in self.options["varyingParameters"]:
            if key == "twist" or key == "shape":
                number = self.options["varyingParameters"][key]["numberOfVariables"]
                dummy = np.append(dummy, [key] * number)
            else:
                dummy = np.append(dummy, key)

        for key in enumerate(self.options["varyingParameters"]):
            self.samples[key] = user_sample[(dummy == key)]

        print(self.samples)

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
        optionalVaryingParameters = ["twist", "shape"]
        varyingParameters.extend(optionalVaryingParameters)
        parameters.extend(varyingParameters)

        userVaryingParameters = options["varyingParameters"].keys()
        userFixedParameters = options["fixedParameters"].keys()

        # Checking if user provided correct verying parameters
        if not set(userVaryingParameters).issubset(set(varyingParameters)):
            self._error("\"varyingParameters\" dictionary contains attribute(s) which are not allowed vary. \
                Only \"aoa\", \"mach\", \"altitude\", \"twist\", \"shape\" are allowed.")

        # Checking for common elements between fixed and varying parameters
        commonElements = set(userVaryingParameters).intersection(set(userFixedParameters))
        if len(commonElements) != 0:
            self._error("\"fixedParameters\" and \"varyingParameters\" dictionary contains common attributes.")

        # Calculating list of requried fixed parameters
        requiredFixedParameters = list(set(parameters) - set(userVaryingParameters) - set(optionalVaryingParameters))

        # Checking if user provided all the requried fixed parameters based on the varying parameters
        if not set(requiredFixedParameters) == set(userFixedParameters):
            temp = set(requiredFixedParameters) - set(userFixedParameters)
            if len(temp) != 0:
                self._error("\"fixedParameters\" dictionary must contain only the required attribute(s).\
                {} attribute(s) is/are missing.".format(temp))
            temp = set(userFixedParameters) - set(requiredFixedParameters)
            if len(temp) != 0:
                self._error("\"fixedParameters\" dictionary contains extra attribute(s).\
                {} attribute(s) is/are extra.".format(temp))
        
        # Checking varyingParameter dictionary
        for key in options["varyingParameters"]:
            if type(options["varyingParameters"][key]) != dict:
                self._error("Value of " + key + " in \"varyingParameters\" dictionary is not a dictionary.")

            if key == "twist" or key == "shape":
                # Checking for FFD file attribute
                if "ffdFile" not in options.keys():
                    self._error("\"ffdFile\" attribute not provided in options dictionary, this is requried for twist or shape.")
                # Checking if FFD file exists or not
                if not os.path.isfile(options["ffdFile"]):
                    self._error("{} file doesn't exists, please check".format(options["ffdFile"]))
                if set(["lowerBound", "upperBound", "numberOfVariables"]) != set(options["varyingParameters"][key]):
                    self._error(key + " dictionary need to have \"lowerBound\", \"upperBound\", and \"numberOfVariables\" attributes only.")
            else:
                if set(["lowerBound", "upperBound"]) != set(options["varyingParameters"][key]):
                    self._error(key + " dictionary need to have \"lowerBound\" and \"upperBound\" attributes only.")

            lbType = type(options["varyingParameters"][key]["lowerBound"])
            if lbType != int and lbType != float:
                print(type(options["varyingParameters"][key]["lowerBound"]))
                self._error("Value of \"lowerBound\" in " + key + " dictionary is not a number.")
                
            ubType = type(options["varyingParameters"][key]["upperBound"])
            if ubType != int and ubType != float:
                self._error("Value of \"upperBound\" in " + key + " dictionary is not a number.")

            if not options["varyingParameters"][key]["upperBound"] > options["varyingParameters"][key]["lowerBound"]:
                self._error("Value of upper bound in " + key + " dictionary is smaller or equal to lower bound.")

        # Checking fixedParameter dictionary
        for key in options["fixedParameters"]:
            valueType = type(options["fixedParameters"][key])
            if valueType != int and valueType != float:
                self._error("Value of " + key + " in \"fixedParameters\" dictionary is not a number.")

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
            if isinstance(options[key], dict):
                if key in self.options.keys():
                    # If the value is dictionary, update the default dictionary.
                    # Otherwise, assign values.
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

        if self.options["type"] == "multi":
            FFD = False
            for key in self.options["varyingParameters"].keys():
                if key == "twist" or key == "shape":
                    FFD = True
                    break

            for sampleNo in range(self.options["numberOfSamples"]):
                os.system("mkdir {}/{}".format(directory,sampleNo))
                pkgdir = sys.modules["datgen"].__path__[0]
                filepath = os.path.join(pkgdir, "runscripts/runscript_aerodynamics.py")
                shutil.copy(filepath, "{}/{}".format(directory,sampleNo))
                os.system("cp -r {} {}/{}/grid.cgns".format(self.options["aeroSolverOptions"]["gridFile"],directory,sampleNo))
                if FFD:
                    pkgdir = sys.modules["datgen"].__path__[0]
                    filepath = os.path.join(pkgdir, "runscripts/deform_mesh.py")
                    shutil.copy(filepath, "{}/{}".format(directory,sampleNo))
                    os.system("cp -r {} {}/{}/ffd.xyz".format(self.options["ffdFile"],directory,sampleNo))

    def _createInputFile(self):
        """
            Method to create an input file for analysis
        """

        directory = self.options["directory"]

        input = {}
        sample = {}
        parameters = self.options["fixedParameters"]
        objectives = self.options["objectives"]

        if self.options["type"] == "multi":
            for sampleNo in range(self.options["numberOfSamples"]):
                os.chdir("{}/{}".format(directory,sampleNo))

                for key in self.samples:
                    if key == "twist" or key == "shape":
                        sample[key] = self.samples[key][sampleNo,:]
                    else:
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
            os.chdir("{}/{}".format(directory, self.sampleNo))

            input = {
                "aeroSolverOptions" : self.options["aeroSolverOptions"],
                "sample" : self.samples,
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
