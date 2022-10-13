# Importing python packages
import os
import shutil
import sys
import numpy as np
from pyDOE2 import lhs
import pickle
from scipy.io import savemat

class DefaultOptions():
    """
        Class creates a default option list (for solvers and other settings)
        which is later edited/appended with user provided options.
    """

    def __init__(self):

        # Solver Options
        self.aeroSolver = "adflow"
        self.structSolver = "tacs"

        # Other options
        self.directory = "output"
        self.noOfProcessors = 4

class AeroStruct():
    """
        Class contains essential methods for generating aero-struct data.
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

        For "single", following is the list of possible attributes:
        
        Optional: (.,.) shows datatype and defualt value respectively.
            "directory" : Folder name where the analysis output will be saved (string, "additional_samples").
            "noOfProcessors" : Number of processors to use (integer, 4).
            "aeroSolver" : Name of the aerodynamics solver (string, "adflow").
            "structSolver" : Name of the structural solver (string, "tacs").
        Compulsory: (.) shows datatype.
            "fixedParameters" : List of all the valid fixed parameters (list of strings).
            "varyingParameters" : List of all the valid varying parameters (list of strings).
            "lowerBound" : List of lower bound values for varying parameters (list of integers).
            "upperBound" : List of upper bound values for varying parameters (list of integers).
            "objectives" : List of desired objectives in y (list of string).
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
        requiredOptions = ["fixedParameters", "varyingParameters", "aeroObjectives", "structObjectives", \
                            "numberOfSamples", "samplingMethod", "aeroSolverOptions", "structGridFile"]
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

        # Checking for common elements between fixed and varying parameters
        commonElements = set(options["aeroSolverOptions"].keys()).intersection(set(["printAllOptions", "printIntro", "outputDirectory"]))
        if len(commonElements) != 0:
            self._error("Please remove {} attribute from the \"aeroSolverOptions\"".format(commonElements))

        ############ Checking structGridFile option
        if type(options["structGridFile"]) == str:
            self._error("\"structGridFile\" option is not a string.")

        if not os.path.isfile(options["structGridFile"]):
            self._error("{} file doesn't exists, please check".format(options["structGridFile"]))

        ############ Checking directory
        if "directory" in userProvidedOptions:
            if type(options["directory"]) is not str:
                self._error("\"directory\" attribute is not string.")

        ############ Checking noOfProcessors
        if "noOfProcessors" in userProvidedOptions:
            if type(options["noOfProcessors"]) is not int:
                self._error("\"noOfProcessors\" attribute is not an integer.")

    # ----------------------------------------------------------------------------
    #               All the methods for multi analysis
    # ----------------------------------------------------------------------------

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
            Method for assigning user provided options.
            This method should be called only after checks.
        """

        # List of allowed design variables
        dvlist = ["aoa", "mach"]

        # Design variables checking and assignment
        if "designVariables" not in options or type(options["designVariables"]) == {}:
            self._error("Design variable option is not provided.")
        
        if type(options["designVariables"]) == dict:
            if not set(options["designVariables"].keys()).issubset(set(dvlist)):
                self._error("One or more design variable(s) is not recognized.")
        else:
            self._error("Design variable option is not a dictionary.")
        
        # Checking whether the provided option is valid
        for key in options.keys():
            # Checking whether the provided option is valid
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
            Method for setting up directories for analysis
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
            filepath = os.path.join(pkgdir, "runscripts/runscript_aerostruct.py")
            shutil.copy(filepath, "{}/{}".format(directory,sampleNo))
            os.system("cp -r {} {}/{}".format(self.options["aeroSolverOptions"]["gridFile"],directory,sampleNo))
            os.system("cp -r {} {}/{}".format(self.options["structSolverOptions"]["gridFile"],directory,sampleNo))
            os.system("cp -r {} {}/{}".format("tacsSetup.py",directory,sampleNo))

    def generateSamples(self):
        """
            Method to generate samples and save the data for further use.
        """

        if self.options["samplingMethod"] == "lhs":
            self._lhs()

        self._createInputFile()

        cl = np.array([])
        cd = np.array([])
        failure = np.array([])

        for sampleNo in range(self.options["numberOfSamples"]):
            os.chdir("{}/{}".format(self.options["directory"],sampleNo))
            print("Running analysis {} of {}".format(sampleNo + 1, self.options["numberOfSamples"]))
            os.system("mpirun -n 10 python runscript_aerostruct.py >> log.txt")
            os.system("f5totec scenario_000.f5")
            os.system("rm -r input.pickle runscript_aerostruct.py reports")

            filehandler = open("output.pickle", 'rb')
            output = pickle.load(filehandler)
            filehandler.close()

            cl = np.append(cl, output["cl"])
            cd = np.append(cd, output["cd"])
            failure = np.append(cd, output["failure"])

            os.system("rm -r output.pickle tacsSetup.py")
            os.system("rm -r {} {}".format(self.options["aeroSolverOptions"]["gridFile"],
                                                  self.options["structSolverOptions"]["gridFile"]))
            os.chdir("../..")

        data = {'cl': cl, 'cd': cd, 'failure': failure}

        for key in self.samples:
                data[key] = self.samples[key]

        os.chdir("{}".format(self.options["directory"]))
        savemat("data.mat", data)
        os.chdir("../")

    def _lhs(self):
        """
            Method for creating a lhs sample.
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

        samples = lowerBound + (upperBound - lowerBound) * samples

        for sampleNo in range(self.options["numberOfSamples"]):
            sample = samples[sampleNo,:]

            for key in self.options["designVariables"]:
                self.samples[key] = np.append(self.samples[key], sample[(dummy == key)])

    def _createInputFile(self):
        """
            Method to create an input file for analysis
        """

        directory = self.options["directory"]

        for sampleNo in range(self.options["numberOfSamples"]):
            os.chdir("{}/{}".format(directory,sampleNo))

            input = {}
            sample = {}

            for key in self.samples:
                sample[key] = self.samples[key][sampleNo]

            input = {
                "aeroSolverOptions" : self.options["aeroSolverOptions"],
                "structSolverOptions" : self.options["structSolverOptions"],
                "sample" : sample
            }

            filehandler = open("input.pickle", "xb")
            pickle.dump(input, filehandler)
            filehandler.close()

            os.chdir("../..")

    # ----------------------------------------------------------------------------
    #          Other required methods, irrespective of type of analysis.
    # ----------------------------------------------------------------------------

    def _checkObjectives(self, options):
        """
            Checking the objectives provided by the user
        """

        allowedAeroObjectives = ["cl", "cd", "lift", "drag"]
        allowedStructObjectives = []

        if type(options["objectives"]) == list:
            if not set(options["objectives"]).issubset(allowedObjectives):
                self._error("One or more objective(s) are not valid. Only these \
                    objectives are supported: {}".format(allowedObjectives))
        else:
            self._error("\"objectives\" option is not a list.")

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
 
        # if comm.rank == 0:
        print(msg, flush=True)

        exit()    
