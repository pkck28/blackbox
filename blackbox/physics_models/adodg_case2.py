import os, sys, shutil, pickle, time
import numpy as np
from pyDOE2 import lhs
from scipy.io import savemat

class DefaultOptions():
    """
        Class creates a default option for physics problems which are later 
        edited/appended with user provided options.
    """

    def __init__(self):
        
        # Aero solver Options
        self.aeroSolver = "adflow"

        # Other options
        self.directory = "output"
        self.noOfProcessors = 4

class ADODGCase2():

    def __init__(self, type="multi", options=None):
        """
            Initialization.
        """

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

        # Setting up default options
        self._getDefaultOptions()

        requiredOptions = ["objectives", "lb", "ub", "numberOfSamples", "samplingMethod", "aeroSolverOptions", "ffd"]

        # Validating user provided options
        self._checkOptions(options, requiredOptions)

        # Setting some default options
        self.options["type"] = "multi"
        self.options["aeroSolverOptions"] = {
                "printAllOptions": False,
                "printIntro": False,
                "outputDirectory": "."
        }

        # Updating/Appending the default option list with user provided options
        self._setOptions(options)

        # Setting up the folders for saving the results
        self._setDirectory()

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

        data = {}
        data['x'] = self.x

        totalTime = 0

        # Running the analysis
        for sampleNo in range(self.options["numberOfSamples"]):

            print("Running analysis {} of {}".format(sampleNo + 1, self.options["numberOfSamples"]))
            os.chdir("{}/{}".format(self.options["directory"],sampleNo))
            
            # creating input file for analysis
            self._creatInputFile(sampleNo)
            
            # Starting time
            t1 = time.time()

            # run the analysis
            os.system("mpirun -n {} --use-hwthread-cpus python runscript_adodg_case2.py >> analysis_log.txt".format(self.options["noOfProcessors"]))

            # Ending time
            t2 = time.time()

            totalTime += (t2-t1)/60

            print("Time taken for analysis: {} min.".format((t2-t1)/60))

            # reading output from the analysis
            filehandler = open("output.pickle", 'rb')
            output = pickle.load(filehandler)
            filehandler.close()

            # cleaning the directory
            os.system("rm -r output.pickle runscript_adodg_case2.py")

            # changing the directory
            os.chdir("../..")

            # storing the results
            for value in output.keys():
                if sampleNo == 0:
                    data[value] = np.array([])
                data[value] = np.append(data[value], output[value])

        print("Total time taken: {} minutes".format(totalTime))

        # saving the results
        os.chdir("{}".format(self.options["directory"]))
        savemat("data.mat", data)
        os.chdir("../")

    # ----------------------------------------------------------------------------
    #               All the methods for single analysis
    # ----------------------------------------------------------------------------

    def _setupSingleAnalysis(self, options):

        # Creating an empty options dictionary
        self.options = {}
        
        # Setting up default options
        self._getDefaultOptions()

        # Changing default directory value for single analysis
        self.options["directory"] = "additional_samples"

        requiredOptions = ["objectives", "aeroSolverOptions", "ffd"]

        # Validating user provided options
        self._checkOptions(options, requiredOptions)

        # Setting some default solver options
        self.options["aeroSolverOptions"] = {
                "printAllOptions": False,
                "printIntro": False,
                "outputDirectory": "."
            }
        self.options["type"] = "single"
        self.sampleNo = 0

        # Updating/Appending the default option list with user provided options
        self._setOptions(options)

        # Setting up the folder for saving the result
        self._setDirectory()

    def getObjectives(self, user_sample):
        """
            Method for running a single sample analysis.
            Input should be a 1D numpy array of correct size.
        """
        
        directory = self.options["directory"]

        # Checking whether the type is consistent
        if self.options["type"] != "single":
            self._error("You can run getObjectives() method only when type is single.")

        if type(user_sample) != np.ndarray:
            self._error("Sample provided by user is not a numpy array")

        # Pasting analysis file
        self._copyAnalysisFile(self.sampleNo)

        # Getting variable value
        alpha = user_sample[0]
        shape = user_sample[1:]

        # Changing directory and running the analysis
        os.chdir("{}/{}".format(directory, self.sampleNo))
        self._creatInputFile(shape=shape, alpha=alpha)
        print("Running analysis {}".format(self.sampleNo+1))
        
        # Starting time
        t1 = time.time()
            
        os.system("mpirun -n {} --use-hwthread-cpus python runscript_adodg_case2.py >> analysis_log.txt".format(self.options["noOfProcessors"]))
        
        # Ending time
        t2 = time.time()

        print("Time taken for analysis: {} min.".format((t2-t1)/60))

        # Reading the output file containing results
        filehandler = open("output.pickle", 'rb')
        result = pickle.load(filehandler)
        filehandler.close()

        os.system("rm -r output.pickle runscript_adodg_case2.py")
        os.chdir("../..")

        self.sampleNo += 1

        return result

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

    def _checkOptions(self, options, requiredOptions):
        """
            This method validates user provided options.
        """

        defaultOptions = list(self.options.keys())
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

        ############ Checking objectives
        self._checkObjectives(options)

        ############ Validating number of samples attribute
        if "numberOfSamples" in userProvidedOptions:
            if type(options["numberOfSamples"]) is not int:
                self._error("\"numberOfSamples\" attribute is not an integer.")
            
            # Setting minimum limit on number of samples
            if options["numberOfSamples"] < 2:
                self._error("Number of samples need to be at-least 2.")

        ############ Validating sampling method
        if "samplingMethod" in userProvidedOptions:
            if options["samplingMethod"] != "lhs":
                self._error("\"samplingMethod\" attribute is not correct. Only \"lhs\" is only allowed.")

        ############ Checking aeroSolverOptions
        if type(options["aeroSolverOptions"]) != dict:
            self._error("\"aeroSolverOptions\" attribute is not a dictionary.")

        # Checking for not allowed aero solver options
        commonElements = set(options["aeroSolverOptions"].keys()).intersection(set(["printAllOptions", "printIntro", "outputDirectory"]))
        if len(commonElements) != 0:
            self._error("Please remove {} attribute from the \"aeroSolverOptions\"".format(commonElements))

        ############ Checking ffdfile
        # To do

        ############ Checking noOfProcessors
        if "noOfProcessors" in userProvidedOptions:
            if type(options["noOfProcessors"]) is not int:
                self._error("\"noOfProcessors\" attribute is not an integer.")

        ############ Validating directory attribute
        if "directory" in userProvidedOptions:
            if type(options["directory"]) is not str:
                self._error("\"directory\" attribute is not string.")

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

    def _checkObjectives(self, options):
        """
            Checking the objectives provided by the user
        """

        allowedObjectives = ["cl", "cd", "cmz", "lift", "drag"]

        if type(options["objectives"]) == list:
            if not set(options["objectives"]).issubset(allowedObjectives):
                self._error("One or more objective(s) are not valid. Only these \
                    objectives are supported: {}".format(allowedObjectives))
        else:
            self._error("\"objectives\" option is not a list.")

    def _setDirectory(self):
        """
            Method for setting up directory.
        """

        directory = self.options["directory"]

        # creating directory for storing the results
        if not os.path.isdir(directory):
            os.system("mkdir {}".format(directory))
        else:
            os.system("rm -r {}".format(directory))
            os.system("mkdir {}".format(directory))

        # creating individual sample directory
        if self.options["type"] == "multi":
            for sampleNo in range(self.options["numberOfSamples"]):
                self._copyAnalysisFile(sampleNo)

    def _lhs(self):
        """
            Method for creating a lhs sample.
            Stores a 2D numpy array of size (samples vs  dimensions).
            Each row represents a new sample and each column corresponds to
            a particular design variable.
        """

        dim = self._numberOfDV()

        lowerBound = np.zeros(dim)
        upperBound = np.zeros(dim)

        # Bounds
        lowerBound[0] = self.options["lb"][0]
        lowerBound[1:] = self.options["lb"][1]
        upperBound[0] = self.options["ub"][0]
        upperBound[1:] = self.options["ub"][1]

        samples = lhs(dim, samples=self.options["numberOfSamples"], criterion='cm', iterations=100)

        self.x = lowerBound + (upperBound - lowerBound) * samples

    def _numberOfDV(self):
        """
            Method to calculate the number of variables (shape + alpha).
        """

        # reading the ffd file
        file = open(self.options["ffd"], "r")
        lines = file.readlines()
        file.close()

        # returning the total DV = ffd points at one section + 1 (alpha)
        return int(lines[1].split()[0]) * int(lines[1].split()[1]) + 1

    def _copyAnalysisFile(self, sampleNo):
        """
            Method to copy analysis files.
        """

        directory = self.options["directory"]
        grid = self.options["aeroSolverOptions"]["gridFile"]
        ffd = self.options["ffd"]

        # creating the sample directory
        os.system("mkdir {}/{}".format(directory,sampleNo))
    
        # copying the runscript
        pkgdir = sys.modules["blackbox"].__path__[0]
        filepath = os.path.join(pkgdir, "runscripts/runscript_adodg_case2.py")
        shutil.copy(filepath, "{}/{}".format(directory,sampleNo))

        # copying the grid file
        os.system("cp -r {} {}/{}/grid.cgns".format(grid,directory,sampleNo))

        # copying the ffd file
        os.system("cp -r {} {}/{}/ffd.xyz".format(ffd,directory,sampleNo))

    def _creatInputFile(self, sampleNo=None, shape=None, alpha=None):
        """
            Method to create an input file for analysis.
        """

        input = {
            "aeroSolverOptions": self.options["aeroSolverOptions"],
            "objectives": self.options["objectives"],
            "sampleNo" : sampleNo
        }

        if self.options["type"] == "multi":
            input["alpha"] = self.x[sampleNo,0]
            input["shape"] = self.x[sampleNo,1:]
        elif self.options["type"] == "single":
            input["alpha"] = alpha
            input["shape"] = shape

        filehandler = open("input.pickle", "xb")
        pickle.dump(input, filehandler)
        filehandler.close()

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
