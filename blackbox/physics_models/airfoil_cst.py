# Imports
import os, sys, shutil, pickle, time
import numpy as np
from scipy.io import savemat
from pyDOE2 import lhs
from mpi4py import MPI
from baseclasses import AeroProblem
from pygeo import DVGeometryCST, DVConstraints
from prefoil.utils import readCoordFile

comm = MPI.COMM_WORLD

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
        self.refine = 0

class AirfoilCST():
    """
        This class provides methods for generating samples for a general airfoil
        using CST parameterization.
    """

    def __init__(self, options):

        # Partial checking of options argument provided by the user.
        if not isinstance(options, dict):
            self._error("The 'options' argument provided is not a dictionary.")
        elif options == {}:
            self._error("The 'options' argument provided is an empty dictionary.")

        # Creating an empty options dictionary
        self.options = {}

        # Setting up default options
        self._getDefaultOptions()

        # Changing to absolute path for airfoil file
        options["airfoilFile"] = os.path.abspath(options["airfoilFile"])

        # Setting up the required options list
        requiredOptions = ["solverOptions", "meshingOptions", "airfoilFile", "aeroProblem", "numCST"]

        # Validating user provided options
        self._checkOptions(options, requiredOptions)

        # Updating/Appending the default option list with user provided options
        self._setOptions(options)

        # Overiding/set some solver options
        self.options["solverOptions"]["printAllOptions"] = False
        self.options["solverOptions"]["printIntro"] = False
        self.options["solverOptions"]["outputDirectory"] = "."
        self.options["solverOptions"]["numberSolutions"] = False

        # Getting abs path for the storage directory
        self.options["directory"] = os.path.abspath(self.options["directory"])

        # Setting up the folder for saving the results
        directory = self.options["directory"]

        # Creating directory for storing the results
        if not os.path.isdir(directory):
            os.system("mkdir {}".format(directory))
        else:
            os.system("rm -r {}".format(directory))
            os.system("mkdir {}".format(directory))

        # Initializing the parametrization object
        self.DVGeo = DVGeometryCST(self.options["airfoilFile"], numCST=self.options["numCST"], comm=comm)

        # Adding pointset to the parametrization
        coords = readCoordFile(self.options["airfoilFile"])
        coords = np.hstack(( coords, np.zeros((coords.shape[0], 1)) ))
        self.DVGeo.addPointSet(coords, "airfoil")

        # Some initializations which will be used later
        self.DV = []
        self.genSamples = 0

    # ----------------------------------------------------------------------------
    #                       Design Variable related methods
    # ----------------------------------------------------------------------------

    def addDV(self, name: str, lowerBound: list, upperBound: list) -> None:
        """
            Method for adding a DV for CST parameterization.
        """

        # Checking
        self._checkDV(name, lowerBound, upperBound)

        if name == "upper" or name == "lower":
            coeff = self.DVGeo.defaultDV[name]
            lowerBound = coeff + np.sign(coeff)*coeff*lowerBound
            upperBound = coeff + np.sign(coeff)*coeff*upperBound
            locator = np.array(["{}".format(name)]*len(lowerBound))
        else:
            locator = np.array(["{}".format(name)])

        # Adding DV into DVGeo
        if name not in ["alpha", "mach", "altitude"]:
            self.DVGeo.addDV("{}".format(name), "{}".format(name))

        # Creating/Appending lower bound, upper bound, and locator
        if len(self.DV) == 0:
            self.upperBound = np.array([upperBound])
            self.lowerBound = np.array([lowerBound])
            self.locator = np.array([locator])
        else:
            self.upperBound = np.append(self.upperBound, upperBound)
            self.lowerBound = np.append(self.lowerBound, lowerBound)
            self.locator = np.append(self.locator, locator)

        # Adding the DV to the list
        self.DV.append(name)

    def removeDV(self, name: str) -> None:
        """
            Method to remove a DV. 
        """

        if name not in self.DV:
            self._error("{} doesn't exists as a DV.".format(name))

        # Finding the indices to be removed
        loc = self.locator == name
        loc = loc.reshape(-1,)

        # Removing the entry in the bounds and locator
        self.lowerBound = np.delete(self.lowerBound, loc)
        self.upperBound = np.delete(self.upperBound, loc)
        self.locator = np.delete(self.locator, loc)

        # Removing the entry from DV list
        self.DV.remove(name)

    # ----------------------------------------------------------------------------
    #                   Methods related to sample generation
    # ----------------------------------------------------------------------------

    def generateSamples(self, numSamples: int) -> None:
        """
            Method for generating samples.
        """

        # Performing checks
        if len(self.DV) == 0:
            self._error("Add design variables before running the analysis.")

        if not isinstance(numSamples, int):
            self._error("Number of samples argument is not an integer.")

        # Number of analysis passed/failed
        noFailed = 0

        # Generating LHS samples
        samples = self._lhs(numSamples)

        # Creating empty dictionary for storing the data
        data = {}

        description = open("{}/description.txt".format(self.options["directory"]), "a")

        description.write("--------------------------------------")
        description.write("\nDescription file for sample generation")
        description.write("\n--------------------------------------")
        description.write("\nDesign variables: {}".format(self.DV))
        description.write("\nLower bound for design variables: {}".format(self.lowerBound))
        description.write("\nUpper bound for design variables: {}".format(self.upperBound))
        description.write("\nTotal number of samples requested: {}".format(numSamples))
        description.write("--------------------------------------")
        description.write("\nAnalysis specific description")
        description.write("\n--------------------------------------")

        # Generate data
        for sampleNo in range(numSamples):

            description.write("\nAnalysis 1: ")

            # Current sample
            x = samples[sampleNo,:]

            # Starting time
            t1 = time.time()

            try:
                # Getting output for specific sample
                output = self.getObjectives(x)

            except:
                print("Error occured during the analysis. Check analysis.log in the respective folder for more details.")
                noFailed += 1

            else:
                if output["fail"] == True: # Check for analysis failure
                    noFailed += 1
                else:
                    # Creating a dictionary of data
                    if self.genSamples - noFailed == 1:
                        data["x"] = np.array(x)
                        for value in output.keys():
                            data[value] = np.array([output[value]])
                    else:
                        data["x"] = np.vstack((data["x"], x))
                        for value in output.keys():
                            data[value] = np.vstack(( data[value], np.array([output[value]]) ))

                    # Saving the results
                    savemat("{}/data.mat".format(self.options["directory"]), data)

            finally:
                # Ending time
                t2 = time.time()

                # Writing time taken to file
                description.write("\nTime taken for analysis: {} min.".format((t2-t1)/60))

        # Closing the description file
        description.close()

    def getObjectives(self, x: np.ndarray) -> dict:
        """
            Method for running a single analysis.
        """

        # Performing checks
        if len(self.DV) == 0:
            self._error("Add design variables before running the analysis.")

        if not isinstance(x, np.ndarray):
            self._error("Input sample is not a numpy array.")

        if x.ndim != 1:
            self._error("Input sample is a single dimensional array.")

        if len(x) != len(self.lowerBound):
            self._error("Input sample is not of correct size.")

        print("Running analysis {}".format(self.genSamples + 1))

        directory = self.options["directory"]

        # Create the folder for saving the results
        os.system("mkdir {}/{}".format(directory, self.genSamples+1))

        # Copying the runscript
        pkgdir = sys.modules["blackbox"].__path__[0]
        filepath = os.path.join(pkgdir, "runscripts/runscript_airfoil_cst.py")
        shutil.copy(filepath, "{}/{}/runscript.py".format(directory, self.genSamples+1))

        # Creating the new design variable dict
        # If there are no shape DV, then DVGeo
        # will not update the airfoil pointset.
        newDV = {}
        for dv in self.DV:
            loc = self.locator == dv
            loc = loc.reshape(-1,)
            newDV[dv] = x[loc]

        # Updating the airfoil pointset based on new DV
        self.DVGeo.setDesignVars(newDV)
        points = self.DVGeo.update("airfoil")[:,0:2]

        # Changing the first and last point for meshing
        # TO DO - This needs to change for more general scenario
        points[0,1] = 0.0
        points[-1,1] = 0.0

        # Changing the directory to analysis folder
        os.chdir("{}/{}".format(directory, self.genSamples+1))

        # Create input file
        self._creatInputFile(x)

        # Writing the surface mesh
        self.DVGeo.foil.writeCoords("surfMesh", points)

        try:
            # Run the runscript
            child_comm = MPI.COMM_WORLD.Spawn(sys.executable, args=["runscript.py"], maxprocs=self.options["noOfProcessors"])
            child_comm.Disconnect()
            time.sleep(0.25) # Very important - do not remove this

        except:
            child_comm.Disconnect()
            time.sleep(0.25) # Very important - do not remove this

            raise Exception

        else:
            if not os.path.exists("output.pickle"):
                raise Exception

            # Reading the output file containing results
            filehandler = open("output.pickle", 'rb')
            output = pickle.load(filehandler)
            filehandler.close()

            # Calculate the volume (here area)
            output = self._calcVol(output)

            return output

        finally:
            # Cleaning the directory
            files = ["surfMesh.xyz", "volMesh.cgns", "input.pickle", "runscript.py", "output.pickle"]
            for file in files:
                if os.path.exists(file):
                    os.system("rm {}".format(file))

            # Changing the directory back to root
            os.chdir("../..")

            # Increase the number of generated samples
            self.genSamples += 1

    # ----------------------------------------------------------------------------
    #                       Methods related to validation
    # ----------------------------------------------------------------------------

    def _checkOptions(self, options: dict, requiredOptions: list) -> None:
        """
            This method validates user provided options.
        """

        defaultOptions = list(self.options.keys())
        allowedUserOptions = defaultOptions
        allowedUserOptions.extend(requiredOptions)
        userProvidedOptions = list(options.keys())

        # Checking if user provided option contains only allowed attributes
        if not set(userProvidedOptions).issubset(allowedUserOptions):
            self._error("Option dictionary contains unrecognized attribute(s): {}"\
                        .format(set(userProvidedOptions) - set(allowedUserOptions)))

        # Checking if user has mentioned all the requried attributes
        if not set(requiredOptions).issubset(userProvidedOptions):
            self._error("Option dictionary doesn't contain following attribute(s): {}"\
                        .format(set(requiredOptions) - set(userProvidedOptions)))

        ############ Validating airfoilFile
        if not os.path.exists(options["airfoilFile"]):
            self._error("\"airfoilFile\" doesn't exists.")            

        ############ Validating numCST
        if not isinstance(options["numCST"], list):
            self._error("\"numCST\" is not a list")
        else:
            if len(options["numCST"]) != 2:
                self._error("\"numCST\" should have only two entries.")

            if not isinstance(options["numCST"][0], int):
                self._error("First entry in \"numCST\" is not an integer.")
            elif options["numCST"][0] <= 0:
                self._error("First entry in \"numCST\" is less than 1.")
                    
            if not isinstance(options["numCST"][1], int):
                self._error("Second entry in \"numCST\" is not an integer.")
            elif options["numCST"][1] <= 0:
                self._error("Second entry in \"numCST\" is less than 1.")

        ############ Validating aeroProblem
        if not isinstance(options["aeroProblem"], AeroProblem):
            self._error("\"aeroProblem\" attribute is not an aeroproblem.")

        ############ Validating solverOptions
        if not isinstance(options["solverOptions"], dict):
            self._error("\"solverOptions\" attribute is not a dictionary.")

        if "gridFile" in options["solverOptions"].keys():
            self._warning("\"gridFile\" attribute in solver options is not required.")

        ############ Validating meshingOptions
        if not isinstance(options["meshingOptions"], dict):
            self._error("\"meshingOptions\" attribute is not a dictionary.")

        if "inputFile" in options["meshingOptions"].keys():
            self._warning("\"inputFile\" attribute in meshing options is not required.")

        ############ Validating noOfProcessors
        if "noOfProcessors" in userProvidedOptions:
            if not isinstance(options["noOfProcessors"], int):
                self._error("\"noOfProcessors\" attribute is not an integer.")

        ############ Validating refine
        if "refine" in userProvidedOptions:
            if not isinstance(options["refine"], int):
                self._error("\"refine\" attribute is not an integer.")

        ############ Validating directory attribute
        if "directory" in userProvidedOptions:
            if not isinstance(options["directory"], str):
                self._error("\"directory\" attribute is not string.")

    def _checkObjectives(self, options: dict) -> None:
        """
            Validating the objectives provided by the user
        """

        allowedObjectives = ["cl", "cd", "cmz", "lift", "drag"]

        if type(options["objectives"]) == list:
            if not set(options["objectives"]).issubset(allowedObjectives):
                self._error("One or more objective(s) are not valid. Only these \
                    objectives are supported: {}".format(allowedObjectives))
        else:
            self._error("\"objectives\" option is not a list.")

    def _checkDV(self, name, lb, ub) -> None:
        """
            Method for validating DV.
        """

        # List of possible DVs
        possibleDVs = ["upper", "lower", "N1", "N2", "N1_upper", "N2_upper", 
                        "N1_lower", "N2_lower", "alpha", "mach", "altitude"]

        # Validating name of the DV
        if not isinstance(name, str):
            self._error("Name argument is not a string.")

        if name not in possibleDVs:
            self._error("{} argument is not a valid DV.".format(name))

        if name in self.DV:
            self._error("{} already exists.".format(name))

        # Validating lb and ub of the DV
        if not isinstance(lb, float):
            self._error("Lower Bound argument is not a float.")

        if not isinstance(ub, float):
            self._error("Upper Bound argument is not a float.")

        if lb >= ub:
            self._error("Lower bound is greater than or equal to upper bound.")

    # ----------------------------------------------------------------------------
    #                               Other methods
    # ----------------------------------------------------------------------------

    def _getDefaultOptions(self) -> None:
        """
            Setting up the initial values of options.
        """

        defaultOptions = DefaultOptions()

        for key in vars(defaultOptions):
            value = getattr(defaultOptions, key)
            self.options[key] = value

    def _setOptions(self, options: dict) -> None:
        """
            Method for assigning user provided options.
        """

        for key in options.keys():
            # update strategy is different for a dictionary
            # and other type of variables
            if isinstance(options[key], dict):
                # if the option is already present in default dictionary, 
                # then append the user provided key-value pairs to the 
                # default dictionary. For example: solverOptions.
                if key in self.options.keys():
                    self.options[key].update(options[key]) 
                else:
                    self.options[key] = options[key]
            else:
                self.options[key] = options[key]

    def _lhs(self, numSamples) -> np.ndarray:
        """
            Method to generate the lhs samples.
        """

        # Number of dimensions
        dim = len(self.lowerBound)

        # Generating normalized lhs samples
        samples = lhs(dim, samples=numSamples, criterion='cm', iterations=100)

        # Scaling the samples
        x = self.lowerBound + (self.upperBound - self.lowerBound) * samples

        return x

    def _creatInputFile(self, x:np.ndarray) -> None:
        """
            Method to create an input file for analysis.
        """

        # Creating input dict
        input = {
            "solverOptions": self.options["solverOptions"],
            "aeroProblem": self.options["aeroProblem"],
            "meshingOptions": self.options["meshingOptions"],
            "refine": self.options["refine"]
        }

        # Adding non-shape DV
        if "alpha" in self.DV:
            loc = self.locator == "alpha"
            loc = loc.reshape(-1,)
            input["alpha"] = x[loc]

        if "mach" in self.DV:
            loc = self.locator == "mach"
            loc = loc.reshape(-1,)
            input["mach"] = x[loc]

        if "altitude" in self.DV:
            loc = self.locator == "altitude"
            loc = loc.reshape(-1,)
            input["altitude"] = x[loc]

        # Saving the input file
        filehandler = open("input.pickle", "xb")
        pickle.dump(input, filehandler)
        filehandler.close()

    def _calcVol(self, output) -> dict:
        """
            Method for calculating the area of the airfoil.
        """

        #### To Do: Entire method can be refactored in a better manner.

        # Setting up constraints
        DVCon = DVConstraints()
        DVCon.setDVGeo(self.DVGeo)
        DVCon.setSurface(output.pop("pts"))

        # Defining the four corners of the 2D plane for volume constraint
        le = 1e-4
        leList = [[le, 0, le], [le, 0, 1.0 - le]]
        teList = [[1.0 - le, 0, le], [1.0 - le, 0, 1.0 - le]]
        
        # Add volume (area in the case of airfoil) constraints
        DVCon.addVolumeConstraint(leList, teList, 2, 100, lower=0.5, upper=3.0, scaled=False, name="area")

        # Calculate the volume
        DVCon.evalFunctions(output)

        return output

    def _warning(self, message: str) -> None:
        """
            Method for printing warnings in nice manner.
        """

        msg = "\n+" + "-" * 78 + "+" + "\n" + "| Blackbox Warning: "
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

    def _error(self, message: str) -> None:
        """
            Method for printing errors in nice manner.
        """

        msg = "\n+" + "-" * 78 + "+" + "\n" + "| Blackbox Error: "
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
