# Imports
import os, sys, shutil, pickle, time, psutil
import numpy as np
from scipy.io import savemat
from scipy import integrate
from pyDOE2 import lhs
from mpi4py import MPI
from baseclasses import AeroProblem
from pygeo import DVGeometry, DVConstraints
from pygeo.geo_utils.polygon import *
from cgnsutilities.cgnsutilities import readGrid
from idwarp import USMesh

# Trying to import matplotlib
try:
    import matplotlib.pyplot as plt
except ImportError:
    msg_matplotlib = "Matplotlib is not installed"
else:
    msg_matplotlib = None

comm = MPI.COMM_WORLD

class DefaultOptions():
    """
        Class creates a default option which are later 
        edited/appended with user provided options.
    """

    def __init__(self):
        
        # Aero solver Options
        self.aeroSolver = "adflow"

        # Other options
        self.directory = "output"
        self.noOfProcessors = 4
        self.sliceLocation = [0.05, 0.25, 0.5, 0.75, 0.95] # defines slice location

class WingFFD():
    """
        This class provides methods for generating samples for a simple wing
        using FFD parameterization.
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

        # Setting up the required options list
        requiredOptions = ["solverOptions", "volumeMesh", "ffdFile", "aeroProblem"]

        # Validating user provided options
        self._checkOptions(options, requiredOptions)

        # Updating/Appending the default option list with user provided options
        self._setOptions(options)

        # Overiding/set some solver options
        self.options["solverOptions"]["printAllOptions"] = False
        self.options["solverOptions"]["printIntro"] = False
        self.options["solverOptions"]["outputDirectory"] = "."
        self.options["solverOptions"]["numberSolutions"] = False
        self.options["solverOptions"]["printTiming"] = False

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

        # Create mesh deformation object
        self.mesh = USMesh(options={"gridFile": self.options["volumeMesh"]})

        # Get the surface mesh coordinates
        surfMesh = self.mesh.getSurfaceCoordinates()

        # Read the volume grid
        grid = readGrid(self.options["volumeMesh"])

        # Extract the surface and write it out for DVCon
        grid.extractSurface("wing_surf.xyz")

        # Creating DVGeometry object
        self.DVGeo = DVGeometry(self.options["ffdFile"])

        # Adding surface mesh co-ordinates as a pointset
        self.DVGeo.addPointSet(surfMesh, "wing_surface_mesh")

        # Creating DVConstraints object
        self.DVCon = DVConstraints()

        # Connecting DVGeo and DVCon
        self.DVCon.setDVGeo(self.DVGeo)

        # Set the surface of DVConstraint
        self.DVCon.setSurface("wing_surf.xyz", name="wing_surface")

        p1, p2, p3 = self.DVCon._getSurfaceVertices("wing_surface")

        vol = volumeTriangulatedMesh(p1, p2, p3)

        # Removing the surface mesh file
        os.system("rm wing_surf.xyz")

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

        if name == "shape":

            # Getting the number of ffd points
            nffd = self.DVGeo.getLocalIndex(0).flatten().shape[0]

            if len(self.DV) == 0:
                self.upperBound = np.array([upperBound]*nffd)
                self.lowerBound = np.array([lowerBound]*nffd)
                self.locator = np.array(["{}".format(name)]*nffd)
            else:
                self.upperBound = np.append(self.upperBound, np.array([upperBound]*nffd))
                self.lowerBound = np.append(self.lowerBound, np.array([lowerBound]*nffd))
                self.locator = np.append(self.locator, np.array(["{}".format(name)]*nffd))

            # Adding ffd as dv
            self.DVGeo.addLocalDV("shape", lower=lowerBound, upper=upperBound, axis="y", scale=1.0)

        else:
            locator = np.array(["{}".format(name)])

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
        failed =[]
        totalTime = 0

        # Generating LHS samples
        samples = self._lhs(numSamples)

        # Creating empty dictionary for storing the data
        data = {}

        # Creating and writing a description file
        description = open("{}/description.txt".format(self.options["directory"]), "a", buffering=1)
        description.write("---------------------------------------------------")
        description.write("\nAirfoil sample generation with FFD parametrization")
        description.write("\n--------------------------------------------------")
        description.write("\nDesign variables: {}".format(self.DV))
        description.write("\nNumber of FFD points: {}".format(samples.shape[1]))
        description.write("\nLower bound for design variables:\n{}".format(self.lowerBound))
        description.write("\nUpper bound for design variables:\n{}".format(self.upperBound))
        description.write("\nTotal number of samples requested: {}".format(numSamples))
        description.write("\n-----------------------------")
        description.write("\nAnalysis specific description")
        description.write("\n-----------------------------")

        # Generate data
        for sampleNo in range(numSamples):

            description.write("\nAnalysis {}: ".format(sampleNo+1))

            # Current sample
            x = samples[sampleNo,:]

            description.write("\nDesign Variable: {}".format(x))

            # Starting time
            t1 = time.time()

            # try:
            # Getting output for specific sample
            output = self.getObjectives(x)

            # finally:
            # Ending time
            t2 = time.time()

            totalTime += (t2-t1)/60

            # Writing time taken to file
            description.write("\nTime taken for analysis: {} min.".format((t2-t1)/60))

        # Making generated samples 0
        self.genSamples = 0

        # Writing final results in the description file
        description.write("\n--------------------------------------")
        description.write("\nTotal time taken for analysis: {} min.".format(totalTime))
        description.write("\nNumber of successful analysis: {}".format(numSamples - len(failed)))
        description.write("\nNumber of failed analysis: {}".format(len(failed)))
        if len(failed) != 0:
            description.write("\nFailed analysis: {}".format(failed))

        # Closing the description file
        description.close()

    def getObjectives(self, x: np.ndarray) -> tuple:
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

        # Getting the directory where package is saved
        pkgdir = sys.modules["blackbox"].__path__[0]

        # Setting filepath based on the how alpha is treated alpha
        filepath = os.path.join(pkgdir, "runscripts/runscript_wing.py")

        # Copy the runscript to analysis directory
        shutil.copy(filepath, "{}/{}/runscript.py".format(directory, self.genSamples+1))

        # Creating the new design variable dict
        # If there are no shape DV, then DVGeo
        # will not update the wing surface.
        newDV = {}
        for dv in self.DV:
            loc = self.locator == dv
            loc = loc.reshape(-1,)
            newDV[dv] = x[loc]

        # Creating the new design variable dict
        self.DVGeo.setDesignVars(newDV)
        newSurfMesh = self.DVGeo.update("wing_surface_mesh")

        # Update the surface mesh in IdWarp
        self.mesh.setSurfaceCoordinates(newSurfMesh)

        # Deform the volume mesh
        self.mesh.warpMesh()

        # Changing the directory to analysis folder
        os.chdir("{}/{}".format(directory, self.genSamples+1))

        # Write the new grid file.
        self.mesh.writeGrid('volMesh.cgns')

        # Create input file
        self._creatInputFile(x)

        # Spawning the runscript on desired number of processors
        child_comm = MPI.COMM_SELF.Spawn(sys.executable, args=["runscript.py"], maxprocs=self.options["noOfProcessors"])

        # Creating empty process id list
        pid_list = []

        # Getting each spawned process
        for processor in range(self.options["noOfProcessors"]):
            pid = child_comm.recv(source=MPI.ANY_SOURCE, tag=processor)
            pid_list.append(psutil.Process(pid))

        # Disconnecting from intercommunicator
        child_comm.Disconnect()

        # Waiting till all the child processors are finished
        while len(pid_list) != 0:
            for pid in pid_list:
                if not pid.is_running():
                    pid_list.remove(pid)

        try:
            # Reading the output file containing results
            filehandler = open("output.pickle", 'rb')

        except:
            raise Exception

        else:
            # Read the output
            output = pickle.load(filehandler)
            filehandler.close()

            # output["volume"] = volume

            return output

        finally:
            # Cleaning the directory
            files = ["surfMesh.xyz", "volMesh.cgns", "input.pickle", "runscript.py",
                    "output.pickle", "fort.6", "opt.hst"]
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

        ############ Validating aeroProblem
        if not isinstance(options["aeroProblem"], AeroProblem):
            self._error("\"aeroProblem\" attribute is not an aeroproblem.")

        ############ Validating solverOptions
        if not isinstance(options["solverOptions"], dict):
            self._error("\"solverOptions\" attribute is not a dictionary.")

        if "gridFile" in options["solverOptions"].keys():
            self._warning("\"gridFile\" attribute in solver options is not required.")

        ############ Validating ffdFile
        if not os.path.exists(os.path.abspath(options["ffdFile"])):
            self._error("Provided FFD file doesn't exists.")
        else:
            options["ffdFile"] = os.path.abspath(options["ffdFile"])

        ############ Validating volumeMesh
        if not os.path.exists(os.path.abspath(options["volumeMesh"])):
            self._error("Provided Volume mesh file doesn't exists.")
        else:
            options["volumeMesh"] = os.path.abspath(options["volumeMesh"])

        ############ Validating noOfProcessors
        if "noOfProcessors" in userProvidedOptions:
            if not isinstance(options["noOfProcessors"], int):
                self._error("\"noOfProcessors\" attribute is not an integer.")

            if psutil.cpu_count(False) < options["noOfProcessors"] + 1:
                self._error("\"noOfProcessors\" requested is more than available processors.")

        ############ Validating sliceLocation
        if "sliceLocation" in userProvidedOptions:
            if not isinstance(options["sliceLocation"], list):
                self._error("\"sliceLocation\" attribute is not a list of relative slice locations on wing.")

        ############ Validating directory attribute
        if "directory" in userProvidedOptions:
            if not isinstance(options["directory"], str):
                self._error("\"directory\" attribute is not string.")

    def _checkDV(self, name: str, lb: float or np.ndarray, ub: float or np.ndarray) -> None:
        """
            Method for validating DV user wants to add.
        """

        # List of possible DVs
        possibleDVs = ["shape", "alpha", "mach", "altitude"]

        # Validating name of the DV
        if not isinstance(name, str):
            self._error("Name argument is not a string.")

        # Checking if the DV is allowed
        if name not in possibleDVs:
            self._error("{} argument is not a valid DV.".format(name))

        # Checking if the DV is already added
        if name in self.DV:
            self._error("{} already exists.".format(name))

        # Checking if these variables are initialized through aero problem or not
        if name in ["mach", "altitude"]:
            if name not in self.options["aeroProblem"].inputs.keys():
                self._error("You need to initialize \"{}\" in the aero problem to set it as design variable.".format(name))

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
        samples = lhs(dim, samples=numSamples, criterion='cm', iterations=100*dim)

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
            "sliceLocation": self.options["sliceLocation"]
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

    def _warning(self, message: str) -> None:
        """
            Method for printing warnings in nice manner.
        """

        ############ To Do: Redundant - error and warning method can be combined.

        # Initial message - total len is 80 characters
        msg = "\n+" + "-" * 78 + "+" + "\n" + "| Blackbox Warning: "

        # Initial number of characters
        i = 16

        for word in message.split():
            if len(word) + i + 1 > 76:  # Finish line and start new one
                msg += " " * (76 - i) + " |\n| " + word + " " # Adding space and word in new line
                i = len(word) + 1 # Setting i value for new line
            else:
                msg += word + " " # Adding the word with a space
                i += len(word) + 1 # Increase the number of characters
        msg += " " * (76 - i) + " |\n" + "+" + "-" * 78 + "+" + "\n" # Adding last line
 
        print(msg, flush=True)

        exit()

    def _error(self, message: str) -> None:
        """
            Method for printing errors in nice manner.
        """

        # Initial message - total len is 80 characters
        msg = "\n+" + "-" * 78 + "+" + "\n" + "| Blackbox Error: "

        # Initial number of characters
        i = 16

        for word in message.split():
            if len(word) + i + 1 > 76:  # Finish line and start new one
                msg += " " * (76 - i) + " |\n| " + word + " " # Adding space and word in new line
                i = len(word) + 1 # Setting i value for new line
            else:
                msg += word + " " # Adding the word with a space
                i += len(word) + 1 # Increase the number of characters
        msg += " " * (76 - i) + " |\n" + "+" + "-" * 78 + "+" + "\n" # Adding last line
 
        print(msg, flush=True)

        exit()
