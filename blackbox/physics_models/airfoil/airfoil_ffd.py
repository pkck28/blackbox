# Imports
import os, sys, shutil, pickle, time, psutil
import numpy as np
from scipy.io import savemat
from scipy import integrate
from pyDOE2 import lhs
from mpi4py import MPI
from baseclasses import AeroProblem
from pygeo import DVGeometry
from prefoil import Airfoil
from prefoil.utils import readCoordFile

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
        self.refine = 0
        self.writeSliceFile = False
        self.writeAirfoilCoordinates = False
        self.plotAirfoil = False
        self.writeDeformedFFD = False

        # FFD related options
        self.fitted = False
        self.xmargin = 0.001
        self.ymarginu = 0.02
        self.ymarginl = 0.02

class AirfoilFFD():
    """
        This class provides methods for generating samples for a general airfoil
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
        requiredOptions = ["solverOptions", "meshingOptions", "airfoilFile", "aeroProblem", "nffd"]

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

        # Raise an error if matplotlib is not installed
        if self.options["plotAirfoil"]:
            if msg_matplotlib != None:
                self._error(msg_matplotlib)

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

        # Read the coordinate file
        self.coords = readCoordFile(self.options["airfoilFile"])

        # Some validation for coordinate file
        if self.coords[0,0] != self.coords[-1,0]:
            self._error("The X coordinate of airfoil doesn't start and end at same point.")
        elif self.coords[0,1] != self.coords[-1,1]:
            self._error("The Y coordinate of airfoil doesn't start and end at same point.")

        # Read the coordinate file
        self.coords = readCoordFile(self.options["airfoilFile"])
        airfoil = Airfoil(self.coords)

        # Creating FFD box
        airfoil.generateFFD(nffd=int(self.options["nffd"]/2), filename=directory + "/ffd", fitted=self.options["fitted"], 
                            xmargin=self.options["xmargin"], ymarginu=self.options["ymarginu"], 
                            ymarginl=self.options["ymarginl"], coords=self.coords)

        # Creating DVGeometry object
        self.DVGeo = DVGeometry(directory + "/ffd.xyz")

        # Adding pointset to the parametrization
        self.coords = np.hstack(( self.coords, np.zeros((self.coords.shape[0], 1)) ))
        self.DVGeo.addPointSet(self.coords, "airfoil")

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
            locator = np.array(["{}".format(name)]*len(lowerBound))

            if len(self.DV) == 0:
                self.upperBound = upperBound
                self.lowerBound = lowerBound
                self.locator = locator
            else:
                self.upperBound = np.append(self.upperBound, upperBound)
                self.lowerBound = np.append(self.lowerBound, lowerBound)
                self.locator = np.append(self.locator, locator)

            # Adding FFD points as a DV
            self.DVGeo.addSpanwiseLocalDV("shape", spanIndex="k", axis="y", lower=lowerBound, upper=upperBound)
            
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
        description.write("\nNumber of FFD points: {}".format(self.options["nffd"]))
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

            try:
                # Getting output for specific sample
                output = self.getObjectives(x)

            except Exception as e:
                print("Error occured during the analysis. Check analysis.log in the respective folder for more details.")
                failed.append(sampleNo + 1)
                description.write("\nAnalysis failed.")

            else:
                # Check for analysis failure
                if output["fail"] == True: # Check for analysis failure
                    failed.append(sampleNo + 1)
                    description.write("\nAnalysis failed.")

                # Creating a dictionary of data
                else:
                    if self.genSamples - len(failed) == 1:
                        data["x"] = np.array(x)
                        for value in output.keys():
                            data[value] = np.array([output[value]])

                    else:
                        # Appending data dictionary created earlier
                        data["x"] = np.vstack((data["x"], x))
                        for value in output.keys():
                            data[value] = np.vstack(( data[value], np.array([output[value]]) ))

                    # Saving the results
                    savemat("{}/data.mat".format(self.options["directory"]), data)

            finally:
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
        filepath = os.path.join(pkgdir, "runscripts/airfoil/runscript_airfoil.py")

        # Copy the runscript to analysis directory
        shutil.copy(filepath, "{}/{}/runscript.py".format(directory, self.genSamples+1))

        # Creating the new design variable dict
        # If there are no shape DV, then DVGeo
        # will not update the airfoil pointset.
        newDV = {}
        for dv in self.DV:
            loc = self.locator == dv
            loc = loc.reshape(-1,)
            newDV[dv] = x[loc]

        # Creating the new design variable dict
        self.DVGeo.setDesignVars(newDV)
        points = self.DVGeo.update("airfoil")[:,0:2]

        # Changing the directory to analysis folder
        os.chdir("{}/{}".format(directory, self.genSamples+1))

        if self.options["writeAirfoilCoordinates"]:
            self._writeCoords(coords=points, filename="deformedAirfoil.dat")

        if self.options["plotAirfoil"]:
            self._plotAirfoil(points)

        if self.options["writeDeformedFFD"]:
            self.DVGeo.writePlot3d("deformedFFD.xyz")

        # Create input file
        self._creatInputFile(x)

        # Writing the surface mesh
        self._writeSurfMesh(coords=points, filename="surfMesh.xyz")

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

            # Calculate the area
            output["area"] = integrate.simpson(points[:,0], points[:,1], even="avg")

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

    def calculateArea(self, x: np.ndarray) -> float:
        """
            Note: This function should not be called in the middle of analysis
            It should ONLY be used from outside. Do not use this method within 
            getObjectives. That method has its own implementation
            of area calculation.

            Function to calculate the area of the airfoil
            based on the value of design variable.

            Input:
            x - 1D numpy array (value of dv).

            Ouput:
            area: area of the airfoil.

            Note: To use this method, shape should be added as DV
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

        if "upper" not in self.DV and "lower" not in self.DV:
            self._error("\"upper\" or \"lower\" surface is not added as design variable.")

        # Creating dictionary from x
        newDV = {}
        for dv in self.DV:
            loc = self.locator == dv
            loc = loc.reshape(-1,)
            newDV[dv] = x[loc]

        # Updating the airfoil pointset based on new DV
        self.DVGeo.setDesignVars(newDV)

        # Getting the updated airfoil points
        points = self.DVGeo.update("airfoil")[:,0:2]
        x = points[:,0]
        y = points[:,1]

        # Calculate the area using simpson's rule
        # Note: x and y are both flipped here
        area = integrate.simpson(x, y, even='avg')

        return area

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
        if not os.path.exists(os.path.abspath(options["airfoilFile"])):
            self._error("\"airfoilFile\" doesn't exists.")
        else:
            options["airfoilFile"] = os.path.abspath(options["airfoilFile"])

        ############ Validating nffd
        if not isinstance(options["nffd"], int):
            self._error("\"nffd\" attribute is not an integer.")

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

            if psutil.cpu_count(False) < options["noOfProcessors"] + 1:
                self._error("\"noOfProcessors\" requested is more than available processors.")

        ############ Validating refine
        if "refine" in userProvidedOptions:
            if not isinstance(options["refine"], int):
                self._error("\"refine\" attribute is not an integer.")

        ############ Validating writeSliceFile
        if "writeSliceFile" in userProvidedOptions:
            if not isinstance(options["writeSliceFile"], bool):
                self._error("\"writeSliceFile\" attribute is not a boolean value.")

        ############ Validating directory attribute
        if "directory" in userProvidedOptions:
            if not isinstance(options["directory"], str):
                self._error("\"directory\" attribute is not string.")

        ############ Validating FFD options
        if "fitted" in userProvidedOptions:
            if not isinstance(options["fitted"], bool):
                self._error("\"fitted\" attribute is not a boolean value.")

        if "xmargin" in userProvidedOptions:
            if not isinstance(options["xmargin"], float):
                self._error("\"xmargin\" attribute is not a float.")

        if "ymarginu" in userProvidedOptions:
            if not isinstance(options["ymarginu"], float):
                self._error("\"ymarginu\" attribute is not a float.")

        if "ymarginl" in userProvidedOptions:
            if not isinstance(options["ymarginl"], float):
                self._error("\"ymarginl\" attribute is not a float.")

        if "writeDeformedFFD" in userProvidedOptions:
            if not isinstance(options["writeDeformedFFD"], bool):
                self._error("\"writeDeformedFFD\" attribute is not a boolean.")

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

        # Validating the bounds for "shape" variable
        if name == "shape":
            if not isinstance(lb, np.ndarray) or lb.ndim != 1:
                self._error("Lower bound for \"shape\" variable should be a 1D numpy array.")

            if not isinstance(ub, np.ndarray) or ub.ndim != 1:
                self._error("Upper bound for \"shape\" variable should be a 1D numpy array.")

            if len(lb) != self.options["nffd"]:
                self._error("Length of lower bound array is not equal to number of FFD points.")

            if len(ub) != self.options["nffd"]:
                self._error("Length of upper bound array is not equal to number of FFD points.")

            if np.any(lb >= ub):
                self._error("Lower bound is greater than or equal to upper bound for atleast one DV.")

            # Checking if the bounds are within the limits
            coeff = self.DVGeo.origFFDCoef
            index = self.DVGeo.getLocalIndex(0)
            dist = coeff[index[:,1,0], 1] - coeff[index[:,0,0], 1]
            allowableLowerBound = np.zeros(self.options["nffd"])
            allowableUpperBound = np.zeros(self.options["nffd"])

            for i in range(dist.shape[0]):
                allowableLowerBound[2*i] = -0.45 * dist[i]
                allowableLowerBound[2*i+1] = -0.45 * dist[i]
                allowableUpperBound[2*i] = 0.45 * dist[i]
                allowableUpperBound[2*i+1] = 0.45 * dist[i]

            if np.any(lb <= allowableLowerBound):
                self._error("Lower bound for some FFD points is greater than or equal to 45% of the \
                            local FFD thickness. Reduce the bound and try again.")
                
            if np.any(ub >= allowableUpperBound):
                self._error("Upper bound for some FFD points is greater than or equal to 45% of the \
                            local FFD thickness. Reduce the bound and try again.")

        else:
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
            "meshingOptions": self.options["meshingOptions"],
            "refine": self.options["refine"],
            "writeSliceFile": self.options["writeSliceFile"]
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

    def _writeCoords(self, coords, filename) -> None:
        """
            Writes out a set of airfoil coordinates in dat format.
        """

        # X and Y ccordinates of the airfoil
        x = coords[:, 0]
        y = coords[:, 1]

        with open(filename, "w") as f:
            for i in range(len(x)):
                f.write(str(round(x[i], 12)) + "\t\t" + str(round(y[i], 12)) + "\n")

        f.close()

    def _writeSurfMesh(self, coords, filename):
        """
            Writes out surface mesh in Plot 3D format (one element in z direction)
        """

        # X and Y ccordinates of the airfoil
        x = coords[:, 0]
        y = coords[:, 1]

        # Writing the file
        with open(filename, "w") as f:
            f.write("1\n")
            f.write("%d %d %d\n" % (len(x), 2, 1))
            for iDim in range(3):
                for j in range(2):
                    for i in range(len(x)):
                        if iDim == 0:
                            f.write("%g\n" % x[i])
                        elif iDim == 1:
                            f.write("%g\n" % y[i])
                        else:
                            f.write("%g\n" % (float(j)))

        f.close()

    def _plotAirfoil(self, points) -> None:
        """
            Method for plotting the base airfoil
            and the deformed airfoil.
        """

        fig, ax = plt.subplots()

        ax.plot(self.coords[:,0], self.coords[:,1], label="Original airfoil")
        ax.plot(points[:,0], points[:,1], label="Deformed airfoil")
        ax.set_xlabel("x/c", fontsize=14)
        ax.set_ylabel("y/c", fontsize=14)
        ax.legend(fontsize=12)

        plt.savefig("airfoil.png", dpi=400)

        plt.close()

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
