# Imports
import os, sys, shutil, pickle, time, psutil
import numpy as np
from scipy.io import savemat
from scipy import integrate
from pyDOE2 import lhs
from mpi4py import MPI
from baseclasses import AeroProblem
from pygeo import DVGeometryCST
from prefoil.utils import readCoordFile

# Trying to import pyvista
try:
    import pyvista
except ImportError:
    msg_pyvista = "pyVista is not installed"
else:
    msg_pyvista = None

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

        # Flow-field related options
        # if getFlowFieldData is false, then all other options are useless
        self.getFlowFieldData = False
        self.region = "surface"

        # Alpha implicit related options
        self.alpha = "explicit"
        self.targetCL = 0.824
        self.targetCLTol = 1e-4

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
        self.options["solverOptions"]["printTiming"] = False

        # Raise an error if pyvista is not installed
        if self.options["getFlowFieldData"]:
            if msg_pyvista != None:
                self._error(msg_pyvista)
            self.options["solverOptions"]["writeSurfaceSolution"] = True

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

        # Initializing the parametrization object
        self.DVGeo = DVGeometryCST(self.options["airfoilFile"], numCST=self.options["numCST"], comm=comm)

        # Adding pointset to the parametrization
        self.coords = np.hstack(( self.coords, np.zeros((self.coords.shape[0], 1)) ))
        self.DVGeo.addPointSet(self.coords, "airfoil")

        # Checking the number of points at trailing edge for blunt TE
        # Only two are allowed for CST. Otherwise, meshing will have problem.
        if not self.DVGeo.sharp:
            if len(np.where(self.coords[1:-1,0] == self.coords[0,0])[0]) > 1:
                self._error("There are more than two points in the trailing edge.")

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
            locator = np.array(["{}".format(name)]*len(lowerBound))

            if len(self.DV) == 0:
                self.upperBound = upperBound
                self.lowerBound = lowerBound
                self.locator = locator
            else:
                self.upperBound = np.append(self.upperBound, upperBound)
                self.lowerBound = np.append(self.lowerBound, lowerBound)
                self.locator = np.append(self.locator, locator)
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

        # Adding DV into DVGeo
        if name not in ["alpha", "mach", "altitude"]:
            self.DVGeo.addDV("{}".format(name), "{}".format(name))

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
        failed =[]
        totalTime = 0

        # Generating LHS samples
        samples = self._lhs(numSamples)

        # Creating empty dictionary for storing the data
        data = {}
        fieldData = {}

        # Creating and writing a description file
        description = open("{}/description.txt".format(self.options["directory"]), "a", buffering=1)
        description.write("---------------------------------------------------")
        description.write("\nAirfoil sample generation with CST parametrization")
        description.write("\n--------------------------------------------------")
        description.write("\nUpper surface CST coefficient: {}".format(self.DVGeo.defaultDV["upper"]))
        description.write("\nLower surface CST coefficient: {}".format(self.DVGeo.defaultDV["lower"]))
        description.write("\nDesign variables: {}".format(self.DV))
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
                output, field = self.getObjectives(x)

            except:
                print("Error occured during the analysis. Check analysis.log in the respective folder for more details.")
                failed.append(sampleNo + 1)
                description.write("\nAnalysis failed.")

            else:
                # Check for analysis failure
                if output["fail"] == True: # Check for analysis failure
                    failed.append(sampleNo + 1)
                    description.write("\nAnalysis failed.")
                    
                # Check for implicit alpha
                elif self.options["alpha"] == "implicit" and abs(output["cl"] - self.options["targetCL"]) > self.options["targetCLTol"]: 
                    failed.append(sampleNo + 1)
                    description.write("\nAnalysis failed.")

                # Creating a dictionary of data
                else:
                    if self.genSamples - len(failed) == 1:
                        data["x"] = np.array(x)
                        for value in output.keys():
                            data[value] = np.array([output[value]])

                        # Creating a dictionary of field data
                        if self.options["getFlowFieldData"]:
                            fieldData["x"] = np.array(x)
                            for value in field.keys():
                                if field[value].ndim == 2:
                                    fieldData[value] = field[value].reshape(1,-1,3)
                                else:
                                    fieldData[value] = field[value].reshape(1,-1)

                    else:
                        # Appending data dictionary created earlier
                        data["x"] = np.vstack((data["x"], x))
                        for value in output.keys():
                            data[value] = np.vstack(( data[value], np.array([output[value]]) ))

                        # Appending field data dictionary created earlier
                        if self.options["getFlowFieldData"]:
                            fieldData["x"] = np.vstack((fieldData["x"], x))
                            for value in field.keys():
                                if field[value].ndim == 2:
                                    fieldData[value] = np.vstack(( fieldData[value], field[value].reshape(1,-1,3) ))
                                else:
                                    fieldData[value] = np.vstack(( fieldData[value], field[value].reshape(1,-1) ))

                    # Saving the results
                    savemat("{}/data.mat".format(self.options["directory"]), data)

                    # Saving the field data
                    if self.options["getFlowFieldData"]:
                        savemat("{}/fieldData.mat".format(self.options["directory"]), fieldData)

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
        if self.options["alpha"] == "explicit":
            filepath = os.path.join(pkgdir, "runscripts/runscript_airfoil.py")
        else:
            # filepath = os.path.join(pkgdir, "runscripts/runscript_airfoil_cst_opt.py")
            filepath = os.path.join(pkgdir, "runscripts/runscript_airfoil_cst_rf.py")

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

        # Updating the airfoil pointset based on new DV
        self.DVGeo.setDesignVars(newDV)
        points = self.DVGeo.update("airfoil")[:,0:2]

        # Changing the directory to analysis folder
        os.chdir("{}/{}".format(directory, self.genSamples+1))

        if self.options["writeAirfoilCoordinates"]:
            self.DVGeo.foil.writeCoords("deformedAirfoil", coords=points, file_format="dat")

        if self.options["plotAirfoil"]:
            self._plotAirfoil(points)

        # Create input file
        self._creatInputFile(x)

        # Writing the surface mesh
        self.DVGeo.foil.writeCoords("surfMesh", points)

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

            if self.options["getFlowFieldData"]:
                # Reading the cgns file
                filename = self.options["aeroProblem"].name + "_surf.cgns"
                reader = pyvista.CGNSReader(filename)
                reader.load_boundary_patch = False

                # Reading the mesh
                mesh = reader.read()

                # Setting region for extraction
                if self.options["region"] == "surface":
                    mesh = mesh[0][0]
                else:
                    mesh = mesh[0][2]

                # Get the values
                fieldData = {}

                for index, var in enumerate(mesh.array_names):
                    # Skipping the first entry in the array
                    if index != 0:
                        # set_active_scalars returns a tuple, and second
                        # entry contains the pyvista numpy array.
                        fieldData[var] = np.asarray(mesh.set_active_scalars(var, "cell")[1])

            else:
                fieldData = None

            return output, fieldData
            
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

            Note: To use this method, atleast lower or upper surface CST
            coefficient should be added as a DV.
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

            if psutil.cpu_count(False) < options["noOfProcessors"] + 1:
                self._error("\"noOfProcessors\" requested is more than available processors.")

        ############ Validating refine
        if "refine" in userProvidedOptions:
            if not isinstance(options["refine"], int):
                self._error("\"refine\" attribute is not an integer.")

        ############ Validating alpha
        if "alpha" in userProvidedOptions:
            if not isinstance(options["alpha"], str):
                self._error("\"alpha\" attribute is not string.")

            if options["alpha"] not in ["explicit", "implicit"]:
                self._error("\"alpha\" attribute is not recognized. It can be either \"explicit\" or \"implicit\".")

            if options["alpha"] == "implicit":
                if "targetCL" in userProvidedOptions:
                    if not isinstance(options["targetCL"], float):
                        self._error("\"targetCL\" option is not float.")

        ############ Validating writeSliceFile
        if "writeSliceFile" in userProvidedOptions:
            if not isinstance(options["writeSliceFile"], bool):
                self._error("\"writeSliceFile\" attribute is not a boolean value.")

        ############ Validating directory attribute
        if "directory" in userProvidedOptions:
            if not isinstance(options["directory"], str):
                self._error("\"directory\" attribute is not string.")

        ############ Validating getFlowFieldData attribute
        if "getFlowFieldData" in userProvidedOptions:
            if not isinstance(options["getFlowFieldData"], bool):
                self._error("\"getFlowFieldData\" attribute is not a boolean value.")

            # Checking the other related options
            if "region" in userProvidedOptions:
                if not isinstance(options["region"], str):
                    self._error("\"region\" attribute is not a string.")

                if options["region"] not in ["surface", "field"]:
                    self._error("\"region\" attribute is not recognized. It can be either \"surface\" or \"field\".")

    def _checkDV(self, name: str, lb: float or np.ndarray, ub: float or np.ndarray) -> None:
        """
            Method for validating DV.
        """

        # List of possible DVs
        possibleDVs = ["upper", "lower", "N1", "N2", "N1_upper", "N2_upper", 
                        "N1_lower", "N2_lower", "alpha", "mach", "altitude"]

        # Validating name of the DV
        if not isinstance(name, str):
            self._error("Name argument is not a string.")

        # Checking if the DV is allowed
        if name not in possibleDVs:
            self._error("{} argument is not a valid DV.".format(name))
        
        # Checking if the DV is already added
        if name in self.DV:
            self._error("{} already exists.".format(name))

        # Checking if alpha can be added as a DV
        if name == "alpha":
            if self.options["alpha"] != "explicit":
                self._error("Alpha cannot be a design variable when \"alpha\" attribute in option is \"implicit\".")

        # Checking if these variables are initialized through aero problem or not
        if name in ["mach", "altitude"]:
            if name not in self.options["aeroProblem"].inputs.keys():
                self._error("You need to initialize \"{}\" in the aero problem to set it as design variable.".format(name))

        # Validating the bounds for "upper" variable
        if name == "upper":
            if not isinstance(lb, np.ndarray) or lb.ndim != 1:
                self._error("Lower bound for \"upper\" variable should be a 1D numpy array.")

            if not isinstance(ub, np.ndarray) or ub.ndim != 1:
                self._error("Upper bound for \"upper\" variable should be a 1D numpy array.")

            if len(lb) != self.options["numCST"][0]:
                self._error("Length of lower bound for \"upper\" variable is not equal to number of CST coeff for upper surface.")

            if len(ub) != self.options["numCST"][0]:
                self._error("Length of upper bound for \"upper\" variable is not equal to number of CST coeff for upper surface.")

        # Validating the bounds for "lower" variable
        elif name == "lower":
            if not isinstance(lb, np.ndarray) or lb.ndim != 1:
                self._error("Lower bound for \"lower\" variable should be a 1D numpy array.")

            if not isinstance(ub, np.ndarray) or ub.ndim != 1:
                self._error("Upper bound for \"lower\" variable should be a 1D numpy array.")

            if len(lb) != self.options["numCST"][1]:
                self._error("Length of lower bound for \"lower\" variable is not equal to number of CST coeff for lower surface.")

            if len(ub) != self.options["numCST"][1]:
                self._error("Length of upper bound for \"lower\" variable is not equal to number of CST coeff for lower surface.")

        # Validating lb and ub of the scalar DVs
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
        samples = lhs(dim, samples=numSamples, criterion='cm', iterations=1000)

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

        # Adding target Cl if alpha is implicit
        if self.options["alpha"] == "implicit":
            input["targetCL"] = self.options["targetCL"]

        # Saving the input file
        filehandler = open("input.pickle", "xb")
        pickle.dump(input, filehandler)
        filehandler.close()

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
