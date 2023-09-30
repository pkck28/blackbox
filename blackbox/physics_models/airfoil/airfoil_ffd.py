# Imports
import os, sys, shutil, pickle, psutil
import numpy as np
from scipy import integrate
from mpi4py import MPI
from pygeo import DVGeometry
from prefoil import Airfoil
from prefoil.utils import readCoordFile
from smt.sampling_methods import LHS
from ...base.airfoilbaseclass import AirfoilBaseClass, DefaultOptions

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

class AirfoilFFDDefaultOptions(DefaultOptions):
    """
        Class creates a default option which are later 
        edited/appended with user provided options.
    """

    def __init__(self):

        # Initializing the parent class
        super().__init__()

        # Other options
        self.writeDeformedFFD = False

        # FFD related options
        self.fitted = False
        self.xmargin = 0.001
        self.ymarginu = 0.02
        self.ymarginl = 0.02

class AirfoilFFD(AirfoilBaseClass):
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
        self._getDefaultOptions(AirfoilFFDDefaultOptions())

        # Setting up the required options list
        defaultOptions = list(self.options.keys())
        requiredOptions = ["solverOptions", "meshingOptions", "airfoilFile", "aeroProblem", "nffd"]

        # Validating user provided options
        self._checkOptions(defaultOptions, requiredOptions, options)

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

        # Read the coordinate file
        self.coords = readCoordFile(self.options["airfoilFile"])
        airfoil = Airfoil(self.coords)

        # Creating FFD box
        airfoil.generateFFD(nffd=int(self.options["nffd"]/2), filename=directory + "/ffd", fitted=self.options["fitted"], 
                            xmargin=self.options["xmargin"], ymarginu=self.options["ymarginu"], 
                            ymarginl=self.options["ymarginl"], coords=self.coords)

        # Creating DVGeometry object
        self.DVGeo = DVGeometry(directory + "/ffd.xyz")
        self.parametrization = "FFD"

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

        # Creating sampler based on internal sampling
        if self.options["sampling"] == "internal":
            
            # Limits for sampler
            xlimits = np.hstack((self.lowerBound.reshape(-1,1), self.upperBound.reshape(-1,1)))

            # Creating the sampler
            self.sampler = LHS(xlimits=xlimits, criterion=self.options["samplingCriterion"], random_state=self.options["randomState"])

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

        # Updating the sampler based on internal sampling
        if self.options["sampling"] == "internal":

            if len(self.DV) == 0:
                delattr(self, "sampler")
            else:
                # Limits for sampler
                xlimits = np.hstack((self.lowerBound.reshape(-1,1), self.upperBound.reshape(-1,1)))

                # Creating the sampler
                self.sampler = LHS(xlimits=xlimits, criterion=self.options["samplingCriterion"], random_state=self.options["randomState"])

    # ----------------------------------------------------------------------------
    #                       Methods related to analysis
    # ----------------------------------------------------------------------------

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
        # Setting filepath based on the how alpha is treated alpha
        if self.options["alpha"] == "explicit":
            filepath = os.path.join(pkgdir, "runscripts/airfoil/runscript_airfoil.py")
        else:
            # filepath = os.path.join(pkgdir, "runscripts/airfoil/runscipt_airfoil_cst_opt.py")
            filepath = os.path.join(pkgdir, "runscripts/airfoil/runscript_airfoil_rf.py")

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
            self._plotAirfoil(plt, self.coords, points)

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

    # ----------------------------------------------------------------------------
    #                       Methods related to validation
    # ----------------------------------------------------------------------------

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

        # Checking if alpha can be added as a DV
        if name == "alpha":
            if self.options["alpha"] != "explicit":
                self._error("Alpha cannot be a design variable when \"alpha\" attribute in option is \"implicit\".")

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
