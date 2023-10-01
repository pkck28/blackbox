# Imports
import os
import numpy as np
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
        self._getDefaultOptions(DefaultOptions())

        # Setting up the required options list
        defaultOptions = list(self.options.keys())
        requiredOptions = ["airfoilFile", "nffd"]

        # Validating user provided options
        self._checkOptions(defaultOptions, requiredOptions, options)

        # Updating/Appending the default option list with user provided options
        self._setOptions(options)

        # Raise an error if pyvista is not installed
        if self.options["getFlowFieldData"]:
            if msg_pyvista != None:
                self._error(msg_pyvista)
            else:
                self.pyvista = pyvista

        # Raise an error if matplotlib is not installed
        if self.options["plotAirfoil"]:
            if msg_matplotlib != None:
                self._error(msg_matplotlib)
            else:
                self.plt = plt

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
