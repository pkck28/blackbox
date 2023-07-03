from blackbox import WingFFD
from baseclasses import AeroProblem
import numpy as np

solverOptions = {
    # Common Parameters
    "gridFile": "crm_volMesh.cgns",
    "monitorvariables": ["cl", "cd", "cmy", "yplus"],
    "writeTecplotSurfaceSolution": True,
    "writeSurfaceSolution": False,
    "writeVolumeSolution": False,
    # Physics Parameters
    "equationType": "RANS",
    "smoother": "DADI",
    "MGCycle": "sg",
    "nsubiterturb": 10,
    "nCycles": 7000,
    "liftIndex": 3, # Very important
    # ANK Solver Parameters
    "useANKSolver": True,
    "ANKSubspaceSize": 400,
    "ANKASMOverlap": 3,
    "ANKPCILUFill": 4,
    "ANKJacobianLag": 5,
    "ANKOuterPreconIts": 3,
    "ANKInnerPreconIts": 3,
    # NK Solver Parameters
    "useNKSolver": True,
    "NKSwitchTol": 1e-6,
    "NKSubspaceSize": 400,
    "NKASMOverlap": 3,
    "NKPCILUFill": 4,
    "NKJacobianLag": 5,
    "NKOuterPreconIts": 3,
    "NKInnerPreconIts": 3,
    # Termination Criteria
    "L2Convergence": 1e-14
}

# Creating aeroproblem for adflow
ap = AeroProblem(
    name="crm", alpha=2.2, mach=0.85, reynolds=5e6, reynoldsLength=1.0, T=298.15,
    areaRef=3.407014, chordRef=1.0, evalFuncs=["cl", "cd", "cmy"], xRef = 1.2077, yRef = 0.0, 
    zRef = 0.007669
)

options = {
    "solverOptions": solverOptions,
    "ffdFile": "crm_ffd.xyz",
    "aeroProblem": ap,
    "noOfProcessors": 2,
    "sliceLocation": [0.883, 1.003, 2.093, 2.612, 3.112, 3.548],
    "writeDeformedFFD": True
}

# Create the wing object
wing = WingFFD(options=options)

# Add alpha as a design variable
wing.addDV("alpha", lowerBound=1.5, upperBound=3.5)

# Add the wing shape as a design variable
lowerBound = np.array([-0.01]*wing.nffd)
upperBound = np.array([0.01]*wing.nffd)
wing.addDV("shape", lowerBound=lowerBound, upperBound=upperBound)

# Add the wing twist as a design variable
lowerBound = np.array([-2.0]*wing.nTwist)
upperBound = np.array([2.0]*wing.nTwist)
wing.addDV("twist", lowerBound=lowerBound, upperBound=upperBound)

# Generate samples
wing.generateSamples(5)
