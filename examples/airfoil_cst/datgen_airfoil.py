from blackbox import AirfoilCST
from baseclasses import AeroProblem
import numpy as np

solverOptions = {
    # Common Parameters
    "monitorvariables": ["cl", "cd", "cmz", "yplus"],
    "writeTecplotSurfaceSolution": True,
    "writeSurfaceSolution": False,
    "writeVolumeSolution": False,
    # Physics Parameters
    "equationType": "RANS",
    "smoother": "DADI",
    "MGCycle": "sg",
    "nsubiterturb": 10,
    "nCycles": 10000,
    # ANK Solver Parameters
    "useANKSolver": True,
    "ANKJacobianLag": 5,
    "ANKPhysicalLSTol": 0.25,
    "ANKOuterPreconIts": 2,
    "ANKInnerPreconIts": 2,
    "ANKASMOverlap": 2,
    "ANKSecondOrdSwitchTol": 1e-3,
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
    "L2Convergence": 1e-14,
    "L2ConvergenceCoarse": 1e-4
}

meshingOptions = {
    # ---------------------------
    #        Input Parameters
    # ---------------------------
    "unattachedEdgesAreSymmetry": False,
    "outerFaceBC": "farfield",
    "autoConnect": True,
    "BC": {1: {"jLow": "zSymm", "jHigh": "zSymm"}},
    "families": "wall",
    # ---------------------------
    #        Grid Parameters
    # ---------------------------
    "N": 129,
    "s0": 1e-6,
    "marchDist": 100.0,
    # ---------------------------
    #   Pseudo Grid Parameters
    # ---------------------------
    "ps0": -1.0,
    "pGridRatio": -1.0,
    "cMax": 3.0,
    # ---------------------------
    #   Smoothing parameters
    # ---------------------------
    "epsE": 1.0,
    "epsI": 2.0,
    "theta": 3.0,
    "volCoef": 0.25,
    "volBlend": 0.0001,
    "volSmoothIter": 100,
}

# Creating aeroproblem for adflow
ap = AeroProblem(
    name="ap", alpha=2.0, mach=0.734, reynolds=6.5e6, reynoldsLength=1.0, T=288.15, 
    areaRef=1.0, chordRef=1.0, evalFuncs=["cl", "cd", "cmz"], xRef = 0.25, yRef = 0.0, zRef = 0.0
)

# Options for blackbox
options = {
    "solverOptions": solverOptions,
    "directory": "multi",
    "noOfProcessors": 11,
    "aeroProblem": ap,
    "airfoilFile": "rae2822.dat",
    "numCST": [6, 6],
    "meshingOptions": meshingOptions,
    "refine": 2
}

# Example for generating samples
airfoil = AirfoilCST(options=options)

######### Multi Analysis

# Adding design variable
airfoil.addDV("alpha", 2.0, 3.0)
airfoil.addDV("lower", -0.3, 0.3)
airfoil.addDV("upper", -0.3, 0.3)

# Generating the samples
airfoil.generateSamples(10)

######### Single Analysis

# upper = np.array([0.12344046, 0.14809657, 0.14858145, 0.2168004, 0.17607825, 0.21018404])

# lower = np.array([-0.13198943, -0.11895939, -0.22056435, -0.12743513, -0.08232715, 0.05055414])

# # Adding design variable

# # airfoil.addDV("upper", -0.3, 0.3)
# # airfoil.addDV("lower", -0.3, 0.3)
# airfoil.addDV("alpha", 2.0, 3.0)

# # x = np.append(upper, lower)

# # x = np.append(x, 2.5)

# x = np.array([2.5])

# output = airfoil.getObjectives(x)

# print(output)
