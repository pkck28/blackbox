from blackbox import AirfoilCST
from baseclasses import AeroProblem

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
    "nCycles": 7000,
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
    "NKSwitchTol": 0.5e-6,
    "NKSubspaceSize": 400,
    "NKASMOverlap": 3,
    "NKPCILUFill": 4,
    "NKJacobianLag": 5,
    "NKOuterPreconIts": 3,
    "NKInnerPreconIts": 3,
    # Adjoint Parameters
    "AdjointSolver": "GMRES",
    "AdjointL2Convergence": 1e-12,
    "ADPC": True,
    "AdjointMaxIter": 5000,
    "AdjointSubspaceSize": 400,
    "ILUFill": 3,
    "ASMOverlap": 3,
    "OuterPreconIts": 3,
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
    "noOfProcessors": 9,
    "aeroProblem": ap,
    "airfoilFile": "rae2822.dat",
    "numCST": [6, 6],
    "meshingOptions": meshingOptions,
    "alpha": "implicit",
    "refine": 1,
    "plotAirfoil": True
}

# Example for generating samples
airfoil = AirfoilCST(options=options)

######### Multi Analysis

# Adding design variable
airfoil.addDV("lower", -0.3, 0.3)
airfoil.addDV("upper", -0.3, 0.3)

# Generating the samples
airfoil.generateSamples(5)
