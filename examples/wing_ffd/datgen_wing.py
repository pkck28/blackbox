from blackbox import WingFFD
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

ap = AeroProblem(name="wing", alpha=2.5, mach=0.85, altitude=10000, areaRef=45.5, chordRef=3.25, evalFuncs=["cl", "cd", "cmz"])

options = {
    "solverOptions": solverOptions,
    "volumeMesh": "wing_vol.cgns",
    "ffdFile": "wing_ffd.xyz",
    "aeroProblem": ap,
    "noOfProcessors": 11,
}

# Create the wing object
wing = WingFFD(options=options)

# Add alpha as a design variable
wing.addDV("alpha", lowerBound=2.0, upperBound=3.0)

# Add the wing shape as a design variable
wing.addDV("shape", lowerBound=-0.03, upperBound=0.03)

wing.generateSamples(2)
