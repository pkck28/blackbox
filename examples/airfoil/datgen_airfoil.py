from datgen import Aerodynamics

aeroSolverOptions = {
    # Common Parameters
    "gridFile": "rae2822.cgns",
    "monitorvariables": ["resrho", "cl", "cd", "mach"],
    "writeTecplotSurfaceSolution": True,
    # Physics Parameters
    "liftindex": 2,  # y is the lift direction
    # Physics Parameters
    "equationType": "RANS",
    # Solver Parameters
    "smoother": "DADI",
    "CFL": 0.5,
    "CFLCoarse": 0.25,
    "MGCycle": "sg",
    "MGStartLevel": -1,
    "nCyclesCoarse": 250,
    # ANK Solver Parameters
    "useNKSolver": True,
    "useanksolver": True,
    "nsubiterturb": 5,
    "anksecondordswitchtol": 1e-4,
    "ankcoupledswitchtol": 1e-6,
    "ankinnerpreconits": 2,
    "ankouterpreconits": 2,
    "anklinresmax": 0.1,
    "infchangecorrection": True,
    # Termination Criteria
    "L2Convergence": 1e-15,
    "L2ConvergenceCoarse": 1e-4,
    "nCycles": 20000
}

varyingParameters = {
    "shape" : {
        "lowerBound": -0.05,
        "upperBound": 0.05,
        "numberOfVariables": 12
    },
    "aoa" : {
        "lowerBound": 0,
        "upperBound": 10
    }
}

fixedParameters = {
    "altitude" : 10000, # in m
    "areaRef" : 1.0, # in sq. m
    "chordRef" : 1.0, # in m
    "mach" : 0.8
}

objectvies = ["cl", "cd"]

options = {
    "aeroSolverOptions": aeroSolverOptions,
    "fixedParameters" : fixedParameters,
    "varyingParameters" : varyingParameters,
    "numberOfSamples": 2,
    "directory" : "multi",
    "noOfProcessors" : 8,
    "objectives" : objectvies,
    "samplingMethod" : "lhs",
    "ffdFile" : "ffd.xyz",
    "shape" : "airfoil"
}

# Example for generating samples
airfoil = Aerodynamics(options=options)

airfoil.generateSamples()
