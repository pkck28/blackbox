from datgen import Aerodynamics
import time

aeroSolverOptions = {
            # I/O Parameters
            "gridFile": "wing_vol.cgns",
            "monitorvariables": ["resrho", "cl", "cd", "mach"],
            "writeTecplotSurfaceSolution": True,
            # Physics Parameters
            "equationType": "RANS",
            ######################################### Very important to change this according to the mesh
            ### Check the mesh and set the value
            "liftindex": 2,  # y is the lift direction
            # Solver Parameters
            "smoother": "DADI",
            "CFL": 0.5,
            "CFLCoarse": 0.25,
            "MGCycle": "sg",
            "MGStartLevel": -1,
            "nCyclesCoarse": 250,
            # ANK Solver Parameters
            "useANKSolver": True,
            "nsubiterturb": 5,
            "anksecondordswitchtol": 1e-4,
            "ankcoupledswitchtol": 1e-6,
            "ankinnerpreconits": 2,
            "ankouterpreconits": 2,
            "anklinresmax": 0.1,
            # Termination Criteria
            "L2Convergence": 1e-12,
            "L2ConvergenceCoarse": 1e-2,
            "nCycles": 1000,
}

varyingParameters = {
    "aoa" : {
        "lowerBound": 0,
        "upperBound": 10
    },
    "mach" : {
        "lowerBound": 0.6,
        "upperBound": 0.8
    }
}

fixedParameters = {
    "altitude" : 10000, # in m
    "areaRef" : 45.5, # in sq. m
    "chordRef" : 3.25, # in m
    # "aoa" : 2 # in deg
}

objectvies = ["cl", "cd"]

options = {
    "aeroSolverOptions": aeroSolverOptions,
    "fixedParameters" : fixedParameters,
    "varyingParameters" : varyingParameters,
    "numberOfSamples": 5,
    "directory" : "aero_training_multi",
    "noOfProcessors" : 10,
    "objectives" : objectvies,
    "samplingMethod" : "lhs"
}

test = Aerodynamics(options=options)

start = time.time()
test.generateSamples()
end = time.time()

print(str(end-start) + " seconds")

# start = time.time()

# sample = [2, 0.8]
# print(test.getObjectives(sample))

# sample = [3, 0.7]
# print(test.getObjectives(sample))

# sample = [0, 0.9]
# print(test.getObjectives(sample))

# end = time.time()

# print(str(end-start) + " seconds")
