from datgen import AeroStruct

aeroOptions = {
            # I/O Parameters
            "gridFile": "wing_vol.cgns",
            "monitorvariables": ["resrho", "mach", "cl", "cd"],
            "writeTecplotSurfaceSolution": True,
            'writevolumesolution':False,
            # 'writesurfacesolution':False,
            # Physics Parameters
            "equationType": "RANS",
            ######################################### Very important to change this according to the mesh
            ### Check the mesh and set the value
            "liftindex": 2,  # y is the lift direction
            # Solver Parameters
            "smoother": "DADI",
            "CFL": 1.5,
            "CFLCoarse": 1.25,
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
            "L2Convergence": 1e-14,
            "L2ConvergenceCoarse": 1e-2,
            "L2ConvergenceRel": 1e-4,
            "nCycles": 5000,
            # force integration
            "forcesAsTractions": False, # Using MELD, If using RLT, then True
}

designVariables = {
    "mach" : {
        "lowerBound": 0.6,
        "upperBound": 0.85
    },
    "aoa" : {
        "lowerBound": 0,
        "upperBound": 5
    }
}

options = {
    "aeroSolverOptions": aeroOptions,
    "designVariables": designVariables,
    "numberOfSamples": 2
}

test = AeroStruct(options=options)

test.generateSamples()
