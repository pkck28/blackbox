from datgen.function import AeroStruct

aeroOptions = {
    # I/O Parameters
    # "gridFile": args.gridFile,
    # "outputDirectory": args.output,
    "monitorvariables": ["resrho", "cl", "cd"],
    "writeTecplotSurfaceSolution": True,
    # Physics Parameters
    "equationType": "RANS",
    # Solver Parameters
    "MGCycle": "sg",
    # ANK Solver Parameters
    "useANKSolver": True,
    # NK Solver Parameters
    "useNKSolver": True,
    "NKSwitchTol": 1e-4,
    # Termination Criteria
    "L2Convergence": 1e-6,
    "nCycles": 1000,
}


options = {
    "folderName":"PK",
}

# options = {}
test = AeroStruct(options=options)

print(test.options)