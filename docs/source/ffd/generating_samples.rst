******************
Generating Samples
******************

This tutorial explains how to use ``AirfoilFFD`` module for generating samples.

Setting up options
------------------

First step involves creating options dictionary which is used for initializating the module. There are various 
options which can be set, please refer to :ref:`options<AirfoilFFD Options>` section for more details. Following snippet 
of the code shows an example::

    # Flow solver options
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

    # Volume meshing options
    meshingOptions = {
        # Input Parameters
        "unattachedEdgesAreSymmetry": False,
        "outerFaceBC": "farfield",
        "autoConnect": True,
        "BC": {1: {"jLow": "zSymm", "jHigh": "zSymm"}},
        "families": "wall",
        # Grid Parameters
        "N": 129,
        "s0": 1e-6,
        "marchDist": 100.0,
    }

    # Creating aero problem
    ap = AeroProblem(
        name="ap", alpha=2.0, mach=0.734, reynolds=6.5e6, reynoldsLength=1.0, T=288.15, 
        areaRef=1.0, chordRef=1.0, evalFuncs=["cl", "cd", "cmz"], xRef = 0.25, yRef = 0.0, zRef = 0.0
    )

    nffd = 20 # Number of FFD points

    # Options for blackbox
    options = {
        "solverOptions": solverOptions,
        "noOfProcessors": 8,
        "aeroProblem": ap,
        "airfoilFile": "rae2822.dat",
        "nffd": 20,
        "meshingOptions": meshingOptions,
        "writeAirfoilCoordinates": True,
        "plotAirfoil": True,
        "writeDeformedFFD": True,
        "fitted": True
    }

Next step is to import the ``AirfoilFFD`` module from Blackbox and initialize it using the options dictionary::

    from blackbox import AirfoilFFD
    airfoil = AirfoilFFD(options=options)
