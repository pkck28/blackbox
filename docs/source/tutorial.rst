.. _adflow_tutorial:

Tutorial
========

Airfoil with CST Parametrization
--------------------------------

First thing to do is to setup the options dictionary. There are two kinds of options - required and optional. 
Required options are:

- ``solveroptions (dict)``: options for the flow solver (only adflow for now) .
- ``meshingOptions (dict)``: options for the volume mesh generator (only pyHyp for now).
- ``airfoilFile (str)``: name of the .dat file with extension.
- ``aeroProblem``: aero-problem from baseclasses defining the problem.
- ``numCST (list)``: two-item list denoting the number of CST coefficients - first entry is for upper surface and second entry is for lower surface.

Optional options are:

- ``directory (str, default="output")``: name of the directory where the results will saved.
- ``noOfProcessors (int, default=4)``: desired number of processors to run the analysis on.
- ``refine (int, default=0)``: value of this options controls how much to refine or coarsen the generated volume mesh.
  When the value is zero, there is no change to the volume mesh. When the value is 1 or 2, the volume mesh is refined
  by one level or two levels respetively. When the value is -1 and -2, the mesh is coarsened by similar levels.

Following snippet of the code shows setting up of options dictionary::

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
        # Pseudo Grid Parameters
        "ps0": -1.0,
        "pGridRatio": -1.0,
        "cMax": 3.0,
        # Smoothing parameters
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
        "noOfProcessors": 8,
        "aeroProblem": ap,
        "airfoilFile": "rae2822.dat",
        "numCST": [6, 6],
        "meshingOptions": meshingOptions,
        "refine": 1
    }

The .dat file 

If you miss any requried option, then Blackbox will throw an error let you know the missed option. Few options
to avoid in solver and meshing options dict:

- ``gridFile`` (for Adflow) or ``inputFile`` (for pyHyp): since these are generated internally. 
- ``printAllOptions``, ``printIntro``, ``outputDirectory`` (for adflow): since this is mentioned through ``directory``` option.

**Note**: Having these options within the options dict will not raise error. These options are anyways over-ridden.