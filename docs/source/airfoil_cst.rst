.. _airfoil_cst_tutorial:

Airfoil with CST Parametrization
================================

The ``AirfoilCST`` module within Blackbox provides various methods for running analysis and 
generating data with CST Parametrization. Before starting to use this module, you will need the .dat file which 
contains the airfoil coordinates. There few important points to note regarding .dat file:

- The .dat file should follow the **selig** format i.e. the points should start from trailing edge and
  go in counter-clockwise direction.
- The first and last point in the .dat file should be same (for both sharp and blunt trailing edge). 
  This ensures that the surface created using those points is closed. The coordinates obtained from
  the UIUC database usually doesn't close the loop, so check the .dat file properly before
  using with Blackbox.
- The leading and trailing edge of the airfoil coordinates should lie on the x-axis (best case scenario)
  or should be very close to it.

The surface and volume mesh is created internally and then passed on to the solver. **Adflow** is the only
supported solver. According, **pyHyp** is the only supported volume mesh generator.

You will also need following packages:

- MDO Lab packages: ADflow, pyHyp, prefoil, pyGeo, baseclasses, cgnsutilities.
- Other packages: pyDOE2, numpy, scipy, mpi4py.

Setting up options
------------------

First thing to do is to setup the options dictionary. There are two kinds of options - required and optional. 
Required options are:

- ``solveroptions (dict)``: options for the flow solver (only adflow for now) .
- ``meshingOptions (dict)``: options for the volume mesh generator (only pyHyp for now).
- ``airfoilFile (str)``: name of the .dat file with extension.
- ``aeroProblem``: aero-problem from baseclasses package defining information related aerodynamic analysis.
- ``numCST (list)``: two-item list denoting the number of CST coefficients - first entry is for upper surface and second entry is for lower surface.

Optional options are:

- ``directory (str, default="output")``: name of the directory where the results will saved.
- ``noOfProcessors (int, default=4)``: desired number of processors to run the analysis on.
- ``slice (bool, default=True)``: adds a slice which will be used to output slice file after the analysis.
- ``refine (int, default=0)``: value of this options controls how much to refine or coarsen the generated volume mesh.
  When the value is zero, there is no change to the volume mesh. When the value is 1 or 2, the volume mesh is refined
  by one level or two levels respetively. When the value is -1 and -2, the mesh is coarsened by similar levels.
- ``getFlowFieldData (bool, default=False)``: flag to specify whether to get the field data or not.
- ``region``: this option only applies when ``getFlowFieldData`` is set to ``True``. This option decides from what
  region to extract the data. There are only two possible values: ``surface`` (will extract the field data at surface) and ``field`` 
  (will extract the entire field).

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

The `rae2822.dat` file used in the tutorial can be found in ``examples/airfoil_cst/`` folder in the 
`repository <https://github.com/ComputationalDesignLab/blackbox>`_. If you miss any requried 
options, then Blackbox will notify regarding missed options. Few options to avoid in solver and meshing options dict:

- ``gridFile`` (for Adflow) or ``inputFile`` (for pyHyp): since these are generated internally. 
- ``printAllOptions``, ``printIntro``, ``outputDirectory`` (for adflow).

**Note**: Having these options within the options dict will not raise error. These options are anyways over-ridden.

Next step is to import the ``AirfoilCST`` module from Blackbox and initialize it using the options dictionary::

    from blackbox import AirfoilCST
    airfoil = AirfoilCST(options=options)

Adding design variables
-----------------------

Now, ``airfoil`` object will be used for adding design variables. The ``addDV`` method needs three arguments:

- ``name``: the design variable to add. The available design variables are: 

    - ``upper``: CST coefficients of upper surface. The number of variables will be equal to first entry 
      in ``numCST`` list in options dictionary.
    - ``lower``: CST coefficients of lower surface. The number of variables will be equal to second entry 
      in ``numCST`` list in options dictionary.
    - ``N1``: First class shape variable for both upper and lower surface. Adds only variable for both surfaces.
    - ``N2``: Second class shape variable for both upper and lower surface. Adds only variable for both surfaces.
    - ``alpha``: Angle of attack for the analysis.
    - ``mach``: Mach number for the analysis.
    - ``altitude``: Altitude for the analysis.

- ``lowerBound``: lower bound for the variable. 
- ``upperBound``: upper bound for the variable.

.. note::
    Only for ``upper`` and ``lower`` variable, the lower and upper bound represent fraction change. For example, 
    if the lower bound for ``lower`` variable is -0.3, then the actual lower bound will be lower surface CST 
    coefficients decreased by 30%. Similarly, if the upper bound for ``upper`` variable is 0.2, then the actual 
    upper bound will be upper surface CST increased by 20%.

In this tutorial, ``alpha``, ``upper`` and ``lower`` are added as the bounds::

    airfoil.addDV("alpha", 2.0, 3.0)
    airfoil.addDV("lower", -0.3, 0.3)
    airfoil.addDV("upper", -0.3, 0.3)

Here, the upper and lower bound for ``lower`` variable is +30% and -30% of the lower surface CST coefficients respectively.
You can also remove a design varialbe using ``removeDV`` method. It takes only one input which is the name of the variable.

Generating samples
------------------

After adding design variables, you can either run a single analysis at a specific value of design variable or generate
data at bunch of design variables. Generating samples using Blackbox is very easy. You just need to use ``generateSamples`` 
method from the initialized object ``airfoil``. This method takes only one integer input which is the number of samples 
to be generated. Following snippet of the code will generate 10 samples::

    airfoil.generateSamples(10)

You can see the following output after completion of smaple generation process:

- A folder is created for each analysis in the specified folder. Each of the folder will contain ``log.txt``.
  There will be other files depending on the options provided to solver and blackbox.

- ``data.mat`` file which contains:

    - **Input variable**: a 2D numpy array ``x`` in which each row represents a specific sample based on which analysis is performed. The number
      of rows will be usually equal to the number of samples argument in the ``generateSamples`` method. But, many times few of the analysis
      fail. It depends a lot on the solver and meshing options, so set those options after some tuning.

      .. note::
          The order of values in each row is based on how you add design variables. In this tutorial, first ``alpha`` is added as
          design variable. Then, lower and upper surface CST coefficients are added. Thus, first value in each row will be alpha, next 6
          values will be upper surface CST coefficients and last 6 will be lower surface CST coefficients.

    - **Output variables**: There are two kinds of output variables - mandatory and user specificed. The ``evalFuncs`` argument in the aero problem
      decides the user desired variables. Along with these variables, `area` of the airfoil is the mandatory objective.

- ``description.txt``: contains various informations about the sample generation such as design variables, bounds, number of failed analysis, etc.

Following snippet shows how to access the data.mat file. In this tutorial, ``evalFuncs`` argument contains 
``cl``, ``cd``, ``cmz``. So, data.mat will contain these variables, along with ``area``::

    from scipy.io import loadmat
    data = loadmat("data.mat") # mention the location of mat file

    x = data["x"]
    cl = data["cl"]
    cd = data["cd"]
    cmz = data["cmz"]
    area = data["area"]

Running single analysis
-----------------------

Along with generating a bunch of samples, you can also just run a single analysis like you would do normally.
You will have to use ``getObjectives`` method from the initialized object ``airfoil``. The methods needs one input
which is the value of design variable as a 1D numpy array. Following snippet shows how to run a single analysis::

    import numpy as np

    # Upper and lower surface CST coefficients
    upper = np.array([0.12344046, 0.14809657, 0.14858145, 0.2168004, 0.17607825, 0.21018404])
    lower = np.array([-0.13198943, -0.11895939, -0.22056435, -0.12743513, -0.08232715, 0.05055414])

    # Creating input x
    x = np.append(np.array([2.5]), upper)
    x = np.append(x, lower)

    # Run a single analysis
    output = airfoil.getObjectives(x)

Note that here ``x`` is 1D numpy array with 13 entires. The values within the array follow the same order in which
design variables are added. ``output`` from the method is a dictionary which contains the same objective as described in
the previous section. Also, anlaysis specific folder will be created in the specificed direcrtory which contains similar output files as 
described in the previous section.

Getting Field data
------------------

You can also get the field data for each gnereated sample. You have to set ``getFlowFieldData`` option as ``True`` in the blackbox options dictionary -
please refer to setting up options section for more details. Following snippet shows how to set blackbox options for extracting field data::

  options = {
        "solverOptions": solverOptions,
        "directory": "multi",
        "noOfProcessors": 8,
        "aeroProblem": ap,
        "airfoilFile": "rae2822.dat",
        "numCST": [6, 6],
        "meshingOptions": meshingOptions,
        "getFlowFieldData": True,
        "region": "field"
  }

There is one more important option associated with extracting the field data - ``region``. This essentially describes what region to extract the data from. 
Refer setting up options section for more details. If the options are set properly, then folder for each analysis will have a file named ``fieldData.mat``.
Once the mat file is load, you will get a dict which contains all the ``surfaceVariables`` you mentioned (or set by default) for the solver. In this tutorial,
``surfaceVariables`` is not set in ``solverOptions``, so by default it contains coeeficient of pressure, mach number, and velocity. Following snippet shows how to 
read from the ``fieldData.mat`` file::

  from scipy.io import loadmat
  data = loadmat("fieldData.mat") # mention the location of mat file

  cp = data["CoefPressure"]
  mach = data["Mach"]
  cmz = data["velocity"]

Here, all the outputs will be 2D numpy array. For scalars values, the first dimension will be number of cells in the grid for field data and
second dimension will be 1. For vector values, first dimension will be same as scalar values, but second dimension will be three which represents x, y, and z
direction.
