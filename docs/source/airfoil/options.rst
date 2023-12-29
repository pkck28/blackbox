******************
Options
******************

Following is the list of options for ``AirfoilCST`` and ``AirfoilFFD`` modules. Options are divided
into three categories: CST related options, FFD related options and general options.

CST related options
--------------------

  ``numCST (list)``
    two-item list denoting the number of CST coefficients - first entry is for upper surface 
    and second entry is for lower surface. This is a required option when using CST parametrization

FFD related options
--------------------

  ``nFFD (int)``
    total number of FFD control points along the chord. This is a required option when using FFD parametrization

  ``fitted (bool, default=False)``
    option to decide whether to use a box FFD or a fitted FFD

  ``xmargin (float, default=0.001)``
    the closest distance of the FFD box to the tip and aft of the airfoil

  ``ymarginu (float, default=0.02)``
    When a box ffd is generated this specifies the top of the box's y values as the maximum y value in the airfoil 
    coordinates plus this margin. When a fitted ffd is generated this is the margin between the FFD point at an xslice 
    location and the upper surface of the airfoil at this location

  ``ymarginl (float, default=0.02)``
    When a box ffd is generated this specifies the bottom of the box's y values as the minimum y value in the airfoil
    coordinates minus this margin. When a fitted ffd is generated this is the margin between the FFD point at an xslice
    location and the lower surface of the airfoil at this location

  ``fixLETE (bool, default=True)``
    option to specify whether to fix the leading edge and trailing edge of the airfoil. This option is only applicable when 
    using FFD parametrization

  ``smoothing (bool, default=False)``
    option to specify whether to use laplacian smoothing or not when generating samples

  ``smoothingTheta (float, default=0.75)``
    option to control the amount of smoothing in an iteration

  ``smoothingMaxIterations (int, default=100)``
    maximum number of iterations for smoothing

  ``smoothingTolerance (float, default=5e-4)``
    stopping criteria for smoothing

  ``writeDeformedFFD (bool, default=False)``
    option to specify whether to write deformed FFD file or not
    
General options
--------------------
  
  ``airfoilFile (str)``
    name of the .dat file with extension. This is a requried option, irrespective of the parametrization

  ``solverOptions (dict, default={})``
    options for the flow solver (only adflow for now)

  ``meshingOptions (dict, default={})``
    options for the volume mesh generator (only pyHyp for now)

  ``aeroProblem (default=None)``
    aero-problem from baseclasses package defining information related aerodynamic analysis

  ``directory (str, default="output")``
    name of the directory where the results will saved

  ``noOfProcessors (int, default=4)``
    desired number of processors to run the analysis on

  ``refine (int, default=0)`` 
    value of this options controls how much to refine or coarsen the generated volume mesh
    When the value is zero, there is no change to the volume mesh. When the value is 1 or 2, the volume mesh is refined
    by one level or two levels respetively. When the value is -1 and -2, the mesh is coarsened by similar levels

  ``writeSliceFile (bool, default=False)``
    option to specify whether to add a slice in the mesh and write the slice file or not

  ``writeAirfoilCoordinates (bool, default=False)``
    option to specify whether to write deformed airfoil coordinates (dat file) or not

  ``plotAirfoil (bool, default=False)``
    option to specify whether to plot deformed airfoil or not

  ``getFlowFieldData (bool, default=False)``
    option to specify whether to get the field data or not

  ``region (str, default="surface")``
    this option only applies when ``getFlowFieldData`` is set to ``True``. This option decides from what
    region to extract the data. There are only two possible values: ``surface`` (will extract the field data at surface) and ``field`` 
    (will extract the entire field)

  ``alpha (str, default="explicit")``
    option to specify whether to consider alpha as an explicit or implicit variable. There are only two possbile values:
    ``explicit`` (normal analysis) and ``implicit`` (internal root finder). When this option is set to implicit, then for each sample secant method
    is used to find alpha such that target CL is achived. So, this option also takes longer to evaluate. **Note**: When this option is set to implicit, then 
    ``alpha`` cannot be added as a DV

  ``targetCL (float, default=0.824)``
    this option only applies when ``alpha`` is set to ``implicit``. 
    This option specifies the value of target CL to be met

  ``targetCLTol (float, default=1e-4)``
    this option only applies when ``alpha`` is set to ``implicit``. 
    This option specifies the required tolerance to be met for target CL

  ``startingAlpha (float, default=2.5)``
    this option only applies when ``alpha`` is set to ``implicit``. 
    This option specifies the initial guess for angle of attack (in degrees) 
    for secant method. Note that this value has a huge impact on the convergence speed of the secant method. 
    So, it is recommended to set this value close to the expected value

  ``samplingCriterion (str, default="cm")``
    this option decides which method to use to generate Latin Hypercube samples. Only four options are available:
    ``c``, ``m``, ``cm``, ``ese``

  ``randomState (int, default=None)``
    this option is used to set the random state while generating Latin Hypercube samples
