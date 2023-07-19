******************
AirfoilFFD Options
******************

There are various options which can be set for this class. Please read the entire list of options carefully:

**Required options**:

- ``solverOptions (dict)``: options for the flow solver (only adflow for now)
- ``meshingOptions (dict)``: options for the volume mesh generator (only pyHyp for now)
- ``airfoilFile (str)``: name of the .dat file with extension
- ``aeroProblem``: aero-problem from baseclasses package defining information related aerodynamic analysis
- ``nffd (int)``: total number of FFD points (upper and lower surface both)

**Optional options**:

- ``directory (str, default="output")``: name of the directory where the results will saved
- ``noOfProcessors (int, default=4)``: desired number of processors to run the analysis on
- ``refine (int, default=0)``: value of this options controls how much to refine or coarsen the generated volume mesh
  When the value is zero, there is no change to the volume mesh. When the value is 1 or 2, the volume mesh is refined
  by one level or two levels respetively. When the value is -1 and -2, the mesh is coarsened by similar levels
- ``writeSliceFile (bool, default=False)``: option to specify whether to add a slice in the mesh and output the slice file or not
- ``writeAirfoilCoordinates (bool, default=False)``: option to specify whether to output deformed airfoil coordinates (dat file) or not. You will find this 
  file in the analysis specific folder
- ``plotAirfoil (bool, default=False)``: option to specify whether to plot deformed airfoil or not. You will find the image in the analysis specific folder.
- ``writeDeformedFFD (boo, default=False)``: option to specify whether to output deformed FFD coordinates (dat file) or not. You will find this 
  file in the analysis specific folder
- ``getFlowFieldData (bool, default=False)``: option to specify whether to get the field data or not
- ``region (str, default="surface")``: this option only applies when ``getFlowFieldData`` is set to ``True``. This option decides from what
  region to extract the data. There are only two possible values: ``surface`` (will extract the field data at surface) and ``field`` 
  (will extract the entire field)
- ``fitted (bool, default=False)``: flag to pick between a fitted FFD (True) and box FFD (False)
- ``xmargin (float, 0.001)``: this option governs closest distance of the FFD box to the tip and aft of the airfoil
- ``ymarginu (float, 0.02)``: When a *box* ffd is generated this specifies the top of the box's y values as the maximum y 
  value in the airfoil coordinates plus this margin. When a *fitted* ffd is generated this is the margin between the FFD point 
  at an xslice location and the upper surface of the airfoil at this location
- ``ymarginl (float, 0.02)``: When a *box* ffd is generated this specifies the bottom of the box's y values as the minimum y 
  value in the airfoil coordinates minus this margin. When a *fitted* ffd is generated this is the margin between the FFD point 
  at an xslice location and the lower surface of the airfoil at this location
- ``alpha (str, default="explicit")``: option to specify whether to consider alpha as an explicit or implicit variable. There are only two possbile values:
  ``explicit`` (normal analysis) and ``implicit`` (internal root finder). When this option is set to implicit, then for each sample secant method
  is used to find alpha such that target CL is achived. So, this option also takes longer to evaluate. **Note**: When this option is set to implicit, then 
  ``alpha`` cannot be added as a DV
- ``targetCL (float, default=0.824)``: this option only applies when ``alpha`` is set to ``implicit``. 
  This option specifies the value of target CL to be met
- ``targetCLTol (float, default=1e-4)`` : this option only applies when ``alpha`` is set to ``implicit``. 
  This option specifies the required tolerance to be met for target CL
- ``startingAlpha (float, default=2.5)`` : this option only applies when ``alpha`` is set to ``implicit``.
  This option specifies the initial guess for angle of attack (in degrees) for secant method. 
  Note that this value has a huge impact on the convergence speed of the secant method. 
  So, it is recommended to set this value close to the expected value