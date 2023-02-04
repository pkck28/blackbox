.. _options_description:

*******
Options
*******

There are various options which can be set for Blackbox. Please read the entire list of options carefully:

**Required options**:

- ``solverOptions (dict)``: options for the flow solver (only adflow for now) .
- ``meshingOptions (dict)``: options for the volume mesh generator (only pyHyp for now).
- ``airfoilFile (str)``: name of the .dat file with extension.
- ``aeroProblem``: aero-problem from baseclasses package defining information related aerodynamic analysis.
- ``numCST (list)``: two-item list denoting the number of CST coefficients - first entry is for upper surface and second entry is for lower surface.

**Optional options**:

- ``directory (str, default="output")``: name of the directory where the results will saved.
- ``noOfProcessors (int, default=4)``: desired number of processors to run the analysis on.
- ``refine (int, default=0)``: value of this options controls how much to refine or coarsen the generated volume mesh.
  When the value is zero, there is no change to the volume mesh. When the value is 1 or 2, the volume mesh is refined
  by one level or two levels respetively. When the value is -1 and -2, the mesh is coarsened by similar levels.
- ``writeSliceFile (bool, default=False)``: option to specify whether to add a slice in the mesh and output the slice file or not.
- ``writeAirfoilCoordinates (bool, default=False)``: option to specify whether to output deformed airfoil coordinates (dat file) or not. You will find this 
  file in the analysis specific folder.
- ``plotAirfoil (bool, default=False)``: option to specify whether to plot deformed airfoil or not. You will find the image in the analysis specific folder.
- ``getFlowFieldData (bool, default=False)``: option to specify whether to get the field data or not.
- ``region (str, default="surface")``: this option only applies when ``getFlowFieldData`` is set to ``True``. This option decides from what
  region to extract the data. There are only two possible values: ``surface`` (will extract the field data at surface) and ``field`` 
  (will extract the entire field).
- ``alpha (str, default="explicit")``: option to specify whether to consider alpha as an explicit or implicit variable. There are only two possbile values:
  ``explicit`` (normal analysis) and ``implicit`` (internal optimization to find target CL). When this option is set to implicit, then for each sample a
  simple one variable optimization is performed to find alpha such that target CL is achived. **Note**: When this option is set to implicit, then ``alpha`` 
  cannot be added as a DV.
- ``targetCL (float, default=0.824)``: this option only applies when ``alpha`` is set to ``implicit``. This option specifies optimizer the value of target CL.
