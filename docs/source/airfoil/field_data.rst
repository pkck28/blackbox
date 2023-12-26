******************
Getting Field Data
******************

This tutorial is very similar to the :ref:`generating samples tutorial<Generating Samples>`, going through that is highly recommended before
proceeding ahead. There are couple of options which need to be set so that field data is also obtained.
Following snippet of the code shows the blackbox options dictionary with those two options::

    # Options for blackbox
    options = {
        "solverOptions": solverOptions,
        "directory": "multi",
        "noOfProcessors": 8,
        "aeroProblem": ap,
        "airfoilFile": "rae2822.dat",
        "numCST": [6, 6],
        "meshingOptions": meshingOptions,
        "getFlowFieldData": True, # Option which enables field data extraction
        "region": "field" # defines the region of extraction
    }

The ``surfaceVariables`` option in solver option dictionary controls which variables to extract. Please refer `adflow options <https://mdolab-adflow.readthedocs-hosted.com/en/latest/options.html#surfaceVariables>`_
for more details. Initializing the ``AirfoilCST`` class and adding design variables is similar to what is described in generating samples tutorial.

Accessing field output
----------------------

Generating the samples and reading all the output data is exactly same as described in generating samples tutorial. There is only one additional file:
``fieldData.mat``. It contains following quatities:

- **Input variable**: a 2D numpy array ``x`` in which each row represents a specific sample based on which analysis is performed. The number
  of rows will be usually equal to the number of samples argument in the ``generateSamples`` method. But, many times few of the analysis
  fail. It depends a lot on the solver and meshing options, so set those options after some tuning.

  .. note::
    The order of values in each row is based on how you add design variables. In this tutorial, first ``alpha`` is added as
    design variable. Then, lower and upper surface CST coefficients are added. Thus, first value in each row will be alpha, next 6
    values will be upper surface CST coefficients and last 6 will be lower surface CST coefficients.

- **Output variables**: it contains all the variables mention in the ``surfaceVariables`` option for solver. When the file is loaded, you will get a dictionary.
  The keys in the dictionary will depend on the entries in the ``surfaceVariables`` option. Value for each key is a numpy array. Refer below table for more
  information. *n* is the number of samples and *k* is the number of points (in surface or entire field).

  .. list-table::
    :widths: 25 50 25
    :header-rows: 1

    * - Possible value in ``surfaceVariables``
      - Matching key in dictionary
      - Size
    * - *rho*
      - *Density*
      - *n x k*
    * - *p*
      - *Pressure*
      - *n x k*
    * - *temp*
      - *Temperature*
      - *n x k*
    * - *cp*
      - *CoefPressure*
      - *n x k*
    * - *ptloss*
      - *RelativePressureStagnationLoss*
      - *n x k*
    * - *mach*
      - *Mach*
      - *n x k*
    * - *cf*
      - *SkinFrictionMagnitude*
      - *n x k*
    * - *ch*
      - *StantonNumber*
      - *n x k*
    * - *yplus*
      - *YPlus*
      - *n x k*
    * - *yplus*
      - *YPlus*
      - *n x k*
    * - *vx*, *vy*, *vz*
      - *Velocity*
      - *n x k x 3*
    * - *cfx*, *cfy*, *cfz*
      - *SkinFriction*
      - *n x k x 3*

  Following snippet shows how to access few of those variables::

    from scipy.io import loadmat
    data = loadmat("fieldData.mat") # mention the location of mat file

    x = data["x"]
    mach = data["Mach"]
    cp = data["CoefPressure"]
    velocity = data["Velocity"]
