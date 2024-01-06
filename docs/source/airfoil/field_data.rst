******************
Getting field data
******************

It is highly recommended to go through sample generation process with CST or FFD before proceeding ahead.
The process to generate field data along with scalar output is independent of parametrization method.
You need to set two options in blackbox options dictionary: ``getFlowFieldData`` and ``region``. The 
``getFlowFieldData`` option is a boolean which enables field data extraction while ``region`` option defines
region of extraction which can be either ``field`` or ``surface``. The ``field`` option extracts 
data from the entire field while ``surface`` option extracts data from the surface of airfoil only.
The ``surfaceVariables`` option in solver option dictionary controls which variables to extract. Please 
refer `ADflow options <https://mdolab-adflow.readthedocs-hosted.com/en/latest/options.html#surfaceVariables>`_
for more details. Rest all process is same as described in sample generation tutorial.

Reading the output data is exactly same as described in sample generation sections. There is only one additional file:
``fieldData.mat`` which contains following quatities:

- **Input variable**: a 2D numpy array ``x`` in which each row represents a specific sample based on which analysis is performed. The number
  of rows will be usually equal to the number of samples argument in the ``generateSamples`` method. But, many times few of the analysis
  fail. It depends a lot on the solver and meshing options, so set those options after some tuning.

  .. note::
    The order of values in each row is based on how you add design variables.

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
