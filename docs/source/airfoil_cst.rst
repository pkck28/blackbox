.. _airfoil_cst_tutorial:

********************************
Airfoil with CST Parametrization
********************************

The ``AirfoilCST`` module within Blackbox provides various methods for running analysis and 
generating data with CST Parametrization. Before starting to use this module, you will need the .dat file which 
contains the airfoil coordinates. There few important points to note regarding .dat file:

- The .dat file should follow the **selig** format i.e. the points should start from trailing edge and
  go in counter-clockwise direction and then back to trailing edge.
- The first and last point in the .dat file should be same (for both sharp and blunt trailing edge). 
  This ensures that the surface created using those points is closed. The coordinates obtained from
  the UIUC database usually doesn't close the loop, so check the .dat file properly before
  using with Blackbox.
- The leading and trailing edge of the airfoil coordinates should lie on the x-axis (best case scenario)
  or should be very close to it.

**Adflow** is the only supported solver. The volume mesh need by the solver is created internally using **pyHyp**. You will need following packages:

- **MDO Lab packages**: `ADflow <https://github.com/mdolab/adflow>`_, `pyHyp <https://github.com/mdolab/pyhyp>`_, `prefoil <https://github.com/mdolab/prefoil>`_, 
  `pyGeo <https://github.com/mdolab/pygeo>`_, `baseclasses <https://github.com/mdolab/baseclasses>`_, `cgnsutilities <https://github.com/mdolab/cgnsutilities>`_.
- **Other packages**: `pyDOE2 <https://github.com/clicumu/pyDOE2>`_, numpy, scipy, `mpi4py <https://github.com/mpi4py/mpi4py>`_.

Each subsection in this page will show how ``AirfoilCST`` module can be used. It is highly recommended to go through options description first
and then proceed to other sub-sections. The `rae2822.dat` file used in the tutorials can be found in ``examples/airfoil_cst/`` folder in the 
`repository <https://github.com/ComputationalDesignLab/blackbox>`_.

.. toctree::
   :maxdepth: 2

   options_description
   generating_samples
   generating_field_data

.. Along with generating a bunch of samples, you can also just run a single analysis like you would do normally.
.. You will have to use ``getObjectives`` method from the initialized object ``airfoil``. The methods needs one input
.. which is the value of design variable as a 1D numpy array. Following snippet shows how to run a single analysis::

..     import numpy as np

..     # Upper and lower surface CST coefficients
..     upper = np.array([0.12344046, 0.14809657, 0.14858145, 0.2168004, 0.17607825, 0.21018404])
..     lower = np.array([-0.13198943, -0.11895939, -0.22056435, -0.12743513, -0.08232715, 0.05055414])

..     # Creating input x
..     x = np.append(np.array([2.5]), lower)
..     x = np.append(x, upper)

..     # Run a single analysis
..     output = airfoil.getObjectives(x)

.. Note that here ``x`` is 1D numpy array with 13 entires. The values within the array follow the same order in which
.. design variables are added. ``output`` from the method is a dictionary which contains the same objective as described in
.. the previous section. Also, anlaysis specific folder will be created in the specificed direcrtory which contains similar output files as 
.. described in the previous section.
