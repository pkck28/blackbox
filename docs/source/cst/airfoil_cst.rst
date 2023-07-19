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

  .. note::
    Whenever the trailing edge is blunt, only two points should be there. More than two points in trailing edge will give an error.

**Adflow** is the only supported solver. The volume mesh need by the solver is created internally using **pyHyp**. You will need following packages:

- **MDO Lab packages**: `ADflow <https://github.com/mdolab/adflow>`_, `pyHyp <https://github.com/mdolab/pyhyp>`_, `prefoil <https://github.com/mdolab/prefoil>`_, 
  `pyGeo <https://github.com/mdolab/pygeo>`_, `baseclasses <https://github.com/mdolab/baseclasses>`_, `cgnsutilities <https://github.com/mdolab/cgnsutilities>`_.
- **Other packages**: `pyDOE2 <https://github.com/clicumu/pyDOE2>`_, numpy, scipy, `mpi4py <https://github.com/mpi4py/mpi4py>`_.

Each subsection in this page will show how ``AirfoilCST`` module can be used. It is highly recommended to go through options description first
and then proceed to other sub-sections. The `rae2822.dat` file used in the tutorials can be found in ``examples/airfoil_cst/`` folder in the 
`repository <https://github.com/ComputationalDesignLab/blackbox>`_.

.. toctree::
   :maxdepth: 2

   options_cst
   generating_samples
   generating_field_data
