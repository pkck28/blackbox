.. _airfoil:

********************************
Airfoil sample generation
********************************

Blackbox provides two modules for generating airfoil samples depending on the parametrization. 
The ``AirfoilCST`` and ``AirfoilFFD`` module provides various methods for running analysis and 
generating data with CST and FFD Parametrization, respectively. Before starting to use these 
modules, a .dat file is needed which contains the airfoil coordinates. There few important 
points to note regarding .dat file:

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
    This requirement is only when using CST parametrization

.. toctree::
    :maxdepth: 3
    :caption: Table of Contents
    
    generating_samples
    field_data
    smoothing
    options
