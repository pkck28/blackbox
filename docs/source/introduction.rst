.. _adflow_introduction:

Introduction
============

High fidelity sample generation requires creating a complex analysis chain
which is time consuming. Blackbox is aimed at easing this process of sample 
generation. **Right now, only airfoil data generation with CST parameterization
is supported**. Work is under progress on various other parameterization techniques
and data generation for wing and aero-structural cases.

Airfoil Data Generation
-----------------------

For generating airfoil data, you will just need the .dat file which contains the airfoil
coordinates. There few important points to note with the .dat file:

- The .dat file should follow the **selig** format i.e. the points should start from trailing edge and
  go in counter-clockwise direction.
- The first and last point in the .dat file should be same (for both sharp and blunt trailing edge). 
  This ensures that the surface created using those points is closed. The coordinates obtained from
  the UIUC database usually doesn't close the loop, so check the .dat file properly before
  using with Blackbox.
- The leading and trailing edge of the airfoil coordinates should lie on the x-axis (best case scenario)
  or should be very close to it.

The surface and volume mesh is created internally and then passed on to the solver. Adflow is the only
supported solver. According, pyHyp is the only supported volume mesh generator.
