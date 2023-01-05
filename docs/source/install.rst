.. _adflow_install:

Installation
============

You will need following packages before intalling Blackbox:

- MDO Lab packages: ADflow, pyHyp, prefoil, pyGeo, baseclasses, cgnsutilities.
- Other packages: pyDOE2, numpy, scipy, mpi4py.

Follow the below steps for installation:

- Download the `repo <https://github.com/ComputationalDesignLab/blackbox>`_.
- Open the terminal and ``cd`` into the downloaded folder and run::

    pip install .

- If you want to make changes to the source code of the package, then install it using::

    pip install -e .

  You don't need to re-install every time you change in the source code.
