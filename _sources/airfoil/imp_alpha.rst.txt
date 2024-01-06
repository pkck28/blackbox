***************************************
Implicit alpha
***************************************

It is highly recommended to go through sample generation process with CST or FFD before proceeding ahead.
Aerodynamic shape optimization problems are often formulated as lift constrained drag minimization problems.
Accordingly, samples can also be generated such that lift constraint is satisfied. Blackbox provides
a way to do this by using an implicit alpha formulation. Under this formulation, analysis for each sample
consists of simple secant method which is used for triming the airfoil to the desired lift coefficient.
Due to this, implicit alpha formulation is also computationally more expensive than standard sample generation process.
ADflow implements this formulation internally, please refer to ADflow's `documentation <https://mdolab-adflow.readthedocs-hosted.com/en/latest/API.html#adflow.pyADflow.ADFLOW.solveCL>`_ 
for more details about the secant method.

The process to generate samples using implicit alpha formulation is similar to the process described in
sample generation section and is independent of parametrization methods. To use implicit alpha formulation, 
the user needs to set ``alpha`` option as ``implicit`` while initializing the Blackbox. There are three more 
options that can be  provided: ``targetCL``, ``targetCLtol`` and ``startingAlpha``, please refer :ref:`options<options>`
for more details.

.. note::
    ``alpha`` cannot be a design variable when using implicit alpha formulation since it is internally
    modified to match desired lift coefficient.
