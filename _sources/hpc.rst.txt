******************
Running on HPC
******************

Blackbox can also run on High Performance Computing (HPC) clusters and process to generate
samples is same as running on local machine. If you encouter a shared memory error while running
on HPC, it is recommended to use following command to launch your python file in slurm job script::

    mpirun -n 1 --mca pml ob1 --mca btl self,tcp python runscript.py

.. note::
    Only one processor is used to launch the python file. The number of processors used for
    meshing and flow solver is specified while initializing Blackbox object, refer sample 
    generation section for more details.

This command changes the mode of communication between the processors to TCP which will avoid the shared memory error.
