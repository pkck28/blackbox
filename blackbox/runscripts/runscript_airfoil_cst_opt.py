############## Script file for running airfoil analysis.
# Imports
import pickle, time, os
from collections import OrderedDict
from mpi4py import MPI
from adflow import ADFLOW
from pyhyp import pyHyp
from cgnsutilities.cgnsutilities import readGrid
from pyoptsparse import Optimization, PSQP, SLSQP
from idwarp import USMesh

# Getting MPI comm
comm = MPI.COMM_WORLD
parent_comm = comm.Get_parent()

# Send the processor
parent_comm.send(os.getpid(), dest=0, tag=comm.rank)

# Redirecting the stdout
stdout = os.dup(1)
log = open("log.txt", "a")
os.dup2(log.fileno(), 1)

# Defining the variables
CL_target = 0.824
alpha = 2.8

try:
    ############## Reading input file for the analysis

    # Reading input file
    filehandler = open("input.pickle", 'rb') 
    input = pickle.load(filehandler)
    filehandler.close()

    # Getting aero problem from input file
    ap = input["aeroProblem"]
    refine = input["refine"]
    slice = input["writeSliceFile"]

    # Assigning non-shape DVs
    if "alpha" in input.keys():
        ap.alpha = input["alpha"][0]

    if "mach" in input.keys():
        ap.mach = input["mach"][0]

    if "altitude" in input.keys():
        ap.altitude = input["altitude"][0]

    # Getting solver and meshing options from input file
    solverOptions = input["solverOptions"]
    solverOptions["gridFile"] = "volMesh.cgns"
    solverOptions["liftindex"] = 2 # Always 2 since meshing is done internally

    meshingOptions = input["meshingOptions"]
    meshingOptions["inputFile"] = "surfMesh.xyz"

    ############## Generating mesh

    if comm.rank == 0:
        print("#" + "-"*129 + "#")
        print(" "*59 + "Meshing Log" + ""*59)
        print("#" + "-"*129 + "#")
        print("")

    hyp = pyHyp(options=meshingOptions, comm=comm)
    hyp.run()
    hyp.writeCGNS("volMesh.cgns")

    ############## Refining the mesh

    if comm.rank == 0:
        # Read the grid
        grid = readGrid("volMesh.cgns")

        # Only one processor has to do this
        if refine == 1:
            grid.refine(['i', 'k'])
        if refine == 2:
            grid.refine(['i', 'k'])
            grid.refine(['i', 'k'])
        if refine == -1:
            grid.coarsen()
        if refine == -2:
            grid.coarsen()
            grid.coarsen()

        grid.writeToCGNS("volMesh.cgns")
        
    else:
        # Other processors need to wait before starting analysis
        time.sleep(0.5)

    ############## Settign up adflow

    if comm.rank == 0:
        print("")
        print("#" + "-"*129 + "#")
        print(" "*59 + "Analysis Log" + ""*59)
        print("#" + "-"*129 + "#")
        print("")

    # Creating adflow object
    CFDSolver = ADFLOW(options=solverOptions, comm=comm)

    # Adding angle of attack as variable
    ap.addDV("alpha", value=alpha, lower=0, upper=4.0, scale=1.0)

    # Adding pressure distribution output
    if slice:
        CFDSolver.addSlices("z", 0.5, sliceType="absolute")

    # Getting triangulated surface mesh, later used in 
    # parent script to calculate volume
    trigSurfMesh = CFDSolver.getTriangulatedMeshSurface()

    ############## Methods for objective and sensitivity

    def Funcs(x):
        """
            Objective Function.
        """

        ap.setDesignVars(x)

        # Run CFD
        CFDSolver(ap)

        # Evaluate functions
        funcs = {}
        CFDSolver.evalFunctions(ap, funcs, ["cl"])
        CFDSolver.checkSolutionFailure(ap, funcs)

        # Calc objective
        funcs["obj"] = (funcs[ap["cl"]] - CL_target)**2
        
        return funcs

    def FuncsSens(x, funcs):
        """
            Function sensitivity.
        """

        # Evaluate sensitivities
        funcsSens = {}
        CFDSolver.evalFunctionsSens(ap, funcsSens, ["cl"])
        CFDSolver.checkAdjointFailure(ap, funcsSens)

        # Calc gradient of objective
        funcsSens["obj"] = OrderedDict()
        funcsSens["obj"]["alpha_ap"] = 2*(funcs["ap_cl"] - CL_target)*funcsSens["ap_cl"]["alpha_ap"]

        return funcsSens

    ############# Optimization

    # Creating optimization problem
    optProb = Optimization("opt", objFun=Funcs, comm=MPI.COMM_WORLD, sens=FuncsSens)

    # Add objective
    optProb.addObj("obj", scale=1e6)

    # Add variables from the AeroProblem
    ap.addVariablesPyOpt(optProb)

    SLSQP_options = {
        "ACC": 1e-4
    }

    PSQP_options = {
        "MET": 2,
        "XMAX": 0.25,
        "MFV": 12,
        "MIT": 12,
        "TOLX": 1e-5
    }

    # Creating optimizer for optimization
    # opt = SLSQP(options=SLSQP_options)
    opt = PSQP(options=PSQP_options)

    # Run Optimization
    sol = opt(optProb, sens=FuncsSens, storeHistory="opt.hst", sensMode="pgc")

    if MPI.COMM_WORLD.rank == 0:
        print(sol)

    ############# Post-processing

    # Run CFD at obtained alpha - it will be quick
    CFDSolver(ap)

    # Evaluating objectives
    funcs = {}
    CFDSolver.evalFunctions(ap, funcs)
    CFDSolver.checkSolutionFailure(ap, funcs)

    # Get the output

    # printing the result
    if MPI.COMM_WORLD.rank == 0:
        print("")
        print("#" + "-"*129 + "#")
        print(" "*59 + "Result" + ""*59)
        print("#" + "-"*129 + "#")
        print("")

        output = {}

        # Printing and storing results based on evalFuncs in aero problem
        for obj in ap.evalFuncs:
            print("{} = ".format(obj), funcs["{}_{}".format(ap.name, obj)])
            output["{}".format(obj)] = funcs["{}_{}".format(ap.name, obj)]

        # Other mandatory outputs
        print("fail = ", funcs["fail"])
        output["fail"] = funcs["fail"]

        # Getting triangulated surface points
        output["pts"] = trigSurfMesh

        # Storing the results in output file
        filehandler = open("output.pickle", "xb")
        pickle.dump(output, filehandler)
        filehandler.close()

except Exception as e:
    if comm.rank == 0:
        print(e)

finally:
    # Redirecting to original stdout
    os.dup2(stdout, 1)

    # close the file and stdout
    log.close()
    os.close(stdout)

    # Getting intercomm and disconnecting
    # Otherwise, program will enter deadlock
    parent_comm.Disconnect()
