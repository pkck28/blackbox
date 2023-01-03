############## Script file for running airfoil analysis.
# Imports
import pickle, os, time
from mpi4py import MPI
from adflow import ADFLOW
from pyhyp import pyHyp

comm = MPI.COMM_WORLD

try:
    ############## Reading input file for the analysis

    # Reading input file
    filehandler = open("input.pickle", 'rb') 
    input = pickle.load(filehandler)
    filehandler.close()

    # Getting aero problem from input file
    ap = input["aeroProblem"]
    refine = input["refine"]

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
    hyp = pyHyp(options=meshingOptions, comm=comm)
    hyp.run()
    hyp.writeCGNS("volMesh.cgns")

    ############## Refining the mesh

    if comm.rank == 0:
        # Only one processor has to do this
        if refine == 1:
            os.system("cgns_utils refine volMesh.cgns volMesh.cgns --axes ['i', 'k']")
        if refine == 2:
            os.system("cgns_utils refine volMesh.cgns volMesh.cgns --axes ['i', 'k']")
            os.system("cgns_utils refine volMesh.cgns volMesh.cgns --axes ['i', 'k']")

        if refine == -1:
            os.system("cgns_utils coarsen volMesh.cgns volMesh.cgns")
        if refine == -2:
            os.system("cgns_utils coarsen volMesh.cgns volMesh.cgns")
            os.system("cgns_utils coarsen volMesh.cgns volMesh.cgns")
        
    else:
        # Other processors need to wait before starting analysis
        time.sleep(0.5)

    ############## Settign up adflow

    # Creating adflow object
    CFDSolver = ADFLOW(options=solverOptions, comm=comm)

    # Adding pressure distribution output
    CFDSolver.addSlices("z", 0.5, sliceType="absolute")

    # Getting triangulated surface mesh, later used in 
    # parent script to calculate volume
    trigSurfMesh = CFDSolver.getTriangulatedMeshSurface()

    ############## Run CFD
    CFDSolver(ap)

    ############## Evaluating objectives
    funcs = {}
    CFDSolver.evalFunctions(ap, funcs)
    CFDSolver.checkSolutionFailure(ap, funcs)

    ############# Post-processing

    # printing the result
    if MPI.COMM_WORLD.rank == 0:

        output = {}

        print("\n------------------- Result -------------------")

        # Printing and storing results based on evalFuncs in aero problem
        for obj in ap.evalFuncs:
            print("{} = ".format(obj), funcs["{}_{}".format(ap.name, obj)])
            output["{}".format(obj)] = funcs["{}_{}".format(ap.name, obj)]

        # Other mandatory outputs
        print("fail = ", funcs["fail"])
        output["fail"] = funcs["fail"]

        output["pts"] = trigSurfMesh

        # Storing the results in output file
        filehandler = open("output.pickle", "xb")
        pickle.dump(output, filehandler)
        filehandler.close()

except:
    pass

finally:
    # Getting intercomm and disconnecting
    # Otherwise, program will enter deadlock
    parent_comm = comm.Get_parent()
    parent_comm.Disconnect()
