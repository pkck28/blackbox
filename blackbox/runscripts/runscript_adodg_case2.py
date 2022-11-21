import pickle
from mpi4py import MPI
from baseclasses import AeroProblem
from adflow import ADFLOW
from idwarp import USMesh
from pygeo import DVGeometry, DVConstraints, geo_utils

############## Reading input file for the analysis
filehandler = open("input.pickle", 'rb') 
input = pickle.load(filehandler)
filehandler.close()

# defining the variables
aero_options = input["aeroSolverOptions"]
aero_options["gridFile"] = "grid.cgns"
evalFuncs = input["objectives"]

############## Settign up adflow
CFDSolver = ADFLOW(options=aero_options, comm=MPI.COMM_WORLD)

# Mesh warping
meshOptions = {"gridFile": aero_options["gridFile"]}
mesh = USMesh(options=meshOptions, comm=MPI.COMM_WORLD)
CFDSolver.setMesh(mesh)

# Adding pressure distribution output
CFDSolver.addSlices("z", 0.5, sliceType="absolute")

# Creating aeroproblem for adflow
ap = AeroProblem(
    name="ap", alpha=2.79, mach=0.734, reynolds=6.5e6, reynoldsLength=1.0, T=288.15, 
    areaRef=1.0, chordRef=1.0, evalFuncs=evalFuncs, xRef = 0.25, yRef = 0.0, zRef = 0.0
)

############## Settign up geometry parametrization
FFDFile = "ffd.xyz"
DVGeo = DVGeometry(FFDFile)
pts = DVGeo.getLocalIndex(0) # location of all the ffd points

# Adding ffd at root as local dv
indexList1 = pts[:, :, 0].flatten()
PS1 = geo_utils.PointSelect("list", indexList1)
DVGeo.addLocalDV("shape0", lower=-0.03, upper=0.03, axis="y", scale=1.0, pointSelect=PS1)

# Adding ffd at tip as local dv
indexList2 = pts[:, :, 1].flatten()
PS2 = geo_utils.PointSelect("list", indexList2)
DVGeo.addLocalDV("shape1", lower=-0.03, upper=0.03, axis="y", scale=1.0, pointSelect=PS2)

# Getting the dv dictionary
designVariables = DVGeo.getValues()

# Add DVGeo object to CFD solver
CFDSolver.setDVGeo(DVGeo)

############## Settign up constraints
DVCon = DVConstraints()
DVCon.setDVGeo(DVGeo)
DVCon.setSurface(CFDSolver.getTriangulatedMeshSurface())

# Defining the four corners of the 2D plane for volume constraint
le = 1e-4
leList = [[le, 0, le], [le, 0, 1.0 - le]]
teList = [[1.0 - le, 0, le], [1.0 - le, 0, 1.0 - le]]
DVCon.addVolumeConstraint(leList, teList, 2, 100, lower=0.5, upper=3.0, scaled=True)

############## Settign up dv
designVariables["shape0"] = input["shape"]
designVariables["shape1"] = input["shape"]

DVGeo.setDesignVars(designVariables)
ap.setDesignVars(designVariables)

designVariables = DVGeo.getValues()

############## Solving for the CL
CFDSolver.solveCL(ap, 0.824, tol=0.0001)

# storing all the requried results
funcs = {}
DVCon.evalFunctions(funcs)
CFDSolver.evalFunctions(ap, funcs)
CFDSolver.checkSolutionFailure(ap, funcs)

# Writing the simulation results
CFDSolver.writeSolution(baseName="sample", number=input["sampleNo"])

# printing the result
if MPI.COMM_WORLD.rank == 0:

    output = {}
    funcs["alpha"] = ap.alpha

    print("\n------------------- Result -------------------")
    print("alpha = ", funcs["alpha"])
    output["alpha"] = funcs["alpha"]

    print("cl = ", funcs["ap_cl"])
    output["cl"] = funcs["ap_cl"]
    
    print("cd = ", funcs["ap_cd"])
    output["cd"] = funcs["ap_cd"]

    print("cmz = ", funcs["ap_cmz"])
    output["cmz"] = funcs["ap_cmz"]

    print("area = ", funcs["DVCon1_volume_constraint_0"])
    output["area"] = funcs["DVCon1_volume_constraint_0"]

    filehandler = open("output.pickle", "xb")
    pickle.dump(output, filehandler)
    filehandler.close()
