############## Script file for generating ADODG case 2 samples.

# General imports
import pickle
from collections import OrderedDict
from mpi4py import MPI

# MDO lab imports
from baseclasses import AeroProblem
from adflow import ADFLOW
from pygeo import DVGeometry, DVConstraints, geo_utils
from pyoptsparse import Optimization, SLSQP, PSQP
from idwarp import USMesh

# Importing warnings package and supressing the warnings
import warnings
warnings.filterwarnings("ignore")

############## Reading input file for the analysis

# Reading input file
filehandler = open("input.pickle", 'rb') 
input = pickle.load(filehandler)
filehandler.close()

# Defining the variables
aero_options = input["aeroSolverOptions"]
evalFuncs = input["objectives"]
CL_target = 0.824
alpha = 2.79

############## Mesh Deformation using PyGeo and IDWarp

idwarp_options = {
    "gridFile": "grid.cgns",
}

# Create the idwarp object
mesh = USMesh(options=idwarp_options)

# Extract all undeformed surface mesh coordinates
coords0 = mesh.getSurfaceCoordinates()

# Initializing the DVGeometry class with the FFD file
DVGeo = DVGeometry("ffd.xyz")

# Getting the ffd points
pts = DVGeo.getLocalIndex(0) # location of all the ffd points

# Adding ffd at root as local dv
indexList1 = pts[:, :, 0].flatten()
PS1 = geo_utils.PointSelect("list", indexList1)
DVGeo.addLocalDV("shape0", lower=-0.03, upper=0.03, axis="y", scale=1.0, pointSelect=PS1)

# Adding ffd at tip as local dv
indexList2 = pts[:, :, 1].flatten()
PS2 = geo_utils.PointSelect("list", indexList2)
DVGeo.addLocalDV("shape1", lower=-0.03, upper=0.03, axis="y", scale=1.0, pointSelect=PS2)

# Adding surface mesh co-ordinates as a pointset
DVGeo.addPointSet(coords0, "airfoil_surface_mesh")

# Get Design Variables
dvDict = DVGeo.getValues()

# Assigning DV
dvDict["shape0"] = input["shape"]
dvDict["shape1"] = input["shape"]

# Set Design Variables
DVGeo.setDesignVars(dvDict)

# Updating suface mesh co-ordinates
newCoords = DVGeo.update("airfoil_surface_mesh")
DVGeo.writePlot3d("deformed_ffd.xyz")

# Assign the newly computed surface coordinates
mesh.setSurfaceCoordinates(newCoords)

# Actually run the mesh warping
mesh.warpMesh()

# Write the new grid file.
mesh.writeGrid("deformed_grid.cgns")

############## Settign up adflow

aero_options["gridFile"] = "deformed_grid.cgns"

# Creating adflow object
CFDSolver = ADFLOW(options=aero_options, comm=MPI.COMM_WORLD)

# Creating aeroproblem for adflow
ap = AeroProblem(
    name="ap", alpha=alpha, mach=0.734, reynolds=6.5e6, reynoldsLength=1.0, T=288.15, 
    areaRef=1.0, chordRef=1.0, evalFuncs=evalFuncs, xRef = 0.25, yRef = 0.0, zRef = 0.0
)

# Adding angle of attack as variable
ap.addDV("alpha", value=alpha, lower=0, upper=15.0, scale=1.0)

############## Settign up constraints

DVCon = DVConstraints()
DVCon.setDVGeo(DVGeo)
DVCon.setSurface(CFDSolver.getTriangulatedMeshSurface())

# Defining the four corners of the 2D plane for volume constraint
le = 1e-4
leList = [[le, 0, le], [le, 0, 1.0 - le]]
teList = [[1.0 - le, 0, le], [1.0 - le, 0, 1.0 - le]]
DVCon.addVolumeConstraint(leList, teList, 2, 100, lower=0.5, upper=3.0, scaled=False)

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

SLSQP_Options = {
    "ACC": 1e-3
}

PSQP_options = {
    "MET": 1,
    "XMAX": 4.0,
    "MFV": 15,
    "MIT": 15
}

# Creating optimizer for optimization
opt = SLSQP(options=SLSQP_Options)
opt = PSQP(options=PSQP_options)

# Run Optimization
sol = opt(optProb, sens=FuncsSens, storeHistory="opt.hst", sensMode="pgc", )

if MPI.COMM_WORLD.rank == 0:
    print(sol)

############# Post-processing

# Writing the surface results in a tecplot file
CFDSolver.writeSurfaceSolutionFileTecplot("sample_{}_surf.plt".format(input["sampleNo"]))

# Writing the slice file
# CFDSolver.writeSlicesFile("sample_{}_slice.dat".format(input["sampleNo"]))

# Evaluating QoI
# Note: This doesn't run the simulation.
funcs = {}
CFDSolver.evalFunctions(ap, funcs, evalFuncs)
DVCon.evalFunctions(funcs)

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
