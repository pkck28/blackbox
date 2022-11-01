from pygeo import DVGeometry
from idwarp import USMesh
import pickle
import os
import numpy as np

# Importing warnings package
# and supressing the warnings
import warnings
warnings.filterwarnings("ignore")

# Reading input file for the analysis
filehandler = open("input.pickle", 'rb') 
input = pickle.load(filehandler)
filehandler.close()

sample = input["sample"]

geoDV = {}

for key in sample:
    if key == "twist" or key == "shape":
        geoDV[key] = sample[key]

print("Design Variables:")
print(geoDV)

# All are default expect useRotation and zeroCornerRotation
options = {
  'gridFile':'grid.cgns',
  'fileType':'CGNS',
  'specifiedSurfaces':None,
  'symmetrySurfaces':None,
  'symmetryPlanes':[],
  'aExp': 3.0,
  'bExp': 5.0,
  'LdefFact': 1.0,
  'alpha': 0.25,
  'errTol': 0.0005,
  'evalMode': 'fast',
  'useRotations': False,
  'zeroCornerRotations': True,
  'cornerAngle': 30.0,
  'bucketSize': 8,
}

# Create the mesh object
mesh = USMesh(options=options)

# Extract all undeformed surface mesh coordinates
coords0 = mesh.getSurfaceCoordinates()

# Initializing the DVGeometry class with the FFD file
DVGeo = DVGeometry("ffd.xyz")

# Adding surface mesh co-ordinates as a pointset
DVGeo.addPointSet(coords0, "wing_surface_mesh")

# Add twist as a DV based on user selection
if "twist" in geoDV.keys():
    # Create reference axis
    nRefAxPts = DVGeo.addRefAxis("wing", xFraction=0.25, alignIndex="j")

    print(nRefAxPts)

    # Set the Twist Variable
    def twist(val, geo):
        for i in range(1, nRefAxPts):
            geo.rot_z["wing"].coef[i] = val[i]

    DVGeo.addGlobalDV(dvName="twist", value=[0] * nRefAxPts, func=twist, lower=-10, upper=10, scale=1.0)

    geoDV["twist"] = np.append(np.array([0]), geoDV["twist"])

# Add shape as a DV based on user selection
if "shape" in geoDV.keys():
    if input["aeroSolverOptions"]["liftindex"] == 3:
        DVGeo.addLocalDV("shape", lower=-0.25, upper=0.25, axis="z", scale=1)
    elif input["aeroSolverOptions"]["liftindex"] == 2:
        DVGeo.addLocalDV("shape", lower=-0.25, upper=0.25, axis="y", scale=1)

# Get Design Variables
dvDict = DVGeo.getValues()

# Set the value of DV from sample
dvDict = geoDV

# Set Design Variables
DVGeo.setDesignVars(dvDict)

# Getting the updating suface mesh co-ordinates
newCoords = DVGeo.update("wing_surface_mesh")
DVGeo.writePlot3d("deformed_ffd.xyz")

# Reset the newly computed surface coordinates
mesh.setSurfaceCoordinates(newCoords)

# Actually run the mesh warping
mesh.warpMesh()

# Write the new grid file.
mesh.writeGrid('deformed_grid.cgns')

# print(input["aeroSolverOptions"]["liftindex"])
# if input["aeroSolverOptions"]["liftindex"] == 3:
#     os.system('cgns_utils symmZero deformed_grid.cgns y')
# elif input["aeroSolverOptions"]["liftindex"] == 2:
#     os.system('cgns_utils symmZero deformed_grid.cgns z')

# Clean up
os.system('rm grid.cgns')
os.system('mv deformed_grid.cgns grid.cgns')
os.system('rm ffd.xyz')
