from pygeo import DVGeometry
from idwarp import USMesh
import pickle
import os

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

geoDV["shape"] = sample["shape"]

print("Design Variables:")
print(geoDV)

# All are default expect useRotation and zeroCornerRotation
options = {
  'gridFile':'grid.cgns',
}

# Create the mesh object
mesh = USMesh(options=options)

# Extract all undeformed surface mesh coordinates
coords0 = mesh.getSurfaceCoordinates()

# Initializing the DVGeometry class with the FFD file
DVGeo = DVGeometry("ffd.xyz")

# Adding surface mesh co-ordinates as a pointset
DVGeo.addPointSet(coords0, "wing_surface_mesh")

# Add shape as a DV based on user selection
if "shape" in geoDV.keys():
    if input["aeroSolverOptions"]["liftindex"] == 3:
        DVGeo.addLocalDV("shape", lower=-0.25, upper=0.25, axis="z", scale=1)
    elif input["aeroSolverOptions"]["liftindex"] == 2:
        DVGeo.addLocalDV("shape", lower=-0.25, upper=0.25, axis="y", scale=1)

# Get Design Variables
dvDict = DVGeo.getValues()

# Set the value of DV from sample
index = DVGeo.getLocalIndex(0)
pts = DVGeo.getLocalIndex(0).shape[0]

# Setting DV values such that same deformation on both sides
for i in range(pts):
    dvDict["shape"][index[i, 0, 0]] = sample["shape"][i]
    dvDict["shape"][index[i, 0, 1]] = sample["shape"][i]

    dvDict["shape"][index[i, 1, 0]] = sample["shape"][i+pts]
    dvDict["shape"][index[i, 1, 1]] = sample["shape"][i+pts]

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

# Clean up
os.system('rm grid.cgns')
os.system('mv deformed_grid.cgns grid.cgns')
os.system('rm ffd.xyz')
