from blackbox import AirfoilFFD
from baseclasses import AeroProblem
import numpy as np
from matplotlib import pyplot as plt

# Creating aeroproblem for adflow
ap = AeroProblem(
    name="ap", alpha=2.0, mach=0.734, reynolds=6.5e6, reynoldsLength=1.0, T=288.15, 
    areaRef=1.0, chordRef=1.0, evalFuncs=["cl", "cd", "cmz"], xRef = 0.25, yRef = 0.0, zRef = 0.0
)

nffd = 20 # Number of FFD points

# Options for blackbox
options = {
    "solverOptions": {},
    "aeroProblem": ap,
    "airfoilFile": "rae2822.dat",
    "nffd": nffd,
    "meshingOptions": {},
    "fitted": True,
    "sampling": "internal",
    "samplingCriterion": "ese",
    "fixLETE": True
}

# Creating blackbox object
airfoil = AirfoilFFD(options=options)

# Lower and upper bounds for shape variables
lower = np.array([-0.01]*(nffd-2))
upper = np.array([0.01]*(nffd-2))

# Add shape as a design variables
airfoil.addDV("shape", lowerBound=lower, upperBound=upper)

# Number of samples
num = 100

# Generate sampling plan
x = airfoil.sampler(num)

# Plotting the airfoils
for i in range(num):
    new_coords = airfoil.getAirfoil(x[i,:])
    fig, ax = plt.subplots()
    ax.plot(airfoil.coords[:,0], airfoil.coords[:,1], label="original")
    ax.plot(new_coords[:,0], new_coords[:,1], label="deformed")
    ax.legend()
    ax.set_title("Sample {}".format(i))
    plt.show()
