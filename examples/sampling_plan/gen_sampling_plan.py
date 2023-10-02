from blackbox import AirfoilFFD
import numpy as np
from matplotlib import pyplot as plt

nffd = 20 # Number of FFD points

# Options for blackbox
options = {
    # Required options
    "airfoilFile": "rae2822.dat",
    "nffd": nffd,
    # FFD Box options
    "fitted": True,
    "ymarginl": 0.015,
    "ymarginu": 0.015,
    # Sampling options
    "sampling": "internal",
    "samplingCriterion": "ese",
    # Fixing the LE/TE
    "fixLETE": True,
    # Smoothing options
    "smoothingTolerance": 5e-4,
    "smoothingTheta": 0.6
}

# Creating blackbox object
airfoil = AirfoilFFD(options=options)

############ Setting bounds for shape design variables ############

coeff = airfoil.DVGeo.origFFDCoef
index = airfoil.DVGeo.getLocalIndex(0)

# Initialize bounds
lowerBound = np.zeros(int(coeff.shape[0]/2))
upperBound = np.zeros(int(coeff.shape[0]/2))

# Lift direction
liftIndex = 1 # y

# Computing the y distance between upper and lower surface FFD points
# 0 at thrid index denotes section location - there are two sections
# (0,1) at second index denotes lower and upper surface of FFD points
dist = coeff[index[:,1,0], liftIndex] - coeff[index[:,0,0], liftIndex]

# Lower and upper surface FFD point index in DV
lowerIndex = np.linspace(0, nffd-2, dist.shape[0], dtype=int) # lower surface ffd point index in DV
upperIndex = np.linspace(1, nffd-1, dist.shape[0], dtype=int) # upper surface ffd point index in DV

# FFD local thickness
delta = 0.2

# Lower bounds
lowerBound[lowerIndex] = -delta * dist
lowerBound[upperIndex] = -delta * dist

# Upper bounds
upperBound[lowerIndex] = delta * dist
upperBound[upperIndex] = delta * dist

# Add shape as a design variables
# First and last entry is dropped for fixing LE/TE
airfoil.addDV("shape", lowerBound=lowerBound[1:-1], upperBound=upperBound[1:-1])

# Number of samples
num = 100

# Generate sampling plan
x = airfoil.sampler(num)

# Plotting the airfoils
for i in range(num):

    # Get airfoil coordinates
    non_smooth = airfoil.getAirfoil(x[i,:])
    y = airfoil.LaplacianSmoothing(x[i,:])
    smooth = airfoil.getAirfoil(y)

    # Plotting
    fig, ax = plt.subplots()
    ax.plot(airfoil.coords[:,0], airfoil.coords[:,1], label="original")
    ax.plot(non_smooth[:,0], non_smooth[:,1], label="non_smooth")
    ax.plot(smooth[:,0], smooth[:,1], label="smooth")
    ax.legend(fontsize=12, loc="lower right")
    ax.set_xlabel("$x/c$", fontsize=14)
    ax.set_ylabel("$y/c$", fontsize=14)
    ax.set_title("Sample {}".format(i+1), fontsize=14)
    plt.show()
