from calendar import c
from scipy.io import loadmat
import matplotlib.pyplot as plt
import numpy as np

data = loadmat("output/data.mat")

x = data["x"]
y = data["y"][:,0]
g = data["y"][:,1]

fig, ax = plt.subplots()

ax.contourf(x[:,0].reshape((10,10)), x[:,1].reshape((10,10)), y.reshape((10,10)), levels = 40)
ax.contour(x[:,0].reshape((10,10)), x[:,1].reshape((10,10)), g.reshape((10,10)), levels = 0, colors="red")
ax.scatter(-np.pi,12.275)
ax.scatter(np.pi,2.275)
ax.scatter(9.42478, 2.475)

# plt.grid()
plt.show()
