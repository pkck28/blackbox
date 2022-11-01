from scipy.io import loadmat
import matplotlib.pyplot as plt

data = loadmat("analytical/output/data.mat")

x = data["x"]
y = data["y"]

print(data["x"])
print(data["y"])

fig, ax = plt.subplots()

ax.scatter(x[:,0], x[:,1])

plt.grid()
plt.show()
