from scipy.io import loadmat

data = loadmat("multi/data.mat")

print(data["x"])

print(data["cl"].shape)

print(data["fail"])
