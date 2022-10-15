from scipy.io import loadmat

data = loadmat("training_single/data.mat")

print(data["x"])
print(data["y"])