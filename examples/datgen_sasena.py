from datgen import Sasena
from scipy.io import loadmat

options = {
    "numberOfSamples" : 34,
    "samplingMethod" : "fullfactorial"
}

test = Sasena(type="single")

print(test.getObjectives([0, 0]))


# test.generateSamples()

# data = loadmat("output/data.mat")

# x = data["x"]
# y = data["y"]

# print(x)
# print(y)

# print(x.shape)
# print(y.shape)

# print(type(x))
# print(type(y))
