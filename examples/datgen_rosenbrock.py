from datgen import Rosenbrock

options = {
    "numberOfSamples" : 100,
    "directory" : "training",
    "samplingMethod" : "lhs"
}

test = Rosenbrock(options)

test.generateSamples()
