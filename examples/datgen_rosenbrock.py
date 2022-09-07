from datgen import Rosenbrock

options = {
    "numberOfSamples" : 100,
    "directory" : "output",
    "samplingMethod" : "lhs"
}

test = Rosenbrock("batch", options)

test.generateSamples()
