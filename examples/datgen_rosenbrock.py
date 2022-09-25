from datgen import Rosenbrock

options = {
    "numberOfSamples" : 5,
    "samplingMethod" : "lhs",
    "lowerBound" : [-2, -2],
    "upperBound" : [2, 2],
    "directory" : "rosenbrock",
    "parameters" : {
        "a" : 1,
        "b" : 100,
        "c" : 1
    }
}

test = Rosenbrock(options=options)

test.generateSamples()
