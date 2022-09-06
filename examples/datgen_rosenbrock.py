from datgen import Rosenbrock

options = {
    "numberOfSamples" : 50,
    "lowerBound": [5, 10],
    "upperBound": [8, 8]
}

test = Rosenbrock(options)

