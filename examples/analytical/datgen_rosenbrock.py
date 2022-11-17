from blackbox import Rosenbrock
import numpy as np

options = {
    "lowerBound" : [-2, -2],
    "upperBound" : [2, 2],
    "numberOfSamples" : 50,
    "samplingMethod" : "lhs"
}

# Example for generating samples
object_1 = Rosenbrock(options=options)
object_1.generateSamples()

# Example for getting value of one sample
object_2 = Rosenbrock(type="single")
print(object_2.getObjectives(np.array([1, 1])))
