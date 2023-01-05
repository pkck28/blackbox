from blackbox import Rosenbrock
import numpy as np

options = {
    "numberOfSamples" : 5,
    "samplingMethod" : "lhs"
}

# Example for generating samples
object_1 = Rosenbrock(options=options)
object_1.generateSamples()

# Example for getting value of one sample
object_2 = Rosenbrock(type="single")
print(object_2.getObjectives(np.array([0.5, 0.5])))
