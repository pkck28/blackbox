from blackbox import Sellar
import numpy as np

options = {
    "numberOfSamples" : 4,
    "samplingMethod" : "lhs"
}

# Example for generating samples
object_1 = Sellar(options=options)
object_1.generateSamples()

# Example for getting value of one sample
object_2 = Sellar(type="single")
x = np.array([3.03, 0, 0])
result, d1, d2 = object_2.getObjectives(x)
print(result)
