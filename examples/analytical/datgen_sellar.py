from datgen import Sellar
import numpy as np

options = {
    "numberOfSamples" : 100,
    "directory" : "forrester",
    "samplingMethod" : "lhs",
    "directory" : "sellar"
}

# # Example for generating samples
# object_1 = Sellar(options=options)
# object_1.generateSamples()

# Example for generating samples
object_1 = Sellar(type="single")
x = np.array([1.977, 0.0, 0.0])
print(object_1.getObjectives(x))
