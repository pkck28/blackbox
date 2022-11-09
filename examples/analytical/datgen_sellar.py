from datgen import Sellar
import numpy as np

options = {
    "numberOfSamples" : 4,
    "directory" : "forrester",
    "samplingMethod" : "lhs",
    "directory" : "sellar"
}

# # Example for generating samples
# object_1 = Sellar(options=options)
# object_1.generateSamples()

# Example for generating samples
object_1 = Sellar(type="single")
x = np.array([-1.0, 7.5, 0.5])
result, d1, d2 = object_1.getObjectives(x)
print(result)
