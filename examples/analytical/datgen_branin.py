from blackbox import Branin
import numpy as np

options = {
    "numberOfSamples" : 400,
    "samplingMethod" : "fullfactorial"
}

# Example for generating samples
object_1 = Branin(options=options)
object_1.generateSamples()

object_2 = Branin(type="single")

x = np.array([-np.pi, 12.275])
print(object_2.getObjectives(x))

x = np.array([np.pi, 2.275])
print(object_2.getObjectives(x))

x = np.array([9.42478, 2.475])
print(object_2.getObjectives(x))
