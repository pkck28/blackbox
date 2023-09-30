from blackbox import Sasena
import numpy as np

options = {
    "numberOfSamples" : 100,
    "samplingMethod" : "lhs"
}

object_1 = Sasena(options=options)
object_1.generateSamples()

object_2 = Sasena(type="single")
print(object_2.getObjectives(np.array([0.5, 0.5])))
