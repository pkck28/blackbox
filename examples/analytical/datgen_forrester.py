from blackbox import Forrester

options = {
    "numberOfSamples" : 25,
}

# Example for generating samples
object_1 = Forrester(options=options)
object_1.generateSamples()

# Example for getting value of one sample
object_2 = Forrester(type="single")
print(object_2.getObjectives(0.75))
