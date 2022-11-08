from datgen import Forrester

options = {
    "lowerBound" : 0,
    "upperBound" : 1,
    "numberOfSamples" : 50,
    "directory" : "forrester"
}

# Example for generating samples
object_1 = Forrester(options=options)
object_1.generateSamples()

# Example for getting value of one sample
# object_2 = Forrester(type="single")
# x = 0.5
# print(object_2.getObjectives(x))
