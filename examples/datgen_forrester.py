import argparse
from datgen import Forrester

# parser = argparse.ArgumentParser(
#             description='Python script for building a surrogate model, \
#                         test and optimize it in a batch.'
#         )
# parser.add_argument('-sample_points', default=[50], type=int, nargs="*")

# args = parser.parse_args()

# options = {
#         # "samplingMethod" : "fullfactorial"
#     }

# for pts in args.sample_points:
    
#     options["numberOfSamples"] = pts
#     # options["directory"] = "output_sample_{}".format(pts)
    
#     test = Forrester("multi", options)

#     test.generateSamples()

options = {
    "lowerBound" : 2,
    "upperBound" : 1,
    "numberOfSamples" : 5
}

test = Forrester(options=options)

# x = 0.5

# print(test.getObjectives(x))