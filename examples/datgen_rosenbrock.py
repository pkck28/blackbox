import argparse
from datgen import Rosenbrock

parser = argparse.ArgumentParser(
            description='Python script for building a surrogate model, \
                        test and optimize it in a batch.'
        )
parser.add_argument('-sample_points', default=[50], type=int, nargs="*")

args = parser.parse_args()

options = {
        "samplingMethod" : "lhs"
    }

for pts in args.sample_points:
    
    options["numberOfSamples"] = pts
    options["directory"] = "output_sample_{}".format(pts)
    
    test = Rosenbrock("multi", options)

    test.generateSamples()
