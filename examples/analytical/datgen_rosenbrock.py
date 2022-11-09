from blackbox import Rosenbrock

options = {
    "parameters" : {
        "a" : 1,
        "b" : 100
    }
}

test = Rosenbrock(type="single")

test.getObjectives([0, 1])
