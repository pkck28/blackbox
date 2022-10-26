# Importing python packages
import os
import shutil
import sys
import numpy as np
from pyDOE2 import lhs, fullfact
import pickle
from scipy.io import savemat
import math

class DefaultOptions():
    """
        Class creates a default option list (for solver and other settings)
        which is later edited/appended with user provided options.
    """

    def __init__(self):

        # Aero solver Options
        self.aeroSolver = "adflow"

        # Other options
        self.directory = "output"
        self.noOfProcessors = 4
        self.ffdFile = None

class Airfoil():
    """
        Class contains essential methods for generating airfoil data.
        There are two values possible for type: "single" and "multi" (default). 
        For "multi", following is the list of possible attributes:
        
        Optional: (.,.) shows datatype and defualt value respectively.
            "directory" : Folder name where the data.mat file will be saved (string, "output").
            "noOfProcessors" : Number of processors to use (integer, 4).
            "aeroSolver" : Name of the aerodynamics solver (string, "adflow").
        Compulsory: (.) shows datatype.
            "numberOfSamples" : number of samples to be generated (integer).
            "fixedParameters" : Dictionary of all the valid fixed parameters (dict).
            "varyingParameters" : Dictionary of all the valid varying parameters (dict).
            "samplingMethod" : name of the sampling method ("lhs" or "fullfactorial") (string).
            "objectives" : List of desired objectives in y (list of string).
            "aeroSolverOptions" : Dictionary of containing various ADflow options (dict).

        For "single", following is the list of possible attributes:
        
        Optional: (.,.) shows datatype and defualt value respectively.
            "directory" : Folder name where the analysis output will be saved (string, "additional_samples").
            "noOfProcessors" : Number of processors to use (integer, 4).
            "aeroSolver" : Name of the aerodynamics solver (string, "adflow").
        Compulsory: (.) shows datatype.
            "fixedParameters" : Dictionary of all the valid fixed parameters (dict).
            "varyingParameters" : Dictionary of all the valid varying parameters (dict).
            "objectives" : List of desired objectives in y (list of string).
            "aeroSolverOptions" : Dictionary of containing various ADflow options (dict).
    """