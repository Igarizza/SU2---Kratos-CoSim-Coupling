from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7

import KratosMultiphysics as Kratos
from KratosMultiphysics.OptimizationApplication.optimization_analysis import OptimizationAnalysis

import sys
import time

if __name__ == "__main__":
    import os
    # os.chdir("/localdata/Su2_FSI/OneraM6_With_Solid")
    with open("optimization_parameters.json",'r') as parameter_file:
        parameters = Kratos.Parameters(parameter_file.read())

    model = Kratos.Model()
    simulation = OptimizationAnalysis(model, parameters)
    simulation.Run()
