from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7

import KratosMultiphysics
from KratosMultiphysics.CoSimulationApplication.co_simulation_analysis import CoSimulationAnalysis
import KratosMultiphysics.StatisticsApplication

import sys
import time

if __name__ == "__main__":

    with open("ProjectParametersCoSimFSI.json",'r') as parameter_file:
        parameters = KratosMultiphysics.Parameters(parameter_file.read())

    simulation = CoSimulationAnalysis(parameters)
    simulation.Run()
