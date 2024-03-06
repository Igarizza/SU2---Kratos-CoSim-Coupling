import KratosMultiphysics
from KratosMultiphysics.CoSimulationApplication.co_simulation_analysis import CoSimulationAnalysis

import sys
import time

if __name__ == "__main__":
    import os
    # os.chdir("/localdata/Su2_FSI/OneraM6_With_Solid")
    with open("ProjectParametersCoSimFSI.json",'r') as parameter_file:
        parameters = KratosMultiphysics.Parameters(parameter_file.read())

    simulation = CoSimulationAnalysis(parameters)
    simulation.Run()
