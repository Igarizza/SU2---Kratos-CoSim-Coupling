import KratosMultiphysics
from KratosMultiphysics.CoSimulationApplication.co_simulation_analysis import CoSimulationAnalysis

if __name__ == "__main__":

    with open("ProjectParametersCoSimFSI.json",'r') as parameter_file:
        parameters = KratosMultiphysics.Parameters(parameter_file.read())

    simulation = CoSimulationAnalysis(parameters)
    simulation.Run()
