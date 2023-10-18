import KratosMultiphysics as Kratos
from KratosMultiphysics.OptimizationApplication.optimization_analysis import OptimizationAnalysis


if __name__ == "__main__":
    with kratos_unittest.WorkFolderScope(".", __file__):
        with open("optimization_parameters.json", "r") as file_input:
            parameters = Kratos.Parameters(file_input.read())

        model = Kratos.Model()
        analysis = OptimizationAnalysis(model, parameters)
        analysis.Run()

