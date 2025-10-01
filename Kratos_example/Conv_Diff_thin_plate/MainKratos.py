import KratosMultiphysics
from KratosMultiphysics.ConvectionDiffusionApplication.convection_diffusion_analysis import ConvectionDiffusionAnalysis
import KratosMultiphysics.MedApplication



"""
For user-scripting it is intended that a new class is derived
from CoSimulationAnalysis to do modifications
Check also "kratos/python_scripts/analysis-stage.py" for available methods that can be overridden
"""

with open("ProjectParameters.json", 'r') as parameter_file:
        parameters = KratosMultiphysics.Parameters(parameter_file.read())

model = KratosMultiphysics.Model()
model_part = model.CreateModelPart("plate_4")

simulation = ConvectionDiffusionAnalysis(model, parameters)

KratosMultiphysics.MedApplication.MedModelPartIO("plate_4.med", KratosMultiphysics.ModelPartIO.READ).ReadModelPart(model_part)
model_part.ProcessInfo[KratosMultiphysics.DOMAIN_SIZE] = 3

print(model_part)

simulation.Run()

KratosMultiphysics.ModelPartIO("plate", KratosMultiphysics.ModelPartIO.WRITE | KratosMultiphysics.ModelPartIO.MESH_ONLY). WriteModelPart(model_part)
