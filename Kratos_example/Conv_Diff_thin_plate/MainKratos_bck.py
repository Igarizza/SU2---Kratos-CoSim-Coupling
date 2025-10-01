import KratosMultiphysics as Kratos
from KratosMultiphysics.ConvectionDiffusionApplication.convection_diffusion_analysis import ConvectionDiffusionAnalysis
import KratosMultiphysics.MedApplication



"""
For user-scripting it is intended that a new class is derived
from CoSimulationAnalysis to do modifications
Check also "kratos/python_scripts/analysis-stage.py" for available methods that can be overridden
"""

with open("ProjectParameters.json", 'r') as parameter_file:
        parameters = KratosMultiphysics.Parameters(parameter_file.read())

model = Kratos.Model()
model_part = model.CreateModelPart("plate")

KratosMultiphysics.MedApplication.MedModelPartIO("plate.med", KratosMultiphysics.ModelPartIO.READ).ReadModelPart(model_part)
model_part.ProcessInfo[KratosMultiphysics.DOMAIN_SIZE] = 3

volume = model_part.CreateSubModelPart("volume")
upper_surface = model_part.CreateSubModelPart("upper_surface")
lower_surface = model_part.CreateSubModelPart("lower_surface")
# sides = model_part.CreateSubModelPart("sides")

# adding nodes

for node in model_part.Nodes:
    # all nodes to volume
    # volume.AddNode(node)

    # nodes to upper surface z == 0
    if node.Z0 == 0.0:
         upper_surface.AddNode(node)

    # nodes to lower z == -2
    if node.Z0 == -2.0:
         lower_surface.AddNode(node)

    # nodes to sides, y=0, y = 1000, x = 2000, x = 3000
    # if node.X0 == 2000.0 or node.X0 == 3000 or node.Y0 == 0.0 :# or node.Y0 == 1000.0:
    #      print(node)
    #      sides.AddNode(node)

params = Kratos.Parameters("""{
    "elements_list" : [
        {
            "model_part_name" : "plate.volume",
            "element_name" : "SmallDisplacementElement3D8N"
        }
    ]
}""")

modeler = Kratos.CreateEntitiesFromGeometriesModeler(model, params)
modeler.SetupModelPart()

print("---- added nodes ----")
print(model_part)


properties = model_part.GetProperties()
for geometry in model_part.Geometries:
    print(geometry)
    model_part.CreateNewElement( properties[0])
    surface = True
    nodes = []
    print(f"number of vertices: {geometry.PointsNumber()}")
    for i in range(geometry.PointsNumber()):
         print(geometry[i])
    upper_surface.AddGeometry(geometry)
    print(upper_surface)
    #print(geometry.GetPoints())
    raise RuntimeError(1)
    for node in element.GetNodes():
        print(node)
        if not node.Z0 == 0.0:
             surface = False
        else:
             nodes.append(node)
    if surface:
        upper_surface.AddElement(element)
        upper_surface.AddNodes(nodes)

print(model_part)

# simulation = ConvectionDiffusionAnalysis(model, parameters)
# simulation.Run()
