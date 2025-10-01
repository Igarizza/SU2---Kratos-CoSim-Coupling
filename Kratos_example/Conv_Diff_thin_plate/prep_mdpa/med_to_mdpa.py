import sys

# sys.path.append("/home/wessels/Software/git/Kratos/bin/Release")
# sys.path.append("/home/wessels/Software/git/Kratos/build/Release")

import KratosMultiphysics as Kratos
import KratosMultiphysics.MedApplication as KratosMED
from KratosMultiphysics.StructuralMechanicsApplication.structural_mechanics_analysis import StructuralMechanicsAnalysis

# Setup basic model to be filled later on
model = Kratos.Model()
plate_4 = model.CreateModelPart("plate_4")

KratosMED.MedModelPartIO("plate_4_m.med", Kratos.ModelPartIO.READ | Kratos.ModelPartIO.MESH_ONLY).ReadModelPart(plate_4)
print(plate_4)

# Use to MED App to read med-file and write it to the previously created model part


# From the med-file only geometries are read. So we need to convert them.
# To tell Kratos how to convert the different geometries we are defining a conversion table in the following:
# Note that you can only have one element type per model part. If you want to have different elements you need to have them grouped in different model parts. 
# The different element_ and condition_names are found in the wide space of Kratos files.
# Use Structure.NAME_OF_GROUP_IN_SALOME.
entity_params = Kratos.Parameters("""{
    "elements_list": [
        {
            "model_part_name": "plate_4",
            "element_name": "Element3D8N"
        }
    ],
    "conditions_list": [
        {
            "model_part_name": "plate_4.upper_surface",
            "condition_name": "SurfaceCondition3D4N"
        },
        {
            "model_part_name": "plate_4.lower_surface",
            "condition_name": "SurfaceCondition3D4N"
        }
    ]    
}""")

# Now, we use the previously defined table to convert our geometries
modeler = Kratos.CreateEntitiesFromGeometriesModeler(model, entity_params)
modeler.SetupModelPart()

# All that's left is to assign a fresh properties container to the model
properties = plate_4.CreateNewProperties(1)
cond: Kratos.Condition
for cond in plate_4.Conditions:
    cond.Properties = properties

for elem in plate_4.Elements:
    elem.Properties = properties

geometry_ids = []
for geom in plate_4.Geometries:
    geometry_ids.append(geom.Id)

for geom_id in geometry_ids:
    plate_4.RemoveGeometryFromAllLevels(geom_id)

# In the end we write the model to a mdpa-file
Kratos.ModelPartIO("plate_m", Kratos.ModelPartIO.MESH_ONLY | Kratos.ModelPartIO.WRITE).WriteModelPart(plate_4)

# ... and clean up unused classes
del model

 
