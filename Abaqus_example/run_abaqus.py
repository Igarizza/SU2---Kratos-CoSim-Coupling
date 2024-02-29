import KratosMultiphysics as Kratos
from KratosMultiphysics.CoSimulationApplication import CoSimIO

from abaqus.ExControl.EX_AbaqusControl import EX_AbaqusControl as abq
import numpy as np
import os

def cosimio_check_equal(a, b):
    assert a == b

s_connection_name = ""

readTemperature = False
readHeatflow = False

# Clean up old files

if "meshFolder" in os.listdir("abaqus/"):
    os.system("rm -r abaqus/meshFolder/")
if "resultFolder" in os.listdir("abaqus/"):
    os.system("mv abaqus/resultFolder abaqus/oldResultFolder")

# Initialize Abaqus Model
abq.createMesh("abaqus/createCommonModel")
# abq.generateMeshOutputFiles()
abq.setMaterial("Aluminum") #Aluminum, Steel, Copper, Titanium, Gold, Silver

def AdvanceInTime(info):

    settings = CoSimIO.Info()
    settings.SetString("identifier", "AdvanceInTime")
    settings.SetString("connection_name", s_connection_name)
    settings.SetDouble("current_time", 1.0)

    CoSimIO.ExportInfo(settings)
    return CoSimIO.Info()

def InitializeSolutionStep(info):

    settings = CoSimIO.Info()
    settings.SetString("connection_name", s_connection_name)
    settings.SetString("identifier", "InitializeSolutionStep")

    CoSimIO.ExportInfo(settings)
    return CoSimIO.Info()

def Predict(info):

    settings = CoSimIO.Info()
    settings.SetString("connection_name", s_connection_name)
    settings.SetString("identifier", "Predict")

    CoSimIO.ExportInfo(settings)
    return CoSimIO.Info()

def SolveSolutionStep(info):
    global readTemperature
    global readHeatflow

    coupling_timestep_size = 1.0

    abq.writeInputFile(coupling_timestep_size)

    abq.solveAbaqus()
    
    readTemperature = False
    readHeatflow = False
    
    return CoSimIO.Info()

def FinalizeSolutionStep(info):
    print("FinalizeSolutionStep")
    # settings = CoSimIO.Info()
    # settings.SetString("connection_name", s_connection_name)
    # settings.SetString("identifier", "info_for_test")
    # settings.SetString("name_for_check", "FinalizeSolutionStep")
    # CoSimIO.ExportInfo(settings)
    return CoSimIO.Info()

def OutputSolutionStep(info):
    print("OutputSolutionStep")
    # settings = CoSimIO.Info()
    # settings.SetString("connection_name", s_connection_name)
    # settings.SetString("identifier", "info_for_test")
    # settings.SetString("name_for_check", "OutputSolutionStep")
    # CoSimIO.ExportInfo(settings)
    return CoSimIO.Info()

def ImportData(info):

    identifier = info.GetString("identifier")

    settings = CoSimIO.Info()
    settings.SetString("connection_name", s_connection_name)
    settings.SetString("identifier", identifier)

    imported_data = CoSimIO.DoubleVector()
    CoSimIO.ImportData(settings, imported_data)

    surfaceID = identifier.split("_")[-1]
    
    BCType = abq.setBCType(surfaceID)
    abq.setBC(BCType, surfaceID, imported_data)

    return CoSimIO.Info()

def ExportData(info):
    global readTemperature
    global readHeatflow

    identifier = info.GetString("identifier")

    settings = CoSimIO.Info()
    settings.SetString("identifier", identifier)
    settings.SetString("connection_name", s_connection_name)

    if not readTemperature:
        abq.readOdbFile("temperature", "abaqus/inputFile.odb", "StructureStep")
        readTemperature = True

    if not readHeatflow:
        abq.readOdbFile("heatflow", "abaqus/inputFile.odb", "StructureStep")
        readHeatflow = True
    
    surfaceData = abq.getSurfaceData(identifier.split("_")[0], identifier.split("_")[-1], "StructureStep")
    data_to_be_send = CoSimIO.DoubleVector(surfaceData)
    CoSimIO.ExportData(settings, data_to_be_send)

    return CoSimIO.Info()

def ImportMesh(info):
    raise RuntimeError("ImportMesh")
    settings = CoSimIO.Info()
    settings.SetString("connection_name", s_connection_name)
    settings.SetString("identifier", "info_for_test")
    settings.SetString("name_for_check", "ImportMesh")
    if (info.Has("identifier")):
        settings.SetString("identifier_control", info.GetString("identifier"))
    CoSimIO.ExportInfo(settings)
    return CoSimIO.Info()

def ExportMesh(info):

    identifier = info.GetString("identifier")
    
    info = CoSimIO.Info()
    info.SetString("identifier", identifier)
    info.SetString("connection_name", s_connection_name)
    
    model = Kratos.Model()
    model_part = model.CreateModelPart(identifier)
    
    with open("abaqus/meshFolder/nodes_" + identifier + "_with_Coords.csv", 'r') as f:
        nodes = f.read().splitlines()
    for i in range(len(nodes)):
        currentNode = nodes[i].split(",")
        model_part.CreateNewNode(int(currentNode[0]),float(currentNode[1]),float(currentNode[2]),float(currentNode[3]))

    CoSimIO.ExportMesh(info, model_part)

    return CoSimIO.Info()

# Connection Settings
settings = CoSimIO.Info()
settings.SetString("my_name", "run_structure")
settings.SetString("connect_to", "structure")
settings.SetInt("echo_level", 1)
settings.SetString("communication_format", "file")
settings.SetString("version", "1.25")

# Connecting
return_info = CoSimIO.Connect(settings)
cosimio_check_equal(return_info.GetInt("connection_status"), CoSimIO.ConnectionStatus.Connected)
s_connection_name = return_info.GetString("connection_name")

# registering the functions
fct_info = CoSimIO.Info()
fct_info.SetString("connection_name", s_connection_name)

fct_info.SetString("function_name", "AdvanceInTime")
CoSimIO.Register(fct_info,           AdvanceInTime)

fct_info.SetString("function_name", "InitializeSolutionStep")
CoSimIO.Register(fct_info,           InitializeSolutionStep)

fct_info.SetString("function_name", "Predict")
CoSimIO.Register(fct_info,           Predict)

fct_info.SetString("function_name", "SolveSolutionStep")
CoSimIO.Register(fct_info,           SolveSolutionStep)

fct_info.SetString("function_name", "FinalizeSolutionStep")
CoSimIO.Register(fct_info,           FinalizeSolutionStep)

fct_info.SetString("function_name", "OutputSolutionStep")
CoSimIO.Register(fct_info,           OutputSolutionStep)

fct_info.SetString("function_name", "ImportData")
CoSimIO.Register(fct_info,           ImportData)

fct_info.SetString("function_name", "ExportData")
CoSimIO.Register(fct_info,           ExportData)

fct_info.SetString("function_name", "ImportMesh")
CoSimIO.Register(fct_info,           ImportMesh)

fct_info.SetString("function_name", "ExportMesh")
CoSimIO.Register(fct_info,           ExportMesh)

# running the simulation
# externally orchestrated
run_info = CoSimIO.Info()
run_info.SetString("connection_name", s_connection_name)
CoSimIO.Run(run_info)

# Disconnecting
disconnect_settings = CoSimIO.Info()
disconnect_settings.SetString("connection_name", s_connection_name)
return_info = CoSimIO.Disconnect(disconnect_settings)
cosimio_check_equal(return_info.GetInt("connection_status"), CoSimIO.ConnectionStatus.Disconnected)
