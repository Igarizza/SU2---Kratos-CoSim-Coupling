import KratosMultiphysics as Kratos
from KratosMultiphysics.CoSimulationApplication import CoSimIO
import numpy as np
from KratosMultiphysics.OptimizationApplication.utilities.logger_utilities import time_decorator


"""
This is just a dummy solver to test coupling.
"""

def cosimio_check_equal(a, b):
    assert a == b

s_connection_name = ""
# constant variables

class Variables:
    temperature_External = [230,240,250,260]
    temperature_Channel0 = [356,358,360,362,364,366,368,370,372,374,376,378,380,382,384,386]
    temperature_Channel1 = [356,358,360,362,364,366,368,370,372,374,376,378,380,382,384,386]
    temperature_Channel2 = [356,358,360,362,364,366,368,370,372,374,376,378,380,382,384,386]
    heatflow_External = []
    heatflow_Channel0 = []
    heatflow_Channel1 = []
    heatflow_Channel2 = []
    """
    heatflow_External = np.empty(len(temperature_External),np.float64)
    heatflow_Channel0 = np.empty(len(temperature_Channel0),np.float64)
    heatflow_Channel1 = np.empty(len(temperature_Channel1),np.float64)
    heatflow_Channel2 = np.empty(len(temperature_Channel2),np.float64)
    """

var = Variables()

@time_decorator()
def AdvanceInTime(info):
    print("AdvanceInTimeInfo: ")
    print(info)

    settings = CoSimIO.Info()
    settings.SetString("identifier", "AdvanceInTime")
    settings.SetString("connection_name", s_connection_name)
    settings.SetDouble("current_time", 1.0)

    CoSimIO.ExportInfo(settings)
    return CoSimIO.Info()

@time_decorator()
def InitializeSolutionStep(info):
    print("InitializeSolutionStepInfo: ")
    print(info)

    settings = CoSimIO.Info()
    settings.SetString("identifier", "InitializeSolutionStep")
    settings.SetString("connection_name", s_connection_name)

    CoSimIO.ExportInfo(settings)
    return CoSimIO.Info()

@time_decorator()
def Predict(info):
    print("PredictInfo: ")
    print(info)

    settings = CoSimIO.Info()
    settings.SetString("identifier", "Predict")
    settings.SetString("connection_name", s_connection_name)

    CoSimIO.ExportInfo(settings)
    return CoSimIO.Info()

@time_decorator()
def SolveSolutionStep(info):
    print("SolveSolutionStepInfo: ")
    print(info)

    for i in range(4):
        var.heatflow_External.append(var.temperature_External[i] * 3000)

    for i in range(16):
        var.heatflow_Channel0.append(var.temperature_Channel0[i] * 3000)
        var.heatflow_Channel1.append(var.temperature_Channel1[i] * 3000)
        var.heatflow_Channel2.append(var.temperature_Channel2[i] * 3000)
    
    settings = CoSimIO.Info()
    settings.SetString("identifier", "SolveSolutionStep")
    settings.SetString("connection_name", s_connection_name)

    return CoSimIO.Info()

@time_decorator()
def FinalizeSolutionStep(info):
    print("FinalizeSolutionStepInfo: ")
    print(info)
    # settings = CoSimIO.Info()
    # settings.SetString("connection_name", s_connection_name)
    # settings.SetString("identifier", "info_for_test")
    # settings.SetString("name_for_check", "FinalizeSolutionStep")
    # CoSimIO.ExportInfo(settings)
    return CoSimIO.Info()

@time_decorator()
def OutputSolutionStep(info):
    print("OutputSolutionStepInfo: ")
    print(info)
    # settings = CoSimIO.Info()
    # settings.SetString("connection_name", s_connection_name)
    # settings.SetString("identifier", "info_for_test")
    # settings.SetString("name_for_check", "OutputSolutionStep")
    # CoSimIO.ExportInfo(settings)
    return CoSimIO.Info()

@time_decorator()
def ImportData(info):
    print("ImportDataInfo: ")
    print(info)

    identifier = info.GetString("identifier")
    
    imported_data = CoSimIO.DoubleVector()

    settings = CoSimIO.Info()
    settings.SetString("connection_name", s_connection_name)
    settings.SetString("identifier", identifier)

    CoSimIO.ImportData(settings, imported_data)
    return CoSimIO.Info()

@time_decorator()
def ExportData(info):
    print("ExportDataInfo: ")
    print(info)
    identifier = info.GetString("identifier")

    if identifier == "heatflow_External":
        data_to_be_send=CoSimIO.DoubleVector(var.heatflow_External)

    elif identifier == "heatflow_Channel0":
        data_to_be_send=CoSimIO.DoubleVector(var.heatflow_Channel0)

    elif identifier == "heatflow_Channel1":
        data_to_be_send=CoSimIO.DoubleVector(var.heatflow_Channel1)

    elif identifier == "heatflow_Channel2":
        data_to_be_send=CoSimIO.DoubleVector(var.heatflow_Channel2)

    elif identifier == "temperature_External":
        data_to_be_send=CoSimIO.DoubleVector(var.temperature_External)

    elif identifier == "temperature_Channel0":
        data_to_be_send=CoSimIO.DoubleVector(var.temperature_Channel0)

    elif identifier == "temperature_Channel1":
        data_to_be_send=CoSimIO.DoubleVector(var.temperature_Channel1)

    elif identifier == "temperature_Channel2":
        data_to_be_send=CoSimIO.DoubleVector(var.temperature_Channel2)

    settings = CoSimIO.Info()
    settings.SetString("connection_name", s_connection_name)
    settings.SetString("identifier", identifier)
        
    CoSimIO.ExportData(settings, data_to_be_send)
    return CoSimIO.Info()

@time_decorator()
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

@time_decorator()
def ExportMesh(info):
    print("ExportMeshInfo: ")
    print(info)
    identifier = info.GetString("identifier")
    
    model = Kratos.Model()
    model_part = model.CreateModelPart(identifier)
    
    if identifier == "External":
        model_part.CreateNewNode(1, 0.00, 0.1, 0.08)
        model_part.CreateNewNode(2, 0.04, 0.2, 0.08)
        model_part.CreateNewNode(3, 0.09, 0.1, 0.08)
        model_part.CreateNewNode(4, 0.14, 0.2, 0.08)
    elif identifier == "Channel0":
        model_part.CreateNewNode(1,  0.03, 0.1, 0.02)
        model_part.CreateNewNode(2,  0.05, 0.1, 0.02)
        model_part.CreateNewNode(3,  0.03, 0.2, 0.02)
        model_part.CreateNewNode(4,  0.05, 0.2, 0.02)
        model_part.CreateNewNode(5,  0.02, 0.1, 0.03)
        model_part.CreateNewNode(6,  0.02, 0.1, 0.05)
        model_part.CreateNewNode(7,  0.02, 0.2, 0.03)
        model_part.CreateNewNode(8,  0.02, 0.2, 0.05)
        model_part.CreateNewNode(9,  0.06, 0.1, 0.03)
        model_part.CreateNewNode(10, 0.06, 0.1, 0.05)
        model_part.CreateNewNode(11, 0.06, 0.2, 0.03)
        model_part.CreateNewNode(12, 0.06, 0.2, 0.05)
        model_part.CreateNewNode(13, 0.03, 0.1, 0.06)
        model_part.CreateNewNode(14, 0.05, 0.1, 0.06)
        model_part.CreateNewNode(15, 0.03, 0.2, 0.06)
        model_part.CreateNewNode(16, 0.05, 0.2, 0.06)
    elif identifier == "Channel1":
        model_part.CreateNewNode(1,  0.08, 0.1, 0.02)
        model_part.CreateNewNode(2,  0.10, 0.1, 0.02)
        model_part.CreateNewNode(3,  0.08, 0.2, 0.02)
        model_part.CreateNewNode(4,  0.10, 0.2, 0.02)
        model_part.CreateNewNode(5,  0.07, 0.1, 0.03)
        model_part.CreateNewNode(6,  0.07, 0.1, 0.05)
        model_part.CreateNewNode(7,  0.07, 0.2, 0.03)
        model_part.CreateNewNode(8,  0.07, 0.2, 0.05)
        model_part.CreateNewNode(9,  0.13, 0.1, 0.03)
        model_part.CreateNewNode(10, 0.13, 0.1, 0.05)
        model_part.CreateNewNode(11, 0.13, 0.2, 0.03)
        model_part.CreateNewNode(12, 0.13, 0.2, 0.05)
        model_part.CreateNewNode(13, 0.08, 0.1, 0.06)
        model_part.CreateNewNode(14, 0.10, 0.1, 0.06)
        model_part.CreateNewNode(15, 0.08, 0.2, 0.06)
        model_part.CreateNewNode(16, 0.10, 0.2, 0.06)
    elif identifier == "Channel2":
        model_part.CreateNewNode(1,  0.13, 0.1, 0.02)
        model_part.CreateNewNode(2,  0.15, 0.1, 0.02)
        model_part.CreateNewNode(3,  0.13, 0.2, 0.02)
        model_part.CreateNewNode(4,  0.15, 0.2, 0.02)
        model_part.CreateNewNode(5,  0.12, 0.1, 0.03)
        model_part.CreateNewNode(6,  0.12, 0.1, 0.05)
        model_part.CreateNewNode(7,  0.12, 0.2, 0.03)
        model_part.CreateNewNode(8,  0.12, 0.2, 0.05)
        model_part.CreateNewNode(9,  0.18, 0.1, 0.03)
        model_part.CreateNewNode(10, 0.18, 0.1, 0.05)
        model_part.CreateNewNode(11, 0.18, 0.2, 0.03)
        model_part.CreateNewNode(12, 0.18, 0.2, 0.05)
        model_part.CreateNewNode(13, 0.13, 0.1, 0.06)
        model_part.CreateNewNode(14, 0.15, 0.1, 0.06)
        model_part.CreateNewNode(15, 0.13, 0.2, 0.06)
        model_part.CreateNewNode(16, 0.15, 0.2, 0.06)
    
    info = CoSimIO.Info()
    info.SetString("identifier", identifier)
    info.SetString("connection_name", s_connection_name)

    CoSimIO.ExportMesh(info, model_part)
    return info

# Connection Settings
settings = CoSimIO.Info()
settings.SetString("my_name", "run_dummy")
settings.SetString("connect_to", "dummy")
settings.SetInt("echo_level", 1)
settings.SetString("version", "1.25")
settings.SetString("communication_format", "file")

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

