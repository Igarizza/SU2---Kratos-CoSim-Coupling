import KratosMultiphysics as Kratos
from KratosMultiphysics.CoSimulationApplication import CoSimIO
import numpy as np
from KratosMultiphysics.OptimizationApplication.utilities.logger_utilities import time_decorator


def cosimio_check_equal(a, b):
    assert a == b

s_connection_name = ""
# constant variables

class Variables:
    B_in_cm: float = 100
    C_in_cm: float = 100
    u: float = 5.111912504285256453e+00
    r: float = 1.905696411667853196e+00
    theta_0_in_rad: float = 0.26
    # input variables
    q_in_N_per_cm_squared = 1
    psi_in_rad = 0.05

    # coupling variables
    L_in_N = 0.0
    phi_in_rad = 0.0

var = Variables()

@time_decorator()
def AdvanceInTime(info):
    settings = CoSimIO.Info()
    settings.SetString("identifier", "AdvanceInTime")
    settings.SetString("connection_name", s_connection_name)
    settings.SetDouble("current_time", 1.0)

    CoSimIO.ExportInfo(settings)
    return CoSimIO.Info()

@time_decorator()
def InitializeSolutionStep(info):
    settings = CoSimIO.Info()
    settings.SetString("connection_name", s_connection_name)
    settings.SetString("identifier", "InitializeSolutionStep")
    CoSimIO.ExportInfo(settings)
    return CoSimIO.Info()

@time_decorator()
def Predict(info):
    settings = CoSimIO.Info()
    settings.SetString("connection_name", s_connection_name)
    settings.SetString("identifier", "Predict")
    CoSimIO.ExportInfo(settings)
    return CoSimIO.Info()

@time_decorator()
def SolveSolutionStep(info):
    S = var.B_in_cm * var.C_in_cm

    q = var.q_in_N_per_cm_squared
    theta_0 = var.theta_0_in_rad

    phi = var.phi_in_rad
    psi = var.psi_in_rad

    u = var.u
    r = var.r

    theta = phi + psi

    C_L = u * theta + r * (1 - np.cos((np.pi / 2) * (theta / theta_0)))
    print(f"C_L = {C_L}, S = {S}")

    L_in_N = q * S * C_L
    var.L_in_N = L_in_N
    print(f"fluid::Computed L_in_N = {var.L_in_N}")
    return CoSimIO.Info()

@time_decorator()
def FinalizeSolutionStep(info):
    print("FinalizeSolutionStep")
    # settings = CoSimIO.Info()
    # settings.SetString("connection_name", s_connection_name)
    # settings.SetString("identifier", "info_for_test")
    # settings.SetString("name_for_check", "FinalizeSolutionStep")
    # CoSimIO.ExportInfo(settings)
    return CoSimIO.Info()

@time_decorator()
def OutputSolutionStep(info):
    print("OutputSolutionStep")
    # settings = CoSimIO.Info()
    # settings.SetString("connection_name", s_connection_name)
    # settings.SetString("identifier", "info_for_test")
    # settings.SetString("name_for_check", "OutputSolutionStep")
    # CoSimIO.ExportInfo(settings)
    return CoSimIO.Info()

@time_decorator()
def ImportData(info):
    settings = CoSimIO.Info()
    settings.SetString("connection_name", s_connection_name)
    settings.SetString("identifier", info.GetString("identifier"))
    imported_data = CoSimIO.DoubleVector()
    CoSimIO.ImportData(settings, imported_data)
    print("import", info.GetString("identifier"), imported_data, len(imported_data))
    var.phi_in_rad = imported_data[0]
    print(f"imported phi_in_rad = {var.phi_in_rad}")
    ### Need to apply the coming displacements
    return CoSimIO.Info()

@time_decorator()
def ExportData(info):
    data_to_be_send = CoSimIO.DoubleVector([var.L_in_N])
    settings = CoSimIO.Info()
    settings.SetString("connection_name", s_connection_name)
    settings.SetString("identifier", info.GetString("identifier"))
    return_info = CoSimIO.ExportData(settings, data_to_be_send)
    print("export", info.GetString("identifier"), data_to_be_send, len(data_to_be_send))
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
    settings = CoSimIO.Info()
    settings.SetString("connection_name", s_connection_name)
    settings.SetString("identifier", "S_Node")
    model = Kratos.Model()
    model_part = model.CreateModelPart("S_Node")
    model_part.CreateNewNode(1, 0.0, 0.0, 0.0)
    return_info = CoSimIO.ExportMesh(settings, model_part)
    return CoSimIO.Info()

# Connection Settings
settings = CoSimIO.Info()
settings.SetString("my_name", "run_fluid")
settings.SetString("connect_to", "fluid")
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

