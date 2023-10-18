import KratosMultiphysics as Kratos
from KratosMultiphysics.CoSimulationApplication import CoSimIO
from KratosMultiphysics.OptimizationApplication.utilities.logger_utilities import time_decorator

def cosimio_check_equal(a, b):
    assert a == b

s_connection_name = ""
# constant variables

class Variables:
    a_overline: float = 0.25
    C_in_cm: float = 80
    B_in_cm: float = 100
    k_1_in_N_per_cm: float = 8000
    k_2_in_N_per_cm: float = 4000
    z_overline_1: float = 0.2
    z_overline_2: float = 0.7
    # input variables
    q_in_N_per_cm_squared = 1
    psi_in_rad = 0.05

    # coupling variables
    L_in_N = 0.0
    phi_in_rad = 0.26

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
    k_1 = var.k_1_in_N_per_cm
    k_2 = var.k_2_in_N_per_cm
    C = var.C_in_cm

    h_1_overline = var.a_overline - var.z_overline_1
    h_2_overline = var.z_overline_2 - var.a_overline
    p = h_1_overline / h_2_overline

    R_1 = var.L_in_N / (1 + p)
    R_2 = var.L_in_N * p / (1 + p)

    d_1 = R_1 / k_1
    d_2 = R_2 / k_2

    phi_in_rad = (d_1 - d_2) / (C * (var.z_overline_2 - var.z_overline_1))
    var.phi_in_rad = phi_in_rad
    print(f"structure::Computed phi_in_rad = {var.phi_in_rad}")
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
    var.L_in_N = imported_data[0]
    print("import", info.GetString("identifier"), imported_data, len(imported_data))
    print(f"imported L_in_N = {var.L_in_N}")
    ### Need to apply the coming displacements
    return CoSimIO.Info()

@time_decorator()
def ExportData(info):
    data_to_be_send = CoSimIO.DoubleVector([var.phi_in_rad])
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
    settings.SetString("identifier", "S_Node_2")
    model = Kratos.Model()
    model_part = model.CreateModelPart("S_Node_2")
    model_part.CreateNewNode(1, 0.0, 0.0, 0.0)
    return_info = CoSimIO.ExportMesh(settings, model_part)
    return CoSimIO.Info()

# Connection Settings
settings = CoSimIO.Info()
settings.SetString("my_name", "run_structure")
settings.SetString("connect_to", "spring")
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
