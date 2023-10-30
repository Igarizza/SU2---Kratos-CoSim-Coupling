from KratosMultiphysics.CoSimulationApplication import CoSimIO
import KratosMultiphysics.CoSimulationApplication as Kratos_CoSim
from KratosMultiphysics.OptimizationApplication.utilities.logger_utilities import time_decorator
import KratosMultiphysics as Kratos
import pysu2
import SU2
from copy import deepcopy
import pylab as plt


from mpi4py import MPI
import numpy as np

MPI = MPI
comm = MPI.COMM_WORLD
myid = comm.Get_rank()
have_MPI = True

### Variables
class Variables:
    C_in_cm: float = 100
    B_in_cm: float = 100
    # coupling variables
    L_in_N = 0.0
    phi_in_rad = 0.0
    phi_in_degrees = 10.0
    psi_in_rad = 0.05
    psi_in_degrees = 0.05 / 3.14 * 180
    initial_run = True

var = Variables()

### SU2 Main start

# Initialize the flow driver of SU2, this includes solver preprocessing

confFile = "inv_ONERAM6.cfg"
new_confFile = "new_inv_ONERAM6.cfg"
FlowMarkerName = "WING"

if have_MPI == True:
    comm.barrier()

model = Kratos.Model()
model_part = model.CreateModelPart(FlowMarkerName)
model_part.CreateNewNode(1, 0.0, 0.0, 0.0)

### CoSIM Register

def cosimio_check_equal(a, b):
    assert a == b

s_connection_name = ""

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
    if myid == 0:
            print("\n------------------------------ Begin Solver -----------------------------\n")

    with open(confFile, 'r') as file:
    # read a list of lines into data
        data = file.readlines()

    line_number = 0
    for number, line in enumerate(data):
        if "AOA" in line:
            line_number = number
    angle = var.phi_in_degrees + var.psi_in_degrees
    data[line_number] = f"AOA = {angle} \n"
    print(data[line_number])

    with open(new_confFile, 'w') as file:
        file.writelines( data )

    send_msg = {'method_name': "SolveSolutionStep"}
    for i in range(1, comm.size):
        comm.send(send_msg, dest=i, tag=12)

    SolveCFD()

    return CoSimIO.Info()

def SolveCFD():
    FlowDriver = pysu2.CSinglezoneDriver(new_confFile, 1, comm)
    FlowDriver.StartSolver()

    FlowMarkerList = FlowDriver.GetAllBoundaryMarkersTag()              # Get all the flow boundary tags
    FlowMarkerIDs = FlowDriver.GetAllBoundaryMarkers()                  # Get all the associated indices to the flow markers
    if FlowMarkerName in FlowMarkerList and FlowMarkerName in FlowMarkerIDs.keys():
        FlowMarkerID = FlowMarkerIDs[FlowMarkerName]                      # Check if the flow FSI marker exists
        nVertex_Marker_Flow = FlowDriver.GetNumberVertices(FlowMarkerID)    # Get the number of vertices of the flow FSI marker
    else:
        nVertex_Marker_Flow = 0

    LocalForce = 0.0
    for iVertex in range(nVertex_Marker_Flow):
        if not FlowDriver.IsAHaloNode(FlowMarkerID, iVertex):
            newF = FlowDriver.GetFlowLoad(FlowMarkerID, iVertex)
            LocalForce += newF[2]

    ForceList=comm.allgather(LocalForce)
    Force = 0.0
    for i in ForceList:
        Force += i

    del LocalForce

    var.L_in_N = Force

@time_decorator()
def FinalizeSolutionStep(info):
    settings = CoSimIO.Info()
    settings.SetString("connection_name", s_connection_name)
    settings.SetString("identifier", "FinalizeSolutionStep")
    CoSimIO.ExportInfo(settings)
    return CoSimIO.Info()

@time_decorator()
def OutputSolutionStep(info):
    settings = CoSimIO.Info()
    settings.SetString("connection_name", s_connection_name)
    settings.SetString("identifier", "OutputSolutionStep")
    CoSimIO.ExportInfo(settings)
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
    var.phi_in_degrees = var.phi_in_rad / 3.14 * 180
    print(f"imported phi_in_degrees = {var.phi_in_degrees}")
    ### Need to apply the coming displacements
    return CoSimIO.Info()

@time_decorator()
def ExportData(info):
    data_to_be_send=CoSimIO.DoubleVector([var.L_in_N])
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

    # Exporting mesh to Kratos
    info = CoSimIO.Info()
    info.SetString("identifier", "airfoil")
    info.SetString("connection_name", s_connection_name)
    return_info = CoSimIO.ExportMesh(info, model_part)
    return info

if myid == 0:
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
    for i in range(1, comm.size):
        send_msg = {'method_name': "Exit"}
        comm.send(send_msg, dest=i, tag=12)
else:
    while True:
        recv_msg = comm.recv(source=0, tag=12)
        if recv_msg["method_name"] == "SolveSolutionStep":
            SolveCFD()
        elif recv_msg["method_name"] == "Exit":
            break
