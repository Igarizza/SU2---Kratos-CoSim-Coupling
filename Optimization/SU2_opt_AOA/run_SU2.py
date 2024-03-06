from KratosMultiphysics.CoSimulationApplication import CoSimIO
# import CoSimIO
from KratosMultiphysics.OptimizationApplication.utilities.logger_utilities import time_decorator
import KratosMultiphysics as Kratos
import pysu2
from mpi4py import MPI
import numpy as np
from ExIO.EX_SU2_IO import EX_SU2_IO


### SU2 Main start

MPI = MPI
comm = MPI.COMM_WORLD
myid = comm.Get_rank()
# Initialize the flow driver of SU2, this includes solver preprocessing

ConfigFile = "inv_ONERAM6.cfg"
new_ConfigFile = "new_inv_ONERAM6.cfg"
FlowMarkerName = "WING"
mesh_file_name = "mesh_ONERAM6_inv_ffd.su2"


### Variables
class Variables:
    phi_in_degrees = 10.0
    EX_SU2_IO = EX_SU2_IO(ConfigFile, FlowMarkerName, mesh_file_name, comm)
    sorted_displcement = None

var = Variables()

### CoSIM Register

s_connection_name = ""

def cosimio_check_equal(a, b):
    assert a == b

@time_decorator()
def AdvanceInTime(info):
    settings = CoSimIO.Info()
    settings.SetString("identifier", "AdvanceInTime")
    settings.SetString("connection_name", s_connection_name)
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

    # Create new conf file
    with open(ConfigFile, 'r') as file:
    # read a list of lines into data
        data = file.readlines()

    line_number = 0
    for number, line in enumerate(data):
        if "AOA" in line:
            line_number = number
    angle = var.phi_in_degrees + var.psi_in_degrees
    data[line_number] = f"AOA = {angle} \n"
    print(data[line_number])

    with open(new_ConfigFile, 'w') as file:
        file.writelines( data )

    # Create new solver
    send_msg = {'method_name': "CreateNewSU2"}
    for i in range(1, comm.size):
        comm.send(send_msg, dest=i, tag=12)
    var.EX_SU2_IO.flow_driver = pysu2.CSinglezoneDriver(new_ConfigFile, 1, comm)

    # Solve CFD simulation using SU2

    send_msg = {'method_name': "SolveSolutionStep"}
    for i in range(1, comm.size):
        comm.send(send_msg, dest=i, tag=12)

    var.EX_SU2_IO.SolveCFD()

    return CoSimIO.Info()

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
    return settings

@time_decorator()
def ImportData(info):
    identifier = info.GetString("identifier")
    settings = CoSimIO.Info()
    settings.SetString("connection_name", s_connection_name)
    settings.SetString("identifier", info.GetString("identifier"))

    if identifier == "alpha":
        imported_data = CoSimIO.DoubleVector()
        CoSimIO.ImportData(settings, imported_data)
        print("import", info.GetString("identifier"), imported_data, len(imported_data))
        var.phi_in_degrees = imported_data[0]
        print(f"imported phi_in_degrees = {var.phi_in_degrees}")
    return settings

@time_decorator()
def ExportData(info):
    raise RuntimeError("ExportData")
    return CoSimIO.Info()

@time_decorator()
def ImportMesh(info):
    raise RuntimeError("ImportMesh")

@time_decorator()
def ExportMesh(info):
    model = Kratos.Model()
    # Check the identifier
    identifier = info.GetString("identifier")
    model_part = model.CreateModelPart(identifier)
    if identifier == "S_Node":
        model_part.CreateNewNode(1, 0.0, 0.0, 0.0)

    info = CoSimIO.Info()
    info.SetString("identifier", identifier)
    info.SetString("connection_name", s_connection_name)
    return_info = CoSimIO.ExportMesh(info, model_part)
    print(f"ExportMesh: {identifier}, number of nodes : {1}")
    return return_info

while True:
    if myid == 0:
        # Connection Settings
        settings = CoSimIO.Info()
        settings.SetString("my_name", "run_SU2")
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
                var.EX_SU2_IO.SolveCFD()
            elif recv_msg["method_name"] == "setFluidDisplacements":
                imported_data = comm.recv(source=0, tag=151)
                var.EX_SU2_IO.setFluidDisplacements(imported_data)
            elif recv_msg["method_name"] == "ExportData":
                var.EX_SU2_IO.get_force()
            elif recv_msg["method_name"] == "CreateNewSU2":
                var.EX_SU2_IO.flow_driver = pysu2.CSinglezoneDriver(new_ConfigFile, 1, comm)
            elif recv_msg["method_name"] == "Exit":
                break