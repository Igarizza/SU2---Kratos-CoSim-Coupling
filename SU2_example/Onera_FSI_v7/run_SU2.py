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
FlowMarkerName = "WING"
mesh_file_name = "mesh_ONERAM6_inv_ffd.su2"
EX_SU2_IO = EX_SU2_IO(ConfigFile, FlowMarkerName, mesh_file_name, comm)


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

    send_msg = {'method_name': "SolveSolutionStep"}
    for i in range(1, comm.size):
        comm.send(send_msg, dest=i, tag=12)

    EX_SU2_IO.SolveCFD()
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
    return CoSimIO.Info()

@time_decorator()
def ImportData(info):
    settings = CoSimIO.Info()
    settings.SetString("connection_name", s_connection_name)
    settings.SetString("identifier", info.GetString("identifier"))
    imported_data = CoSimIO.DoubleVector()
    CoSimIO.ImportData(settings, imported_data)
    data_size = len(imported_data)
    imported_data_numpy = np.arange(data_size, dtype=np.float64)
    for i, data in enumerate(imported_data):
        imported_data_numpy[i] = data

    # Kratos mesh returns displacement in sorted order. To apply to SU2 we need to reorder the displacement vector to match SU2 order
    node_ids_sort = np.copy(EX_SU2_IO.nodeIDs)
    node_ids_sort = np.sort(node_ids_sort)
    sorted_displcement = np.empty(EX_SU2_IO.nFluidInterfacePhysicalNodes * 3)

    for node_index in range(EX_SU2_IO.nFluidInterfacePhysicalNodes):
        node_id = node_ids_sort[node_index]
        old_index = np.where(EX_SU2_IO.nodeIDs==node_id)[0][0]
        sorted_displcement[old_index * 3] = imported_data_numpy[node_index * 3]
        sorted_displcement[old_index * 3 + 1] = imported_data_numpy[node_index * 3 + 1]
        sorted_displcement[old_index * 3 + 2] = imported_data_numpy[node_index * 3 + 2]

    send_msg = {'method_name': "setFluidDisplacements"}
    for i in range(1, comm.size):
        comm.send(send_msg, dest=i, tag=12)
    for i in range(1, comm.size):
        comm.send(sorted_displcement, dest=i, tag=151)

    EX_SU2_IO.setFluidDisplacements(sorted_displcement)
    ### Need to apply the coming displacements
    return CoSimIO.Info()

@time_decorator()
def ExportData(info):
    settings = CoSimIO.Info()
    settings.SetString("connection_name", s_connection_name)
    settings.SetString("identifier", info.GetString("identifier"))

    send_msg = {'method_name': "ExportData"}
    for i in range(1, comm.size):
        comm.send(send_msg, dest=i, tag=12)

    EX_SU2_IO.get_force()
    fluid_force = []

    node_ids_sort = np.copy(EX_SU2_IO.nodeIDs)
    node_ids_sort = np.sort(node_ids_sort)


    for node_index in range(EX_SU2_IO.nFluidInterfacePhysicalNodes):
        node_id = node_ids_sort[node_index]
        # We need to find correct value of force using node Id
        node_force = EX_SU2_IO.searchDisplacement(node_id, 3)
        fluid_force.append(node_force[0])
        fluid_force.append(node_force[1])
        fluid_force.append(node_force[2])

    data_to_be_send = CoSimIO.DoubleVector(fluid_force)
    return_info = CoSimIO.ExportData(settings, data_to_be_send)
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
    info.SetString("identifier", FlowMarkerName)
    info.SetString("connection_name", s_connection_name)

    # To speed up the process of creating ModelPart in Kratos, we need to sort the ideas. Otherwise Kratos will sort it on his side.
    node_ids_sort = np.copy(EX_SU2_IO.nodeIDs)
    node_ids_sort = np.sort(node_ids_sort)

    model = Kratos.Model()
    model_part = model.CreateModelPart(FlowMarkerName)

    for node_index in range(EX_SU2_IO.nFluidInterfacePhysicalNodes):
        node_id = node_ids_sort[node_index]
        coords = EX_SU2_IO.searchDisplacement(node_id, 2)
        model_part.CreateNewNode(node_id, coords[0], coords[1], coords[2])

    properties = model_part.CreateNewProperties(1)
    for element_index in range(EX_SU2_IO.nFluidInterfaceElements):
        elem_nodes = [EX_SU2_IO.FluidInterfaceElementsNodes[element_index*3], EX_SU2_IO.FluidInterfaceElementsNodes[element_index*3 + 1], EX_SU2_IO.FluidInterfaceElementsNodes[element_index*3 + 2]]
        model_part.CreateNewElement("Element2D3N", element_index + 1, elem_nodes, properties) # trinagle mesh
    return_info = CoSimIO.ExportMesh(info, model_part)
    return info

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
            EX_SU2_IO.SolveCFD()
        elif recv_msg["method_name"] == "setFluidDisplacements":
            imported_data = comm.recv(source=0, tag=151)
            EX_SU2_IO.setFluidDisplacements(imported_data)
        elif recv_msg["method_name"] == "ExportData":
            EX_SU2_IO.get_force()
        elif recv_msg["method_name"] == "Exit":
            break