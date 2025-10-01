from KratosMultiphysics.CoSimulationApplication import CoSimIO
import KratosMultiphysics as Kratos


def SolveSolutionStep(info):

    print(info)

    print("Solving Solution Step")

    eff_value = 1.0

    info_to_export = CoSimIO.Info()
    info_to_export.SetString("connection_name", s_connection_name)
    info_to_export.SetString("identifier", "info_exchange")
    info_to_export.SetDouble("value", eff_value)
    return_info = CoSimIO.ExportInfo(info_to_export)

    return CoSimIO.Info()

settings = CoSimIO.Info()
settings.SetString("my_name", "a")
settings.SetString("connect_to", "b")
settings.SetInt("echo_level", 1)
settings.SetString("version", "1.25")
settings.SetString("communication_format", "file")

# Connecting
return_info = CoSimIO.Connect(settings)
s_connection_name = return_info.GetString("connection_name")


fct_info = CoSimIO.Info()
fct_info.SetString("connection_name", s_connection_name)
fct_info.SetString("function_name", "SolveSolutionStep")
CoSimIO.Register(fct_info,           SolveSolutionStep)

# externally orchestrated
run_info = CoSimIO.Info()
run_info.SetString("connection_name", s_connection_name)
CoSimIO.Run(run_info)


# Disconnecting
disconnect_settings = CoSimIO.Info()
disconnect_settings.SetString("connection_name", s_connection_name)
return_info = CoSimIO.Disconnect(disconnect_settings)