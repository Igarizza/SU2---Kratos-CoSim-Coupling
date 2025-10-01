from KratosMultiphysics.CoSimulationApplication import CoSimIO
import KratosMultiphysics as Kratos


settings = CoSimIO.Info()
settings.SetString("my_name", "b")
settings.SetString("connect_to", "a")
settings.SetInt("echo_level", 1)
settings.SetString("version", "1.25")
settings.SetString("communication_format", "file")

# Connecting
return_info = CoSimIO.Connect(settings)
s_connection_name = return_info.GetString("connection_name")


settings = Kratos.Parameters("""
    {"adjoint": false}
                                """) 
info = CoSimIO.Info()
info.SetString("connection_name", s_connection_name)
info.SetString("identifier", "run_control")
info.SetString("control_signal", "SolveSolutionStep")
if settings:
    info.SetInfo("info", CoSimIO.InfoFromParameters(settings))
CoSimIO.ExportInfo(info)


info = CoSimIO.Info()
info.SetString("connection_name", s_connection_name)
info.SetString("identifier", "info_exchange")
imported_info = CoSimIO.ImportInfo(info)
data = imported_info.GetDouble("value")
print(f"{s_connection_name}: {data}")


# Disconnecting


disconnect_settings = CoSimIO.Info()
disconnect_settings.SetString("connection_name", s_connection_name)
return_info = CoSimIO.Disconnect(disconnect_settings)