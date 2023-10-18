import pysu2
import KratosMultiphysics as Kratos
from mpi4py import MPI
import numpy as np

class EX_SU2_IO():
    def __init__(self, config_file: str, flow_marker_name: str, mesh_file: str, comm: MPI.COMM_WORLD) -> None:
        self.config_file = config_file
        self.flow_driver = pysu2.CSinglezoneDriver(config_file, 1, comm)
        self.flow_marker_name = flow_marker_name
        self.mesh_file = mesh_file
        self.comm = comm
        self.rootProcess = 0
        self.FlowMarkerID = -1
        self.nLocalFluidInterfaceNodes = 0
        self.nLocalFluidInterfaceHaloNode = 0
        self.nLocalFluidInterfacePhysicalNodes = 0
        self.fluidInterfaceProcessors = list()          #list of partitions where there are fluid interface nodes
        self.fluidSolverProcessors = list()		        #list of partitions where the fluid solver is initialized
        self.fluidGlobalIndexRange = {}			        #contains the global FSI indexing of each fluid interface node for all partitions
        self.FluidInterfaceNodeCoord=[]                 # Nodes Interface Coordinates x,y,z
        self.FluidInterfaceNodeCoordX=[]
        self.FluidInterfaceNodeCoordY=[]
        self.FluidInterfaceNodeCoordZ=[]
        self.FluidForceList={}
        self.nNodesInElementList = []                   # Array of numner of nodes in the Element for CoSim
        self.mesh = []                                  # Mesh info for CoSim
        self.FluidInterfaceElementsIds = []
        self.index_Global_ID_physical_nodes_map = {}    # key - Id, value local index
        self.haveFluidInterface = False
        self.FluidLocalHaloNodeList = {}
        self.FluidHaloNodeList = {}
        self.set_up_variables()
        self.time_iteration = 0.0

    def set_up_variables(self):
        self.MPIPrint("Start set_up_variables ...")
        myid = self.comm.Get_rank()
        MPIsize = self.comm.Get_size()
        FlowMarkerIDs = self.flow_driver.GetAllBoundaryMarkers()
        if self.flow_marker_name in FlowMarkerIDs.keys():
            self.FlowMarkerID = FlowMarkerIDs[self.flow_marker_name]                      # Check if the flow FSI marker exists
            self.nLocalFluidInterfaceNodes = self.flow_driver.GetNumberVertices(self.FlowMarkerID)    # Get the number of vertices of the flow FSI marker

        if self.nLocalFluidInterfaceNodes != 0:
            self.haveFluidInterface = True
        print(f"myid = {myid}, self.nLocalFluidInterfaceNodes = {self.nLocalFluidInterfaceNodes}, self.haveFluidInterface = {self.haveFluidInterface}")

        if self.flow_driver != None:
            self.haveFluidSolver = True
        # --- Exchange information about processors on which the solvers are defined and where the interface nodes are lying ---

        if self.haveFluidSolver == True:
            sendBufFluid = np.array(int(1))
        else:
            sendBufFluid = np.array(int(0))
        if self.haveFluidInterface == True:
            sendBufFluidInterface = np.array(int(1))
        else:
            sendBufFluidInterface = np.array(int(0))
        rcvBufFluid = np.zeros(MPIsize, dtype = int)
        rcvBufFluidInterface = np.zeros(MPIsize, dtype = int)
        self.comm.Allgather(sendBufFluid, rcvBufFluid)
        self.comm.Allgather(sendBufFluidInterface, rcvBufFluidInterface)
        for iProc in range(MPIsize):
            if rcvBufFluid[iProc] == 1:
                self.fluidSolverProcessors.append(iProc)
            if rcvBufFluidInterface[iProc] == 1:
                self.fluidInterfaceProcessors.append(iProc)
        del sendBufFluid, rcvBufFluid, sendBufFluidInterface, rcvBufFluidInterface

        self.MPIBarrier()
        # --- Calculate the total number of nodes at the fluid interface (sum over all the partitions) ---
        # Calculate the number of halo nodes on each partition

        for iVertex in range(self.nLocalFluidInterfaceNodes):
            if self.flow_driver.IsAHaloNode(self.FlowMarkerID, iVertex) == True:
                GlobalIndex = self.flow_driver.GetVertexGlobalIndex(self.FlowMarkerID, iVertex)
                self.FluidLocalHaloNodeList[GlobalIndex] = iVertex
                self.nLocalFluidInterfaceHaloNode += 1
         # Calculate the number of physical (= not halo) nodes on each partition
        self.nLocalFluidInterfacePhysicalNodes = self.nLocalFluidInterfaceNodes - self.nLocalFluidInterfaceHaloNode
        # print(f"myid = {myid}, number of nodes = {self.nLocalFluidInterfacePhysicalNodes}")
        self.FluidHaloNodeList = self.comm.allgather(self.FluidLocalHaloNodeList)

         # --- Calculate the total number of nodes (with and without halo) at the fluid interface (sum over all the partitions) and broadcast the number accross all processors ---
        sendBuffHalo = np.array(int(self.nLocalFluidInterfaceNodes))
        sendBuffPhysical = np.array(int(self.nLocalFluidInterfacePhysicalNodes))
        rcvBuffHalo = np.zeros(1, dtype=int)
        rcvBuffPhysical = np.zeros(1, dtype=int)

        self.comm.barrier()
        self.comm.Allreduce(sendBuffHalo,rcvBuffHalo,op=MPI.SUM)
        self.comm.Allreduce(sendBuffPhysical,rcvBuffPhysical,op=MPI.SUM)
        self.nFluidInterfaceNodes = rcvBuffHalo[0]
        self.nFluidInterfacePhysicalNodes = rcvBuffPhysical[0]
        del sendBuffHalo, rcvBuffHalo, sendBuffPhysical, rcvBuffPhysical

        # --- Store the number of physical interface nodes on each processor and allgather the information ---
        self.fluidPhysicalInterfaceNodesDistribution = np.zeros(MPIsize, dtype=int)

        sendBuffPhysical = np.array(int(self.nLocalFluidInterfacePhysicalNodes))
        self.comm.Allgather(sendBuffPhysical,self.fluidPhysicalInterfaceNodesDistribution)
        del sendBuffPhysical

         # --- Calculate and store the global indexing of interface physical nodes on each processor and allgather the information ---

        if myid in self.fluidInterfaceProcessors:
            globalIndexStart = 0
            for iProc in range(myid):
                globalIndexStart += self.fluidPhysicalInterfaceNodesDistribution[iProc]
            globalIndexStop = globalIndexStart + self.nLocalFluidInterfacePhysicalNodes-1
        else:
            globalIndexStart = 0
            globalIndexStop = 0
        self.fluidGlobalIndexRange[myid] = [globalIndexStart,globalIndexStop]
        self.fluidGlobalIndexRange = self.comm.allgather(self.fluidGlobalIndexRange)
        # print(self.fluidGlobalIndexRange)

        self.MPIPrint('Total number of fluid interface nodes (halo nodes included) : {}'.format(self.nFluidInterfaceNodes))
        self.MPIPrint('Total number of fluid interface nodes : {}'.format(self.nFluidInterfacePhysicalNodes))

        ## --- Get and store the coordinates of the verticies
        FluidLocalInterfaceCoord=[]
        FluidLocalInterfaceCoordX=[]
        FluidLocalInterfaceCoordY=[]
        FluidLocalInterfaceCoordZ=[]
        for iVertex in range(self.nLocalFluidInterfaceNodes):
            if self.flow_driver.IsAHaloNode(self.FlowMarkerID, iVertex) == False:
                coords = self.flow_driver.GetInitialMeshCoord(self.FlowMarkerID, iVertex)
                FluidLocalInterfaceCoord.append(coords[0])
                FluidLocalInterfaceCoord.append(coords[1])
                FluidLocalInterfaceCoord.append(coords[2])
                FluidLocalInterfaceCoordX.append(coords[0])
                FluidLocalInterfaceCoordY.append(coords[1])
                FluidLocalInterfaceCoordZ.append(coords[2])

        FluidInterfaceCoordList = self.comm.allgather(FluidLocalInterfaceCoord)
        FluidInterfaceCoordListX = self.comm.allgather(FluidLocalInterfaceCoordX)
        FluidInterfaceCoordListY = self.comm.allgather(FluidLocalInterfaceCoordY)
        FluidInterfaceCoordListZ = self.comm.allgather(FluidLocalInterfaceCoordZ)
        for i in range (MPIsize):
            for j in FluidInterfaceCoordList[i]:
                self.FluidInterfaceNodeCoord.append(j)
            for j in FluidInterfaceCoordListX[i]:
                self.FluidInterfaceNodeCoordX.append(j)
            for j in FluidInterfaceCoordListY[i]:
                self.FluidInterfaceNodeCoordY.append(j)
            for j in FluidInterfaceCoordListZ[i]:
                self.FluidInterfaceNodeCoordZ.append(j)

        self.FluidInterfaceNodeCoord = np.array (self.FluidInterfaceNodeCoord)
        self.FluidInterfaceNodeCoordX = np.array (self.FluidInterfaceNodeCoordX)
        self.FluidInterfaceNodeCoordY = np.array (self.FluidInterfaceNodeCoordY)
        self.FluidInterfaceNodeCoordZ = np.array (self.FluidInterfaceNodeCoordZ)
        del FluidLocalInterfaceCoord

    ## --- Get Global indexes of the nodes (Wet-Interface)
        nodeIDs=[]
        self.localNodesIDs = []
        for iVertex in range (self.nLocalFluidInterfaceNodes):
            if self.flow_driver.IsAHaloNode(self.FlowMarkerID, iVertex) == False:
                node_id = int(self.flow_driver.GetVertexGlobalIndex(self.FlowMarkerID,iVertex))
                self.localNodesIDs.append(node_id)

        NodeIDsList = self.comm.allgather(self.localNodesIDs)
        for i in range (MPIsize):
            for j in NodeIDsList[i]:
                # nodeIDs.append(j+1)
                nodeIDs.append(j)
        self.nodeIDs = np.array (nodeIDs)

        del nodeIDs

        ### Get Elements from a mesh file
        self.get_Elements(self.mesh_file, self.flow_marker_name)

        if myid == 0:
            for node_id in self.nodeIDs:
                pass

    # --- set the Initial Values of the coordinates of the mesh (halo is includied)
        self.Initil_particle_displacements(self.FluidInterfaceNodeCoord)

    def get_Elements (self, mesh_file, wet_interface):
        lookup = 'MARKER_TAG= '
        for i in wet_interface:
            if i!=')' and i!='(' and i!=' ':
                lookup = lookup + i
        self.nNodesInElement = 3
        filename = mesh_file
        localElements = []
        with open(filename) as currentLoad:
            for num, line in enumerate(currentLoad, 1):
                if lookup in line:
                    num=num+1
                    break


            for num, line in enumerate(currentLoad, 1):
                if num == 1:
                    line_string = line.split(' ')
                    self.nFluidInterfaceElements = int(line_string[1])
                if num < self.nFluidInterfaceElements+2 and num > 1:
                    line_string = line.split('\t')

                    if line_string [0] == 5:
                        self.nNodesInElement = 3
                    if line_string [0] == 9:
                        self.nNodesInElement = 4
                    for i,line in enumerate(line_string,1):
                        if i>1 and i < self.nNodesInElement + 2:
                            localElements.append(int(line))


        self.FluidInterfaceElementsNodes = np.array(localElements)
        for i in range (self.nFluidInterfaceElements):
            self.nNodesInElementList.append(self.nNodesInElement)
        self.nNodesInElementList = np.array(self.nNodesInElementList)

    def particle_displacements (self,FluidSolver,newDisplacements):

        myid = self.comm.Get_rank()

        self.DisplacementListX = []
        self.DisplacementListY = []
        self.DisplacementListZ = []
        localDisplacementArrayX = []
        localDisplacementArrayY = []
        localDisplacementArrayZ = []

        localDisplacementListX = {}
        localDisplacementListY = {}
        localDisplacementListZ = {}
        if self.haveFluidInterface:
            myLocalIndexStart = self.fluidGlobalIndexRange[myid].get(myid)[0]*3
            for i in range (self.nLocalFluidInterfacePhysicalNodes):
                GlobalIndex = self.localNodesIDs[i]
                localDisplacementListX[GlobalIndex] =  newDisplacements[myLocalIndexStart + i*3]
                localDisplacementListY[GlobalIndex] =  newDisplacements[myLocalIndexStart + i*3+1]
                localDisplacementListZ[GlobalIndex] =  newDisplacements[myLocalIndexStart + i*3+2]

        self.DisplacementListX = self.comm.allgather(localDisplacementListX)
        self.DisplacementListY = self.comm.allgather(localDisplacementListY)
        self.DisplacementListZ = self.comm.allgather(localDisplacementListZ)
        del localDisplacementArrayX, localDisplacementArrayY, localDisplacementArrayZ, localDisplacementListX, localDisplacementListY, localDisplacementListZ

    def Initil_particle_displacements (self,newDisplacements):

        myid = self.comm.Get_rank()

        localInitDisplacementArrayX = []
        localInitDisplacementArrayY = []
        localInitDisplacementArrayZ = []

        localInitDisplacementListX = {}
        localInitDisplacementListY = {}
        localInitDisplacementListZ = {}
        if self.haveFluidInterface:
            myLocalIndexStart = self.fluidGlobalIndexRange[myid].get(myid)[0]*3
            for i in range (0,self.nLocalFluidInterfacePhysicalNodes*3,3):
                localInitDisplacementArrayX.append(newDisplacements[myLocalIndexStart+i])
                localInitDisplacementArrayY.append(newDisplacements[myLocalIndexStart+i+1])
                localInitDisplacementArrayZ.append(newDisplacements[myLocalIndexStart+i+2])

        if self.haveFluidInterface:
            myLocalIndex = 0
            for iVertex in range (self.nLocalFluidInterfaceNodes):
                if self.flow_driver.IsAHaloNode(self.FlowMarkerID, iVertex) == False:
                    GlobalIndex = self.flow_driver.GetVertexGlobalIndex(self.FlowMarkerID, iVertex)
                    localInitDisplacementListX[GlobalIndex] =  localInitDisplacementArrayX[myLocalIndex]
                    localInitDisplacementListY[GlobalIndex] =  localInitDisplacementArrayY[myLocalIndex]
                    localInitDisplacementListZ[GlobalIndex] =  localInitDisplacementArrayZ[myLocalIndex]
                    myLocalIndex += 1

        self.InitCoordListX = self.comm.allgather(localInitDisplacementListX)
        self.InitCoordListY = self.comm.allgather(localInitDisplacementListY)
        self.InitCoordListZ = self.comm.allgather(localInitDisplacementListZ)
        del localInitDisplacementArrayX, localInitDisplacementArrayY, localInitDisplacementArrayZ, localInitDisplacementListX, localInitDisplacementListY, localInitDisplacementListZ

    def get_mesh(self):
        for i in range (self.nFluidInterfaceElements):
            self.FluidInterfaceElementsIds.append(i)
        self.mesh.append(self.nFluidInterfacePhysicalNodes)
        self.mesh.append(self.nFluidInterfaceElements)
        self.mesh.append(self.FluidInterfaceNodeCoord)
        self.mesh.append(self.nodeIDs)
        self.mesh.append(self.nNodesInElementList)
        self.mesh.append(self.FluidInterfaceElementsNodes)
        self.mesh.append(np.array(self.FluidInterfaceElementsIds))
        return self.mesh


    def get_force(self):
        LocalForce = []
        self.Force = []
        SingleLocalForce = 0.0
        self.FluidForceLocal = {}
        # --- Get the fluid interface loads from the fluid solver and directly fill the corresponding PETSc vector ---
        for iVertex in range(self.nLocalFluidInterfaceNodes):
            if not self.flow_driver.IsAHaloNode(self.FlowMarkerID, iVertex):
                nodeId = self.flow_driver.GetVertexGlobalIndex(self.FlowMarkerID, iVertex)
                newF = self.flow_driver.GetFlowLoad(self.FlowMarkerID, iVertex)
                self.FluidForceLocal[nodeId] = newF
                LocalForce.append(newF[0])
                LocalForce.append(newF[1])
                LocalForce.append(newF[2])
                SingleLocalForce += newF[2]

        self.FluidForceList = self.comm.allgather(self.FluidForceLocal)
        ForceList=self.comm.allgather(LocalForce)
        self.MPIBarrier()
        for i in ForceList:
            for j in i:
                self.Force.append(j)

        SingleForceList=self.comm.allgather(SingleLocalForce)
        SingleForce = 0.0
        for i in SingleForceList:
            SingleForce += i
        self.Force = np.array(self.Force)
        del LocalForce
        return SingleForce


    def setFluidDisplacements (self,newDisplacements):
        # --- Send the new fluid interface position to the fluid solver (on each partition, halo nodes included) ---
        myid = self.comm.Get_rank()

        self.particle_displacements (self.flow_driver, newDisplacements)
        number = 0
        if self.haveFluidInterface:
            for iVertex in range(self.nLocalFluidInterfaceNodes):
                GlobalIndex = self.flow_driver.GetVertexGlobalIndex(self.FlowMarkerID, iVertex)
                disp = self.searchDisplacement(GlobalIndex,1)
                init_coord = self.searchDisplacement(GlobalIndex,2)
                if disp[0]==0.0 and disp[1]==0.0 and disp[2]==0.0:
                    # print(GlobalIndex)
                    number +=1
                self.flow_driver.SetMeshDisplacement(self.FlowMarkerID, iVertex, disp[0], disp[1], disp[2])
                # print(f"GlobalIndex = {GlobalIndex}, coord_x = {init_coord[0] + disp[0]}")
        self.flow_driver.CommunicateMeshDisplacement()

    def searchDisplacement(self,GlobalIndex,flag=1):
        disp = []
        myid = 0
        if flag == 1:
            while True:
                if GlobalIndex in self.DisplacementListX[myid]:
                    disp.append(self.DisplacementListX[myid].get(GlobalIndex))
                    disp.append(self.DisplacementListY[myid].get(GlobalIndex))
                    disp.append(self.DisplacementListZ[myid].get(GlobalIndex))
                    break
                myid += 1
                if myid == self.comm.Get_size():
                    raise RuntimeError(f"Error code 875, GlobalIndex = {GlobalIndex}")
        if flag == 2:
            while True:
                if self.InitCoordListX[myid].get(GlobalIndex) != None :
                    disp.append(self.InitCoordListX[myid].get(GlobalIndex))
                    disp.append(self.InitCoordListY[myid].get(GlobalIndex))
                    disp.append(self.InitCoordListZ[myid].get(GlobalIndex))
                    break
                myid += 1
                if myid == self.comm.Get_size():
                    raise RuntimeError(f"Error code 876, GlobalIndex = {GlobalIndex}")
        if flag == 3:
            while True:
                if GlobalIndex in self.FluidForceList[myid]:
                    disp = self.FluidForceList[myid].get(GlobalIndex)
                    break
                myid +=1
                if myid == self.comm.Get_size():
                    raise RuntimeError(f"Error code 877, GlobalIndex = {GlobalIndex}")
        return disp


    def MPIPrint(self, message):
        """
        Print a message on screen only from the master process.
        """
        myid = self.comm.Get_rank()

        if myid == self.rootProcess:
            print(message)


    def MPIBarrier(self):
        """
        Perform a synchronization barrier in case of parallel run with MPI.
        """
        self.comm.barrier()

    def SolveCFD(self):
        self.flow_driver.StartSolver()