import os, sys, shutil, copy
import numpy as np
import scipy as sp
import scipy.spatial.distance as spdist
from math import *

# ----------------------------------------------------------------------
#  FSI Interface Class
# ----------------------------------------------------------------------

class EX_Interface:

    """
    FSI interface class that handles fluid/solid solvers synchronisation and communication
    """

    def __init__(self, FluidSolver, have_MPI):
	"""
	Class constructor. Declare some variables and do some screen outputs.
	"""
        if have_MPI == True:
            from mpi4py import MPI
            self.MPI = MPI
            self.comm = MPI.COMM_WORLD			#MPI World communicator
            self.have_MPI = True
            self.MPIsize = self.comm.Get_size()
            myid = self.comm.Get_rank()
        else:
            self.comm = 0
            self.have_MPI = False
            myid = 0
        self.rootProcess = 0
        """
        Add to common config file
        """
        self.nDim = 3

        self.haveFluidSolver = False  #True if the fluid solver is initialized on the current rank
        self.haveFluidInterface = False			#True if the current rank owns at least one fluid interface node
        self.fluidSolverProcessors = list()		#list of partitions where the fluid solver is initialized
        self.fluidInterfaceProcessors = list()          #list of partitions where there are fluid interface nodes

        self.fluidInterfaceIdentifier = None		#object that can identify the f/s interface within the fluid solver

        self.fluidGlobalIndexRange = {}			#contains the global FSI indexing of each fluid interface node for all partitions

        self.FluidHaloNodeList = {}
        self.FluidLocalHaloNodeList = {}			#contains the the indices (fluid solver indexing) of the halo nodes for each partition
        self.fluidIndexing = {}				#links between the fluid solver indexing and the FSI indexing for the interface nodes

        self.nLocalFluidInterfaceNodes = 0		#number of nodes (halo nodes included) on the fluid interface, on each partition
        self.nLocalFluidInterfaceHaloNode = 0		#number of halo nodes on the fluid intrface, on each partition
        self.nLocalFluidInterfacePhysicalNodes = 0	#number of physical (= non halo) nodes on the fluid interface, on each partition
        self.nFluidInterfaceNodes = 0			#number of nodes on the fluid interface, sum over all the partitions
        self.nFluidInterfacePhysicalNodes = 0		#number of physical nodes on the fluid interface, sum over all partitions
        self.nodeIDs = None                     #Global indexes (Ids) of the wet-Interface
        self.mesh = []                        # Mesh info for a Empire
        self.nFluidInterfaceElements = 0      # Number of elements on the wet interface
        self.nNodesInElement = 0            # Number of nodes in the elements
        self.nNodesInElementList = []         # Array of numner of nodes in the Element for Empire
        self.FluidInterfaceElementsNodes = None    # Array of the nodes Ids for elements
        self.FluidInterfaceNodeCoord=[]                # Nodes Interface Coordinates x,y,z
        self.FluidInterfaceElementsIds = []
        self.FluidInterfaceInitialCoordXNodes = []       #Store initial coordinates X of the nodes in the wet-interface
        self.FluidInterfaceInitialCoordYNodes = []
        self.FluidInterfaceInitialCoordZNodes = []


        self.InitCoordListX = {}
        self.InitCoordListY = {}
        self.InitCoordListZ = {}

        self.DisplacementListX = {}
        self.DisplacementListY = {}
        self.DisplacementListZ = {}




        self.localDisplacementListX = {}
        self.localDisplacementListY = {}
        self.localDisplacementListZ = {}   #local list of new displacements
        self.localDisplacementArrayX = []
        self.localDisplacementArrayY = []
        self.localDisplacementArrayZ = []
        self.fluidCoordNodes = None  #coordinates of the nodes in the wet interface
        self.localFluidInterface_array_X_init = None	#initial fluid interface position on each partition (used for the meshes mapping)
        self.localFluidInterface_array_Y_init = None
        self.localFluidInterface_array_Z_init = None

        self.haloNodesPositionsInit = {}		#initial position of the halo nodes (fluid side only)

        self.fluidInterface_array_DispX = None		#fluid interface displacement
        self.fluidInterface_array_DispY = None
        self.fluidInterface_array_DispZ = None

        self.fluidLoads_array_X = None			#loads on the fluid side of the f/s interface
        self.fluidLoads_array_Y = None
        self.fluidLoads_array_Z = None

        self.FSIIter = 0				#current FSI iteration
        self.unsteady = False		#flag for steady or unsteady simulation (default is steady)

        self.Normal = []                                #nodal normals on wet-interface
        self.CurrCoords = []                            #current nodal coords on wet-interface

    # ---Some screen output ---
        self.MPIPrint('Fluid solver : SU2_CFD')





    def MPIPrint(self, message):
        """
        Print a message on screen only from the master process.
        """
        if self.have_MPI == True:
            myid = self.comm.Get_rank()
        else:
            myid = 0

        if myid == self.rootProcess:
            print(message)


    def MPIBarrier(self):
        """
        Perform a synchronization barrier in case of parallel run with MPI.
        """
        if self.have_MPI == True:
            self.comm.barrier()

    def set_up_variables(self,FluidSolver,mesh_file,wet_interface):
        """
	    Connection between solvers.
	    Creates the communication support between the two solvers.
	    Gets information about f/s interfaces from the two solvers.
	    """
        if self.have_MPI == True:
            myid = self.comm.Get_rank()
	    MPIsize = self.comm.Get_size()
        else:
            myid = 0
            MPIsize = 1
	# --- Identify the fluid and solid interfaces and store the number of nodes on both sides (and for each partition) ---
        if FluidSolver != None:
            self.haveFluidSolver = True
	    self.fluidInterfaceIdentifier = FluidSolver.GetMovingMarker()
	    self.nLocalFluidInterfaceNodes = FluidSolver.GetNumberVertices(self.fluidInterfaceIdentifier)
	    if self.nLocalFluidInterfaceNodes != 0:
                self.haveFluidInterface = True
	    else:
	        pass
# --- Exchange information about processors on which the solvers are defined and where the interface nodes are lying ---
        if self.have_MPI == True:
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
        else:
            self.fluidSolverProcessors.append(0)
            self.fluidInterfaceProcessors.append(0)
        self.MPIBarrier()
 	# --- Calculate the total number of nodes at the fluid interface (sum over all the partitions) ---
        # Calculate the number of halo nodes on each partition
        self.nLocalFluidInterfaceHaloNode = 0
        for iVertex in range(self.nLocalFluidInterfaceNodes):
            if FluidSolver.IsAHaloNode(self.fluidInterfaceIdentifier, iVertex) == True:
                GlobalIndex = FluidSolver.GetVertexGlobalIndex(self.fluidInterfaceIdentifier, iVertex)
                self.FluidLocalHaloNodeList[GlobalIndex] = iVertex
                self.nLocalFluidInterfaceHaloNode += 1
 # Calculate the number of physical (= not halo) nodes on each partition
            self.nLocalFluidInterfacePhysicalNodes = self.nLocalFluidInterfaceNodes - self.nLocalFluidInterfaceHaloNode
        if self.have_MPI == True:
            self.FluidHaloNodeList = self.comm.allgather(self.FluidLocalHaloNodeList)
        else:
            self.FluidHaloNodeList = self.FluidLocalHaloNodeList
 # --- Calculate the total number of nodes (with and without halo) at the fluid interface (sum over all the partitions) and broadcast the number accross all processors ---
        sendBuffHalo = np.array(int(self.nLocalFluidInterfaceNodes))
        sendBuffPhysical = np.array(int(self.nLocalFluidInterfacePhysicalNodes))
        rcvBuffHalo = np.zeros(1, dtype=int)
        rcvBuffPhysical = np.zeros(1, dtype=int)
        if self.have_MPI == True:
            self.comm.barrier()
            self.comm.Allreduce(sendBuffHalo,rcvBuffHalo,op=self.MPI.SUM)
            self.comm.Allreduce(sendBuffPhysical,rcvBuffPhysical,op=self.MPI.SUM)
            self.nFluidInterfaceNodes = rcvBuffHalo[0]
            self.nFluidInterfacePhysicalNodes = rcvBuffPhysical[0]
        else:
            self.nFluidInterfaceNodes = np.copy(sendBuffHalo)
            self.nFluidInterfacePhysicalNodes = np.copy(sendBuffPhysical)
        del sendBuffHalo, rcvBuffHalo, sendBuffPhysical, rcvBuffPhysical

# --- Store the number of physical interface nodes on each processor and allgather the information ---
        self.fluidPhysicalInterfaceNodesDistribution = np.zeros(MPIsize, dtype=int)
        if self.have_MPI == True:
            sendBuffPhysical = np.array(int(self.nLocalFluidInterfacePhysicalNodes))
            self.comm.Allgather(sendBuffPhysical,self.fluidPhysicalInterfaceNodesDistribution)
            del sendBuffPhysical
        else:
            self.fluidPhysicalInterfaceNodesDistribution[0] = self.nFluidInterfacePhysicalNodes
         # --- Calculate and store the global indexing of interface physical nodes on each processor and allgather the information ---
        if self.have_MPI == True:
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
        else:
            temp = {}
            temp[0] = [0,self.nLocalFluidInterfacePhysicalNodes-1]
            self.fluidGlobalIndexRange = list()
            self.fluidGlobalIndexRange.append(temp)
        self.MPIPrint('Total number of fluid interface nodes (halo nodes included) : {}'.format(self.nFluidInterfaceNodes))
        self.MPIPrint('Total number of fluid interface nodes : {}'.format(self.nFluidInterfacePhysicalNodes))


## --- Get and store the coordinates of the verticies
        FluidLocalInterfaceCoord=[]
        for iVertex in range(self.nLocalFluidInterfaceNodes):
            if FluidSolver.IsAHaloNode(self.fluidInterfaceIdentifier, iVertex) == False:
                FluidLocalInterfaceCoord.append(FluidSolver.GetVertexCoordX(self.fluidInterfaceIdentifier,iVertex))
                FluidLocalInterfaceCoord.append(FluidSolver.GetVertexCoordY(self.fluidInterfaceIdentifier,iVertex))
                FluidLocalInterfaceCoord.append(FluidSolver.GetVertexCoordZ(self.fluidInterfaceIdentifier,iVertex))

        if self.have_MPI == True:
            FluidInterfaceCoordList = self.comm.allgather(FluidLocalInterfaceCoord)
            for i in range (MPIsize):
                for j in FluidInterfaceCoordList[i]:
                    self.FluidInterfaceNodeCoord.append(j)
        else:
            self.FluidInterfaceNodeCoord = FluidLocalInterfaceCoord




        self.FluidInterfaceNodeCoord = np.array (self.FluidInterfaceNodeCoord)
        del FluidLocalInterfaceCoord

## --- Get Global indexes of the nodes (Wet-Interface)
        nodeIDs=[]
        localNodesIDs = []
        for iVertex in range (self.nLocalFluidInterfaceNodes):
            if FluidSolver.IsAHaloNode(self.fluidInterfaceIdentifier, iVertex) == False:
                localNodesIDs.append(int(FluidSolver.GetVertexGlobalIndex(self.fluidInterfaceIdentifier,iVertex)))
        if self.have_MPI == True:
            NodeIDsList = self.comm.allgather(localNodesIDs)
        for i in range (self.MPIsize):
            for j in NodeIDsList[i]:
                nodeIDs.append(j+1)
        self.nodeIDs = np.array (nodeIDs)

        del nodeIDs, localNodesIDs

        ### Get Elements from a mesh file
        self.get_Elements(mesh_file, wet_interface)


# --- set the Initial Values of the coordinates of the mesh (halo is includied)
        self.Initil_particle_displacements(FluidSolver,self.FluidInterfaceNodeCoord)

    def get_Elements (self,mesh_file,wet_interface):


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
                            localElements.append(int(line)+1)


        self.FluidInterfaceElementsNodes = np.array(localElements)
        for i in range (self.nFluidInterfaceElements):
            self.nNodesInElementList.append(self.nNodesInElement)
        self.nNodesInElementList = np.array(self.nNodesInElementList)

    def get_mesh(self,FluidSolver):
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


    def get_force(self,FluidSolver):
        if self.have_MPI == True:
            myid = self.comm.Get_rank()
        else:
            myid = 0
        LocalForce = []
        self.Force = []
# --- Get the fluid interface loads from the fluid solver and directly fill the corresponding PETSc vector ---
        for iVertex in range(self.nLocalFluidInterfaceNodes):
            halo = FluidSolver.ComputeVertexForces(self.fluidInterfaceIdentifier, iVertex) # !!we have to ignore halo node coming from mesh partitioning because they introduice non-physical forces
            if halo==False:
                newFx = FluidSolver.GetVertexForceX(self.fluidInterfaceIdentifier, iVertex)
                newFy = FluidSolver.GetVertexForceY(self.fluidInterfaceIdentifier, iVertex)
                newFz = FluidSolver.GetVertexForceZ(self.fluidInterfaceIdentifier, iVertex)
                LocalForce.append(newFx)
                LocalForce.append(newFy)
                LocalForce.append(newFz)
        if self.have_MPI == True:
            ForceList=self.comm.allgather(LocalForce)
        self.MPIBarrier()
        for i in ForceList:
            for j in i:
                self.Force.append(j)
        self.Force = np.array(self.Force)
        del LocalForce
        return self.Force


    def setFluidDisplacements (self,FluidSolver,newDisplacements):

        if self.have_MPI == True:
    	    myid = self.comm.Get_rank()
        else:
            myid = 0
# --- Send the new fluid interface position to the fluid solver (on each partition, halo nodes included) ---
        localIndex = 0
        self.particle_displacements (FluidSolver, newDisplacements)
        if self.haveFluidInterface:
            for iVertex in range(self.nLocalFluidInterfaceNodes):
                GlobalIndex = FluidSolver.GetVertexGlobalIndex(self.fluidInterfaceIdentifier, iVertex)
                disp = self.searchDisplacement(GlobalIndex,1)
                Initdisp = self.searchDisplacement(GlobalIndex,2)
                #print ('X',GlobalIndex, FluidSolver.GetVertexCoordX(self.fluidInterfaceIdentifier, iVertex), Initdisp[0], disp [0])
                #print ('Y',GlobalIndex, FluidSolver.GetVertexCoordY(self.fluidInterfaceIdentifier, iVertex), Initdisp[1], disp [1])
                #print ('Z',GlobalIndex, FluidSolver.GetVertexCoordZ(self.fluidInterfaceIdentifier, iVertex), Initdisp[2], disp [2])
                FluidSolver.SetVertexCoordX(self.fluidInterfaceIdentifier, iVertex, Initdisp[0] + disp[0])
                FluidSolver.SetVertexCoordY(self.fluidInterfaceIdentifier, iVertex, Initdisp[1] + disp[1])
                FluidSolver.SetVertexCoordZ(self.fluidInterfaceIdentifier, iVertex, Initdisp[2] + disp[2])
                FluidSolver.SetVertexVarCoord(self.fluidInterfaceIdentifier, iVertex)




    def particle_displacements (self,FluidSolver,newDisplacements):
        if self.have_MPI == True:
	        myid = self.comm.Get_rank()
        else:
            myid = 0
        if self.have_MPI == True:
    	        myid = self.comm.Get_rank()
        else:
            myid = 0

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
            for i in range (0,self.nLocalFluidInterfacePhysicalNodes*3,3):
                localDisplacementArrayX.append(newDisplacements[myLocalIndexStart+i])
                localDisplacementArrayY.append(newDisplacements[myLocalIndexStart+i+1])
                localDisplacementArrayZ.append(newDisplacements[myLocalIndexStart+i+2])

        if self.haveFluidInterface:
            myLocalIndex = 0
            for iVertex in range (self.nLocalFluidInterfaceNodes):
                if FluidSolver.IsAHaloNode(self.fluidInterfaceIdentifier, iVertex) == False:
                    GlobalIndex = FluidSolver.GetVertexGlobalIndex(self.fluidInterfaceIdentifier, iVertex)
                    localDisplacementListX[GlobalIndex] =  localDisplacementArrayX[myLocalIndex]
                    localDisplacementListY[GlobalIndex] =  localDisplacementArrayY[myLocalIndex]
                    localDisplacementListZ[GlobalIndex] =  localDisplacementArrayZ[myLocalIndex]
                    myLocalIndex += 1

        if self.have_MPI:
            self.DisplacementListX = self.comm.allgather(localDisplacementListX)
            self.DisplacementListY = self.comm.allgather(localDisplacementListY)
            self.DisplacementListZ = self.comm.allgather(localDisplacementListZ)
        del localDisplacementArrayX, localDisplacementArrayY, localDisplacementArrayZ, localDisplacementListX, localDisplacementListY, localDisplacementListZ

    def Initil_particle_displacements (self,FluidSolver,newDisplacements):
        if self.have_MPI == True:
    	        myid = self.comm.Get_rank()
        else:
            myid = 0
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
                if FluidSolver.IsAHaloNode(self.fluidInterfaceIdentifier, iVertex) == False:
                    GlobalIndex = FluidSolver.GetVertexGlobalIndex(self.fluidInterfaceIdentifier, iVertex)
                    localInitDisplacementListX[GlobalIndex] =  localInitDisplacementArrayX[myLocalIndex]
                    localInitDisplacementListY[GlobalIndex] =  localInitDisplacementArrayY[myLocalIndex]
                    localInitDisplacementListZ[GlobalIndex] =  localInitDisplacementArrayZ[myLocalIndex]
                    myLocalIndex += 1

        if self.have_MPI:
            self.InitCoordListX = self.comm.allgather(localInitDisplacementListX)
            self.InitCoordListY = self.comm.allgather(localInitDisplacementListY)
            self.InitCoordListZ = self.comm.allgather(localInitDisplacementListZ)
        del localInitDisplacementArrayX, localInitDisplacementArrayY, localInitDisplacementArrayZ, localInitDisplacementListX, localInitDisplacementListY, localInitDisplacementListZ

    def searchDisplacement(self,GlobalIndex,flag=1):
        disp = []
        myid = 0
        if flag == 1:
            while True:

                if self.DisplacementListX[myid].get(GlobalIndex) != None :
                    disp.append(self.DisplacementListX[myid].get(GlobalIndex))
                    disp.append(self.DisplacementListY[myid].get(GlobalIndex))
                    disp.append(self.DisplacementListZ[myid].get(GlobalIndex))
                    break
                myid += 1
        if flag == 2:
            while True:

                if self.InitCoordListX[myid].get(GlobalIndex) != None :
                    disp.append(self.InitCoordListX[myid].get(GlobalIndex))
                    disp.append(self.InitCoordListY[myid].get(GlobalIndex))
                    disp.append(self.InitCoordListZ[myid].get(GlobalIndex))
                    break
                myid += 1
        return disp

        # print ('my id is {} and my local Disp X is {} and my physical nodes are {} '. format (
        #     myid, len(self.localDisplacementListX),self.nLocalFluidInterfacePhysicalNodes
        # ))


    def __getGlobalIndex(self, physics, iProc, iLocalVertex):
        """
        Calculate the global indexing of interface nodes accross all the partitions. This does not include halo nodes.
        """

        if physics == 'fluid':
            globalStartIndex = self.fluidGlobalIndexRange[iProc][iProc][0]

        globalIndex = globalStartIndex + iLocalVertex

        return globalIndex

    def getSensitivity (self,FluidSolver):
      LocalSensitivity = []
      self.Sensitivity = []
      maxSens = []
      maxSens1 = []
      if self.have_MPI == True:
        myid = self.comm.Get_rank()
      else:
        myid = 0
      if self.have_MPI == True:
        myid = self.comm.Get_rank()
      else:
        myid = 0
      if self.haveFluidInterface:
        for iVertex in range(self.nLocalFluidInterfaceNodes):
          GlobalIndex = FluidSolver.GetVertexGlobalIndex(self.fluidInterfaceIdentifier, iVertex)
          halo = FluidSolver.ComputeVertexForces(self.fluidInterfaceIdentifier, iVertex) # !!we have to ignore halo node coming from mesh partitioning because they introduice non-physical forces
          if halo==False:
            newSens = FluidSolver.GetSensitivity(self.fluidInterfaceIdentifier, iVertex)
            #Normal Vectors
            Nx = FluidSolver.GetNormalX(self.fluidInterfaceIdentifier, iVertex);
            Ny = FluidSolver.GetNormalY(self.fluidInterfaceIdentifier, iVertex);
            Nz = FluidSolver.GetNormalZ(self.fluidInterfaceIdentifier, iVertex);
            Nx = -1 * Nx;
            Ny = -1 * Ny;
            Nz = -1 * Nz;

            MagnitudeNormal = sqrt(Nx*Nx+Ny*Ny+Nz*Nz)
            newSensX = (newSens/MagnitudeNormal) * (Nx)
            newSensY = (newSens/MagnitudeNormal) * (Ny)
            newSensZ = (newSens/MagnitudeNormal) * (Nz)
            LocalSensitivity.append(newSensX)
            LocalSensitivity.append(newSensY)
            LocalSensitivity.append(newSensZ)
      if self.have_MPI == True:
        SensitivityList=self.comm.allgather(LocalSensitivity)
        maxSensList = self.comm.allgather(maxSens)
        maxSensList1 = self.comm.allgather(maxSens1)
      self.MPIBarrier()
      for i in SensitivityList:
        for j in i:
          self.Sensitivity.append(j)
      self.Sensitivity = np.array(self.Sensitivity)

      del LocalSensitivity, maxSens, maxSens1
      return self.Sensitivity


    def getLBterm (self,FluidSolver):
      LocalSensitivity = []
      self.Sensitivity = []
      maxSens = []
      maxSens1 = []
      if self.have_MPI == True:
        myid = self.comm.Get_rank()
      else:
        myid = 0
      if self.have_MPI == True:
        myid = self.comm.Get_rank()
      else:
        myid = 0
      if self.haveFluidInterface:
        for iVertex in range(self.nLocalFluidInterfaceNodes):
          GlobalIndex = FluidSolver.GetVertexGlobalIndex(self.fluidInterfaceIdentifier, iVertex)
          halo = FluidSolver.ComputeVertexForces(self.fluidInterfaceIdentifier, iVertex) # !!we have to ignore halo node coming from mesh partitioning because they introduice non-physical forces
          if halo==False:

            newSensX = FluidSolver.GetLBtermX(self.fluidInterfaceIdentifier, iVertex)
            newSensY = FluidSolver.GetLBtermY(self.fluidInterfaceIdentifier, iVertex)
            newSensZ = FluidSolver.GetLBtermZ(self.fluidInterfaceIdentifier, iVertex)
            LocalSensitivity.append(newSensX)
            LocalSensitivity.append(newSensY)
            LocalSensitivity.append(newSensZ)
      if self.have_MPI == True:
        SensitivityList=self.comm.allgather(LocalSensitivity)
        maxSensList = self.comm.allgather(maxSens)
        maxSensList1 = self.comm.allgather(maxSens1)
      self.MPIBarrier()
      for i in SensitivityList:
        for j in i:
          self.Sensitivity.append(j)
      self.Sensitivity = np.array(self.Sensitivity)

      del LocalSensitivity, maxSens, maxSens1
      return self.Sensitivity


    def Get_local_sensitivity(self,FluidSolver):
      LocalSensitivity = []
      self.Sensitivity = []
      maxSens = []
      maxSens1 = []
      if self.have_MPI == True:
        myid = self.comm.Get_rank()
      else:
        myid = 0
      if self.have_MPI == True:
        myid = self.comm.Get_rank()
      else:
        myid = 0
      if self.haveFluidInterface:
        for iVertex in range(self.nLocalFluidInterfaceNodes):
          GlobalIndex = FluidSolver.GetVertexGlobalIndex(self.fluidInterfaceIdentifier, iVertex)
          halo = FluidSolver.ComputeVertexForces(self.fluidInterfaceIdentifier, iVertex) # !!we have to ignore halo node coming from mesh partitioning because they introduice non-physical forces
          if halo==False:

            newSensX = FluidSolver.Get_LS_X(self.fluidInterfaceIdentifier, iVertex)
            newSensY = FluidSolver.Get_LS_Y(self.fluidInterfaceIdentifier, iVertex)
            newSensZ = FluidSolver.Get_LS_Z(self.fluidInterfaceIdentifier, iVertex)
            LocalSensitivity.append(newSensX)
            LocalSensitivity.append(newSensY)
            LocalSensitivity.append(newSensZ)
      if self.have_MPI == True:
        SensitivityList=self.comm.allgather(LocalSensitivity)
        maxSensList = self.comm.allgather(maxSens)
        maxSensList1 = self.comm.allgather(maxSens1)
      self.MPIBarrier()
      for i in SensitivityList:
        for j in i:
          self.Sensitivity.append(j)
      self.Sensitivity = np.array(self.Sensitivity)
      del LocalSensitivity, maxSens, maxSens1, SensitivityList
      return self.Sensitivity

    def Get_state_sensitivity(self,FluidSolver):
      LocalSensitivity = []
      self.Sensitivity = []
      maxSens = []
      maxSens1 = []
      if self.have_MPI == True:
        myid = self.comm.Get_rank()
      else:
        myid = 0
      if self.have_MPI == True:
        myid = self.comm.Get_rank()
      else:
        myid = 0
      if self.haveFluidInterface:
        for iVertex in range(self.nLocalFluidInterfaceNodes):
          GlobalIndex = FluidSolver.GetVertexGlobalIndex(self.fluidInterfaceIdentifier, iVertex)
          halo = FluidSolver.ComputeVertexForces(self.fluidInterfaceIdentifier, iVertex) # !!we have to ignore halo node coming from mesh partitioning because they introduice non-physical forces
          if halo==False:

            newSensX = FluidSolver.Get_SS_X(self.fluidInterfaceIdentifier, iVertex)
            newSensY = FluidSolver.Get_SS_Y(self.fluidInterfaceIdentifier, iVertex)
            newSensZ = FluidSolver.Get_SS_Z(self.fluidInterfaceIdentifier, iVertex)
            LocalSensitivity.append(newSensX)
            LocalSensitivity.append(newSensY)
            LocalSensitivity.append(newSensZ)
      if self.have_MPI == True:
        SensitivityList=self.comm.allgather(LocalSensitivity)
        maxSensList = self.comm.allgather(maxSens)
        maxSensList1 = self.comm.allgather(maxSens1)
      self.MPIBarrier()
      for i in SensitivityList:
        for j in i:
          self.Sensitivity.append(j)
      self.Sensitivity = np.array(self.Sensitivity)
      del LocalSensitivity, maxSens, maxSens1, SensitivityList
      return self.Sensitivity




    def getAdjSolution (self,FluidSolver):
      if self.have_MPI == True:
        myid = self.comm.Get_rank()
      else:
        myid = 0
      if self.have_MPI == True:
        myid = self.comm.Get_rank()
      else:
        myid = 0
      if self.haveFluidInterface:
        for iVertex in range(self.nLocalFluidInterfaceNodes):
          GlobalIndex = FluidSolver.GetVertexGlobalIndex(self.fluidInterfaceIdentifier, iVertex)
          halo = FluidSolver.ComputeVertexForces(self.fluidInterfaceIdentifier, iVertex) # !!we have to ignore halo node coming from mesh partitioning because they introduice non-physical forces
          if halo==False:
            PfiX = FluidSolver.GetAdjSolution(self.fluidInterfaceIdentifier, iVertex, 1)
            PfiY = FluidSolver.GetAdjSolution(self.fluidInterfaceIdentifier, iVertex, 2)
            PfiZ = FluidSolver.GetAdjSolution(self.fluidInterfaceIdentifier, iVertex, 3)




    def setNewAdjointDisplacments (self, FluidSolver, newAdjVelocity):
        if self.have_MPI == True:
    	    myid = self.comm.Get_rank()
        else:
            myid = 0
        # --- Send the new fluid interface position to the fluid solver (on each partition, halo nodes included) ---
        localIndex = 0
        self.particle_velocity (FluidSolver, newAdjVelocity)
        if self.haveFluidInterface:
            for iVertex in range(self.nLocalFluidInterfaceNodes):
                GlobalIndex = FluidSolver.GetVertexGlobalIndex(self.fluidInterfaceIdentifier, iVertex)
                Adjdisp = self.searchAdjDispl(GlobalIndex)
                FluidSolver.SetAdjSolutionToFinestMesh(self.fluidInterfaceIdentifier, iVertex, Adjdisp[0], Adjdisp[1], Adjdisp[2])
        FluidSolver.SetAdjSolutionToMultiGrids(self.fluidInterfaceIdentifier)




    def particle_velocity (self,FluidSolver,newAdjVelocity):
        if self.have_MPI == True:
	        myid = self.comm.Get_rank()
        else:
            myid = 0
        if self.have_MPI == True:
    	        myid = self.comm.Get_rank()
        else:
            myid = 0

        self.AdjVelocityListX = []
        self.AdjVelocityListY = []
        self.AdjVelocityListZ = []
        localVelocityArrayX = []
        localVelocityArrayY = []
        localVelocityArrayZ = []

        localVelocityListX = {}
        localVelocityListY = {}
        localVelocityListZ = {}
        if self.haveFluidInterface:
            myLocalIndexStart = self.fluidGlobalIndexRange[myid].get(myid)[0]*3
            for i in range (0,self.nLocalFluidInterfacePhysicalNodes*3,3):
                localVelocityArrayX.append(newAdjVelocity[myLocalIndexStart+i])
                localVelocityArrayY.append(newAdjVelocity[myLocalIndexStart+i+1])
                localVelocityArrayZ.append(newAdjVelocity[myLocalIndexStart+i+2])

        if self.haveFluidInterface:
            myLocalIndex = 0
            for iVertex in range (self.nLocalFluidInterfaceNodes):
                if FluidSolver.IsAHaloNode(self.fluidInterfaceIdentifier, iVertex) == False:
                    GlobalIndex = FluidSolver.GetVertexGlobalIndex(self.fluidInterfaceIdentifier, iVertex)
                    localVelocityListX[GlobalIndex] =  localVelocityArrayX[myLocalIndex]
                    localVelocityListY[GlobalIndex] =  localVelocityArrayY[myLocalIndex]
                    localVelocityListZ[GlobalIndex] =  localVelocityArrayZ[myLocalIndex]
                    myLocalIndex += 1

        if self.have_MPI:
            self.AdjVelocityListX = self.comm.allgather(localVelocityListX)
            self.AdjVelocityListY = self.comm.allgather(localVelocityListY)
            self.AdjVelocityListZ = self.comm.allgather(localVelocityListZ)
        del localVelocityArrayX, localVelocityArrayY, localVelocityArrayZ, localVelocityListX, localVelocityListY, localVelocityListZ

    def searchAdjDispl(self,GlobalIndex):
        disp = []
        myid = 0
        while True:
            if self.AdjVelocityListX[myid].get(GlobalIndex) != None :
                disp.append(self.AdjVelocityListX[myid].get(GlobalIndex))
                disp.append(self.AdjVelocityListY[myid].get(GlobalIndex))
                disp.append(self.AdjVelocityListZ[myid].get(GlobalIndex))
                break
            myid += 1
        return disp
        if self.haveFluidInterface:
            for iVertex in range(self.nLocalFluidInterfaceNodes):
                GlobalIndex = FluidSolver.GetVertexGlobalIndex(self.fluidInterfaceIdentifier, iVertex)
                disp = self.searchDisplacement(GlobalIndex,1)
                Initdisp = self.searchDisplacement(GlobalIndex,2)
                #print ('X',GlobalIndex, FluidSolver.GetVertexCoordX(self.fluidInterfaceIdentifier, iVertex), Initdisp[0], disp [0])
                #print ('Y',GlobalIndex, FluidSolver.GetVertexCoordY(self.fluidInterfaceIdentifier, iVertex), Initdisp[1], disp [1])
                #print ('Z',GlobalIndex, FluidSolver.GetVertexCoordZ(self.fluidInterfaceIdentifier, iVertex), Initdisp[2], disp [2])
                FluidSolver.SetVertexCoordX(self.fluidInterfaceIdentifier, iVertex, Initdisp[0] + disp[0])
                FluidSolver.SetVertexCoordY(self.fluidInterfaceIdentifier, iVertex, Initdisp[1] + disp[1])
                FluidSolver.SetVertexCoordZ(self.fluidInterfaceIdentifier, iVertex, Initdisp[2] + disp[2])
                FluidSolver.SetVertexVarCoord(self.fluidInterfaceIdentifier, iVertex)


    def getSensitivityCoupleLoop (self,FluidSolver):
      LocalSensitivity = []
      Sensitivity = []
      maxSens = []
      maxSens1 = []
      if self.have_MPI == True:
        myid = self.comm.Get_rank()
      else:
        myid = 0
      if self.have_MPI == True:
        myid = self.comm.Get_rank()
      else:
        myid = 0
      if self.haveFluidInterface:
        for iVertex in range(self.nLocalFluidInterfaceNodes):
          GlobalIndex = FluidSolver.GetVertexGlobalIndex(self.fluidInterfaceIdentifier, iVertex)
          halo = FluidSolver.ComputeVertexForces(self.fluidInterfaceIdentifier, iVertex) # !!we have to ignore halo node coming from mesh partitioning because they introduice non-physical forces
          if halo==False:
            newSens = FluidSolver.GetSensitivity(self.fluidInterfaceIdentifier, iVertex)
            #Normal Vectors
            Nx = FluidSolver.GetNormalX(self.fluidInterfaceIdentifier, iVertex);
            Ny = FluidSolver.GetNormalY(self.fluidInterfaceIdentifier, iVertex);
            Nz = FluidSolver.GetNormalZ(self.fluidInterfaceIdentifier, iVertex);
            Nx = -1 * Nx;
            Ny = -1 * Ny;
            Nz = -1 * Nz;

            Adjdisp = self.searchAdjDispl(GlobalIndex)
            #Pressure Gradients
            gradCpX = FluidSolver.GetPressureGradientX(self.fluidInterfaceIdentifier, iVertex)
            gradCpY = FluidSolver.GetPressureGradientY(self.fluidInterfaceIdentifier, iVertex)
            gradCpZ = FluidSolver.GetPressureGradientZ(self.fluidInterfaceIdentifier, iVertex)
            MagnitudeNormal = sqrt(Nx*Nx+Ny*Ny+Nz*Nz)
            DotProdDispGradPressure = gradCpX*(0.998574-Adjdisp[0]) + gradCpY*(-Adjdisp[1]) +gradCpZ*(0.0533817 - Adjdisp[2])
            newSensX = ((newSens/MagnitudeNormal) * (Nx))
            newSensY = ((newSens/MagnitudeNormal) * (Ny))
            newSensZ = ((newSens/MagnitudeNormal) * (Nz))
            LocalSensitivity.append(newSensX)
            LocalSensitivity.append(newSensY)
            LocalSensitivity.append(newSensZ)
      if self.have_MPI == True:
        SensitivityList=self.comm.allgather(LocalSensitivity)
        maxSensList = self.comm.allgather(maxSens)
        maxSensList1 = self.comm.allgather(maxSens1)
      self.MPIBarrier()
      for i in SensitivityList:
        for j in i:
          Sensitivity.append(j)
      Sensitivity = np.array(Sensitivity)

      del LocalSensitivity, maxSens, maxSens1
      return Sensitivity



        # --- apply mesh perturbation on node by its index
    def SetPerturbationComponent(self,FluidSolver, pointIndex, component, perturbation):
        if self.have_MPI == True:
        	    myid = self.comm.Get_rank()
        else:
            myid = 0
# --- Send the new fluid interface position to the fluid solver (on each partition, halo nodes included) ---
        localIndex = 0
        if self.haveFluidInterface:
            for iVertex in range(self.nLocalFluidInterfaceNodes):
                GlobalIndex = FluidSolver.GetVertexGlobalIndex(self.fluidInterfaceIdentifier, iVertex)
                if GlobalIndex == pointIndex:
                    if (component==0):
                        FluidSolver.SetPerturbation(self.fluidInterfaceIdentifier, iVertex, FluidSolver.GetVertexCoordX(self.fluidInterfaceIdentifier, iVertex) + perturbation, FluidSolver.GetVertexCoordY(self.fluidInterfaceIdentifier, iVertex), FluidSolver.GetVertexCoordZ(self.fluidInterfaceIdentifier, iVertex))

                    if (component==1):
                        FluidSolver.SetPerturbation(self.fluidInterfaceIdentifier, iVertex, FluidSolver.GetVertexCoordX(self.fluidInterfaceIdentifier, iVertex), FluidSolver.GetVertexCoordY(self.fluidInterfaceIdentifier, iVertex) + perturbation, FluidSolver.GetVertexCoordZ(self.fluidInterfaceIdentifier, iVertex))

                    if (component==2):
                        FluidSolver.SetPerturbation(self.fluidInterfaceIdentifier, iVertex, FluidSolver.GetVertexCoordX(self.fluidInterfaceIdentifier, iVertex), FluidSolver.GetVertexCoordY(self.fluidInterfaceIdentifier, iVertex), FluidSolver.GetVertexCoordZ(self.fluidInterfaceIdentifier, iVertex) + perturbation)

        FluidSolver.StaticMeshPertirbation()

    def calcNormals(self,globalLocalMap,CurrentNodalCoords):

      #nodalArea = np.zeros (self.mesh[0])
      self.nodalNormal = np.zeros (3*self.mesh[0])

      for i in range (0,len(self.mesh[5]),3):
        nodeGID1 = self.mesh[5][i]
        nodeGID2 = self.mesh[5][i+1]
        nodeGID3 = self.mesh[5][i+2]

        nodeLID1 = globalLocalMap[nodeGID1]
        nodeLID2 = globalLocalMap[nodeGID2]
        nodeLID3 = globalLocalMap[nodeGID3]


        e10 = CurrentNodalCoords[3*nodeLID1  ] - CurrentNodalCoords[3*nodeLID2  ]
        e11 = CurrentNodalCoords[3*nodeLID1+1] - CurrentNodalCoords[3*nodeLID2+1]
        e12 = CurrentNodalCoords[3*nodeLID1+2] - CurrentNodalCoords[3*nodeLID2+2]

        e20 = CurrentNodalCoords[3*nodeLID3  ] - CurrentNodalCoords[3*nodeLID2  ]
        e21 = CurrentNodalCoords[3*nodeLID3+1] - CurrentNodalCoords[3*nodeLID2+1]
        e22 = CurrentNodalCoords[3*nodeLID3+2] - CurrentNodalCoords[3*nodeLID2+2]

        normalx = e11 * e22 - e21 * e12;
        normaly = e12 * e20 - e22 * e10;
        normalz = e10 * e21 - e20 * e11;

        norm2 = normalx*normalx+normaly*normaly+normalz*normalz;

        norm2 = sqrt( norm2 )

        #nodalArea[nodeLID1] = nodalArea[nodeLID1] + norm2/6.0
        #nodalArea[nodeLID2] = nodalArea[nodeLID2] + norm2/6.0
        #nodalArea[nodeLID3] = nodalArea[nodeLID3] + norm2/6.0

        self.nodalNormal[3*nodeLID1+0] += normalx/6
        self.nodalNormal[3*nodeLID1+1] += normaly/6
        self.nodalNormal[3*nodeLID1+2] += normalz/6

        self.nodalNormal[3*nodeLID2+0] += normalx/6
        self.nodalNormal[3*nodeLID2+1] += normaly/6
        self.nodalNormal[3*nodeLID2+2] += normalz/6

        self.nodalNormal[3*nodeLID3+0] += normalx/6
        self.nodalNormal[3*nodeLID3+1] += normaly/6
        self.nodalNormal[3*nodeLID3+2] += normalz/6


      return  self.nodalNormal

      #for i in range (len(nodalArea)):
          #Nx = nodalNormal[3*i]
          #Ny = nodalNormal[3*i+1]
          #Nz = nodalNormal[3*i+2]
          ##MagnitudeNormal = sqrt(Nx*Nx+Ny*Ny+Nz*Nz)
          ##Nx = Nx/MagnitudeNormal
          ##Ny = Ny/MagnitudeNormal
          ##Nz = Nz/MagnitudeNormal
          #print("i = ",i, "SU2 x are = ",SU2normals[3*i],"myArea x = ",Nx)
          #print("i = ",i, "SU2 y are = ",SU2normals[3*i+1],"myArea y = ",Ny)
          #print("i = ",i, "SU2 z are = ",SU2normals[3*i+2],"myArea z = ",Nz)
          #print("***********************")


    def getNodalAreaNormal (self,FluidSolver):
      LocalNormal = []
      self.Normal = []
      if self.have_MPI == True:
        myid = self.comm.Get_rank()
      else:
        myid = 0
      if self.have_MPI == True:
        myid = self.comm.Get_rank()
      else:
        myid = 0
      if self.haveFluidInterface:
        for iVertex in range(self.nLocalFluidInterfaceNodes):
          GlobalIndex = FluidSolver.GetVertexGlobalIndex(self.fluidInterfaceIdentifier, iVertex)
          halo = FluidSolver.ComputeVertexForces(self.fluidInterfaceIdentifier, iVertex) # !!we have to ignore halo node coming from mesh partitioning because they introduice non-physical sensitivities
          if halo==False:
            #Normal Vectors
            Nx = FluidSolver.GetNormalX(self.fluidInterfaceIdentifier, iVertex);
            Ny = FluidSolver.GetNormalY(self.fluidInterfaceIdentifier, iVertex);
            Nz = FluidSolver.GetNormalZ(self.fluidInterfaceIdentifier, iVertex);
            LocalNormal.append(-Nx)
            LocalNormal.append(-Ny)
            LocalNormal.append(-Nz)
      if self.have_MPI == True:
        NormalList=self.comm.allgather(LocalNormal)
      self.MPIBarrier()
      for i in NormalList:
        for j in i:
          self.Normal.append(j)
      del LocalNormal
      return np.array(self.Normal)

    def getCurrentNodalCoords (self,FluidSolver):
      LocalCoords = []
      self.CurrCoords = []
      if self.have_MPI == True:
        myid = self.comm.Get_rank()
      else:
        myid = 0
      if self.have_MPI == True:
        myid = self.comm.Get_rank()
      else:
        myid = 0
      if self.haveFluidInterface:
        for iVertex in range(self.nLocalFluidInterfaceNodes):
          GlobalIndex = FluidSolver.GetVertexGlobalIndex(self.fluidInterfaceIdentifier, iVertex)
          halo = FluidSolver.ComputeVertexForces(self.fluidInterfaceIdentifier, iVertex) # !!we have to ignore halo node coming from mesh partitioning because they introduice non-physical sensitivities
          if halo==False:
            #Normal Vectors
            x = FluidSolver.GetVertexCoordX(self.fluidInterfaceIdentifier, iVertex);
            y = FluidSolver.GetVertexCoordY(self.fluidInterfaceIdentifier, iVertex);
            z = FluidSolver.GetVertexCoordZ(self.fluidInterfaceIdentifier, iVertex);
            LocalCoords.append(x)
            LocalCoords.append(y)
            LocalCoords.append(z)
      if self.have_MPI == True:
        NormalList=self.comm.allgather(LocalCoords)
      self.MPIBarrier()
      for i in NormalList:
        for j in i:
          self.CurrCoords.append(j)
      del LocalCoords
      return np.array(self.CurrCoords)

    def calcPsiGardForce(self,globalLocalMap,SU2normals,CurrentNodalCoords,pressureForce,adjDisp):

        #check wether we have only triangles
        for i in range (len(self.mesh[4])):
            if (self.mesh[4][i]!=3):
                print("not valid element type")
                exit()


        self.result = np.zeros (3*self.mesh[0])
        nodalPressure = np.zeros (self.mesh[0])

        #calculate Pressure out of force
        for i in range (len(nodalPressure)):
            nodalPressure[i] = pressureForce[3*i] * SU2normals[3*i] + pressureForce[3*i+1] * SU2normals[3*i+1] + pressureForce[3*i+2] * SU2normals[3*i+2]
            norm2 = sqrt( SU2normals[3*i] * SU2normals[3*i] + SU2normals[3*i+1] * SU2normals[3*i+1] + SU2normals[3*i+2] * SU2normals[3*i+2] )
            nodalPressure[i] /=(norm2*norm2)

        #localPressForceTest = np.zeros (3*self.mesh[0])


        for i in range (0,len(self.mesh[5]),3):
            nodeGID1 = self.mesh[5][i]
            nodeGID2 = self.mesh[5][i+1]
            nodeGID3 = self.mesh[5][i+2]

            nodeLID1 = globalLocalMap[nodeGID1]
            nodeLID2 = globalLocalMap[nodeGID2]
            nodeLID3 = globalLocalMap[nodeGID3]


            e10 = CurrentNodalCoords[3*nodeLID1+0] - CurrentNodalCoords[3*nodeLID2+0]
            e11 = CurrentNodalCoords[3*nodeLID1+1] - CurrentNodalCoords[3*nodeLID2+1]
            e12 = CurrentNodalCoords[3*nodeLID1+2] - CurrentNodalCoords[3*nodeLID2+2]

            e20 = CurrentNodalCoords[3*nodeLID3+0] - CurrentNodalCoords[3*nodeLID2+0]
            e21 = CurrentNodalCoords[3*nodeLID3+1] - CurrentNodalCoords[3*nodeLID2+1]
            e22 = CurrentNodalCoords[3*nodeLID3+2] - CurrentNodalCoords[3*nodeLID2+2]

            areaNormalX = (e11 * e22 - e21 * e12)/6;
            areaNormalY = (e12 * e20 - e22 * e10)/6;
            areaNormalZ = (e10 * e21 - e20 * e11)/6;

            #localPressForceTest[3*nodeLID1+0] += areaNormalX * nodalPressure[nodeLID1]
            #localPressForceTest[3*nodeLID1+1] += areaNormalY * nodalPressure[nodeLID1]
            #localPressForceTest[3*nodeLID1+2] += areaNormalZ * nodalPressure[nodeLID1]

            #localPressForceTest[3*nodeLID2+0] += areaNormalX * nodalPressure[nodeLID2]
            #localPressForceTest[3*nodeLID2+1] += areaNormalY * nodalPressure[nodeLID2]
            #localPressForceTest[3*nodeLID2+2] += areaNormalZ * nodalPressure[nodeLID2]

            #localPressForceTest[3*nodeLID3+0] += areaNormalX * nodalPressure[nodeLID3]
            #localPressForceTest[3*nodeLID3+1] += areaNormalY * nodalPressure[nodeLID3]
            #localPressForceTest[3*nodeLID3+2] += areaNormalZ * nodalPressure[nodeLID3]


            self.result[3*nodeLID1+0] +=  (-e22/6) * (nodalPressure[nodeLID1] * adjDisp[3*nodeLID1+1]) + \
                                     (-e22/6) * (nodalPressure[nodeLID2] * adjDisp[3*nodeLID2+1]) + \
                                     (-e22/6) * (nodalPressure[nodeLID3] * adjDisp[3*nodeLID3+1]) + \
                                     ( e21/6) * (nodalPressure[nodeLID1] * adjDisp[3*nodeLID1+2]) + \
                                     ( e21/6) * (nodalPressure[nodeLID2] * adjDisp[3*nodeLID2+2]) + \
                                     ( e21/6) * (nodalPressure[nodeLID3] * adjDisp[3*nodeLID3+2])


            self.result[3*nodeLID1+1] +=  ( e22/6) * (nodalPressure[nodeLID1] * adjDisp[3*nodeLID1])   + \
                                     ( e22/6) * (nodalPressure[nodeLID2] * adjDisp[3*nodeLID2])   + \
                                     ( e22/6) * (nodalPressure[nodeLID3] * adjDisp[3*nodeLID3])   + \
                                     (-e20/6) * (nodalPressure[nodeLID1] * adjDisp[3*nodeLID1+2]) + \
                                     (-e20/6) * (nodalPressure[nodeLID2] * adjDisp[3*nodeLID2+2]) + \
                                     (-e20/6) * (nodalPressure[nodeLID3] * adjDisp[3*nodeLID3+2])


            self.result[3*nodeLID1+2] +=  (-e21/6) * (nodalPressure[nodeLID1] * adjDisp[3*nodeLID1])   + \
                                     (-e21/6) * (nodalPressure[nodeLID2] * adjDisp[3*nodeLID2])   + \
                                     (-e21/6) * (nodalPressure[nodeLID3] * adjDisp[3*nodeLID3])   + \
                                     ( e20/6) * (nodalPressure[nodeLID1] * adjDisp[3*nodeLID1+1]) + \
                                     ( e20/6) * (nodalPressure[nodeLID2] * adjDisp[3*nodeLID2+1]) + \
                                     ( e20/6) * (nodalPressure[nodeLID3] * adjDisp[3*nodeLID3+1])


            self.result[3*nodeLID2+0] +=  ( e22/6) * (nodalPressure[nodeLID1] * adjDisp[3*nodeLID1+1]) + \
                                     ( e22/6) * (nodalPressure[nodeLID2] * adjDisp[3*nodeLID2+1]) + \
                                     ( e22/6) * (nodalPressure[nodeLID3] * adjDisp[3*nodeLID3+1]) + \
                                     (-e21/6) * (nodalPressure[nodeLID1] * adjDisp[3*nodeLID1+2]) + \
                                     (-e21/6) * (nodalPressure[nodeLID2] * adjDisp[3*nodeLID2+2]) + \
                                     (-e21/6) * (nodalPressure[nodeLID3] * adjDisp[3*nodeLID3+2])


            self.result[3*nodeLID2+0] +=  (-e12/6) * (nodalPressure[nodeLID1] * adjDisp[3*nodeLID1+1]) + \
                                     (-e12/6) * (nodalPressure[nodeLID2] * adjDisp[3*nodeLID2+1]) + \
                                     (-e12/6) * (nodalPressure[nodeLID3] * adjDisp[3*nodeLID3+1]) + \
                                     ( e11/6) * (nodalPressure[nodeLID1] * adjDisp[3*nodeLID1+2]) + \
                                     ( e11/6) * (nodalPressure[nodeLID2] * adjDisp[3*nodeLID2+2]) + \
                                     ( e11/6) * (nodalPressure[nodeLID3] * adjDisp[3*nodeLID3+2])



            self.result[3*nodeLID2+1] +=  (-e22/6) * (nodalPressure[nodeLID1] * adjDisp[3*nodeLID1])   + \
                                     (-e22/6) * (nodalPressure[nodeLID2] * adjDisp[3*nodeLID2])   + \
                                     (-e22/6) * (nodalPressure[nodeLID3] * adjDisp[3*nodeLID3])   + \
                                     ( e20/6) * (nodalPressure[nodeLID1] * adjDisp[3*nodeLID1+2]) + \
                                     ( e20/6) * (nodalPressure[nodeLID2] * adjDisp[3*nodeLID2+2]) + \
                                     ( e20/6) * (nodalPressure[nodeLID3] * adjDisp[3*nodeLID3+2])


            self.result[3*nodeLID2+1] +=  ( e12/6) * (nodalPressure[nodeLID1] * adjDisp[3*nodeLID1])   + \
                                     ( e12/6) * (nodalPressure[nodeLID2] * adjDisp[3*nodeLID2])   + \
                                     ( e12/6) * (nodalPressure[nodeLID3] * adjDisp[3*nodeLID3])   + \
                                     (-e10/6) * (nodalPressure[nodeLID1] * adjDisp[3*nodeLID1+2]) + \
                                     (-e10/6) * (nodalPressure[nodeLID2] * adjDisp[3*nodeLID2+2]) + \
                                     (-e10/6) * (nodalPressure[nodeLID3] * adjDisp[3*nodeLID3+2])


            self.result[3*nodeLID2+2] +=  ( e21/6) * (nodalPressure[nodeLID1] * adjDisp[3*nodeLID1])   + \
                                     ( e21/6) * (nodalPressure[nodeLID2] * adjDisp[3*nodeLID2])   + \
                                     ( e21/6) * (nodalPressure[nodeLID3] * adjDisp[3*nodeLID3])   + \
                                     (-e20/6) * (nodalPressure[nodeLID1] * adjDisp[3*nodeLID1+1]) + \
                                     (-e20/6) * (nodalPressure[nodeLID2] * adjDisp[3*nodeLID2+1]) + \
                                     (-e20/6) * (nodalPressure[nodeLID3] * adjDisp[3*nodeLID3+1])


            self.result[3*nodeLID2+2] +=  (-e11/6) * (nodalPressure[nodeLID1] * adjDisp[3*nodeLID1])   + \
                                     (-e11/6) * (nodalPressure[nodeLID2] * adjDisp[3*nodeLID2])   + \
                                     (-e11/6) * (nodalPressure[nodeLID3] * adjDisp[3*nodeLID3])   + \
                                     ( e10/6) * (nodalPressure[nodeLID1] * adjDisp[3*nodeLID1+1]) + \
                                     ( e10/6) * (nodalPressure[nodeLID2] * adjDisp[3*nodeLID2+1]) + \
                                     ( e10/6) * (nodalPressure[nodeLID3] * adjDisp[3*nodeLID3+1])


            self.result[3*nodeLID3+0] +=  ( e12/6) * (nodalPressure[nodeLID1] * adjDisp[3*nodeLID1+1]) + \
                                     ( e12/6) * (nodalPressure[nodeLID2] * adjDisp[3*nodeLID2+1]) + \
                                     ( e12/6) * (nodalPressure[nodeLID3] * adjDisp[3*nodeLID3+1]) + \
                                     (-e11/6) * (nodalPressure[nodeLID1] * adjDisp[3*nodeLID1+2]) + \
                                     (-e11/6) * (nodalPressure[nodeLID2] * adjDisp[3*nodeLID2+2]) + \
                                     (-e11/6) * (nodalPressure[nodeLID3] * adjDisp[3*nodeLID3+2])


            self.result[3*nodeLID3+1] +=  (-e12/6) * (nodalPressure[nodeLID1] * adjDisp[3*nodeLID1])   + \
                                     (-e12/6) * (nodalPressure[nodeLID2] * adjDisp[3*nodeLID2])   + \
                                     (-e12/6) * (nodalPressure[nodeLID3] * adjDisp[3*nodeLID3])   + \
                                     ( e10/6) * (nodalPressure[nodeLID1] * adjDisp[3*nodeLID1+2]) + \
                                     ( e10/6) * (nodalPressure[nodeLID2] * adjDisp[3*nodeLID2+2]) + \
                                     ( e10/6) * (nodalPressure[nodeLID3] * adjDisp[3*nodeLID3+2])


            self.result[3*nodeLID3+2] +=  ( e11/6) * (nodalPressure[nodeLID1] * adjDisp[3*nodeLID1])   + \
                                     ( e11/6) * (nodalPressure[nodeLID2] * adjDisp[3*nodeLID2])   + \
                                     ( e11/6) * (nodalPressure[nodeLID3] * adjDisp[3*nodeLID3])   + \
                                     (-e10/6) * (nodalPressure[nodeLID1] * adjDisp[3*nodeLID1+1]) + \
                                     (-e10/6) * (nodalPressure[nodeLID2] * adjDisp[3*nodeLID2+1]) + \
                                     (-e10/6) * (nodalPressure[nodeLID3] * adjDisp[3*nodeLID3+1])


        #for i in range (len(localPressForceTest)/3):
            #print("i = ",i, "SU2 x are = ",pressureForce[3*i],"myArea x = ",localPressForceTest[3*i])
            #print("i = ",i, "SU2 y are = ",pressureForce[3*i+1],"myArea y = ",localPressForceTest[3*i+1])
            #print("i = ",i, "SU2 z are = ",pressureForce[3*i+2],"myArea z = ",localPressForceTest[3*i+2])
            #print("***********************")

        del nodalPressure
        return self.result







































































































