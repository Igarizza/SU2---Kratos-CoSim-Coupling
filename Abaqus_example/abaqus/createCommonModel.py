# -*- coding: mbcs -*-
# Do not delete the following import lines
import os
import shutil

from abaqus import *
from abaqusConstants import *
import section
import regionToolset
import displayGroupMdbToolset as dgm
import part
import material
import assembly
import step
import interaction
import load
import mesh
import optimization
import job
import sketch
import visualization
import xyPlot
import displayGroupOdbToolset as dgo
import connectorBehavior

def createPart(numberOfChannels, xChannel, yChannel, zChannel, xSide, xBridge, yBelow, yAbove):
    
    overallX = 2*xSide + numberOfChannels*xChannel + (numberOfChannels-1)*xBridge
    overallY = yBelow + yChannel + yAbove

    s = mdb.models['Model-1'].ConstrainedSketch(name='__profile__', 
        sheetSize=200.0)
    g, v, d, c = s.geometry, s.vertices, s.dimensions, s.constraints
    s.setPrimaryObject(option=STANDALONE)
    s.rectangle(point1=(0, 0), point2=(overallX, overallY))
    for i in range(numberOfChannels):
        s.rectangle(point1=(xSide + i*(xChannel + xBridge), yBelow),
                    point2=(xSide + i*(xChannel + xBridge) + xChannel, yBelow + yChannel))
    p = mdb.models['Model-1'].Part(name='Panel', dimensionality=THREE_D, 
        type=DEFORMABLE_BODY)
    p.BaseSolidExtrude(sketch=s, depth=zChannel)
    s.unsetPrimaryObject()
    session.viewports['Viewport: 1'].setValues(displayedObject=p)
    del mdb.models['Model-1'].sketches['__profile__']

def createFillet(filletRadius, numberOfChannels, xChannel, yChannel, zChannel, xSide, xBridge, yBelow):
    if filletRadius <= 0.0:
        pass
    else:
        p = mdb.models['Model-1'].parts['Panel']
        e = p.edges

        coordList = []
        for i in range(numberOfChannels):
            coordList.append([xSide + i*(xChannel + xBridge), yBelow + yChannel, zChannel/2])
            coordList.append([xSide + i*(xChannel + xBridge), yBelow, zChannel/2])
            coordList.append([xSide + i*(xChannel + xBridge) + xChannel , yBelow + yChannel, zChannel/2])
            coordList.append([xSide + i*(xChannel + xBridge) + xChannel, yBelow, zChannel/2])
     
        edgeList = []
        for i in range(numberOfChannels*4):
            edgeList.append(e.findAt(coordinates=(coordList[i])))
        edgeTuple = tuple(edgeList)
        p.Round(radius=min(filletRadius, min(xChannel/2,yChannel/2)), edgeList=(edgeTuple))

def createSurfacesAndSets(numberOfChannels, xChannel, yChannel, zChannel, xSide, xBridge, yBelow, yAbove):
    boundaryList = []
    p = mdb.models['Model-1'].parts['Panel']
    f = p.faces
    faces = f.findAt(((xSide, yBelow + yChannel + yAbove, zChannel/2),))
    p.Set(faces=faces, name='External')
    p.Surface(side1Faces=faces, name='Surface_External')
    boundaryList.append('Surface_External')
    for i in range(numberOfChannels):
        faces = f.getByBoundingBox(xMin=xSide + i*(xChannel + xBridge), yMin=yBelow, zMin=0,
                                   xMax=xSide + xChannel + i*(xChannel + xBridge), yMax=yBelow+yChannel, zMax=zChannel)
        p.Set(faces=faces, name='Channel' + str(i))
        p.Surface(side1Faces=faces, name='Surface_Channel' + str(i))
        boundaryList.append('Surface_Channel' + str(i))

    return boundaryList

def meshPart(generateMesh, meshSize, meshDeviationFactor, meshMinSizeFactor):
    if generateMesh:
        p = mdb.models['Model-1'].parts['Panel']

        cells = p.cells.getSequenceFromMask(mask=('[#1 ]', ), )
        p.setMeshControls(regions=cells, elemShape=TET, technique=FREE)
        
        print("\nMeshing...\n")
        cellList =(p.cells.getSequenceFromMask(mask=('[#1 ]', ), ), )
        elemType1 = mesh.ElemType(elemCode=DC3D20, elemLibrary=STANDARD)
        elemType2 = mesh.ElemType(elemCode=DC3D15, elemLibrary=STANDARD)
        elemType3 = mesh.ElemType(elemCode=DC3D10, elemLibrary=STANDARD)
        p.setElementType(regions=cellList, elemTypes=(elemType1, elemType2, elemType3))

        p.seedPart(size=meshSize, deviationFactor=meshDeviationFactor, minSizeFactor=meshMinSizeFactor)
        p.generateMesh()
    else:
        pass

def createAssembly():
    a = mdb.models['Model-1'].rootAssembly
    a.DatumCsysByDefault(CARTESIAN)
    p = mdb.models['Model-1'].parts['Panel']
    a.Instance(name='Structure_Instance', part=p, dependent=ON)
    a.rotate(instanceList=('Structure_Instance', ), axisPoint=(0.0, 0.0, 0.0), 
        axisDirection=(1.0, 0.0, 0.0), angle=90.0)

def exportInputFile(generateInputFile):
    if generateInputFile:
        mdb.Job(name='setupCase', model='Model-1', description='', type=ANALYSIS, 
            atTime=None, waitMinutes=0, waitHours=0, queue=None, memory=90, 
            memoryUnits=PERCENTAGE, getMemoryFromAnalysis=True, 
            explicitPrecision=SINGLE, nodalOutputPrecision=SINGLE, echoPrint=OFF, 
            modelPrint=OFF, contactPrint=OFF, historyPrint=OFF, userSubroutine='', 
            scratch='', resultsFormat=ODB, numDomains=1, 
            activateLoadBalancing=False, numThreadsPerMpiProcess=1, 
            multiprocessingMode=DEFAULT, numCpus=1, numGPUs=0)
        mdb.jobs['setupCase'].writeInput(consistencyChecking=OFF)

    else:
        pass

def createBoundaryInputFile(createSurfaceCouplingFiles,couplingBoundaryList):
    if createSurfaceCouplingFiles == False:
        print("BoundaryInputFile was NOT created.")
        return
    mdb.models['Model-1'].HeatTransferStep(name='temporaryStep', previous='Initial', 
        deltmx=10000.0)
    mdb.models['Model-1'].ExpressionField(name='temporaryExpressionField', localCsys=None, 
        description='', expression='X+Y+Z')
    a = mdb.models['Model-1'].rootAssembly
    for i in range(len(couplingBoundaryList)):
        region = a.instances['Structure_Instance'].surfaces[couplingBoundaryList[i]]
        mdb.models['Model-1'].SurfaceHeatFlux(name='temporaryLoad_' + str(i), createStepName='temporaryStep', 
            region=region, magnitude=1.0, distributionType=FIELD, field='temporaryExpressionField')
        
    mdb.Job(name='withDummyLoads', model='Model-1', description='', type=ANALYSIS, 
        atTime=None, waitMinutes=0, waitHours=0, queue=None, memory=90, 
        memoryUnits=PERCENTAGE, getMemoryFromAnalysis=True, 
        explicitPrecision=SINGLE, nodalOutputPrecision=SINGLE, echoPrint=OFF, 
        modelPrint=OFF, contactPrint=OFF, historyPrint=OFF, userSubroutine='', 
        scratch='', resultsFormat=ODB, numThreadsPerMpiProcess=1, 
        multiprocessingMode=DEFAULT, numCpus=1, numGPUs=0)
    mdb.jobs['withDummyLoads'].writeInput(consistencyChecking=OFF)
    del mdb.jobs['setupCase']
    del mdb.jobs['withDummyLoads']
    del mdb.models['Model-1'].steps['temporaryStep']
    del mdb.models['Model-1'].analyticalFields['temporaryExpressionField']


path = os.getcwd() + "/abaqus/"
os.chdir(path)
generateMesh = True
generateInputFile = True # Can only be True if generateMesh is True
createSurfaceCouplingFiles = True

numberOfChannels = 3
xSide = 0.02
xChannel = 0.04
xBridge = 0.01
yBelow = 0.02
yChannel = 0.04
yAbove = 0.02
zChannel = 0.32
filletRadius = 0.0

meshSize = 0.005
meshDeviationFactor = 0.1
meshMinSizeFactor = 0.1

mdb.models['Model-1'].setValues(absoluteZero=0, stefanBoltzmann=5.670374419E-8, universalGas=8.314) #Temperature in Kelvin!

if generateMesh == False and generateInputFile == True:
    generateInputFile = False
    print("Input File can not be written without a mesh!\n")
    print("Enable the generateMesh-Flag first!\n")

createPart(numberOfChannels, xChannel, yChannel, zChannel, xSide, xBridge, yBelow, yAbove)
createFillet(filletRadius, numberOfChannels, xChannel, yChannel, zChannel, xSide, xBridge, yBelow)
boundaryList = createSurfacesAndSets(numberOfChannels, xChannel, yChannel, zChannel, xSide, xBridge, yBelow, yAbove)
meshPart(generateMesh, meshSize, meshDeviationFactor, meshMinSizeFactor)
createAssembly()
exportInputFile(generateInputFile)
createBoundaryInputFile(createSurfaceCouplingFiles,boundaryList)
