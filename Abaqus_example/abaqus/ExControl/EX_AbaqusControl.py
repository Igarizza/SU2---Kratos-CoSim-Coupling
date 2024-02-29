import KratosMultiphysics as Kratos
import numpy as np
import json
import os
import re
import time


def abqLog(msg):
    print("[ABAQUS:] " + msg)

def errorLog(msg):
    print("[ERROR:] " + msg)

BCString = ""

class EX_AbaqusControl():

    def checkSurfaceID(surfaceID):
        with open("abaqus/meshFolder/couplingPatchNames.csv",'r') as f:
            patchNames = f.read().splitlines()
        if surfaceID in patchNames:
            return
        else:
            errorLog("Nodeset \'" + surfaceID + "\' does not exist.")
            exit() 

    def getNodesOfSurface(surfaceID):
        EX_AbaqusControl.checkSurfaceID(surfaceID)
        with open("abaqus/meshFolder/nodes_" + surfaceID + ".csv",'r') as f:
            nodeIDs = f.read().splitlines()
        return nodeIDs

    def getSurfaceIDs():
        with open("abaqus/meshFolder/couplingPatchNames.csv",'r') as f:
            patchNames = f.read().splitlines()
        return patchNames
    
    def createMesh(script_file):
        abqLog("Generating Mesh...")
        os.system("abaqus cae noGUI=" + script_file)
        EX_AbaqusControl.generateMeshOutputFiles()

    def generateMeshOutputFiles():
        abqLog("Reading Mesh...")
        with open("abaqus/withDummyLoads.inp",'r') as f:
            lines = f.readlines()

        counter = 0
        flagString = ''
        nodes = []
        elements = []
        patchNames = []
        nodeSets = []
        elementSets = []
        surfaceSets = []
        sortedList = []
        surfaceCounter = 0

        for line in lines:
            counter += 1
            line = line.strip()
            if(line[0:5].upper()=='*NODE'):
                flagString = 'nodes'
                continue
            if(line[0:8].upper()=='*ELEMENT'):
                elTyp = line.split('=')[-1]
                flagString = 'elements'
                continue
            if(line[0:6].upper()=='*ELSET'):
                flagString = ''
                continue
            if(line[0:5].upper()=='*NSET'):
                if '_PickedSet' in line:
                    flagString = ''
                    continue
                flagString = 'nodeSets'
                patchNames.append(line.split("=")[-1])
                nodeSets.append([])
                continue
            if(line[0:6].upper()=='*DFLUX'):
                flagString = 'elementSets'
                surfaceName = patchNames[surfaceCounter]
                surfaceCounter += 1
                elementSets.append([])
                continue
            if(line[0:2].upper()=='**'):
                flagString = ''
                continue
            if(line[0:8].upper()=='*SURFACE'):
                flagString = ''
                continue

            if flagString == 'nodes':
                nodes.append(line.replace(" ","") + "\n")
            if flagString == 'elements':
                elements.append(line.replace(" ","") + "," + elTyp + "\n")
            if flagString == 'nodeSets':
                lineList = line.replace(" ","").split(",")
                for i in range(len(lineList)):
                    nodeSets[-1].append(lineList[i])
            if flagString == 'elementSets':
                lineList = line.replace(" ","").split(",")
                elementSets[surfaceCounter-1].append(lineList[0].split('.')[-1] + "," + lineList[1] + "\n")


        ### Output of CSV-Files
        ## Creating the meshFolder
                
        os.system("mkdir abaqus/meshFolder")
        os.system("mkdir abaqus/resultFolder")
        ## Output of all Nodes and all elements
        
        ## Output all the nodes
        with open("abaqus/meshFolder/nodes.csv", 'w') as f:
            f.writelines(nodes)
        
        ## Output all the elements
        with open("abaqus/meshFolder/elements.csv", 'w') as f:
            f.writelines(elements)

        ## Output all the coupling nodes
        for i in range(len(patchNames)):
            with open("abaqus/meshFolder/nodes_" + patchNames[i] + ".csv", 'w') as f:
                for j in range(len(nodeSets[i])):
                    f.writelines(str(nodeSets[i][j]) + "\n")


        ## Output all the coupling elements and respective coupling face
        for i in range(len(patchNames)):
            with open("abaqus/meshFolder/elements_" + patchNames[i] + ".csv", 'w') as f:
                sortedList = sorted(elementSets[i], key=lambda s: int(re.search(r'\d+', s).group()))
                for j in range(len(elementSets[i])):
                    f.writelines(sortedList[j])

        ## Output all the coupling elements and nodes of the respective coupling face
        with open("abaqus/meshFolder/elements.csv",'r') as f:
            elementList = f.readlines()
            np.array(elementList)
        with open("abaqus/ExControl/elementTypes.json",'r') as f:
            elementDict = json.load(f)
        for i in range(len(patchNames)):
            with open("abaqus/meshFolder/elements_" + patchNames[i] + ".csv",'r') as f:
                elementPatchList = f.readlines()
            with open("abaqus/meshFolder/elements_" + patchNames[i] + "_with_nodes.csv",'w') as f:
                for j in range(len(elementPatchList)):
                    elementType = elementList[int(elementPatchList[j].split(",")[0])-1].split(",")[-1].strip()
                    for k in range(len(elementDict.keys())):
                        if elementType not in elementDict[list(elementDict.keys())[k]]["ElementTypes"]:
                            pass
                        else:
                            lineString = ""
                            lineString += elementPatchList[j].split(",")[0]
                            nodeAssignmentList = elementDict[list(elementDict.keys())[k]]["FaceAssignment"][0][elementPatchList[j].split(",")[-1].split("\n")[0]]
                            for l in range(len(nodeAssignmentList)):
                                node = elementList[int(elementPatchList[j].split(",")[0])-1].split(",")[nodeAssignmentList[l]+1]
                                lineString += "," + node
                            f.write(lineString + "\n")

        ## Output all the coupling nodes with coordinates
        for i in range(len(patchNames)):
            lastRun = False
            lineCounter = 0
            outputList = []
            for j in range(len(nodes)):
                if nodeSets[i][lineCounter] == nodeSets[i][-1]:
                    lastRun = True
                if nodes[j].split(",")[0] == nodeSets[i][lineCounter]:
                    outputList.append(nodes[j])
                    lineCounter += 1
                    if lastRun:
                        break
            with open("abaqus/meshFolder/nodes_" + patchNames[i] + "_with_Coords.csv", "w") as f:
                f.writelines(outputList)

        ## Output the boundary names
        with open("abaqus/meshFolder/couplingPatchNames.csv","w") as f:
            for patchName in patchNames:
                f.writelines(patchName + '\n')

        ## Checking the files
        files = os.listdir("abaqus/meshFolder")
        for fileName in files:
            if "couplingPatchNames.csv" == fileName:
                pass
            else:
                with open("abaqus/meshFolder/" + fileName,'r') as f:
                    lines = f.readlines()
                for j in range(len(lines)-1):
                    if int(lines[j].split(",")[0]) > int(lines[j+1].split(",")[0]):
                        errorLog("The nodes/elements of file \'" + fileName + "\' are not in ascending order!")
                        exit()
                    else:
                        pass

        #os.system("mv setupCase.inp abaqus/")
        #os.system("mv withDummyLoads.inp abaqus/")

    def setMaterial(materialName):
        abqLog("Setting material...")
        with open("abaqus/materials.json",'r') as material_file:
            materialDict = json.loads(material_file.read())
        
        if materialName in list(materialDict.keys()):
            inputFileString = \
            "**\n" + \
            "** MATERIALS\n" + \
            "**\n" + \
            "*Material, name=" + materialName + "\n" + \
            "*Density\n" + \
            str(materialDict[materialName]["DENSITY"]) + "\n" + \
            "*Elastic\n" + \
            str(materialDict[materialName]["YOUNG_MODULUS"]) + ", " + str(materialDict[materialName]["POISSON_RATIO"]) + "\n"
            
            inputFileString += "*Conductivity\n"
            for i in range(len(list(materialDict[materialName]["HEAT_CONDUCTIVITY"]["VALUES"]))):
                inputFileString += str(materialDict[materialName]["HEAT_CONDUCTIVITY"]["VALUES"][i]) + ", " + str(materialDict[materialName]["HEAT_CONDUCTIVITY"]["TEMPERATURES"][i]) + "\n"
            
            inputFileString += "*Specific Heat\n"
            for i in range(len(list(materialDict[materialName]["SPECIFIC_HEAT_CAPACITY"]["VALUES"]))):
                inputFileString += str(materialDict[materialName]["SPECIFIC_HEAT_CAPACITY"]["VALUES"][i]) + ", " + str(materialDict[materialName]["SPECIFIC_HEAT_CAPACITY"]["TEMPERATURES"][i]) + "\n"
                        
        else:
            print("Material \'" + materialName + "\' is not yet implemented. Available materials are:")
            for i in range(len(list(materialDict.keys()))):
                print("\t" + list(materialDict.keys())[i])
            exit()
        with open("abaqus/materialProperties.inp",'w') as f:
            f.writelines(inputFileString)

    def setBCType(identifier):
        with open("BCTypes.json",'r') as bc_file:
            bc_dict = json.loads(bc_file.read())
        
        BCType = bc_dict["Abaqus"][identifier]
        
        return BCType

    def setBC(BCType, surfaceID, values):
        EX_AbaqusControl.checkSurfaceID(surfaceID)
        global BCString

        if BCType == "Dirichlet":

            BCString += "**\n"
            BCString += "** Name: mappedTempDirichlet, surfaceName:" + surfaceID + "\n"
            BCString += "**\n"
            BCString += "*Boundary\n"

            nodeIDs = EX_AbaqusControl.getNodesOfSurface(surfaceID)
            for i in range(len(nodeIDs)):
                BCString += "Structure_Instance." + str(nodeIDs[i]) + ", 11, 11, " + str(values[i]) + "\n"

        elif BCType == "Neumann":

            BCString += "**\n"
            BCString += "** Name: mappedFluxNeumann, surfaceName:" + surfaceID + "\n"
            BCString += "**\n"
            BCString += "*CFlux\n"

            nodeIDs = EX_AbaqusControl.getNodesOfSurface(surfaceID)
            for i in range(len(nodeIDs)):
                BCString += "Structure_Instance." + str(nodeIDs[i]) + ", 11, " + str(values[i]) + "\n"

        elif BCType == "Robin":

            BCString += "**\n"
            BCString += "** Name: mappedFluxRobin, surfaceName:" + surfaceID + "\n"
            BCString += "**\n"
            BCString += "*Cfilm\n"

            nodeIDs = EX_AbaqusControl.getNodesOfSurface(surfaceID)
            for i in range(len(nodeIDs)):
                BCString += "Structure_Instance." + str(nodeIDs[i]) + ", 1, " + str(values[0][i]) + ", " + str(values[1][i]) + "\n"

        else:
            errorLog("BCType \'" + BCType + "\' is not available.")
            return
        
        abqLog("Set " + BCType + " BC on surface \'" + surfaceID + "\'")

    def writeInputFile(timePeriod, **kwargs):

        global BCString

        initialIncrementSize = kwargs.get('initialIncrementSize', timePeriod)
        minimumIncrementSize = kwargs.get('minimumIncrementSize', timePeriod / 1000)
        maximumIncrementSize = kwargs.get('maximumIncrementSize', timePeriod)
        maximumNumberOfIncrements = kwargs.get('maximumNumberOfIncrements', 10000)

        stepString = "*INCLUDE, INPUT=setupCase.inp \n" + \
        "*INCLUDE, INPUT=materialProperties.inp \n" + \
        "** ----------------------------------------------------------------\n" + \
        "**\n" + \
        "** STEP: StructureStep\n" + \
        "**\n" + \
        "*Step, name=StructureStep, nlgeom=NO, inc=" + str(maximumNumberOfIncrements) + "\n" + \
        "*Heat Transfer, end=PERIOD, deltmx=1000\n" + \
        str(initialIncrementSize) + ", " + str(timePeriod) + ", " + str(minimumIncrementSize) + ", " + str(maximumIncrementSize) + "\n" + \
        "**\n"
        
        endString = "**\n" + \
        "** OUTPUT REQUESTS\n" + \
        "**\n" + \
        "*Restart, write, frequency=1\n" + \
        "**\n" + \
        "** FIELD OUTPUT: F-Output-1\n" + \
        "**\n" + \
        "*Output, field, variable=PRESELECT\n" + \
        "*Output, history, frequency=0\n" + \
        "*End Step\n"

        inputFileString = stepString + BCString + endString
        with open("inputFile.inp",'w') as f:
            f.writelines(inputFileString)

        BCString = ""

    def readOdbFile(dataType, odbFileName, stepName):
        abqLog("Reading  field \'" + dataType + "\'...")
        os.system("abaqus cae noGUI='abaqus/ExControl/EX_OdbReader.py' -- " + dataType + " " + odbFileName + " " + stepName)

    def getSurfaceData(dataType, surfaceID, stepName):

        abqLog("Reading \'" + dataType + "\' on surface \'" + surfaceID + "\'")

        nodeIDs = np.array(EX_AbaqusControl.getNodesOfSurface(surfaceID))
        data = np.empty(len(nodeIDs),np.float64)

        with open("abaqus/resultFolder/allNodes_" + dataType.upper() + "_" + stepName + ".csv",'r') as result_file:
            all_data = result_file.read().splitlines()

        counter = 0
        for i in range(len(all_data)):
            if i == int(nodeIDs[counter]):
                data[counter] = all_data[i]
                counter += 1
                if counter == len(data):
                    break

        return data

    def solveAbaqus():
        abqLog("Starting Abaqus Analysis...")
        
        os.system("mkdir abaqus/kreuzeberg_abaqus")

        start_time = time.time()
        os.system("abaqus job=abaqus/inputFile scratch=abaqus/kreuzeberg_abaqus")
        while "inputFile.lck" in os.listdir():
            time.sleep(1)
        end_time = time.time() - 2
        abqLog("Finished Analysis in " + str(end_time - start_time) + " s.")
