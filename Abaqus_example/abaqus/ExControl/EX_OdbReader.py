from odbAccess import openOdb
import numpy
import sys

dataType = sys.argv[-3]
odbFileName = sys.argv[-2]
stepName = sys.argv[-1]

dataTypes = {	'TEMPERATURE':'NT11',   # node based result
                'DEFORMATION':'U',      # node based result
                'HEATFLOW':'RFL11',     # node based result
                'VELOCITY':'V',         # node based result
                'ACCELERATION':'A',     # node based result
                'FORCE':'RF'}           # node based result --> reaction force

if dataType.upper() in dataTypes:
    pass
else:
    print("There is no data type called \'" + dataType + "\'.")
    print("Aborted due to reading error.")
    exit()

odb=openOdb(odbFileName)
values = odb.steps[stepName].frames[-1].fieldOutputs[dataTypes[dataType.upper()]].values
counter = 0

numberOfValues = len(values)

if dataType.upper() in list(dataTypes.keys()):
    data = numpy.empty(numberOfValues,numpy.float64)

if str(values[0].precision) == 'SINGLE_PRECISION':
    if dataType.upper() in ['TEMPERATURE','HEATFLOW']:
        for i in range(numberOfValues):
            data[i] = values[i].data
    elif dataType.upper() in ['DEFORMATION','VELOCITY','ACCELERATION','FORCE']:
        for i in range(numberOfValues):
            data[i] = values[i].data[:3]
else:
    if dataType.upper() in ['TEMPERATURE','HEATFLOW']:
        for i in range(numberOfValues):
            data[i] = values[i].dataDouble
    elif dataType.upper() in ['DEFORMATION','VELOCITY','ACCELERATION','FORCE']:
        for i in range(numberOfValues):
            data[i] = values[i].dataDouble[:3]
odb.close()

numpy.savetxt("abaqus/resultFolder/allNodes_" + dataType.upper() + "_" + stepName + ".csv", data)