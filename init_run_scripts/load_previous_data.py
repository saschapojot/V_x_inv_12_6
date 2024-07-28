import sys
import glob
import re
import json
from decimal import Decimal
import pandas as pd
import numpy as np
import subprocess


#this script loads previous data
numArgErr=4
if (len(sys.argv)!=3):
    print("wrong number of arguments.")
    exit(numArgErr)


jsonDataFromConf =json.loads(sys.argv[1])
jsonFromSummary=json.loads(sys.argv[2])

potential_function_name=jsonDataFromConf["potential_function_name"]
U_dist_dataDir=jsonFromSummary["U_dist_dataDir"]
startingFileInd=jsonFromSummary["startingFileInd"]
startingVecPosition=jsonFromSummary["startingVecPosition"]


#search and read U_dist files
#give arbitrary values to L,x0A, x0B,x1A,x1B without reading data
UInit=6

x0AInit=1
x0BInit=2
x1AInit=3
x1BInit=4
LInit=x1BInit+0.7

loopLastFile=-1

#search csv files
csvFileList=[]
loopEndAll=[]
for file in glob.glob(U_dist_dataDir+"/*.csv"):
    csvFileList.append(file)
    matchEnd=re.search(r"loopEnd(\d+)",file)
    if matchEnd:
        loopEndAll.append(int(matchEnd.group(1)))


def create_loadedJsonData(UVal,LVal,x0AVal,x0BVal,x1AVal,x1BVal,loopLastFileVal):
    """

    :param UVal:
    :param LVal:
    :param x0AVal:
    :param x0BVal:
    :param x1AVal:
    :param x1BVal:
    :param loopLastFileVal:
    :return:
    """
    initDataDict={
        "U":str(UVal),
        "L":str(LVal),
        "x0A":str(x0AVal),
        "x0B":str(x0BVal),
        "x1A":str(x1AVal),
        "x1B":str(x1BVal),
        "loopLastFile":str(loopLastFileVal)
    }
    return json.dumps(initDataDict)

#if no data found, return the arbitrary values
if len(csvFileList)==0:
    loadedJsonDataStr=create_loadedJsonData(UInit,LInit,x0AInit,x0BInit,x1AInit,x1BInit,loopLastFile)
    loadedJsonData_stdout="loadedJsonData="+loadedJsonDataStr
    print(loadedJsonData_stdout)
    exit(0)


#if found csv data
sortedEndInds=np.argsort(loopEndAll)
sortedLoopEnd=[loopEndAll[ind] for ind in sortedEndInds]
sortedCsvFileNames=[csvFileList[ind] for ind in sortedEndInds]
loopLastFile=sortedLoopEnd[-1]

lastFileName=sortedCsvFileNames[-1]

def get_last_row(csv_file):
    result = subprocess.run(['tail', '-n', '1', csv_file], stdout=subprocess.PIPE)
    last_row = result.stdout.decode('utf-8').strip()
    return last_row


csvLastRowStr=get_last_row(lastFileName)
valsInLastRow = [float(value) for value in csvLastRowStr.split(',')]

UInit=valsInLastRow[0]
LInit=valsInLastRow[1]
x0AInit=valsInLastRow[2]
x0BInit=valsInLastRow[3]
x1AInit=valsInLastRow[4]
x1BInit=valsInLastRow[5]


loadedJsonDataStr=create_loadedJsonData(UInit,LInit,x0AInit,x0BInit,x1AInit,x1BInit,loopLastFile)
loadedJsonData_stdout="loadedJsonData="+loadedJsonDataStr
print(loadedJsonData_stdout)
exit(0)