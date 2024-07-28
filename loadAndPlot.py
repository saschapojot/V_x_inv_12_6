import numpy as np

import glob

import re
import matplotlib.pyplot as plt
import pandas as pd

dataPath="./dataAll/row0/T1/U_dist_dataFiles/"

csvFileNamesAll=[]
loopEndAll=[]
for file in glob.glob(dataPath+"/*.csv"):
    matchEnd=re.search(r"loopEnd(\d+)",file)
    if matchEnd:
        loopEndAll.append(int(matchEnd.group(1)))
        csvFileNamesAll.append(file)

sortedLpEndInds=np.argsort(loopEndAll)

sortedCsvFilesAll=[csvFileNamesAll[ind] for ind in sortedLpEndInds]
# print(sortedCsvFilesAll)

plt_name="x0A"

vec=np.array([])

for file in sortedCsvFilesAll:
    dfTmp=pd.read_csv(file)
    colValsTmp=np.array(dfTmp[plt_name])
    vec=np.r_[vec,colValsTmp]



startingInd=int(len(vec)*1/4)
plt_vec=vec[startingInd:]
plt.figure()
plt.plot(plt_vec,color="black")
plt.title(plt_name)
plt.savefig("tmp.png")
plt.close()

