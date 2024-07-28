import pickle
import numpy as np
from datetime import datetime

import pandas as pd
import statsmodels.api as sm
import sys
import re
import warnings
from scipy.stats import ks_2samp
import glob
from pathlib import Path
import os
import json

#This script checks if U, L, x0A,x0B,x1A, x1B values reach equilibrium and writes summary file of dist


argErrCode=2
sameErrCode=3
if (len(sys.argv)!=3):
    print("wrong number of arguments")
    exit(argErrCode)


jsonFromSummaryLast=json.loads(sys.argv[1])
jsonDataFromConf=json.loads(sys.argv[2])

TDirRoot=jsonFromSummaryLast["TDirRoot"]
U_dist_dataDir=jsonFromSummaryLast["U_dist_dataDir"]
effective_data_num_required=int(jsonDataFromConf["effective_data_num_required"])


summary_U_distFile=TDirRoot+"/summary_U_dist.txt"


def sort_data_files_by_loopEnd(oneDir):
    dataFilesAll=[]
    loopEndAll=[]
    for oneDataFile in glob.glob(oneDir+"/*.csv"):
        # print(oneDataFile)
        dataFilesAll.append(oneDataFile)
        matchEnd=re.search(r"loopEnd(\d+)",oneDataFile)
        if matchEnd:
            indTmp=int(matchEnd.group(1))
            loopEndAll.append(indTmp)
    endInds=np.argsort(loopEndAll)
    sortedDataFiles=[dataFilesAll[i] for i in endInds]
    return sortedDataFiles


def parseSummaryU_Dist():
    startingFileInd=-1
    startingVecPosition=-1

    summaryFileExists=os.path.isfile(summary_U_distFile)
    if summaryFileExists==False:
        return startingFileInd,startingVecPosition

    with open(summary_U_distFile,"r") as fptr:
        lines=fptr.readlines()
    for oneLine in lines:
        #match startingFileInd
        matchStartingFileInd=re.search(r"startingFileInd=(\d+)",oneLine)
        if matchStartingFileInd:
            startingFileInd=int(matchStartingFileInd.group(1))

        #match startingVecPosition
        matchStartingVecPosition=re.search(r"startingVecPosition=(\d+)",oneLine)
        if matchStartingVecPosition:
            startingVecPosition=int(matchStartingVecPosition.group(1))

    return startingFileInd, startingVecPosition



def checkU_distDataFilesForOneT(U_dist_csv_dir):
    """

    :param dist_csv_dir:
    :return:
    """
    U_dist_sortedDataFilesToRead=sort_data_files_by_loopEnd(U_dist_csv_dir)
    if len(U_dist_sortedDataFilesToRead)==0:
        print("no data for U_dist.")
        exit(0)

    startingFileInd,startingVecPosition=parseSummaryU_Dist()

    if startingFileInd<0:
        #we guess that the equilibrium starts at this file
        startingFileInd=int(len(U_dist_sortedDataFilesToRead)*5/7)
    startingFileName=U_dist_sortedDataFilesToRead[startingFileInd]
    # print(startingFileName)
    #read the starting U_dist csv file
    in_dfStart=pd.read_csv(startingFileName,dtype={"U":float,"L":float,"x0A":float,"x0B":float,"x1A":float,"x1B":float})

    in_nRowStart,in_nColStart=in_dfStart.shape
    if startingVecPosition<0:
        #we guess equilibrium starts at this position
        startingVecPosition=int(in_nRowStart/2)

    UVec=np.array(in_dfStart.loc[startingVecPosition:,"U"])
    LVec=np.array(in_dfStart.loc[startingVecPosition:,"L"])
    x0AVec=np.array(in_dfStart.loc[startingVecPosition:,"x0A"])
    x0BVec=np.array(in_dfStart.loc[startingVecPosition:,"x0B"])
    x1AVec=np.array(in_dfStart.loc[startingVecPosition:,"x1A"])
    x1BVec=np.array(in_dfStart.loc[startingVecPosition:,"x1B"])

    #read the rest of the csv files
    for csv_file in U_dist_sortedDataFilesToRead[(startingFileInd+1):]:
        in_df=pd.read_csv(csv_file,dtype={"U":float,"L":float,"x0A":float,"x0B":float,"x1A":float,"x1B":float})
        in_U=np.array(in_df.loc[:,"U"])
        in_L=np.array(in_df.loc[:,"L"])
        in_x0A=np.array(in_df.loc[:,"x0A"])
        in_x0B=np.array(in_df.loc[:,"x0B"])
        in_x1A=np.array(in_df.loc[:,"x1A"])
        in_x1B=np.array(in_df.loc[:,"x1B"])

        UVec=np.r_[UVec,in_U]
        LVec=np.r_[LVec,in_L]
        x0AVec=np.r_[x0AVec,in_x0A]
        x0BVec=np.r_[x0BVec,in_x0B]
        x1AVec=np.r_[x1AVec,in_x1A]
        x1BVec=np.r_[x1BVec,in_x1B]


    NLags=int(len(LVec)*3/4)
    eps=1e-3
    lagVal=-1
    pThreshHold=0.1
    same=False
    with warnings.catch_warnings():
        warnings.filterwarnings("error")
    try:
        acfOfVecU=sm.tsa.acf(UVec,nlags=NLags)
    except Warning as w:
        same=True

    with warnings.catch_warnings():
        warnings.filterwarnings("error")
    try:
        acfOfVecL=sm.tsa.acf(LVec,nlags=NLags)
    except Warning as w:
        same=True

    with warnings.catch_warnings():
        warnings.filterwarnings("error")
    try:
        acfOfVecx0A=sm.tsa.acf(x0AVec,nlags=NLags)
    except Warning as w:
        same=True

    with warnings.catch_warnings():
        warnings.filterwarnings("error")
    try:
        acfOfVecx0B=sm.tsa.acf(x0BVec,nlags=NLags)
    except Warning as w:
        same=True

    with warnings.catch_warnings():
        warnings.filterwarnings("error")
    try:
        acfOfVecx1A=sm.tsa.acf(x1AVec,nlags=NLags)
    except Warning as w:
        same=True


    with warnings.catch_warnings():
        warnings.filterwarnings("error")
    try:
        acfOfVecx1B=sm.tsa.acf(x1BVec,nlags=NLags)
    except Warning as w:
        same=True

    #all values are the same, exit with err code
    if same==True:
        with open(summary_U_distFile,"w+") as fptr:
            msg="error: same\n"
            fptr.writelines(msg)

    acfOfVecLAbs=np.abs(acfOfVecL)
    acfOfVecx0AAbs=np.abs(acfOfVecx0A)
    acfOfVecx0BAbs=np.abs(acfOfVecx0B)
    acfOfVecx1AAbs=np.abs(acfOfVecx1A)
    acfOfVecx1BAbs=np.abs(acfOfVecx1B)

    acfOfVecUAbs=np.abs(acfOfVecU)

    minAutcL=np.min(acfOfVecLAbs)
    minAutcx0A=np.min(acfOfVecx0AAbs)
    minAutcx0B=np.min(acfOfVecx0BAbs)
    minAutcx1A=np.min(acfOfVecx1AAbs)
    minAutcx1B=np.min(acfOfVecx1BAbs)
    minAutcU=np.min(acfOfVecUAbs)


    if minAutcL<eps and minAutcx0A< eps \
            and minAutcx0B<eps and minAutcx1A<eps and minAutcx1B<eps\
            and minAutcU<eps:
        lagL=np.where(acfOfVecLAbs<=eps)[0][0]

        lagx0A=np.where(acfOfVecx0AAbs<=eps)[0][0]
        lagx0B=np.where(acfOfVecx0BAbs<=eps)[0][0]
        lagx1A=np.where(acfOfVecx1AAbs<=eps)[0][0]
        lagx1B=np.where(acfOfVecx1BAbs<=eps)[0][0]
        lagU=np.where(acfOfVecUAbs<=eps)[0][0]
        print("lagL="+str(lagL))
        print("lagx0A="+str(lagx0A))
        print("lagx0B="+str(lagx0B))
        print("lagx1A="+str(lagx1A))
        print("lagx1B="+str(lagx1B))
        print("lagU="+str(lagU))

        lagVal=np.max([lagL,lagx0A,lagx0B,lagx1A,lagx1B,lagU])

        # #     #select values by lagVal
        LVecSelected=LVec[::lagVal]
        x0AVecSelected=x0AVec[::lagVal]
        x0BVecSelected=x0BVec[::lagVal]
        x1AVecSelected=x1AVec[::lagVal]
        x1BVecSelected=x1BVec[::lagVal]
        UVecSelected=UVec[::lagVal]

        lengthTmp=len(LVecSelected)
        if lengthTmp%2==1:
            lengthTmp-=1
        lenPart=int(lengthTmp/2)

        LVecValsToCompute=LVecSelected[-lengthTmp:]
        x0AVecValsToCompute=x0AVecSelected[-lengthTmp:]
        x0BVecValsToCompute=x0BVecSelected[-lengthTmp:]
        x1AVecValsToCompute=x1AVecSelected[-lengthTmp:]
        x1BVecValsToCompute=x1BVecSelected[-lengthTmp:]
        UVecValsToCompute=UVecSelected[-lengthTmp:]

        #     # #test distributions

        #test L
        selectedLVecPart0=LVecValsToCompute[:lenPart]
        selectedLVecPart1=LVecValsToCompute[lenPart:]
        resultL=ks_2samp(selectedLVecPart0,selectedLVecPart1)
        # print(resultL)

        #test x0A
        selectedx0AVecPart0=x0AVecValsToCompute[:lenPart]
        selectedx0AVecPart1=x0AVecValsToCompute[lenPart:]
        resultx0A=ks_2samp(selectedx0AVecPart0,selectedx0AVecPart1)

        #test x0B
        selectedx0BVecPart0=x0BVecValsToCompute[:lenPart]
        selectedx0BVecPart1=x0BVecValsToCompute[lenPart:]
        resultx0B=ks_2samp(selectedx0BVecPart0,selectedx0BVecPart1)


        #test x1A
        selectedx1AVecPart0=x1AVecValsToCompute[:lenPart]
        selectedx1AVecPart1=x1AVecValsToCompute[lenPart:]
        resultx1A=ks_2samp(selectedx1AVecPart0,selectedx1AVecPart1)

        #test x1B
        selectedx1BVecPart0=x1BVecValsToCompute[:lenPart]
        selectedx1BVecPart1=x1BVecValsToCompute[lenPart:]
        resultx1B=ks_2samp(selectedx1BVecPart0,selectedx1BVecPart1)

        #test U
        selectedUVecPart0=UVecValsToCompute[:lenPart]
        selectedUVecPart1=UVecValsToCompute[lenPart:]
        resultU=ks_2samp(selectedUVecPart0,selectedUVecPart1)

        numDataPoints=len(selectedLVecPart0)+len(selectedLVecPart1)
        pValsAll=np.array([resultL.pvalue,resultx0A.pvalue,resultx0B.pvalue,resultx1A.pvalue,resultx1B.pvalue,resultU.pvalue])

        if np.min(pValsAll)>=pThreshHold and numDataPoints>=200:
            if numDataPoints>=effective_data_num_required:
                newDataPointNum=0
            else:
                newDataPointNum=effective_data_num_required-numDataPoints
            msg="equilibrium\n" \
                +"lag="+str(lagVal)+"\n" \
                +"numDataPoints="+str(numDataPoints)+"\n" \
                +"startingFileInd="+str(startingFileInd)+"\n" \
                +"startingVecPosition="+str(startingVecPosition)+"\n" \
                +"newDataPointNum="+str(newDataPointNum)+"\n"

            with open(summary_U_distFile,"w+") as fptr:
                fptr.writelines(msg)
            exit(0)

        continueMsg="continue\n"
        if np.min(pValsAll)<pThreshHold:
            #not the same distribution
            continueMsg+="p value: "+str(np.min(pValsAll))+"\n"
        if numDataPoints<200:
            #not enough data number

            continueMsg+="numDataPoints="+str(numDataPoints)+" too low\n"
            continueMsg+="lag="+str(lagVal)+"\n"
        with open(summary_U_distFile,"w+") as fptr:
            fptr.writelines(continueMsg)
        exit(0)

    else:
        highMsg="high correlation"
        with open(summary_U_distFile,"w+") as fptr:
            fptr.writelines(msg)
        exit(0)


checkU_distDataFilesForOneT(U_dist_dataDir)