import re
import subprocess
import sys

import json
argErrCode=2
if (len(sys.argv)!=2):
    print("wrong number of arguments")
    print("example: python launch_one_run.py /path/to/mc.conf")
    exit(argErrCode)


confFileName=str(sys.argv[1])
invalidValueErrCode=1
summaryErrCode=2
loadErrCode=3
confErrCode=4

#################################################
#parse conf, get jsonDataFromConf
confResult=subprocess.run(["python3", "./init_run_scripts/parseConf.py", confFileName], capture_output=True, text=True)
confJsonStr2stdout=confResult.stdout
# print(confJsonStr2stdout)
if confResult.returncode !=0:
    print("Error running parseConf.py with code "+str(confResult.returncode))
    # print(confResult.stderr)
    exit(confErrCode)



match_confJson=re.match(r"jsonDataFromConf=(.+)$",confJsonStr2stdout)
if match_confJson:
    jsonDataFromConf=json.loads(match_confJson.group(1))
else:
    print("jsonDataFromConf missing.")
    exit(confErrCode)
# print(jsonDataFromConf)
################################################

##################################################
#read summary file, get jsonFromSummary
parseSummaryResult=subprocess.run(["python3","./init_run_scripts/search_and_read_summary.py", json.dumps(jsonDataFromConf)],capture_output=True, text=True)
# print(parseSummaryResult.stdout)
if parseSummaryResult.returncode!=0:
    print("Error in parsing summary with code "+str(parseSummaryResult.returncode))
    # print(parseSummaryResult.stdout)
    # print(parseSummaryResult.stderr)
    exit(summaryErrCode)

match_summaryJson=re.match(r"jsonFromSummary=(.+)$",parseSummaryResult.stdout)
if match_summaryJson:
    jsonFromSummary=json.loads(match_summaryJson.group(1))
# print(jsonFromSummary)
##################################################

###############################################
#load previous data, to get L, x0A,x0B,x1A,x1B,
#get loadedJsonData
loadResult=subprocess.run(["python3","./init_run_scripts/load_previous_data.py", json.dumps(jsonDataFromConf), json.dumps(jsonFromSummary)],capture_output=True, text=True)
# print(loadResult.stdout)
if loadResult.returncode!=0:
    print("Error in loading with code "+str(loadResult.returncode))
    exit(loadErrCode)

match_loadJson=re.match(r"loadedJsonData=(.+)$",loadResult.stdout)
if match_loadJson:
    loadedJsonData=json.loads(match_loadJson.group(1))
else:
    print("loadedJsonData missing.")
    exit(loadErrCode)
###############################################



###########################################################



###########################################################

##########################################################


##########################################################
#statistics
checkU_distErrCode=5
checkU_distResult=subprocess.run(["python3","./oneTCheckObservables/check_U_distOneT.py",json.dumps(jsonFromSummary),json.dumps(jsonDataFromConf)],capture_output=True, text=True)
print(checkU_distResult.stdout)
if checkU_distResult.returncode!=0:
    print("Error in checking dist with code "+str(checkU_distResult.returncode))
    exit(checkU_distErrCode)

##########################################################



