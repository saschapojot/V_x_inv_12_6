from pathlib import Path
from decimal import Decimal
import pandas as pd


#This script creates directories and conf files for mc
TVals=[0.5,1]
rowNum=0
inParamFileName="./V_inv_12_6Params.csv"

inDf=pd.read_csv(inParamFileName)
oneRow=inDf.iloc[rowNum,:]
a1=float(oneRow.loc["a1"])
b1=float(oneRow.loc["b1"])
a2=float(oneRow.loc["a2"])
b2=float(oneRow.loc["b2"])


dataRoot="./dataAll/row"+str(rowNum)+"/"

def format_using_decimal(value):
    # Convert the float to a Decimal
    decimal_value = Decimal(value)
    # Remove trailing zeros and ensure fixed-point notation
    formatted_value = decimal_value.quantize(Decimal(1)) if decimal_value == decimal_value.to_integral() else decimal_value.normalize()
    return str(formatted_value)

TDirsAll=[]
TStrAll=[]
# print(TDirsAll)
for k in range(0,len(TVals)):
    T=TVals[k]
    # print(T)
    TStr=format_using_decimal(T)
    TStrAll.append(TStr)
    TDir=dataRoot+"/T"+TStr+"/"
    TDirsAll.append(TDir)
    Path(TDir).mkdir(exist_ok=True,parents=True)



def contents_to_conf(k):
    """

    :param k: index of T
    :return:
    """

    contents=[
        "#This is the configuration file for mc computations\n",
        "\n"
        "potential_function_name=V_inv_12_6\n",
        "\n" ,
        "#parameters of coefficients\n",
        "#parameter_row=row0\n",
        "#the following row is the values of a1, b1, a2, b2\n"
        "coefs=["+format_using_decimal(a1)+","+format_using_decimal(b1)+","\
        +format_using_decimal(a2)+", "+format_using_decimal(b2)+"]\n",
        "\n",
        "#Temperature\n",
        "T="+TStrAll[k]+"\n"
        "\n",
        "erase_data_if_exist=False\n",
        "\n",
        "search_and_read_summary_file=True\n"
        "\n",
        "#For the observable name, only digits 0-9, letters a-zA-Z, underscore _ are allowed\n",
        "\n",
        "observable_name=U_dist\n",
        "\n",
        "effective_data_num_required=2000\n",
        "\n",
        "loop_to_write=10000\n",
        "\n",
        "#within each flush,  loop_to_write mc computations are executed\n",
        "\n",
        "default_flush_num=10\n"



    ]

    outConfName=TDirsAll[k]+"/run_T"+TStrAll[k]+".mc.conf"
    with open(outConfName,"w+") as fptr:
        fptr.writelines(contents)



for k in range(0,len(TDirsAll)):
    contents_to_conf(k)