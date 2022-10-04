#!/usr/bin/python

import argparse, os, sys
import numpy as np
import pandas as pd

"""
#python3 ./normalize_column.py -i TEST.BETA -o OUT
usage: normalize_column.py [-h] -i INFILE -o OUTFILE [-c COLUMN]

Please, normalize a column

optional arguments:
  -h, --help            show this help message and exit
  -i INFILE, --infile INFILE
                        (required) The name of the input file
  -o OUTFILE, --outfile OUTFILE
                        (required) The name of the input file
  -c COLUMN, --column COLUMN
                        The number of the column to normalize (default is 4)
"""

# get arguments
parser = argparse.ArgumentParser(description="Please, normalize a column")
parser.add_argument("-i","--infile", help="(required) The name of the input file", type=str, required=True)
parser.add_argument("-o","--outfile", help="(required) The name of the input file", type=str, required=True)
parser.add_argument("-c","--column", help="The number of the column to normalize (default is 4)", type=int, required=False, default=4)
args = parser.parse_args()


#########################################################

#almost identical to normalizeFeatureVec from diploSHIC

def normalize(statVec):
    minVal = min(statVec)
    if minVal < 0:
        statVec = [x-minVal for x in statVec]
    normStatVec = []
    statSum = float(sum(statVec))
    if statSum == 0 or any(np.isinf(statVec)) or any(np.isnan(statVec)):
        normStatVec = [1.0/len(statVec)]*len(statVec)
    else:
        for k in range(len(statVec)):
            normStatVec.append(statVec[k]/statSum)
    return normStatVec



##########################################################

#load a file
data = pd.read_csv(args.infile,sep="\t",header=None)

out = normalize(data[data.columns[args.column-1]])
out=pd.DataFrame(out)
#or reshape it to 4 lines (25 subwindows)
out = out.values.reshape((3, 25)).astype(float)
out = pd.DataFrame(out)
out.columns = ["beta_win" + str(x) for x in range(0,25)]

out.to_csv(args.outfile, header=True, index=None, sep='\t', mode='w')
