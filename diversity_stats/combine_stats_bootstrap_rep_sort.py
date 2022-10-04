#!/usr/bin/python

import argparse, os, sys
import numpy as np
import pandas as pd


"""

python3 ./combine_stats_bootstrap_rep_sorted.py --help
usage: combine_stats_bootstrap_rep_sorted.py [-h] -i INFILE -b BETA -d DOMAINS -o
                                      OUTFILE [--repl REPL]

Please, normalize a column

optional arguments:
  -h, --help            show this help message and exit
  -i INFILE, --infile INFILE
                        (required) The name of the input stats file
  -b BETA, --beta BETA  (required) The name of the input beta file
  -d DOMAINS, --domains DOMAINS
                        (required) The domain file
  -o OUTFILE, --outfile OUTFILE
                        (required) The name of the input file
  --repl REPL           The number of replicates (default is 50)

"""

# get arguments
parser = argparse.ArgumentParser(description="Please, normalize a column")
parser.add_argument("-i","--infile", help="(required) The name of the input stats file", type=str, required=True)
parser.add_argument("-b","--beta", help="(required) The name of the input beta file", type=str, required=True)
parser.add_argument("-d","--domains", help="(required) The domain file", type=str, required=True)
parser.add_argument("-o","--outfile", help="(required) The name of the input file", type=str, required=True)
parser.add_argument("--repl", help="The number of replicates (default is 50)", type=int, required=False, default=50)

args = parser.parse_args()


#########################################################

#identical to normalizeFeatureVec from diploSHIC
def normalizeColumn(statVec):
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


def Reshape(data):
	out = data.values.reshape((3, 25)).astype(float)
	out = pd.DataFrame(out)
	return out


##########################################################

#read the files
stats = pd.read_csv(args.infile,sep="\t")
beta = pd.read_csv(args.beta,sep="\t", names=["chrom","classifiedWinStart", "classifiedWinEnd","beta"])
domains = pd.read_csv(args.domains,sep="\t")
pseudodoord=pd.read_csv("/projects/phillipslab/ateterina/CR_popgen/scripts/4columns.txt",sep="\t")


#change columns in stats
#2nd
stats["classifiedWinStart"] = stats["classifiedWinStart"].apply(lambda x: x+49999)
#3rd +5000
stats["classifiedWinEnd"] = stats["classifiedWinEnd"].apply(lambda x: x+50000)


#merge diploSHIC' stats and beta
stats = stats.merge(beta, how='left', on=["chrom","classifiedWinStart", "classifiedWinEnd"])
stats.columns = stats.columns.str.replace("_win0", "")


#left and right boundaries of the Central domain
left = pd.to_numeric(domains[domains['CHR'] == stats['chrom'].values[0]].DOML.values[0])
right = pd.to_numeric(domains[domains['CHR'] == stats['chrom'].values[0]].DOMR.values[0])


#subset domains
LEFTARM = stats[stats["classifiedWinEnd"] < left]
CENTER = stats[(stats["classifiedWinEnd"] > left) & (stats["classifiedWinStart"] < right)]
RIGHTARM = stats[stats["classifiedWinStart"] > right]


print(stats['chrom'].values[0])
print("Left arm:  "+ str(LEFTARM.shape) + "  Center:  " + str(CENTER.shape) + "  Right arm:  " + str(RIGHTARM.shape))
print("Replicates "+ str(args.repl))


##########################################################################

#bootstrap
for repl in range(0,args.repl):
	OUT = pd.concat([LEFTARM.sample(n = 25, replace = True), CENTER.sample(n = 25, replace = True),RIGHTARM.sample(n = 25, replace = True)], axis=0,ignore_index=True)

	OUT=OUT.sort_values('classifiedWinStart')

	OUT.drop(labels=["chrom","classifiedWinStart","classifiedWinEnd","bigWinRange"],axis = 1, inplace = True)

	#normalize
	OUT=OUT.apply(normalizeColumn, axis=0)

	#reshape and concatenate
	OUT2=pd.DataFrame()

	for name in OUT.columns:
		OUT2 = pd.concat([OUT2,Reshape(OUT[name])], axis=1, ignore_index=True)


	OUT2.columns = [str(name) + "_win_" + str(x) for name in OUT.columns for x in range(0,25)]

	#add pseudocoordinates to have the same 1-4 columns as in SLiM files
	OUT2=pd.concat([pseudodoord,OUT2],axis=1,sort=False)
	filename = args.outfile + "_sorted_" +  str(repl) + ".txt"

	##save into a file
	OUT2.to_csv(filename, header=True, index=None, sep='\t', mode='w')

##########################################################################

print("Done!")
