#!/usr/bin/python

import random
import pandas as pd


names = pd.read_csv("/projects/phillipslab/ateterina/CR_popgen/data/reads/BAMS/EXONS_INTRONS/LIST_WILD_IND.txt",names=['IND'])


indlist=[]
n=0

#generate 100 random combinations of 8 individuals splitted by comma
while n<101:
    A=names.sample(n=8)
    A=[','.join(val) for val in A.values.tolist()]

    indlist.append(A)
    n=n+1



with open("/projects/phillipslab/ateterina/CR_popgen/data/reads/BAMS/EXONS_INTRONS/SMCPP/NAMES8.txt","w") as txt:
	txt.write('\n'.join(','.join(elems) for elems in indlist))
