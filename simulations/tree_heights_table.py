#!/usr/bin/python

import argparse
import os, subprocess, msprime, statistics, pyslim
import numpy as np
#import numpy
#import matplotlib.pyplot as plt

"""
This script converts a SLiM tree file to the VCF format, adds neutral mutations and extract a subset of individuals.
usage: tree2vcf.py [-h] -t TREE -v VCF [-s SUB]
optional arguments:
  -h, --help            show this help message and exit
  -t TREE, --tree TREE    (required) The name of the tree file from SLiM
  -v VCF, --vcf VCF     (required) The name of the output VCF
  -s SUB, --subset SUB  The number of individuals to subset (default 100)
"""

# get arguments
parser = argparse.ArgumentParser(description="Convert a SLiM tree file to the VCF format")
parser.add_argument("-t","--tree", help="(required) (required) The name of the tree file from SLiM", type=str, required=True)
parser.add_argument("-v","--vcf", help="(required) The name of the output VCF file", type=str, required=False)
parser.add_argument("-s","--subset", help="The number of individuals to subset (default 100)", type=int, required=False, default=100)
args = parser.parse_args()

#get a sample of individuals

# Run the SLiM model and load the resulting .trees file  subprocess.check_output(["slim", "-m", "-s", "0", "ex2_TS.slim"])

ts = pyslim.load(args.tree).simplify()

# Measure the tree height at each base position
height_for_pos = np.zeros(int(ts.sequence_length))
for tree in ts.trees():
	mean_height = statistics.mean([tree.time(root) for root in tree.roots])
	left, right = map(int, tree.interval)
#	print(left,"_",right)
	height_for_pos[left: right] = mean_height





def sliding_window(data, size, stepsize=1, padded=False, axis=-1, copy=True):
    if axis >= data.ndim:
        raise ValueError(
            "Axis value out of range"
        )

    if stepsize < 1:
        raise ValueError(
            "Stepsize may not be zero or negative"
        )

    if size > data.shape[axis]:
        raise ValueError(
            "Sliding window size may not exceed size of selected axis"
        )

    shape = list(data.shape)
    shape[axis] = numpy.floor(data.shape[axis] / stepsize - size / stepsize + 1).astype(int)
    shape.append(size)

    strides = list(data.strides)
    strides[axis] *= stepsize
    strides.append(data.strides[axis])

    strided = numpy.lib.stride_tricks.as_strided(
        data, shape=shape, strides=strides
    )

    if copy:
        return strided.copy()
    else:
        return strided



a=height_for_pos.astype(int)
i=1
for win in sliding_window(a, size=40000,stepsize=40000):
	print((i*40000+(i-1)*40000)/2,np.mean(win))
	i=i+1
