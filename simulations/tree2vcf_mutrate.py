#!/usr/bin/python

import argparse, os, sys, re, msprime, pyslim
import numpy as np
import pandas as pd
import random
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
parser.add_argument("-v","--vcf", help="(required) The name of the output VCF file", type=str, required=True)
parser.add_argument("-s","--subset", help="The number of individuals to subset (default 100)", type=int, required=False, default=100)
parser.add_argument("-m","--mu", help="Mutation rate", type=float, required=False, default=2e-8)
parser.add_argument("-d","--diff", help="Mutation rate difference between domains", type=float, required=False, default=1)
parser.add_argument("-n","--Ne", help="Population size", type=int, required=False, default=5000)
parser.add_argument("-f","--frac", help="Fraction of non neutral variation", type=float, required=False, default=0)


args = parser.parse_args()

#get a sample of individuals
def treeSubset(ts,sample_size):
    ts = ts.simplify()
    inds = [ind.id for ind in ts.individuals()]
    sample_subset = np.sort(np.random.choice(inds,sample_size,replace=False))
    sample_nodes = []
    for i in sample_subset:
        ind = ts.individual(i)
        sample_nodes.append(ind.nodes[0])
        sample_nodes.append(ind.nodes[1])

    ts = ts.simplify(sample_nodes)
    return ts


#difference in the mu between domians
prop=1-1/args.diff

#neutral variation
mscale=(1-args.frac)*args.diff*args.mu

#load the trees and the recombination map
tr = pyslim.load(args.tree)
recomb_map = msprime.RecombinationMap.read_hapmap("/projects/phillipslab/ateterina/slim/worms_snakemake/worms.hapmap")

#recapitate & add neutral variation - warning, mostly it will be used in the outcrossing populations, and that's good, because recapitation for selfers is not available
recapped = tr.recapitate(recombination_map=recomb_map, Ne=args.Ne) #I used it for all data
ts = pyslim.SlimTreeSequence(msprime.mutate(recapped, rate=mscale, keep=True))


if args.diff>1:
	print("Removing "+str(prop*100)+"% of mutations within the central region")
	pos = ts.tables.sites.position
	is_post_recap = np.repeat(False, ts.num_sites)
	temp = ts.tables.nodes.time[ts.tables.mutations.node] < ts.slim_generation
	is_post_recap[ts.tables.mutations.site] = temp
	in_regions=np.logical_and(pos>999999, pos<1999999)
	is_msp = (np.diff(ts.tables.mutations.metadata_offset) == 0)
	is_msp_site = np.repeat(False, ts.num_sites)
	is_msp_site[ts.tables.mutations.site] = is_msp
	removable_sites = np.where(np.logical_and(np.logical_and(in_regions, is_msp_site), is_post_recap))[0]
	remove = np.array(random.sample(list(removable_sites),int(len(removable_sites)*prop)))
	new_table = ts.tables
	new_table.delete_sites(remove)
	filtered=new_table.tree_sequence()
	subtree=treeSubset(filtered,args.subset)
else:
	print("No need to remove extra mutations")
	subtree=treeSubset(ts,args.subset)


subtree=pyslim.SlimTreeSequence(subtree)


#output all indivudualt to a vcf
with open(args.vcf, "w") as vcf_file:
	subtree.write_vcf(vcf_file)
