#!/usr/bin/env Rscript

library("optparse")

option_list = list(
				make_option(c("-i", "--infile"), type="character", default=NULL,
							help="input file", metavar="character"),
				make_option(c("-o", "--outfile"), type="character", default=NULL,
							help="output file", metavar="character"),
				make_option(c("-r", "--rectable"), type="character", default=NULL,
							help="recombination table", metavar="character"));


opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);




#read the map file
map <- read.csv(file=opt$infile, header=FALSE, sep=" ", stringsAsFactors=FALSE)

#read the table file
rec <- read.csv(file=opt$rectable, header=FALSE, sep=" ", stringsAsFactors=FALSE)
colnames(rec)<-c("chr","start","end","rate")

for (i in 1:nrow(map)){

A <- map[i,];
val <- (A[1,4]*rec[rec$chr== A[1,1] & rec$start<= A[1,4] & rec$end> A[1,4],"rate"])/1000000;


cat(c(A[1,1], " ", as.character(A[1,2]), " ", val," ",A[1,4],"\n"));
