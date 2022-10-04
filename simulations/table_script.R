#a script to combine SLiM stats

args = commandArgs(trailingOnly=TRUE)

if (length(args)==0) {
    args[1] = "SLiM_worms.table.out"
}

my_files<-list.files(pattern="*12stats.txt")
my_files
my_names<-gsub(".txt","", perl=T, my_files)
STATS<-c();for (i in 1:(length(my_files)) ) {
  B=read.csv(file=my_files[i], header=T, sep="\t");
 STATS<-rbind(STATS,data.frame(B,sample=my_names[i]));
}

write.table(STATS,file=args[1],quote = FALSE, sep = "\t",row.names = FALSE, col.names = TRUE)
