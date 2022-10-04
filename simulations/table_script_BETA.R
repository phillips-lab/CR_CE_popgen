#a script to combine SLiM stats

args = commandArgs(trailingOnly=TRUE)

if (length(args)==0) {
    args[1] = "SLiM_worms.BETA.table.out"
}

my_files<-list.files(pattern="*kb.BETA")
my_files
my_names<-gsub(".BETA","", perl=T, my_files)
STATS<-c();for (i in 1:(length(my_files)) ) { B=read.csv(file=my_files[i], header=F, sep="\t"); colnames(B)<-c("CHR","start","end","BETA"); STATS<-rbind(STATS,data.frame(B,sample=my_names[i]));}

STATS<-as.data.frame(STATS)

write.table(STATS,file=args[1],quote = FALSE, sep = "\t",row.names = FALSE, col.names = TRUE)
