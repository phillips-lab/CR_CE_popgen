#####################################################
########### VCF Filters #############################
#####################################################


setwd("/Users/anastasia/Documents/Phillips_lab/drafts/CR_popgen/TEST_VCF_FILTER")

library(ggplot2)
library(gridExtra)
library(magick)
library(dichromat)


#diploSHIC statistics
my_files<-list.files(pattern="CR_WILD_population")
my_names<-gsub("CR_WILD_population","", perl=T, my_files)
my_names<-gsub(".100K.0.1.stats","", perl=T, my_names)
my_names<-gsub(".NOPHASE","", perl=T, my_names)
my_names<-gsub("_.T","T", perl=F, my_names)


STATS<-c();for (i in 1:(length(my_files)) ) { B=read.csv(file=my_files[i], header=T, sep="\t"); STATS<-rbind(STATS,data.frame(B,sample=my_names[i]));}

colnames(STATS)<-gsub("_win0","",colnames(STATS))
STATS$Type<-gsub(".[IVX]+$","",perl=T,STATS$sample)
colnames(STATS)<-gsub("diplo_","",colnames(STATS))
colnames(STATS)[16]<-"omega"
colnames(STATS)[6]<-"theta"
colnames(STATS)[7] <-"TajimaD"


chrom <-c("I","II","III","IV","V","X")
###chromosome sizes
endsCR<-c(17247545,19935723,17877849,25790997,22502457,21501900)


# the central domains
domains<-data.frame(classifiedWinEnd=c(5803000,13007000,6828000, 14292000,5849000,11888000,8639000,17108000,6767000,13852000,9040000,16070000), chrom=c("I","I","II","II","III","III","IV","IV","V","V","X","X"))


colorsMUT<-c("#183843", "#196770", "#73A8BD", "#C0BDBC", "#f2c76f", "#916b20")


colfunc<-colorRampPalette(colorsMUT)

STATS2<-STATS[STATS$Type %in% c("5-100_0.25","5-100_0.5","5-100_0.75","10-100_0.25","10-100_0.5","5-200_0.25","5-200_0.5","5-200_0.75","10-200_0.25","10-200_0.5"),]
#STATS3<-STATS[STATS$Type %in% c("14_filt_snps_5-100_0.5_fin","14_filt_snps_5-100_0.75_fin","TESTFILT.5-100_0.5","TESTFILT.5-100_0.75","TESTFILT.3-100_0.5","TESTFILT.3-100_0.75","TESTFILT.10-100_0.5"),]

# "14iltnps-100.75in","TESTFILT.5-100_0.5","TESTFILT.5-100_0.75","TESTFILT.5-200_0.5","TESTFILT.5-200_0.75"),]

STATS$Type<-gsub("14_filt_snps","14_ind",STATS$Type)
STATS$Type<-gsub("_fin","",STATS$Type)
STATS$Type<-gsub("TESTFILT.","17_ind_",STATS$Type)



STATS3<-STATS[STATS$Type %in% c("14_ind_5-100_0.5","14_ind_5-100_0.75","17_ind_5-100_0.5","17_ind_5-100_0.75","17_ind_3-100_0.5","17_ind_3-100_0.75","17_ind_10-100_0.5"),]
filtcolors<-colfunc(7)
filtcolors[1]<-"black"
for (stat in 5:16){
  #  geom_point(alpha=0.2,size=0.3) +
  pl<-ggplot(STATS3, aes(x = classifiedWinEnd, y = STATS3[,stat], ordered = FALSE,color=Type))  +
    labs(title ="", x = "Genome position",y = ggplot2:::parse_safe(colnames(STATS)[stat])) +  scale_color_manual(values=filtcolors,name="Filter") +
    scale_y_continuous(labels = function(x) format(x, scientific = FALSE)) + scale_x_continuous(labels = function(x) format(x/1000000)) +
    theme_bw(base_size = 12) +
    theme(axis.text.x = element_text(angle=0,vjust=0.5, hjust = 0.5), plot.title = element_text(hjust = 0, face = "bold"),strip.background = element_blank(),panel.grid.major = element_blank(), panel.grid.minor = element_blank(), strip.text = element_text(size = rel(1))) +
    theme(legend.position = "right") +
    theme(strip.text.y = element_text(face = "italic",angle=0))+
    theme(strip.background = element_blank(),panel.grid.major = element_blank(), panel.grid.minor = element_blank(), strip.text = element_text(size = rel(1))) +
    facet_grid(. ~ chrom, space = "free", scales="free" ) +
    geom_vline(data = domains,aes(xintercept = classifiedWinEnd),colour = "#bfbfbf",size = 0.75,linetype = "dashed") +
    geom_point(alpha=0.2,size=0.2) +
    geom_smooth(alpha=0.5,size=0.5,se =FALSE,method = 'loess',span=0.4) + theme(plot.margin = unit(c(1, 3, 1.5, 3), "pt"))

  tiff(paste0(colnames(STATS3)[stat],".filters.tiff"),width = 7.2,height = 2.5,res=300, units="in")
  plot(pl)
  dev.off()

}







PIC <- file.info(list.files(pattern="*.filters.tiff"))
PIC <- PIC[with(PIC, order(as.POSIXct(mtime))), ]
# just 5 of them pi,  theta, Tajima'sD, ZnS, and omega
PIC<-rownames(PIC)[c(1,2,3,11,12)]
combo<-image_append(do.call("c", lapply(PIC, function(h){image_read(h)})), stack = TRUE)
image_write(combo, path = "A1.7.tiff", format = "tiff", quality = 100)
