#######################################################
######### Plot diversity statistics ###################
#########  C. elegans & C.remanei   ###################
#######################################################
###########    Fig 2, S2          ######################


setwd("~/Documents/Phillips_lab/drafts/CR_popgen/CR_CE_pi_draft/NEW_HARD_MASK_CE_CR/")

library(ggplot2)
library(gridExtra)
library(lsr)
library(coin)
library(magick)
library(boot)


#diploSHIC statistics
my_files<-list.files(pattern=".100K.0.1.stats")
my_names<-gsub(".100K.0.1.stats","", perl=T, my_files)

STATS<-c();for (i in 1:(length(my_files)) ) { B=read.csv(file=my_files[i], header=T, sep="\t"); STATS<-rbind(STATS,data.frame(B,sample=my_names[i]));}

STATS$Species<-gsub("^([CER]+)_.*","\\1",perl=T,STATS$sample)
STATS$Species<-gsub("CE","C. elegans", STATS$Species)
STATS$Species<-gsub("CR","C. remanei", STATS$Species)
colnames(STATS)<-gsub("_win0","",colnames(STATS))
colnames(STATS)<-gsub("diplo_","",colnames(STATS))
colnames(STATS)[16]<-"omega"
colnames(STATS)[6]<-"theta"
colnames(STATS)[7] <-"TajimaD"

###STATS[is.nan(STATS$tajD),]$tajD <- 0

chrom <-c("I","II","III","IV","V","X")
###chromosome sizes
endsCR<-c(17247545,19935723,17877849,25790997,22502457,21501900)
endsCE<-c(15331301,15525148,14108536,17759200,21243235,18110855)

STATS$NORMPOS<-888

for (i in 1:6){

  STATS[STATS$Species=="C. elegans" & STATS$chrom==chrom[i],]$NORMPOS <- STATS[STATS$Species=="C. elegans" & STATS$chrom==chrom[i],]$classifiedWinEnd/endsCE[i]
  STATS[STATS$Species=="C. remanei" & STATS$chrom==chrom[i],]$NORMPOS <- STATS[STATS$Species=="C. remanei" & STATS$chrom==chrom[i],]$classifiedWinEnd/endsCR[i]

}

# the central domains
domains<-data.frame(classifiedWinEnd=c(3858000,11040000,4879000,12020000,3722000,10340000,3896000,12970000,5897000,16550000,6137000,12480000,5803000,13007000,6828000, 14292000,5849000,11888000,8639000,17108000,6767000,13852000,9040000,16070000), chrom=rep(c("I","I","II","II","III","III","IV","IV","V","V","X","X"),2), Species=c(rep("C. elegans",12),rep("C. remanei",12)))


domains$NORMPOS<-888

for (i in 1:6){
  domains[domains$Species=="C. elegans" & domains$chrom==chrom[i],]$NORMPOS <- domains[domains$Species=="C. elegans" & domains$chrom==chrom[i],]$classifiedWinEnd/endsCE[i]
  domains[domains$Species=="C. remanei" & domains$chrom==chrom[i],]$NORMPOS <- domains[domains$Species=="C. remanei" & domains$chrom==chrom[i],]$classifiedWinEnd/endsCR[i]

}


ggplot(STATS, aes(x = NORMPOS, y = pi, ordered = FALSE))  +
  labs(title ="A", x = "Normalized genome position", y = bquote("Nucleotide diversity ("~pi~")")) +
  scale_y_continuous(labels = function(x) format(x, scientific = FALSE)) +
  scale_x_continuous(breaks=c(0,0.5,1),labels=c(0,0.5,1)) +
  theme_bw(base_size = 12) +
  theme(axis.text.x = element_text(angle=0,vjust=0.5, hjust = 0.5), plot.title = element_text(hjust = 0, face = "bold"),strip.background = element_blank(),panel.grid.major = element_blank(), panel.grid.minor = element_blank(), strip.text = element_text(size = rel(1))) +
  theme(legend.position = "none") +
  theme(strip.text.y = element_text(face = "italic"))+
  theme(strip.background = element_blank(),panel.grid.major = element_blank(), panel.grid.minor = element_blank(), strip.text = element_text(size = rel(1))) +
  facet_grid(Species ~ chrom ) +
  geom_vline(data = domains,aes(xintercept = NORMPOS),colour = "#bfbfbf",size = 0.75,linetype = "dashed") +
  geom_point(alpha=0.2,size=0.6, col="#196770") +
  geom_smooth(col="#196770",alpha=1,size=1.1,se =FALSE,method = 'loess',span=0.4) + theme(plot.margin = unit(c(1, 3, 1.5, 3), "pt")) -> F1A


##############exons & introns #######################

CEEX<-read.csv("CE_exons_pi_100kb_FIN.bed",header=F,sep="\t")
CEIN<-read.csv("CE_introns_pi_100kb_FIN.bed",header=F,sep="\t")

CEEXC<-read.csv("CE_Exons.COVERAGE_100Kb.bed",header=F,sep="\t")
CEINC<-read.csv("CE_Introns.COVERAGE_100Kb.bed",header=F,sep="\t")

CREX<-read.csv("CR_exons_pi_100kb_FIN.bed",header=F,sep="\t")
CRIN<-read.csv("CR_introns_pi_100kb_FIN.bed",header=F,sep="\t")

CREXC<-read.csv("CR_Exons.COVERAGE_100Kb.bed",header=F,sep="\t")
CRINC<-read.csv("CR_Introns.COVERAGE_100Kb.bed",header=F,sep="\t")

CEEX<-cbind(CEEX,CEEXC)
CEIN<-cbind(CEIN,CEINC)

CREX<-cbind(CREX,CREXC)
CRIN<-cbind(CRIN,CRINC)

CEEX$Species<-"C. elegans"
CEEX$Type<-"Exon"

CEIN$Species<-"C. elegans"
CEIN$Type<-"Intron"

CREX$Species<-"C. remanei"
CREX$Type<-"Exon"

CRIN$Species<-"C. remanei"
CRIN$Type<-"Intron"


COMBO<-rbind(CEEX,CEIN)
COMBO<-rbind(COMBO,CREX)
COMBO<-rbind(COMBO,CRIN)

COMBO<-COMBO[,c(1:4,9,11,12,13)]
colnames(COMBO)<-c("CHR","POS","POS2","PISUM","COVERED","FRAC_COV","Species","Type")
#dots!
COMBO$PISUM<-gsub('\\.$', '0', perl=TRUE,COMBO$PISUM)

COMBO$PISUM<-as.numeric(as.character(COMBO$PISUM))
COMBO$PI<-COMBO$PISUM/COMBO$COVERED

colors<-c("orange","#008080", "#181C43", "#0E5EBC", "#73A8BD", "#C0BDBC", "#CF8971", "#A52224")
COMBO<-COMBO[COMBO$CHR!="MtDNA",]

COMBO$Type <- factor(COMBO$Type, levels = c("Exon","Intron"))

COMBO$NORMPOS<-888
colnames(COMBO)[1] <-"chrom"
for (i in 1:6){
  COMBO[COMBO$Species=="C. elegans" & COMBO$chrom==chrom[i],]$NORMPOS <- COMBO[COMBO$Species=="C. elegans" & COMBO$chrom==chrom[i],]$POS/endsCE[i]
  COMBO[COMBO$Species=="C. remanei" & COMBO$chrom==chrom[i],]$NORMPOS <- COMBO[COMBO$Species=="C. remanei" & COMBO$chrom==chrom[i],]$POS/endsCR[i]

}

###figure 1 B
colors <- c("#E6B243","#808080")

COMBO<-as.data.frame(COMBO)
ggplot(COMBO[COMBO$FRAC_COV>0.05,], aes(x = NORMPOS, y = PI, color = Type, ordered = FALSE))  +
  labs(title ="B", x = "Normalized genome position", y = bquote("Nucleotide diversity ("~pi~")")) +
  scale_y_continuous(labels = function(x) format(x, scientific = FALSE),breaks=c(0,0.02,0.04)) +
  scale_x_continuous(breaks=c(0,0.5,1),labels=c(0,0.5,1)) +
  theme_bw(base_size = 12) +
  theme(axis.text.x = element_text(angle=0,vjust=0.5, hjust = 0.5), plot.title = element_text(hjust = 0, face = "bold"),strip.background = element_blank(),panel.grid.major = element_blank(), panel.grid.minor = element_blank(), strip.text = element_text(size = rel(1))) +
  theme(legend.position = c(0.89,1.3)) +
  theme(legend.key =element_blank(),legend.background=element_blank()) +
  theme(strip.text.y = element_text(face = "italic"))+
  theme(strip.background = element_blank(),panel.grid.major = element_blank(), panel.grid.minor = element_blank(), strip.text = element_text(size = rel(1))) +
  facet_grid(Species ~ chrom ) +
  geom_vline(data = domains,aes(xintercept = NORMPOS),colour = "#bfbfbf",size = 0.75,linetype = "dashed") +
  geom_point(alpha=0.2,size=0.6) +
  geom_smooth(alpha=0.35,size=1.1,se =FALSE,method = 'loess',span=0.4) +
  scale_color_manual(values = colors, name="") +
  guides(colour = guide_legend(nrow = 1)) + theme(plot.margin = unit(c(1, 3, 1.5, 3), "pt")) -> F1B


grid.arrange(F1A,F1B, nrow=2) ->F1

#ggsave(filename = "CR_CE_diversity_fig1.png",F1,width = 7.2,height = 5.8, units = "in", dpi=300,scale=1)
tiff(filename = "Fig2.tiff",width = 7.2,height = 5.8, res=300, units = "in")

plot(F1)
dev.off()



#############
#add  Arms/Centers

COMBO$CHRTYPE<-"Arm"
for (chr in 1:length(levels(as.factor(as.character(COMBO$chrom))))){
  COMBO[COMBO$Species=="C. elegans" & COMBO$chrom==levels(as.factor(as.character(COMBO$chrom)))[chr] & COMBO$POS>domains$classifiedWinEnd[chr*2-1] & COMBO$POS<domains$classifiedWinEnd[chr*2],]$CHRTYPE<-"Center"
  COMBO[COMBO$Species=="C. remanei" & COMBO$chrom==levels(as.factor(as.character(COMBO$chrom)))[chr] & COMBO$POS>domains$classifiedWinEnd[12+chr*2-1] & COMBO$POS<domains$classifiedWinEnd[12+chr*2],]$CHRTYPE<-"Center"

}



STATS$CHRTYPE<-"Arm"
for (chr in 1:length(levels(as.factor(as.character(STATS$chrom))))){
  STATS[STATS$Species=="C. elegans" & STATS$chrom==levels(as.factor(as.character(STATS$chrom)))[chr] & STATS$classifiedWinEnd>domains$classifiedWinEnd[chr*2-1] & STATS$classifiedWinEnd<domains$classifiedWinEnd[chr*2],]$CHRTYPE<-"Center"
  STATS[STATS$Species=="C. remanei" & STATS$chrom==levels(as.factor(as.character(STATS$chrom)))[chr] & STATS$classifiedWinEnd>domains$classifiedWinEnd[12+chr*2-1] & STATS$classifiedWinEnd<domains$classifiedWinEnd[12+chr*2],]$CHRTYPE<-"Center"

}


########################################
######## Estimate statistics ###########
########################################


STATS<-STATS[complete.cases(STATS),]
STATS$CHRTYPEN<-as.factor(STATS$CHRTYPE)
COMBO<-COMBO[complete.cases(COMBO),]
COMBO$CHRTYPEN<-as.factor(COMBO$CHRTYPE)
COMBO2<-COMBO[COMBO$FRAC_COV>0.05,]

###difference between domains in total pi

cohensD(pi ~CHRTYPE, data=STATS[STATS$Species=="C. elegans",], method="corrected")
#[1] 0.619952

cohensD(pi ~CHRTYPE, data=STATS[STATS$Species=="C. remanei",], method="corrected")
#[1] 1.083254

oneway_test(pi ~ CHRTYPEN, data=STATS[STATS$Species=="C. elegans",],alternative="greater",distribution=approximate(nresample=10000))

#Approximative Two-Sample Fisher-Pitman Permutation Test
#
#data:  pi by CHRTYPEN (Arm, Center)
#Z = 9.2777, p-value < 1e-04
#alternative hypothesis: true mu is greater than 0

oneway_test(pi ~ CHRTYPEN, data=STATS[STATS$Species=="C. remanei",],alternative="greater",distribution=approximate(nresample=10000))

#Approximative Two-Sample Fisher-Pitman Permutation Test
#
#data:  pi by CHRTYPEN (Arm, Center)
#Z = 13.958, p-value < 1e-04
#alternative hypothesis: true mu is greater than 0


#####difference between gene features
cohensD(PI ~Type, data=COMBO2[COMBO2$Species=="C. elegans" & COMBO2$CHRTYPE=="Arm",], method="corrected")
#[1] 0.06346452

cohensD(PI ~Type, data=COMBO2[COMBO2$Species=="C. elegans" & COMBO2$CHRTYPE=="Center",], method="corrected")
#[1] 0.03664263

cohensD(PI ~Type, data=COMBO2[COMBO2$Species=="C. remanei" & COMBO2$CHRTYPE=="Arm",], method="corrected")
#[1] 0.4320014

cohensD(PI ~Type, data=COMBO2[COMBO2$Species=="C. remanei" & COMBO2$CHRTYPE=="Center",], method="corrected")
#[1] 1.358204

cohensD(PI ~CHRTYPE, data=COMBO2[COMBO2$Species=="C. elegans" & COMBO2$Type=="Exon",], method="corrected")
#[1] 0.5182527

cohensD(PI ~CHRTYPE, data=COMBO2[COMBO2$Species=="C. elegans" & COMBO2$Type=="Intron",], method="corrected")
#[1] 0.6496155

cohensD(PI ~CHRTYPE, data=COMBO2[COMBO2$Species=="C. remanei" & COMBO2$Type=="Exon",], method="corrected")
#[1] 1.362185

cohensD(PI ~CHRTYPE, data=COMBO2[COMBO2$Species=="C. remanei" & COMBO2$Type=="Intron",], method="corrected")
#[1] 0.700596

oneway_test(PI ~ CHRTYPEN, data=COMBO2[COMBO2$Species=="C. elegans" & COMBO2$Type=="Exon",],alternative="greater",distribution=approximate(nresample=10000))

#Approximative Two-Sample Fisher-Pitman Permutation Test
#
#data:  PI by CHRTYPEN (Arm, Center)
#Z = 7.9395, p-value < 1e-04
#alternative hypothesis: true mu is greater than 0

oneway_test(PI ~ CHRTYPEN, data=COMBO2[COMBO2$Species=="C. elegans" & COMBO2$Type=="Intron",],alternative="greater",distribution=approximate(nresample=10000))

#Approximative Two-Sample Fisher-Pitman Permutation Test
#
#data:  PI by CHRTYPEN (Arm, Center)
#Z = 9.7283, p-value < 1e-04
#alternative hypothesis: true mu is greater than 0


oneway_test(PI ~ CHRTYPEN, data=COMBO2[COMBO2$Species=="C. remanei" & COMBO2$Type=="Exon",],alternative="greater",distribution=approximate(nresample=10000))

#Approximative Two-Sample Fisher-Pitman Permutation Test
#
#data:  PI by CHRTYPEN (Arm, Center)
#Z = 16.162, p-value < 1e-04
#alternative hypothesis: true mu is greater than 0


oneway_test(PI ~ CHRTYPEN, data=COMBO2[COMBO2$Species=="C. remanei" & COMBO2$Type=="Intron",],alternative="greater",distribution=approximate(nresample=10000))

#Approximative Two-Sample Fisher-Pitman Permutation Test
#
#data:  PI by CHRTYPEN (Arm, Center)
#Z = 8.3312, p-value < 1e-04
#alternative hypothesis: true mu is greater than 0

oneway_test(PI ~ Type, data=COMBO2[COMBO2$Species=="C. remanei" & COMBO2$CHRTYPE=="Center",],alternative="less",distribution=approximate(nresample=10000))

#Approximative Two-Sample Fisher-Pitman Permutation Test
#
#data:  PI by Type (Exon, Intron)
#Z = -15.862, p-value < 1e-04
#alternative hypothesis: true mu is less than 0


oneway_test(PI ~ Type, data=COMBO2[COMBO2$Species=="C. remanei" & COMBO2$CHRTYPE=="Arm",],alternative="less",distribution=approximate(nresample=10000))

#Approximative Two-Sample Fisher-Pitman Permutation Test
#
#data:  PI by Type (Exon, Intron)
#Z = -5.3878, p-value < 1e-04
#alternative hypothesis: true mu is less than 0


oneway_test(PI ~ Type, data=COMBO2[COMBO2$Species=="C. elegans" & COMBO2$CHRTYPE=="Arm",],alternative="less",distribution=approximate(nresample=10000))

#Approximative Two-Sample Fisher-Pitman Permutation Test
#
#data:  PI by Type (Exon, Intron)
#Z = 1.0309, p-value = 0.8456
#alternative hypothesis: true mu is less than 0


oneway_test(PI ~ Type, data=COMBO2[COMBO2$Species=="C. elegans" & COMBO2$CHRTYPE=="Center",],alternative="less",distribution=approximate(nresample=10000))

#Approximative Two-Sample Fisher-Pitman Permutation Test
#
#data:  PI by Type (Exon, Intron)
#Z = -0.56267, p-value = 0.2927
#alternative hypothesis: true mu is less than 0




#################################################
########### plot all other stats ################
#################################################

for (stat in 6:16){
  pl<- ggplot(STATS, aes(x = NORMPOS, y = STATS[,stat], ordered = FALSE))  +
    labs(title =LETTERS[stat-4], x = "Normalized genome position", y = ggplot2:::parse_safe(colnames(STATS)[stat])) +
    scale_y_continuous(labels = function(x) format(x, scientific = FALSE)) +
    scale_x_continuous(breaks=c(0,0.5,1),labels=c(0,0.5,1)) +
    theme_bw(base_size = 12) +
    theme(axis.text.x = element_text(angle=0,vjust=0.5, hjust = 0.5), plot.title = element_text(hjust = 0, face = "bold"),strip.background = element_blank(),panel.grid.major = element_blank(), panel.grid.minor = element_blank(), strip.text = element_text(size = rel(1))) +
    theme(legend.position = "none") +
    theme(strip.text.y = element_text(face = "italic"))+
    theme(strip.background = element_blank(),panel.grid.major = element_blank(), panel.grid.minor = element_blank(), strip.text = element_text(size = rel(1))) +
    facet_grid(Species ~ chrom ) +
    geom_vline(data = domains,aes(xintercept = NORMPOS),colour = "#bfbfbf",size = 0.75,linetype = "dashed") +
    geom_point(alpha=0.2,size=0.6, col="#196770") +
    geom_smooth(col="#196770",alpha=1,size=1.1,se =FALSE,method = 'loess',span=0.4) + theme(plot.margin = unit(c(1, 3, 1.5, 3), "pt"))

  #ggsave(paste0(colnames(STATS)[stat],".CE_CR_100Kb_Supl_norm.png"),pl,width = 7.2,height = 3.77,dpi=300,scale=1)
  tiff(paste0(colnames(STATS)[stat],".CE_CR_100Kb_Supl_norm.tiff"),width = 7.2,height = 3.77,res=300, units="in")
  plot(pl)
  dev.off()
}


###############################
#also BETA
#load
my_files<-list.files(pattern="*BETA")
my_names<-gsub("_5-100_0.5.100kb.BETA","", perl=T, my_files)

BETA<-c();for (i in 1:(length(my_files)) ) { B=read.csv(file=my_files[i], header=F, sep="\t"); BETA<-rbind(BETA,data.frame(B,sample=my_names[i]));}

BETA$Species<-gsub("CE","C. elegans", BETA$sample)
BETA$Species<-gsub("CR","C. remanei", BETA$Species)

colnames(BETA)<-c("chrom", "classifiedWinStart", "classifiedWinEnd", "beta","sample","Species")
BETA<-BETA[BETA$chrom!="MtDNA",]
BETA$classifiedWinEnd <- BETA$classifiedWinEnd-50000
BETA<-merge(BETA,STATS,by=c("chrom","classifiedWinEnd","Species"),all.y=TRUE)

pl<- ggplot(BETA, aes(x = NORMPOS, y = beta, ordered = FALSE))  +
  labs(title =LETTERS[13], x = "Normalized genome position", y = ggplot2:::parse_safe("beta")) +
  scale_y_continuous(labels = function(x) format(x, scientific = FALSE)) +
  scale_x_continuous(breaks=c(0,0.5,1),labels=c(0,0.5,1)) +
  theme_bw(base_size = 12) +
  theme(axis.text.x = element_text(angle=0,vjust=0.5, hjust = 0.5), plot.title = element_text(hjust = 0, face = "bold"),strip.background = element_blank(),panel.grid.major = element_blank(), panel.grid.minor = element_blank(), strip.text = element_text(size = rel(1))) +
  theme(legend.position = "none") +
  theme(strip.text.y = element_text(face = "italic"))+
  theme(strip.background = element_blank(),panel.grid.major = element_blank(), panel.grid.minor = element_blank(), strip.text = element_text(size = rel(1))) +
  facet_grid(Species ~ chrom ) +
  geom_vline(data = domains,aes(xintercept = NORMPOS),colour = "#bfbfbf",size = 0.75,linetype = "dashed") +
  geom_point(alpha=0.2,size=0.6, col="#196770") +
  geom_smooth(col="#196770",alpha=1,size=1.1,se =FALSE,method = 'loess',span=0.4) + theme(plot.margin = unit(c(1, 3, 1.5, 3), "pt"))


#ggsave("beta.CE_CR_100Kb_Supl_norm.png",pl,width = 7.2,height = 3.77,dpi=300,scale=1)

tiff("beta.CE_CR_100Kb_Supl_norm.tiff",width = 7.2,height = 3.77,res=300,units="in")
plot(pl)
dev.off()



############Fis
my_files<-list.files(pattern="WILD_FIS_MAF_100_FILTER.txt")
my_names<-gsub("_WILD_FIS_MAF_100_FILTER.txt","", perl=T, my_files)

FIS<-c();
for (i in 1:(length(my_files)) ) {
  B=read.csv(file=my_files[i], header=F, sep="\t");
  FIS<-rbind(FIS,data.frame(B,sample=my_names[i]));
}

colnames(FIS)<-c("chrom","start","end","Fis","Species")
FIS$Species<-gsub("CE","C. elegans", FIS$Species)
FIS$Species<-gsub("CR","C. remanei", FIS$Species)

FIS<-FIS[complete.cases(FIS),]
FIS$NORMPOS<-888

for (i in 1:6){

  FIS[FIS$Species=="C. elegans" & FIS$chrom==chrom[i],]$NORMPOS <- FIS[FIS$Species=="C. elegans" & FIS$chrom==chrom[i],]$start/endsCE[i]
  FIS[FIS$Species=="C. remanei" & FIS$chrom==chrom[i],]$NORMPOS <- FIS[FIS$Species=="C. remanei" & FIS$chrom==chrom[i],]$start/endsCR[i]

}

setwd("~/Documents/Phillips_lab/drafts/CR_popgen/CR_CE_pi_draft/NEW_HARD_MASK_CE_CR/")


pl<- ggplot(FIS, aes(x = NORMPOS, y = Fis, ordered = FALSE))  +
  labs(title =LETTERS[14], x = "Normalized genome position", y = ggplot2:::parse_safe("F[IS]")) +
  scale_y_continuous(labels = function(x) format(x, scientific = FALSE)) +
  scale_x_continuous(breaks=c(0,0.5,1),labels=c(0,0.5,1)) +
  theme_bw(base_size = 12) +
  theme(axis.text.x = element_text(angle=0,vjust=0.5, hjust = 0.5), plot.title = element_text(hjust = 0, face = "bold"),strip.background = element_blank(),panel.grid.major = element_blank(), panel.grid.minor = element_blank(), strip.text = element_text(size = rel(1))) +
  theme(legend.position = "none") +
  theme(strip.text.y = element_text(face = "italic"))+
  theme(strip.background = element_blank(),panel.grid.major = element_blank(), panel.grid.minor = element_blank(), strip.text = element_text(size = rel(1))) +
  facet_grid(Species ~ chrom ) +
  geom_vline(data = domains,aes(xintercept = NORMPOS),colour = "#bfbfbf",size = 0.75,linetype = "dashed") +
  geom_point(alpha=0.2,size=0.6, col="#196770") +
  geom_smooth(col="#196770",alpha=1,size=1.1,se =FALSE,method = 'loess',span=0.4) + theme(plot.margin = unit(c(1, 3, 1.5, 3), "pt"))

#ggsave("FIS_CR_CE_100Kb_norm.png",pl,width = 7.2,height = 3.77,dpi=300,scale=1)

tiff("FIS_CR_CE_100Kb_Supl_norm.tiff",width = 7.2,height = 3.77,res=300,units="in")
plot(pl)
dev.off()

######## Combine stats to 1 pdf


TIF <- file.info(list.files(pattern="*_Supl_norm.tiff"))
TIF <- TIF[with(TIF, order(as.POSIXct(mtime))), ]
#rownames(TIF)
TIF<-rownames(TIF)
#[-1]

#!!!
combo<-image_append(do.call("c", lapply(TIF, function(h){image_read(h)})), stack = TRUE)
image_write(combo, path = "S2.tiff", format = "tiff", quality = 100)


########################################
mean(FIS[FIS$Species=="C. elegans",]$Fis)
#[1] 0.8582874

sd(FIS[FIS$Species=="C. elegans",]$Fis)
#[1] 0.2247299

mean(FIS[FIS$Species=="C. remanei",]$Fis)
#[1] 0.3785114

sd(FIS[FIS$Species=="C. remanei",]$Fis)
#[1] 0.1517498


########## mean and sd of pi
###########################
###C. elegans
mean(STATS[STATS$Species=="C. elegans" & STATS$CHRTYPE=="Center",]$pi)
#[1] 0.000541659

sd(STATS[STATS$Species=="C. elegans" & STATS$CHRTYPE=="Center",]$pi)
#[1] 0.001224117

mean(STATS[STATS$Species=="C. elegans" & STATS$CHRTYPE=="Arm",]$pi)
#[1] 0.001763041

sd(STATS[STATS$Species=="C. elegans" & STATS$CHRTYPE=="Arm",]$pi)
#[1] 0.00243815

########################
###C. remanei
mean(STATS[STATS$Species=="C. remanei" & STATS$CHRTYPE=="Center",]$pi)
#[1] 0.0123735


#sd(STATS[STATS$Species=="C. remanei" & STATS$CHRTYPE=="Center",]$pi)
#[1] 0.004230264

mean(STATS[STATS$Species=="C. remanei" & STATS$CHRTYPE=="Arm",]$pi)
#[1] 0.01732536

sd(STATS[STATS$Species=="C. remanei" & STATS$CHRTYPE=="Arm",]$pi)
#[1] 0.004860457

###total
mean(STATS[STATS$Species=="C. elegans",]$pi)
#[1] 0.001194552

sd(STATS[STATS$Species=="C. elegans",]$pi)
#[1] 0.002059868


mean(STATS[STATS$Species=="C. remanei",]$pi)
#[1] 0.014936

sd(STATS[STATS$Species=="C. remanei",]$pi)
#[1] 0.005192832

###total
mean(STATS[STATS$Species=="C. elegans",]$theta)
#[1] 0.0006267073

sd(STATS[STATS$Species=="C. elegans",]$theta)
#[1] 0.001070122


mean(STATS[STATS$Species=="C. remanei",]$theta)
#[1] 0.01331987

sd(STATS[STATS$Species=="C. remanei",]$theta)
#[1] 0.00458332

###################
#Population size
mean(sapply(STATS[STATS$Species=="C. elegans",]$theta, function(th) { return(th/(2*4*2.3e-09))}))
#[1] 34060.18
sd(sapply(STATS[STATS$Species=="C. elegans",]$theta, function(th) { return(th/(2*4*2.3e-09))}))
#[1] 58158.8

#library(boot)
popsize <- function(data, indices) {
  th <- mean(data[indices])
  return(th/(4*2.3e-09))
}
bootNeCe <- boot(STATS[STATS$Species=="C. elegans",]$theta*0.5,popsize, R = 1000)
bciNeCe<-boot.ci(boot.out = bootNeCe, conf = 0.95, type = "basic")
bciNeCe
#BOOTSTRAP CONFIDENCE INTERVAL CALCULATIONS
#Based on 1000 bootstrap replicates
#
#CALL :
#  boot.ci(boot.out = bootNeCe, conf = 0.95, type = "basic")
#
#Intervals :
#  Level      Basic
#95%   (30361, 37458 )
#Calculations and Intervals on Original Scale


mean(sapply(STATS[STATS$Species=="C. remanei",]$theta, function(th) { return(th/(4*2.3e-09))}))
#[1] 1447812
sd(sapply(STATS[STATS$Species=="C. remanei",]$theta, function(th) { return(th/(4*2.3e-09))}))
#[1] 498187

bootNeCr <- boot(STATS[STATS$Species=="C. remanei",]$theta,popsize, R = 1000)
bciNeCr<-boot.ci(boot.out = bootNeCr, conf = 0.95, type = "basic")
bciNeCr
#BOOTSTRAP CONFIDENCE INTERVAL CALCULATIONS
#Based on 1000 bootstrap replicates
#
#CALL :
#  boot.ci(boot.out = bootNeCr, conf = 0.95, type = "basic")
#
#Intervals :
#  Level      Basic
#95%   (1415907, 1481330 )
#Calculations and Intervals on Original Scale

######
