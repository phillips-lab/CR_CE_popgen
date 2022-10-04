#######################################################
####### Divergence, diversity, substitutions  #########
##################  in C.remanei   ####################
#######################################################


#https://www.megasoftware.net/mega1_manual/Distance.html
TAMURA92<-function(TS,TV,N,GC){

  P=TS/N;
  Q=TV/N;
  d = -2*GC*(1 - GC)*log(1 - P/(2*GC*(1 - GC)) - Q) -(1 -2*GC*(1 -GC))*log(1 - 2*Q)/2;

  return(d);

}

PDIST<-function(TS,TV,N){

  d=(TS+TV)/N;
  return(d);

}


setwd("/Users/anastasia/Documents/Phillips_lab/drafts/CR_popgen/Ancestral/dir_100kb/")

library(ggplot2)
#devtools::install_github("kwstat/pals")
library(pals)
library(dplyr)
library(ggpubr)
library(gridExtra)

library(data.table)
library(lsr)
library(coin)

data<-read.csv("CHR_POS_PX506_ANC0_100kb.txt",header=F,sep="\t")
colnames(data)<-c("Count","CHR","POS","PX506","Anc0")
data$CHR<-gsub("CM021144.1","I",data$CHR)
data$CHR<-gsub("CM021145.1","II",data$CHR)
data$CHR<-gsub("CM021146.1","III",data$CHR)
data$CHR<-gsub("CM021147.1","IV",data$CHR)
data$CHR<-gsub("CM021148.1","V",data$CHR)
data$CHR<-gsub("CM021149.1","X",data$CHR)


sum(data$Count)
#[1] 3553584 #number of polymorphic positions

data$replacement<-paste0(data$Anc0,"->",data$PX506)
data$group<-data$replacement
data$group<-gsub("T->A","A->T",data$group)
data$group<-gsub("T->G","A->C",data$group)
data$group<-gsub("T->C","A->G",data$group)
data$group<-gsub("G->T","C->A",data$group)
data$group<-gsub("G->A","C->T",data$group)
data$group<-gsub("G->C","C->G",data$group)
data$group<-gsub("A->T","A→T|T→A",data$group)
data$group<-gsub("A->C","A→C|T→G",data$group)
data$group<-gsub("A->G","A→G|T→C",data$group)
data$group<-gsub("C->A","C→A|G→T",data$group)
data$group<-gsub("C->T","C→T|G→A",data$group)
data$group<-gsub("C->G","C→G|G→C",data$group)

data$TYPE<-data$group
data$TYPE<-gsub("A→T\\|T→A|A→C\\|T→G|C→A\\|G→T|C→G\\|G→C|C→A\\|G→T","TV",perl=TRUE,data$TYPE)
data$TYPE<-gsub("A→G\\|T→C|C→T\\|G→A","TS",perl=TRUE,data$TYPE)

bed<-read.csv("Ancestral_states_coverage_100kb.bed",sep="\t",header=F)

#### some stats on coverage
sum(bed$V5)/sum(bed$V6)*100
#[1] 56.34172  #the number of position with ancestral states

for(ch in levels(as.factor(bed$V1))){print(sum(bed[bed$V1==ch,]$V5)/sum(bed[bed$V1==ch,]$V6)*100);}
#[1] 59.88874
#[1] 55.70638
#[1] 57.41857
#[1] 41.89103
#[1] 58.98964
#[1] 67.75234



bed$V1<-gsub("CM021144.1","I",bed$V1)
bed$V1<-gsub("CM021145.1","II",bed$V1)
bed$V1<-gsub("CM021146.1","III",bed$V1)
bed$V1<-gsub("CM021147.1","IV",bed$V1)
bed$V1<-gsub("CM021148.1","V",bed$V1)
bed$V1<-gsub("CM021149.1","X",bed$V1)
bed$POS<-as.integer(bed$V2/100000)
BED<-bed[,c("V1","POS","V5")]
colnames(BED)<-c("CHR","POS","COVERAGE")


#percent of polymorphic positions from covered
var<-c(); for(ch in levels(as.factor(data$CHR))){
  print(sum(data[data$CHR==ch,]$Count)/sum(bed[bed$V1==ch,]$V5)*100);
  var<-c(var,(sum(data[data$CHR==ch,]$Count)/sum(bed[bed$V1==ch,]$V5)*100));
}
#[1] 4.968421
#[1] 5.014295
#[1] 5.114169
#[1] 5.445749
#[1] 4.702497
#[1] 5.12051

mean(var)
#[1] 5.06094
sd(var)
#[1] 0.2423538


GC<-read.csv("ANCESTRAL_GC.txt",header=F,sep="\t")
colnames(GC)<-c("CHR","POS","XXX","GC")
GC<-GC[,-3]
GC$POS<-GC$POS/100000


GC$CHR<-gsub("CM021144.1","I",GC$CHR)
GC$CHR<-gsub("CM021145.1","II",GC$CHR)
GC$CHR<-gsub("CM021146.1","III",GC$CHR)
GC$CHR<-gsub("CM021147.1","IV",GC$CHR)
GC$CHR<-gsub("CM021148.1","V",GC$CHR)
GC$CHR<-gsub("CM021149.1","X",GC$CHR)

DATA<-merge(data,BED,by=c("CHR","POS"),all.x=TRUE,ordered=FALSE)
DATA$Fraction<-DATA$Count/DATA$COVERAGE*100
DATA<-merge(DATA,GC,by=c("CHR","POS"),all.x=TRUE,ordered=FALSE)

DATA$FractionGC<-4

DATA[DATA$Anc0 %in% c("G", "C"), ]$FractionGC <-
  (DATA[DATA$Anc0 %in% c("G", "C"), ]$Count / (DATA[DATA$Anc0 %in% c("G", "C"), ]$COVERAGE *
                                                 DATA[DATA$Anc0 %in% c("G", "C"), ]$GC)) * 100
DATA[DATA$Anc0 %in% c("A", "T"), ]$FractionGC <-
  (DATA[DATA$Anc0 %in% c("A", "T"), ]$Count / (DATA[DATA$Anc0 %in% c("A", "T"), ]$COVERAGE *(1 - DATA[DATA$Anc0 %in% c("A", "T"), ]$GC))) * 100


DATA %>%
  group_by(group,CHR,POS)  %>%
  summarise(Fraction = sum(FractionGC))-> A

df<-as.data.frame(A)
#add asterics to the groups that show significant changes between domains
df$group <-gsub("C→T.G→A","C→T|G→A*", perl=T, df$group)
df$group <-gsub("A→T.T→A","A→T|T→A*", perl=T, df$group)
df$group <-gsub("C→G.G→C","C→G|G→C*", perl=T,df$group)
df$group <- factor(df$group, levels = c("C→T|G→A*", "A→G|T→C","C→A|G→T", "A→T|T→A*", "C→G|G→C*", "A→C|T→G"))

#domainsCR<-data.frame(classifiedWinEnd=c(5382370,12864400,5953500,14377500,5403040,11895800,8164350,20180000,5528390,15095000,8016160,15396200), chrom=c("I","I","II","II","III","III","IV","IV","V","V","X","X"), Species=rep("C. remanei",12))
domainsCR<-data.frame(classifiedWinEnd=c(5803000,13007000,6828000, 14292000,5849000,11888000,8639000,17108000,6767000,13852000,9040000,16070000), chrom=c("I","I","II","II","III","III","IV","IV","V","V","X","X"), Species=rep("C. remanei",12))

domainsCR$POS<-as.integer(domainsCR$classifiedWinEnd/100000)

domainsCR2<-domainsCR
colnames(domainsCR2)[2]<-"CHR"

#FigS3
colorsMUT<-c("#183843", "#196770", "#73A8BD", "#C0BDBC", "#f2c76f", "#916b20")


ggplot(df, aes(x = POS/10, y = Fraction, color= group,ordered = FALSE))  +
  labs(title ="B", x = "Genome position (Mb)", y = "Substitutions\n from ancestral (%)") +
  scale_y_continuous(labels = function(x) format(x, scientific = FALSE))  +
  scale_x_continuous(breaks = seq(0, 100, by = 5)) +
  theme_bw(base_size = 12) +
  theme(axis.text.x = element_text(angle=0,vjust=0.5, hjust = 0.5), plot.title = element_text(hjust = 0, face = "bold"),strip.background = element_blank(),panel.grid.major = element_blank(), panel.grid.minor = element_blank(), strip.text = element_text(size = rel(1))) +
  theme(legend.position = "right") +
  theme(strip.text.y = element_text(face = "italic"))+
  theme(strip.background = element_blank(),panel.grid.major = element_blank(), panel.grid.minor = element_blank(), strip.text = element_text(size = rel(1))) +
  facet_grid(.~ CHR, scale="free",space="free" ) +
  geom_vline(data = domainsCR2,aes(xintercept = POS/10),colour = "#bfbfbf",size = 0.75,linetype = "dashed") +
  geom_smooth(alpha=1,size=1.1,se =FALSE,method = 'loess',span=0.4) + theme(plot.margin = unit(c(1, 3, 1.5, 6), "pt")) + scale_color_manual(values=colorsMUT, name="") ->F3


ggplot(df, aes(x = POS/10, y = Fraction, color= group,fill=group,ordered = FALSE))  +
  labs(title ="A", x = "Genome position (Mb)", y = "Fraction of substitutions\n from ancestral (%)") +
  scale_y_continuous(labels = function(x) format(x*100, scientific = FALSE))  +
  scale_x_continuous(breaks = seq(0, 100, by = 5)) +
  theme_bw(base_size = 12) +
  theme(axis.text.x = element_text(angle=0,vjust=0.5, hjust = 0.5), plot.title = element_text(hjust = 0, face = "bold"),strip.background = element_blank(),panel.grid.major = element_blank(), panel.grid.minor = element_blank(), strip.text = element_text(size = rel(1))) +
  theme(legend.position = "right") +
  theme(strip.text.y = element_text(face = "italic"))+
  theme(strip.background = element_blank(),panel.grid.major = element_blank(), panel.grid.minor = element_blank(), strip.text = element_text(size = rel(1))) +
  facet_grid(.~ CHR, scale="free",space="free" ) +
  theme(plot.margin = unit(c(1, 3, 1.5, 3), "pt")) +
  scale_color_manual(values=colorsMUT, name="") +
  scale_fill_manual(values=colorsMUT, name="") +
  geom_bar(stat="identity",position="fill") +
  geom_vline(data = domainsCR2,aes(xintercept = POS/10),colour = "#bfbfbf",size = 0.75,linetype = "dashed") ->F3S




ggarrange(F3S, F3, ncol=1, nrow=2, common.legend = TRUE, legend="right") ->FIGS3


#tiff(filename = "S3.tiff",width = 7.2,height = 5,res=300, units="in")
#plot(FIGS3)
#dev.off()

df$group<-gsub("[*]","", df$group)
df$group <- factor(df$group, levels = c("C→T|G→A", "A→G|T→C","C→A|G→T", "A→T|T→A", "C→G|G→C", "A→C|T→G"))



########TAMURA##########


DATA %>%
  group_by(CHR,POS,TYPE,GC,COVERAGE)  %>%
  summarise(count = sum(Count))-> B

df2<-as.data.frame(B)



TMP1<-df2[df2$TYPE=="TS",]
TMP2<-df2[df2$TYPE=="TV",]
TMP1<-TMP1[,-3]
TMP2<-TMP2[,-3]


colnames(TMP1)[5]<-"TS"
colnames(TMP2)[5]<-"TV"

sum(TMP1[TMP1$COVERAGE>30000,]$TS)/sum(TMP2[TMP2$COVERAGE>30000,]$TV)
#[1] 1.1665


DIST<-merge(TMP1,TMP2,by=c("CHR","POS","GC","COVERAGE"),all.x=TRUE,ordered=FALSE)
colnames(DIST)[4]<-"N"
RES<-c();for(i in 1:nrow(DIST)){RES<-rbind(RES,data.frame(CHR=DIST$CHR[i],POS=DIST$POS[i],TAMURA=TAMURA92(DIST$TS[i],DIST$TV[i],DIST$N[i],DIST$GC[i])));}


###########################
####### add pi ############
###########################


setwd("/Users/anastasia/Documents/Phillips_lab/drafts/CR_popgen/CR_CE_pi_draft/NEW_HARD_MASK_CE_CR/")


pi_files<-list.files(pattern="CR_WILD_population14_filt_snps_5-100_0.5_fin.NOPHASE*")
pi_names<-gsub(".100K.0.1.stats","", perl=T, pi_files)

STATS<-c();for (i in 1:(length(pi_files)) ) { B=read.csv(file=pi_files[i], header=T, sep="\t"); STATS<-rbind(STATS,data.frame(B,sample=pi_names[i]));}
STATS$Species<-gsub("^([CER]+)_.*","\\1",perl=T,STATS$sample)

STATS$Species<-gsub("CR","C. remanei", STATS$Species)
colnames(STATS)<-gsub("_win0","",colnames(STATS))

STATS$CHR<-STATS$chrom

STATS$POS<-as.integer(STATS$classifiedWinEnd/100000)

STATS<-merge(STATS,RES,by=c("CHR","POS"), all.x=TRUE,ordered=FALSE)
STATS<-merge(STATS,DATA,by=c("CHR","POS"),all.x=TRUE,ordered=FALSE)
STATS2<-STATS[STATS$COVERAGE>30000,]  #########!

STATS2$CHRTYPE<-"Arm"
for (chr in 1:length(levels(as.factor(as.character(STATS2$chrom))))){
  STATS2[STATS2$chrom==levels(as.factor(as.character(STATS2$chrom)))[chr] & STATS2$POS>domainsCR$classifiedWinEnd[chr*2-1]/100000 & STATS2$POS<domainsCR$classifiedWinEnd[chr*2]/100000,]$CHRTYPE<-"Center"

}

#difference in divergence
mean(STATS2[STATS2$CHRTYPE=="Arm", ]$TAMURA)/mean(STATS2[STATS2$CHRTYPE=="Center", ]$TAMURA)
#[1] 1.554929




##########################
##### final figure 4 #####
##########################

ggplot(STATS2, aes(x = POS/10, y = TAMURA, ordered = FALSE))  +
  labs(title ="A", x = "Genome position (Mb)", y = "Divergence") +
  scale_y_continuous(labels = function(x) format(x, scientific = FALSE)) +
  scale_x_continuous(breaks = seq(0, 100, by = 5)) +
  theme_bw(base_size = 12) +
  theme(axis.text.x = element_text(angle=0,vjust=0.5, hjust = 0.5), plot.title = element_text(hjust = 0, face = "bold"),strip.background = element_blank(),panel.grid.major = element_blank(), panel.grid.minor = element_blank(), strip.text = element_text(size = rel(1))) +
  theme(legend.position = "none") +
  theme(strip.text.y = element_text(face = "italic"))+
  theme(strip.background = element_blank(),panel.grid.major = element_blank(), panel.grid.minor = element_blank(), strip.text = element_text(size = rel(1))) +
  facet_grid(.~ chrom, scale="free",space="free" ) +
  geom_vline(data = domainsCR,aes(xintercept = POS/10),colour = "#bfbfbf",size = 0.75,linetype = "dashed") +
  geom_point(alpha=0.10,size=0.6, col="#183843") +
  geom_smooth(col="#183843",alpha=1,size=1.1,se =FALSE,method = 'loess',span=0.4) + theme(plot.margin = unit(c(1, 3, 1.5, 3), "pt")) ->F4A


#  coord_cartesian(ylim=c(0,3)) +
ggplot(STATS2, aes(x = POS/10, y = pi/TAMURA, ordered = FALSE))  +
  labs(title ="B", x = "Genome position (Mb)", y = "Diversity(theta)/Divergence") +
  scale_y_continuous(labels = function(x) format(x, scientific = FALSE)) +
  scale_x_continuous(breaks = seq(0, 100, by = 5)) +
  theme_bw(base_size = 12) +
  theme(axis.text.x = element_text(angle=0,vjust=0.5, hjust = 0.5), plot.title = element_text(hjust = 0, face = "bold"),strip.background = element_blank(),panel.grid.major = element_blank(), panel.grid.minor = element_blank(), strip.text = element_text(size = rel(1))) +
  theme(legend.position = "none") +
  theme(strip.text.y = element_text(face = "italic"))+
  theme(strip.background = element_blank(),panel.grid.major = element_blank(), panel.grid.minor = element_blank(), strip.text = element_text(size = rel(1))) +
  facet_grid(.~ chrom, scale="free",space="free" ) +
  geom_vline(data = domainsCR,aes(xintercept = POS/10),colour = "#bfbfbf",size = 0.75,linetype = "dashed") +
  geom_point(alpha=0.10,size=0.6, col="#73A8BD") +
  geom_smooth(col="#73A8BD",alpha=1,size=1.1,se =FALSE,method = 'loess',span=0.4) + theme(plot.margin = unit(c(1, 3, 1.5, 3), "pt")) ->F4B




grid.arrange(F4A,F4B, nrow=2) ->F4

setwd("/Users/anastasia/Documents/Phillips_lab/drafts/CR_popgen/Ancestral/dir_100kb/")
#ggsave(filename = "CR_CE_divergence_fig4.png",F4,width = 7.2,height = 5,dpi=300,scale=1)
tiff(filename = "Fig3.tiff",width = 7.2,height = 5,res=300,units="in")
plot(F4)
dev.off()

STATS2$ratio<-STATS2$TAMURA/STATS2$pi
STATS2$CHRTYPE<-factor(STATS2$CHRTYPE)
independence_test(ratio ~ CHRTYPE, data=STATS2, alternative="two.sided",distribution=approximate(nresample=10000))


#Approximative General Independence Test
#
#data:  ratio by CHRTYPE (Arm, Center)
#Z = 13.546, p-value < 1e-04
#alternative hypothesis: two.sided

cohensD(ratio ~ CHRTYPE, data=STATS2, method="corrected")
#[1] 0.2770436

####correlations


cor(STATS2[STATS2$CHRTYPE=="Center",]$TAMURA, STATS2[STATS2$CHRTYPE=="Center",]$pi)
#[1] 0.6059364
cor(STATS2[STATS2$CHRTYPE=="Arm",]$TAMURA, STATS2[STATS2$CHRTYPE=="Arm",]$pi)
#[1] 0.3818657

#####stats

#arm/center
df$CHRTYPE<-"Arm"
for (chr in 1:length(levels(as.factor(as.character(df$CHR))))){
  df[df$CHR==levels(as.factor(as.character(df$CHR)))[chr] & df$POS>domainsCR$classifiedWinEnd[chr*2-1]/1000000 & df$POS<domainsCR$classifiedWinEnd[chr*2]/1000000,]$CHRTYPE<-"Center"

}



df3<-df
setDT(df3)[, frac := Fraction / sum(Fraction), by=c("CHR","POS")]

df3<-data.frame(df3,stringsAsFactors = TRUE)
df3$CHRTYPE<-factor(df3$CHRTYPE)

independence_test(frac ~ group + CHRTYPE, data=df3, alternative="two.sided",distribution=approximate(nresample=10000))


#Approximative General Independence Test
#
#data:  frac by group, CHRTYPE
#maxT = 79.219, p-value < 1e-04
#alternative hypothesis: two.sided

mutstat<-c();
for (i in levels(as.factor(df3$group))){
  a<-df3[df3$group==i,];

  D<-cohensD(frac ~ CHRTYPE, data=a, method="corrected");

  it<-independence_test(frac ~ CHRTYPE, data=a, alternative="two.sided",distribution=approximate(nresample=10000));
  mutstat<-rbind(mutstat,data.frame(group=i, D=D,stat=statistic(it)[1], pval=pvalue(it)));
}

mutstat$adj<-p.adjust(mutstat$pval)

mutstat
#      group         D      stat   pval    adj
#1 C→T|G→A 0.5542107  3.634444 0.0004 0.0024
#2 A→G|T→C 0.3418894  2.249454 0.0272 0.0816
#3 C→A|G→T 0.3099842 -2.040271 0.0437 0.0816
#4 A→T|T→A 0.4218895 -2.772876 0.0069 0.0276
#5 C→G|G→C 0.4361383 -2.865920 0.0044 0.0220
#6 A→C|T→G 0.3286937 -2.162965 0.0319 0.0816

setDT(STATS2)[, frac := Fraction / sum(Fraction), by=c("CHR","POS")]


for (i in levels(as.factor(STATS2$group))){
  a<-STATS2[STATS2$group==i,];
  print(i);
  print(cor(STATS2[STATS2$group==i,]$TAMURA,STATS2[STATS2$group==i,]$frac));
  print(cor(STATS2[STATS2$group==i,]$pi,STATS2[STATS2$group==i,]$frac));
}

#[1] "A→C|T→G"
#[1] 0.2339641
#[1] 0.2920816
#[1] "A→G|T→C"
#[1] -0.3846189 ##!
#[1] -0.2968253
#[1] "A→T|T→A"
#[1] 0.1940253
#[1] 0.09296976
#[1] "C→A|G→T"
#[1] 0.1278696
#[1] 0.195708
#[1] "C→G|G→C"
#[1] 0.6152506 ###Wow!
#[1] 0.4428357
#[1] "C→T|G→A"
#[1] -0.3011542
#[1] -0.2919256

dataSUB<-read.csv("CR_polim_counts.txt",header=F,sep="\t")

head(dataSUB)
colnames(dataSUB)<-c("Count","CHR","POS","REF","ALT")

dataSUB$group<-paste0(dataSUB$REF,"/",dataSUB$ALT)
dataSUB$group<-gsub("^A\\/C$","A\\/C\\|T\\/G",perl=T,dataSUB$group)
dataSUB$group<-gsub("^A\\/G$","A\\/G\\|T\\/C",perl=T,dataSUB$group)
dataSUB$group<-gsub("^A\\/T$","A\\/T\\|T\\/A",perl=T,dataSUB$group)
dataSUB$group<-gsub("^C\\/A$","A\\/C\\|T\\/G",perl=T,dataSUB$group)
dataSUB$group<-gsub("^C\\/G$","C\\/G\\|G\\/C",perl=T,dataSUB$group)
dataSUB$group<-gsub("^C\\/T$","A\\/G\\|T\\/C",perl=T,dataSUB$group)
dataSUB$group<-gsub("^G\\/A$","A\\/G\\|T\\/C",perl=T,dataSUB$group)
dataSUB$group<-gsub("^G\\/C$","C\\/G\\|G\\/C",perl=T,dataSUB$group)
dataSUB$group<-gsub("^G\\/T$","A\\/C\\|T\\/G",perl=T,dataSUB$group)
dataSUB$group<-gsub("^T\\/A$","A\\/T\\|T\\/A",perl=T,dataSUB$group)
dataSUB$group<-gsub("^T\\/C$","A\\/G\\|T\\/C",perl=T,dataSUB$group)
dataSUB$group<-gsub("^T\\/G$","A\\/C\\|T\\/G",perl=T,dataSUB$group)


#"A/C" "A/G" "A/T" "C/A" "C/G" "C/T" "G/A" "G/C" "G/T" "T/A" "T/C" "T/G"

dataCOV<-read.csv("CR_fraction_of_masked.bed",header=F,sep="\t")

head(dataCOV)

dataCOV$Covered<-dataCOV$V6*(1-dataCOV$V7)
dataCOV$POS<- dataCOV$V2/100000
colnames(dataCOV)[1]<-"CHR"

dataSUB2<-merge(dataSUB,dataCOV[,c(1,8,9)],by=c("CHR","POS"),all.x=TRUE,ordered=FALSE)

dataSUB2$Freq<-dataSUB2$Count/dataSUB2$Covered

dataSUB2$group <- factor(dataSUB2$group, levels = c(
  "A/G|T/C", "A/T|T/A", "A/C|T/G", "C/G|G/C"))



colorsMUT3<-c("#183843", "#C0BDBC", "#73A8BD", "#f2c76f")

#expression(italic(E.coli)~Area~(mm^2)))
#bquote("Density of substitutions\n in the "~italic("C. remanei")~" population") )+

ggplot(dataSUB2, aes(x = POS/10, y = Freq*100, color= group, ordered = FALSE))  +
  labs(title ="C", x = "Genome position (Mb)", y = expression("Density of polymorphic\n sites in the population (%)"))+
  scale_y_continuous(labels = function(x) format(x, scientific = FALSE))  +
  scale_x_continuous(breaks = seq(0, 100, by = 5)) +
  theme_bw(base_size = 12) +
  theme(axis.text.x = element_text(angle=0,vjust=0.5, hjust = 0.5), plot.title = element_text(hjust = 0, face = "bold"),strip.background = element_blank(),panel.grid.major = element_blank(), panel.grid.minor = element_blank(), strip.text = element_text(size = rel(1))) +
  theme(legend.position = "right") +
  theme(strip.text.y = element_text(face = "italic"))+
  theme(strip.background = element_blank(),panel.grid.major = element_blank(), panel.grid.minor = element_blank(), strip.text = element_text(size = rel(1))) +
  facet_grid(.~ CHR, scale="free",space="free" ) +
  geom_vline(data = domainsCR2,aes(xintercept = POS/10),colour = "#bfbfbf",size = 0.75,linetype = "dashed") +
  geom_smooth(alpha=1,size=1.1,se =FALSE,method = 'loess',span=0.4) + theme(plot.margin = unit(c(1, 3, 1.5, 14), "pt")) + scale_color_manual(values=colorsMUT3, name="") ->F3P




ggarrange(FIGS3, F3P, ncol=1, nrow=2, heights=c(2,1)) ->FIGS3ABC

FIGS3ABC


tiff(filename = "S3.tiff",width = 7.2,height = 8,res=300, units="in")
plot(FIGS3ABC)
dev.off()
