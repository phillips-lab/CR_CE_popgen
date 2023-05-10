######################################################################################
################ Distribution of dominance and selection coefficients ################
######################################################################################
############### Fig S14,


setwd("/Users/anastasia/Documents/Phillips_lab/drafts/CR_popgen/CR_CE_pi_draft/SLiM/100kb_data/BEN/")

library(ggplot2)
library(dplyr)
library(coin)
library(ggpubr)


#expected number of mutations and their effects
Ne=5000
Gen=Ne*10
L=3000000
mu=2*1e-8

nummut=mu*L*Gen*Ne
#[1] 1.5e+23

#hom many deleterius
ndel<-nummut*0.1
#hom many dbeneficial
nben<-nummut*0.01

#S=3,7.5,15
weakdel<- -rgamma(ndel, 0.3, scale =3/Ne/0.3)*Ne
strdel<- -rgamma(ndel, 0.3, scale =15/Ne/0.3)*Ne
middel<- -rgamma(ndel, 0.3, scale =7.5/Ne/0.3)*Ne
midben<-rgamma(nben, 0.3, scale =7.5/Ne/0.3)*Ne
strben<-rgamma(nben, 0.3, scale =15/Ne/0.3)*Ne


### distibution of effect and dominance

domdel<-c(rbeta(3*ndel/4, 2, 6),runif(ndel/4, min = 0, max = 1))

domben<-c(rbeta(3*nben/4, 5, 5),runif(nben/4, min = 0, max = 1))


COMBO<- c()
COMBO<-data.frame(Dominance=c(sample(domdel,10000),sample(domdel,10000),sample(domdel,10000),sample(domben,10000),sample(domben,10000)), Selection=c(sample(weakdel,10000),sample(middel,10000),sample(strdel,10000),sample(midben,10000),sample(strben,10000)), Type=c(rep("weakly\ndeleterious",10000),rep("moderately\ndeleterious",10000),rep("strongly\ndeleteriuos",10000),rep("moderately\nbeneficial",10000),rep("strongly\nbeneficial",10000)), stringsAsFactors = FALSE)
COMBO$Type <- factor(COMBO$Type, levels = c("weakly\ndeleterious","moderately\ndeleterious", "strongly\ndeleteriuos","moderately\nbeneficial","strongly\nbeneficial"))

colorsMUT<-c("#183843", "#196770", "#73A8BD", "#C0BDBC", "#f2c76f", "#916b20")



COMBO$DomClass<-"OOO"

COMBO[COMBO$Dominance<=0.2,]$DomClass<- "h <= 0.2"
COMBO[COMBO$Dominance>0.2 & COMBO$Dominance <= 0.4,]$DomClass<- "0.2 < h <= 0.4"
COMBO[COMBO$Dominance>0.4 & COMBO$Dominance <= 0.6,]$DomClass<- "0.4 < h <= 0.6"
COMBO[COMBO$Dominance>0.6 & COMBO$Dominance <= 0.8,]$DomClass<- "0.6 < h <= 0.8"
COMBO[COMBO$Dominance>0.8 & COMBO$Dominance <= 1,]$DomClass<- "0.8 < h <= 1"

COMBO$SelClass<-"OOO"
COMBO$Selection<-as.numeric(as.character(COMBO$Selection))
COMBO[abs(COMBO$Selection)<=1,]$SelClass<- "|Ns| <= 1"
COMBO[abs(COMBO$Selection)>1 & abs(COMBO$Selection) <=10,]$SelClass<- "1 < |Ns| <= 10"
COMBO[abs(COMBO$Selection) >10,]$SelClass<- "|Ns| > 10"



COMBO %>%
  group_by(Type) %>%
  mutate(group_size = n()) %>%
  group_by(Selection, Dominance,SelClass,DomClass, Type, group_size) %>%
  summarise(perc = n()/max(group_size)*100) -> pld

pld$SelClass <- factor(pld$SelClass, levels = c("|Ns| <= 1","1 < |Ns| <= 10", "|Ns| > 10"))
pld$DomClass <- factor(pld$DomClass, levels = c("0.8 < h <= 1","0.6 < h <= 0.8", "0.4 < h <= 0.6", "0.2 < h <= 0.4","h <= 0.2"))


domselhist<- ggplot(pld, aes(x = SelClass,y=perc,fill=DomClass,ordered = FALSE))  +
  labs(title ="", x = "Selection", y='%') +
  theme_bw(base_size = 12) +
  theme(axis.text.x = element_text(angle=0,vjust=0.5, hjust = 0.5), plot.title = element_text(hjust = 0, face = "bold"),strip.background = element_blank(),panel.grid.major = element_blank(), panel.grid.minor = element_blank(), strip.text = element_text(size = rel(1))) +
  theme(legend.position = "right") +
  theme(strip.background = element_blank(),panel.grid.major = element_blank(), panel.grid.minor = element_blank(), strip.text = element_text(size = rel(1)))+
  theme(plot.margin = unit(c(1, 3, 1.5, 3), "pt")) +
  scale_color_manual(values=colorsMUT, name="Dominance") +
  scale_fill_manual(values=colorsMUT, name="Dominance") +
  geom_bar(stat="identity") +
  facet_grid(Type ~.,space="free" ) +
  theme(axis.text.x = element_text(angle = 0, vjust = 0.5, hjust=0.5))

#ggsave("SLiM_dom_sel_distributions_hist.png",domselhist, width =8,height = 6, scale=1)


tiff(filename = "S14.tiff", width =8,height = 6,res=300,units="in")
plot(domselhist)
dev.off()





#####################################
################ After selection ####
#####################################

#only the uniform mutation landscape
BEN<-read.csv("/Users/anastasia/Documents/Phillips_lab/drafts/CR_popgen/CR_CE_pi_draft/SLiM/100kb_data/BEN/SLIM_SEL_BENEFICIAL_PART.out",header=T,sep="\t")

colnames(BEN)<-c("Type", "Position","Dominance","AF","TT","PP","Age","Selection","sample")



BEN$scenario<-BEN$sample
BEN$scenario<-gsub("sim_[0-9]+_", "", perl=T, BEN$scenario)
BEN$scenario<-gsub(".ALLMUT.txt", "", perl=T, BEN$scenario)
BEN$Class<-gsub("2", "del_arm",BEN$Type)
BEN$Class<-gsub("3", "ben_arm",BEN$Class)
BEN$Class<-gsub("4", "del_cent",BEN$Class)
BEN$Class<-gsub("5", "ben_cent",BEN$Class)

BEN$mut<-gsub(".*_Mut_([12.15]+)_FrD.*", "\\1", perl=T, BEN$scenario)
BEN$sel<-gsub(".*(FrD_.*_FrB.*_.*_SDA_.*_SDC_.*_SB.*A_.*_SB.*C_.*)","\\1",perl=T,BEN$scenario)
BEN$self<-gsub(".*Self_(.*)_Mut_.*","\\1",perl=T, BEN$scenario)
BEN$self<-as.numeric(as.character(BEN$self))
BEN$self<-BEN$self*100

BEN2<-BEN[BEN$self==0,]

#colorsMUT<-c("#181C43", "#0E5EBC", "#73A8BD", "#C0BDBC", "#CF8971", "#A52224")
BEN2$Dominance<-as.numeric(as.character(BEN2$Dominance))
BEN2$Selection<-as.numeric(as.character(BEN2$Selection))
BEN2$Age<-as.numeric(as.character(BEN2$Age))
BEN2$Position<-as.numeric(as.character(BEN2$Position))
BEN2$AF<-as.numeric(as.character(BEN2$AF))


#only 1 scenario with the same ("Strongly *") distribution of selection coefficients for deleterious and beneficial mutations

BEN3<-BEN2[BEN2$scenario=="Ne_5000_Self_0_Mut_1_FrD_0.1_FrB_0.01_SDA_15_SDC_15_SBA_15_SBC_15",]


BEN3[BEN3$Class %in% c("ben_arm","ben_cent") & BEN3$AF>=0.5 & BEN3$Position<2000000 & BEN3$Age >40000,] %>% group_by(sample) %>% count(Class) %>% filter(Class=="ben_arm") -> A
BEN3[BEN3$Class %in% c("ben_arm","ben_cent") & BEN3$AF>=0.5 & BEN3$Position<2000000 & BEN3$Age >40000,] %>% group_by(sample) %>% count(Class) %>% filter(Class=="ben_cent") -> B
BEN3[BEN3$Class %in% c("ben_arm","ben_cent") & BEN3$AF>=0.5 & BEN3$Position<2000000 & BEN3$Age >40000,] %>% group_by(sample) %>% count(Class)  -> AL
AL<-as.data.frame(AL)
A<-merge(A,B, by="sample")
hist(A$n.x/A$n.y,breaks = 30)
table(A$n.x/A$n.y>1)
cohensD(x=A$n.x, y=A$n.y, method ="corrected")
#[1] 4.75894

AL$Class<-as.factor(AL$Class)
AL$n<-as.numeric(as.character(AL$n))
oneway_test(n ~ Class, data=AL,alternative="greater",distribution=approximate(nresample=10000))

#Approximative Two-Sample Fisher-Pitman Permutation Test
#
#data:  n by Class (ben_arm, ben_cent)
#Z = 13.022, p-value < 1e-04
#alternative hypothesis: true mu is greater than 0


#### AF < 0.5
BEN3[BEN3$Class %in% c("ben_arm","ben_cent") & BEN3$AF<0.5 & BEN3$Position<2000000 & BEN3$Age >40000,] %>% group_by(sample) %>% count(Class) %>% filter(Class=="ben_arm") -> A2
BEN3[BEN3$Class %in% c("ben_arm","ben_cent") & BEN3$AF<0.5 & BEN3$Position<2000000 & BEN3$Age >40000,] %>% group_by(sample) %>% count(Class) %>% filter(Class=="ben_cent") -> B2
BEN3[BEN3$Class %in% c("ben_arm","ben_cent") & BEN3$AF<0.5 & BEN3$Position<2000000 & BEN3$Age >40000,] %>% group_by(sample) %>% count(Class)  -> AL2
AL2<-as.data.frame(AL2)
A2<-merge(A2,B2, by="sample")
hist(A2$n.x/A2$n.y,breaks = 30)
table(A2$n.x/A2$n.y>1)
cohensD(x=A2$n.x, y=A2$n.y, method ="corrected")
#[1] 1.342378

AL2$Class<-as.factor(AL2$Class)
AL2$n<-as.numeric(as.character(AL2$n))
oneway_test(n ~ Class, data=AL2,alternative="greater",distribution=approximate(nresample=10000))

#Approximative Two-Sample Fisher-Pitman Permutation Test

#data:  n by Class (ben_arm, ben_cent)
#Z = 7.9095, p-value < 1e-04
#alternative hypothesis: true mu is greater than 0


BEN4<-BEN3[BEN3$Position<2000000 & BEN3$Age >40000,]

###NOW LET'S COMPARE WITH THE INITIAL DISTRIBUTION OF DOMINANCE AND SELECTION:
Ne=5000

NbenA<-nrow(BEN3[BEN3$Class=="ben_arm" & BEN3$AF>=0.5 & BEN3$Position<2000000 & BEN3$Age >40000,])
NbenC<-nrow(BEN3[BEN3$Class=="ben_cent" & BEN3$AF>=0.5 & BEN3$Position<2000000 & BEN3$Age >40000,])
NbenAm<-nrow(BEN3[BEN3$Class=="ben_arm" & BEN3$AF<0.5 & BEN3$Position<2000000 & BEN3$Age >40000,])
NbenCm<-nrow(BEN3[BEN3$Class=="ben_cent" & BEN3$AF<0.5 & BEN3$Position<2000000 & BEN3$Age >40000,])

NdelA<-nrow(BEN3[BEN3$Class=="del_arm" & BEN3$AF>=0.5 & BEN3$Position<2000000 & BEN3$Age >40000,])
NdelC<-nrow(BEN3[BEN3$Class=="del_cent" & BEN3$AF>=0.5 & BEN3$Position<2000000 & BEN3$Age >40000,])
NdelAm<-nrow(BEN3[BEN3$Class=="del_arm" & BEN3$AF<0.5 & BEN3$Position<2000000 & BEN3$Age >40000,])
NdelCm<-nrow(BEN3[BEN3$Class=="del_cent" & BEN3$AF<0.5 & BEN3$Position<2000000 & BEN3$Age >40000,])

strbenA<-rgamma(NbenA, 0.3, scale =15/Ne/0.3)
strbenC<-rgamma(NbenC, 0.3, scale =15/Ne/0.3)
dombenA<-c(rbeta(3*NbenA/4, 5, 5),runif(NbenA/4+1, min = 0, max = 1))
dombenC<-c(rbeta(3*NbenC/4, 5, 5),runif(NbenC/4+1, min = 0, max = 1))

strbenAm<-rgamma(NbenAm, 0.3, scale =15/Ne/0.3)
strbenCm<-rgamma(NbenCm, 0.3, scale =15/Ne/0.3)
dombenAm<-c(rbeta(3*NbenAm/4, 5, 5),runif(NbenAm/4, min = 0, max = 1))
dombenCm<-c(rbeta(3*NbenCm/4, 5, 5),runif(NbenCm/4+1, min = 0, max = 1))

#different dominance distribution
strdelA<- -rgamma(NdelA, 0.3, scale =15/Ne/0.3)
strdelC<- -rgamma(NdelC, 0.3, scale =15/Ne/0.3)
domdelA<-c(rbeta(3*NdelA/4, 2, 6),runif(NdelA/4, min = 0, max = 1))
domdelC<-c(rbeta(3*NdelC/4, 2, 6),runif(NdelC/4, min = 0, max = 1))

strdelAm<- -rgamma(NdelAm, 0.3, scale =15/Ne/0.3)
strdelCm<- -rgamma(NdelCm, 0.3, scale =15/Ne/0.3)
domdelAm<-c(rbeta(3*NdelAm/4, 2, 6),runif(NdelAm/4, min = 0, max = 1))
domdelCm<-c(rbeta(3*NdelCm/4, 2, 6),runif(NdelCm/4+1, min = 0, max = 1))

#combine all things
BEN4$AFClass<-"AA"
BEN4[BEN4$AF>=0.5,]$AFClass<-"AF >= 0.5"
BEN4[BEN4$AF<0.5,]$AFClass<-"AF < 0.5"
BEN4$SelectionType<-gsub("ben_.*","Beneficial",perl=TRUE,BEN4$Class)
BEN4$SelectionType<-gsub("del_.*","Deleterious",perl=TRUE,BEN4$SelectionType)
BEN4$Recomb<-gsub(".*_arm","Arm",perl=TRUE,BEN4$Class)
BEN4$Recomb<-gsub(".*_cent","Center",perl=TRUE,BEN4$Recomb)



COMBOEMP<-c()
COMBOEMP<-data.frame(Selection=c(strbenA,strbenC,strbenAm,strbenCm,strdelA,strdelC,strdelAm,strdelCm),
                  Dominance=c(dombenA,dombenC,dombenAm,dombenCm,domdelA,domdelC,domdelAm,domdelCm),
                  AFClass=rep("Initial",nrow(BEN4)),
                  Recomb=rep("Initial",nrow(BEN4)),
                  SelectionType=c(rep("Beneficial",nrow(BEN4[BEN4$Class %in% c("ben_arm","ben_cent"),])),rep("Deleterious",nrow(BEN4[BEN4$Class %in% c("del_arm","del_cent"),]))),
                  stringsAsFactors = FALSE)
COMBOEMP<-rbind(COMBOEMP,BEN4[,c("Selection", "Dominance", "AFClass",  "Recomb", "SelectionType")])


COMBOEMP$DomClass<-"OOO"

COMBOEMP[COMBOEMP$Dominance<=0.2,]$DomClass<- "h <= 0.2"
COMBOEMP[COMBOEMP$Dominance>0.2 & COMBOEMP$Dominance <= 0.4,]$DomClass<- "0.2 < h <= 0.4"
COMBOEMP[COMBOEMP$Dominance>0.4 & COMBOEMP$Dominance <= 0.6,]$DomClass<- "0.4 < h <= 0.6"
COMBOEMP[COMBOEMP$Dominance>0.6 & COMBOEMP$Dominance <= 0.8,]$DomClass<- "0.6 < h <= 0.8"
COMBOEMP[COMBOEMP$Dominance>0.8 & COMBOEMP$Dominance <= 1,]$DomClass<- "0.8 < h <= 1"

COMBOEMP$SelClass<-"OOO"
COMBOEMP$Selection<-as.numeric(as.character(COMBOEMP$Selection))
COMBOEMP[abs(COMBOEMP$Selection*Ne)<=1,]$SelClass<- "|Ns| <= 1"
COMBOEMP[abs(COMBOEMP$Selection*Ne)>1 & abs(COMBOEMP$Selection*Ne) <=10,]$SelClass<- "1 < |Ns| <= 10"
COMBOEMP[abs(COMBOEMP$Selection*Ne) >10,]$SelClass<- "|Ns| > 10"



COMBOEMP[COMBOEMP$AFClass!="AF < 0.5",] %>%
  group_by(SelectionType, Recomb) %>%
  mutate(group_size = n()) %>%
  group_by(SelectionType, Recomb,SelClass,DomClass,group_size) %>%
  summarise(perc = n()/max(group_size)*100) -> plot_data


as.data.frame(plot_data) %>% select (-c(SelClass, DomClass, perc)) %>%
  distinct() %>% filter(Recomb!= 'Initial') ->texttable

texttable$group_size<-gsub("^","n = ",perl=T, texttable$group_size)
texttable$Recomb <- factor(texttable$Recomb, levels=c("Initial","Arm","Center"))




plot_data$SelClass <- factor(plot_data$SelClass, levels = c("|Ns| <= 1","1 < |Ns| <= 10", "|Ns| > 10"))
plot_data$DomClass <- factor(plot_data$DomClass, levels = c("0.8 < h <= 1","0.6 < h <= 0.8", "0.4 < h <= 0.6", "0.2 < h <= 0.4","h <= 0.2"))


plot_data$Recomb <- factor(plot_data$Recomb, levels=c("Initial","Arm","Center"))


pl<-
  ggplot(plot_data, aes(x = SelClass,y=perc,fill=DomClass,ordered = FALSE))  +
  labs(title ="B. Allele Frequency ≥ 0.5", x = "Selection", y='%') +
  theme_bw(base_size = 12) +
  theme(axis.text.x = element_text(angle=0,vjust=0.5, hjust = 0.5), plot.title = element_text(hjust = 0, face = "bold"),strip.background = element_blank(),panel.grid.major = element_blank(), panel.grid.minor = element_blank(), strip.text = element_text(size = rel(1))) +
  theme(legend.position = "right") +
  theme(strip.background = element_blank(),panel.grid.major = element_blank(), panel.grid.minor = element_blank(), strip.text = element_text(size = rel(1)))+
  theme(plot.margin = unit(c(1, 3, 1.5, 3), "pt")) +
  scale_color_manual(values=colorsMUT, name="Dominance") +
  scale_fill_manual(values=colorsMUT, name="Dominance") +
  geom_bar(stat="identity") +
  facet_grid(SelectionType~ Recomb, scale="free",space="free" ) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  geom_text(data = texttable, aes(label = group_size, x = -Inf, y = Inf, fill=NULL, hjust=-1.25, vjust = 1.1),colour = "#6e6e6e")




COMBOEMP[COMBOEMP$AFClass!="AF >= 0.5",] %>%
  group_by(SelectionType, Recomb) %>%
  mutate(group_size = n()) %>%
  group_by(SelectionType, Recomb,SelClass,DomClass,group_size) %>%
  summarise(perc = n()/max(group_size)*100) -> plot_data


as.data.frame(plot_data) %>% select (-c(SelClass, DomClass, perc)) %>%
  distinct() %>% filter(Recomb!= 'Initial') ->texttable

texttable$group_size<-gsub("^","n = ",perl=T, texttable$group_size)
texttable$Recomb <- factor(texttable$Recomb, levels=c("Initial","Arm","Center"))





plot_data$SelClass <- factor(plot_data$SelClass, levels = c("|Ns| <= 1","1 < |Ns| <= 10", "|Ns| > 10"))
plot_data$DomClass <- factor(plot_data$DomClass, levels = c("0.8 < h <= 1","0.6 < h <= 0.8", "0.4 < h <= 0.6", "0.2 < h <= 0.4","h <= 0.2"))


plot_data$Recomb <- factor(plot_data$Recomb, levels=c("Initial","Arm","Center"))


pl1<-
  ggplot(plot_data, aes(x = SelClass,y=perc,fill=DomClass,ordered = FALSE))  +
  labs(title ="A. Allele Frequency < 0.5", x = "Selection", y='%') +
  theme_bw(base_size = 12) +
  theme(axis.text.x = element_text(angle=0,vjust=0.5, hjust = 0.5), plot.title = element_text(hjust = 0, face = "bold"),strip.background = element_blank(),panel.grid.major = element_blank(), panel.grid.minor = element_blank(), strip.text = element_text(size = rel(1))) +
  theme(legend.position = "right") +
  theme(strip.background = element_blank(),panel.grid.major = element_blank(), panel.grid.minor = element_blank(), strip.text = element_text(size = rel(1)))+
  theme(plot.margin = unit(c(1, 3, 1.5, 3), "pt")) +
  scale_color_manual(values=colorsMUT, name="Dominance") +
  scale_fill_manual(values=colorsMUT, name="Dominance") +
  geom_bar(stat="identity") +
  facet_grid(SelectionType~ Recomb, scale="free",space="free" ) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  geom_text(data = texttable, aes(label = group_size, x = -Inf, y = Inf, fill=NULL, hjust=-1.25, vjust = 1.1),colour = "#6e6e6e")




###combine them

ggarrange(pl1, pl, ncol=1, nrow=2, common.legend = TRUE, legend="right") ->FIG
tiff("S9.tiff",width = 8.2,height = 11,res=300,units="in")
plot(FIG)
dev.off()


######## stats

###difference between domains after selection
independence_test(Dominance ~ Recomb, data=COMBOEMPTEST[COMBOEMPTEST$SelectionType=="Deleterious" & COMBOEMPTEST$AFClass =="AF < 0.5",],distribution = approximate(nresample = 10000))

#Approximative General Independence Test
#
#data:  Dominance by Recomb (Arm, Center)
#Z = -2.0314, p-value = 0.0446
#alternative hypothesis: two.sided


independence_test(Dominance ~ Recomb, data=COMBOEMPTEST[COMBOEMPTEST$SelectionType=="Deleterious" & COMBOEMPTEST$AFClass =="AF >= 0.5",],distribution = approximate(nresample = 10000))

#Approximative General Independence Test

#data:  Dominance by Recomb (Arm, Center)
#Z = -0.98086, p-value = 0.3277
#alternative hypothesis: two.sided

independence_test(Selection ~ Recomb, data=COMBOEMPTEST[COMBOEMPTEST$SelectionType=="Deleterious" & COMBOEMPTEST$AFClass =="AF >= 0.5",],distribution = approximate(nresample = 10000))

#Approximative General Independence Test
#
#data:  Selection by Recomb (Arm, Center)
#Z = 9.5556, p-value < 1e-04
#alternative hypothesis: two.sided

independence_test(Selection ~ Recomb, data=COMBOEMPTEST[COMBOEMPTEST$SelectionType=="Deleterious" & COMBOEMPTEST$AFClass =="AF < 0.5",],distribution = approximate(nresample = 10000))

#Approximative General Independence Test
#
#data:  Selection by Recomb (Arm, Center)
#Z = 7.7892, p-value < 1e-04
#alternative hypothesis: two.sided



independence_test(Dominance ~ Recomb, data=COMBOEMPTEST[COMBOEMPTEST$SelectionType=="Beneficial" & COMBOEMPTEST$AFClass =="AF >= 0.5",],distribution = approximate(nresample = 10000))

#Approximative General Independence Test
#
#data:  Dominance by Recomb (Arm, Center)
#Z = -3.7106, p-value = 5e-04
#alternative hypothesis: two.sided

independence_test(Dominance ~ Recomb, data=COMBOEMPTEST[COMBOEMPTEST$SelectionType=="Beneficial" & COMBOEMPTEST$AFClass =="AF < 0.5",],distribution = approximate(nresample = 10000))

#Approximative General Independence Test
#
#data:  Dominance by Recomb (Arm, Center)
#Z = 2.0827, p-value = 0.0385
#alternative hypothesis: two.sided


independence_test(Selection ~ Recomb, data=COMBOEMPTEST[COMBOEMPTEST$SelectionType=="Beneficial" & COMBOEMPTEST$AFClass =="AF >= 0.5",],distribution = approximate(nresample = 10000))

#Approximative General Independence Test
#
#data:  Selection by Recomb (Arm, Center)
#Z = -16.775, p-value < 1e-04
#alternative hypothesis: two.sided

independence_test(Selection ~ Recomb, data=COMBOEMPTEST[COMBOEMPTEST$SelectionType=="Beneficial" & COMBOEMPTEST$AFClass =="AF < 0.5",],distribution = approximate(nresample = 10000))

#Approximative General Independence Test
#
#data:  Selection by Recomb (Arm, Center)
#Z = -4.7806, p-value < 1e-04
#alternative hypothesis: two.sided


####differencee from the initial distribution

independence_test(Selection + Dominance ~ AFClass, data=COMBOEMPTEST[COMBOEMPTEST$SelectionType=="Beneficial" & COMBOEMPTEST$Recomb !="Arm" & COMBOEMPTEST$AFClass !="AF >= 0.5",],distribution = approximate(nresample = 10000))
#Approximative General Independence Test
#
#data:  Selection, Dominance by AFClass (AF < 0.5, Initial)
#maxT = 8.1201, p-value < 1e-04
#alternative hypothesis: two.sided

independence_test(Selection + Dominance ~ AFClass, data=COMBOEMPTEST[COMBOEMPTEST$SelectionType=="Beneficial" & COMBOEMPTEST$Recomb !="Arm" & COMBOEMPTEST$AFClass !="AF < 0.5",],distribution = approximate(nresample = 10000))

#Approximative General Independence Test
#
#data:  Selection, Dominance by AFClass (AF >= 0.5, Initial)
#maxT = 76.122, p-value < 1e-04
#alternative hypothesis: two.sided

independence_test(Selection + Dominance ~ AFClass, data=COMBOEMPTEST[COMBOEMPTEST$SelectionType=="Beneficial" & COMBOEMPTEST$Recomb !="Center" & COMBOEMPTEST$AFClass !="AF >= 0.5",],distribution = approximate(nresample = 10000))

#Approximative General Independence Test
#
#data:  Selection, Dominance by AFClass (AF < 0.5, Initial)
#maxT = 4.125, p-value = 1e-04
#alternative hypothesis: two.sided


independence_test(Selection + Dominance ~ AFClass, data=COMBOEMPTEST[COMBOEMPTEST$SelectionType=="Beneficial" & COMBOEMPTEST$Recomb !="Center" & COMBOEMPTEST$AFClass !="AF < 0.5",],distribution = approximate(nresample = 10000))

#Approximative General Independence Test
#
#data:  Selection, Dominance by AFClass (AF >= 0.5, Initial)
#maxT = 71.424, p-value < 1e-04
#alternative hypothesis: two.sided

independence_test(Selection + Dominance ~ AFClass, data=COMBOEMPTEST[COMBOEMPTEST$SelectionType=="Deleterious" & COMBOEMPTEST$Recomb !="Arm" & COMBOEMPTEST$AFClass !="AF >= 0.5",],distribution = approximate(nresample = 10000))

#Approximative General Independence Test
#
#data:  Selection, Dominance by AFClass (AF < 0.5, Initial)
#maxT = 7.2004, p-value < 1e-04
#alternative hypothesis: two.sided

independence_test(Selection + Dominance ~ AFClass, data=COMBOEMPTEST[COMBOEMPTEST$SelectionType=="Deleterious" & COMBOEMPTEST$Recomb !="Arm" & COMBOEMPTEST$AFClass !="AF < 0.5",],distribution = approximate(nresample = 10000))

#Approximative General Independence Test
#
#data:  Selection, Dominance by AFClass (AF >= 0.5, Initial)
#maxT = 13.748, p-value < 1e-04
#alternative hypothesis: two.sided

independence_test(Selection + Dominance ~ AFClass, data=COMBOEMPTEST[COMBOEMPTEST$SelectionType=="Deleterious" & COMBOEMPTEST$Recomb !="Center" & COMBOEMPTEST$AFClass !="AF >= 0.5",],distribution = approximate(nresample = 10000))

#Approximative General Independence Test
#
#data:  Selection, Dominance by AFClass (AF < 0.5, Initial)
#maxT = 16.72, p-value < 1e-04
#alternative hypothesis: two.sided

independence_test(Selection + Dominance ~ AFClass, data=COMBOEMPTEST[COMBOEMPTEST$SelectionType=="Deleterious" & COMBOEMPTEST$Recomb !="Center" & COMBOEMPTEST$AFClass !="AF < 0.5",],distribution = approximate(nresample = 10000))

#Approximative General Independence Test
#
#data:  Selection, Dominance by AFClass (AF >= 0.5, Initial)
#maxT = 11.209, p-value < 1e-04
#alternative hypothesis: two.sided
