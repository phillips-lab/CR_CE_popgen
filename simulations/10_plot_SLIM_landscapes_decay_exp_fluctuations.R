#####################################################
######### Plot SLiM simulation results ##############
#####################################################


setwd("/Users/anastasia/Documents/Phillips_lab/drafts/CR_popgen/CR_CE_pi_draft/SLiM/100kb_data/")



library(ggplot2)
library(dplyr)
library(magick)

LAND<-read.table("SLIM_SEL_TABLE_40KB_FIN.out",header=TRUE,sep="\t")

LAND$scenario<-LAND$sample
LAND$scenario<-gsub("sim_[0-9]+_", "", perl=T, LAND$scenario)
LAND$scenario<-gsub(".12stats_40kb", "", perl=T, LAND$scenario)
LAND$scenario<-gsub("Ne_5000_", "", perl=T, LAND$scenario)
colnames(LAND)<-gsub("_win0","",colnames(LAND))
LAND$mut<-gsub(".*_Mut_([12.15]+)_FrD.*", "\\1", perl=T, LAND$scenario)
LAND$sel<-gsub(".*(FrD_.*_FrB.*_.*_SDA_.*_SDC_.*_SB.*A_.*_SB.*C_.*)","\\1",perl=T,LAND$scenario)
LAND$sel<-gsub(".bottleneck","",perl=T,LAND$sel)

LAND$class<-gsub("FrD_0_FrB_0_SDA_0_SDC_0_SBA_0_SBC_0","Neutral",perl=T,LAND$sel)
LAND$class<-gsub(".*FrBa_0.01_.*","Neutral_Del_Bal",perl=T,LAND$class)
LAND$class<-gsub(".*FrB_0.01_.*","Neutral_Del_Ben",perl=T,LAND$class)
LAND$class<-gsub(".*FrB.*","Neutral_Del",perl=T,LAND$class)


LAND$self<-gsub("Self_(.*)_Mut_.*","\\1",perl=T, LAND$scenario)
LAND$self<-as.numeric(as.character(LAND$self))
LAND$self<-LAND$self*100

LAND$self<-gsub("$","%\nselfing",perl=T,LAND$self)
LAND[grep("bottleneck",LAND$scenario,fixed=TRUE),]$self<-"outcrossing\ninbreeding"
LAND$self<-gsub("^0%\nselfing","outcrossing",perl=T,LAND$self)


LAND$self <- factor(LAND$self , levels = c("outcrossing","outcrossing\ninbreeding", "90%\nselfing", "98%\nselfing", "99.9%\nselfing", "100%\nselfing"))

LAND$mut <-gsub("1.15","1.15-1-1.15", perl=FALSE, LAND$mut)
LAND$mut <-gsub("1.5","1.5-1-1.5", perl=FALSE, LAND$mut)
LAND$mut <-gsub("^1$","relative\nmutation rate\n(arm-center-arm)\n1-1-1",perl=TRUE,LAND$mut)
LAND$mut <-gsub("2","2-1-2",LAND$mut)
LAND$mut <- factor(LAND$mut , levels = c("relative\nmutation rate\n(arm-center-arm)\n1-1-1","1.15-1-1.15","1.5-1-1.5","2-1-2"))
colnames(LAND)<-gsub("diplo_","",colnames(LAND))
colnames(LAND)[16]<-"omega"
colnames(LAND)[6]<-"theta"
colnames(LAND)[7] <-"TajimaD"


colors<-c("#b3b1b1", "#E6B243", "#73A8BD", "#196770")



LAND2<-LAND[LAND$sel %in% c("FrD_0_FrB_0_SDA_0_SDC_0_SBA_0_SBC_0","FrD_0.1_FrB_0_SDA_15_SDC_15_SBA_0_SBC_0","FrD_0.1_FrBa_0.01_SDA_15_SDC_15_SBaA_15_SBaC_15","FrD_0.1_FrB_0.01_SDA_15_SDC_15_SBA_15_SBC_15"),]

PILINES<-ggplot(LAND2[LAND2$mut!="1.15-1-1.15",]
                , aes(x = classifiedWinEnd, y = LAND2[LAND2$mut!="1.15-1-1.15",5], color= class,fill=class,ordered = FALSE))  +
  labs(title ="A", x = "Genome position (Mb)", y = bquote("Nucleotide diversity ("~pi~")")) +
  scale_y_continuous(labels = function(x) format(x, scientific = FALSE))  +
  scale_x_continuous(labels = function(x) x / 1000000) +
  theme_bw(base_size = 12) +
  theme(axis.text.x = element_text(angle=0,vjust=0.5, hjust = 0.5), plot.title = element_text(hjust = 0, face = "bold"),strip.background = element_blank(),panel.grid.major = element_blank(), panel.grid.minor = element_blank(), strip.text = element_text(size = rel(0.9)), strip.text.y = element_text(angle =0)) +
  theme(legend.position="top",legend.background=element_blank()) +
  theme(legend.position = c(0.65,1.21)) +
  facet_grid(mut ~ self , space = "free") +
  geom_vline(xintercept = c(1000000,2000000),colour = "#bfbfbf",size = 0.75,linetype = "dashed") +
  stat_summary(fun = mean, geom="line",size=0.75) +
  theme(plot.margin = unit(c(0.6, 0.46, 0.05, 0.05), "cm")) +
  scale_color_manual(values = colors,name = "",labels = c("Neutral","Neutral & Deleterious","Neutral & Deleterious & Balancing","Neutral & Deleterious & Beneficial")) +
  scale_fill_manual(values = colors,name = "",labels = c("Neutral","Neutral & Deleterious","Neutral & Deleterious & Balancing","Neutral & Deleterious & Beneficial")) +
  guides(col=guide_legend(ncol=2,byrow=FALSE))






########### fig6 #https://stackoverflow.com/questions/21239418/use-stat-summary-in-ggplot2-to-calculate-the-mean-and-sd-then-connect-mean-poin

#thanks to https://stackoverflow.com/questions/3472980/how-to-change-facet-labels
class_names <- list(
  'Neutral'="Neutral",
  'Neutral_Del'="Neutral &\n Deleterious",
  'Neutral_Del_Bal'="Neutral &\n Deleterious &\n Balancing",
  'Neutral_Del_Ben'="Neutral &\n Deleterious &\n Beneficial"
)

plot_labeller <- function(variable,value){
  if (variable=='class') {
    return(class_names[value])
  } else {
    return(as.character(value))
  }
}




colors3<- c("#b3b1b1", "#E6B243", "#196770")

#stat_summary(fun.y = mean,
#             fun.ymin = function(x) mean(x) - sd(x),
#             fun.ymax = function(x) mean(x) + sd(x),
#             geom = "ribbon", color=NA)

#scale_y_continuous(breaks = seq(0, 1, by = 0.2)) +

PILINESSD<- ggplot(LAND2[LAND2$mut=="relative\nmutation rate\n(arm-center-arm)\n1-1-1",], aes(x = classifiedWinEnd, y = LAND2[LAND2$mut=="relative\nmutation rate\n(arm-center-arm)\n1-1-1", 5], color= class,fill=class,ordered = FALSE))  +
  labs(title ="B", x = "Genome position (Mb)", y = expression(SD(pi))) +
  scale_y_continuous(labels = function(x) format(x, scientific = FALSE),breaks = seq(0, 0.1, by = 0.0001))  +
  scale_x_continuous(labels = function(x) x / 1000000) +
  theme_bw(base_size = 12) +
  theme(axis.text.x = element_text(angle=0,vjust=0.5, hjust = 0.5), plot.title = element_text(hjust = 0, face = "bold"),strip.background = element_blank(),panel.grid.major = element_blank(), panel.grid.minor = element_blank(), strip.text = element_text(size = rel(0.9)), strip.text.y = element_text(angle =0)) +
  theme(legend.position="none") +
  facet_grid(class  ~ self , space = "free", labeller=plot_labeller) +
  stat_summary(fun.y = "sd",geom="line",size=0.75) +
  geom_vline(xintercept = c(1000000,2000000),colour = "#bfbfbf",size = 0.75,linetype = "dashed") +
  theme(plot.margin = unit(c(0.05, 0.05, 0.05, 0.05), "cm")) +
  scale_color_manual(values = colors,name = "",labels = c("Neutral","Neutral & Deleterious","Neutral & Deleterious & Balancing","Neutral & Deleterious & Beneficial")) +
  scale_fill_manual(values = colors,name = "",labels = c("Neutral","Neutral & Deleterious","Neutral & Deleterious & Balancing", "Neutral & Deleterious & Beneficial"))




library(gridExtra)

grid.arrange(PILINES,PILINESSD, nrow=2,heights=c(1.4,1)) ->F8

tiff(filename = "Fig6.tiff",width = 7.4,height = 8.3,res=300,units="in")
plot(F8)
dev.off()


#####################################################################################
############################### Supplementary figures ###############################
#####################################################################################



LAND$CLASS<-LAND$sel
LAND$CLASS<- gsub("FrD_0_FrB_0_SDA_0_SDC_0_SBA_0_SBC_0","N",LAND$CLASS)
LAND$CLASS<- gsub("FrD_0.1_FrB_0_SDA_15_SDC_15_SBA_0_SBC_0","ND",LAND$CLASS)
LAND$CLASS<- gsub("FrD_0.1_FrB_0_SDA_3_SDC_3_SBA_0_SBC_0","ND",LAND$CLASS)
LAND$CLASS<- gsub("FrD_0.1_FrB_0_SDA_7.5_SDC_15_SBA_0_SBC_0","ND",LAND$CLASS)
LAND$CLASS<- gsub("FrD_0.1_FrB_0.01_SDA_3_SDC_3_SBA_15_SBC_15","NDB",LAND$CLASS)
LAND$CLASS<- gsub("FrD_0.1_FrB_0.01_SDA_15_SDC_15_SBA_15_SBC_15","NDB",LAND$CLASS)
LAND$CLASS<- gsub("FrD_0.1_FrB_0.01_SDA_7.5_SDC_15_SBA_7.5_SBC_15","NDB",LAND$CLASS)
LAND$CLASS<- gsub("FrD_0.1_FrBa_0.01_SDA_15_SDC_15_SBaA_15_SBaC_15","NDB",LAND$CLASS)


LAND$MOD<-LAND$sel
LAND$MOD<- gsub("FrD_0_FrB_0_SDA_0_SDC_0_SBA_0_SBC_0","N",LAND$MOD)
LAND$MOD<- gsub("FrD_0.1_FrB_0_SDA_15_SDC_15_SBA_0_SBC_0","SD-SD",LAND$MOD)
LAND$MOD<- gsub("FrD_0.1_FrB_0_SDA_3_SDC_3_SBA_0_SBC_0","WD-WD",LAND$MOD)
LAND$MOD<- gsub("FrD_0.1_FrB_0_SDA_7.5_SDC_15_SBA_0_SBC_0","MD-SD",LAND$MOD)
LAND$MOD<- gsub("FrD_0.1_FrB_0.01_SDA_3_SDC_3_SBA_15_SBC_15","WD-WD-SB-SB",LAND$MOD)
LAND$MOD<- gsub("FrD_0.1_FrB_0.01_SDA_15_SDC_15_SBA_15_SBC_15","SD-SD-SB-SB",LAND$MOD)
LAND$MOD<- gsub("FrD_0.1_FrB_0.01_SDA_7.5_SDC_15_SBA_7.5_SBC_15","MD-SD-MB-SB",LAND$MOD)
LAND$MOD<- gsub("FrD_0.1_FrBa_0.01_SDA_15_SDC_15_SBaA_15_SBaC_15","SD-SD-SBa-SBa",LAND$MOD)



### Neutral + Deleterious + Beneficial

colorsBEN<-c("#5de3da","#73A8BD","#266e69","#343635")


for (stat in 5:16) {
  #stat_summary(fun = mean, geom="line",size=1.05)
  #geom_smooth(alpha=1,size=1.1,se =FALSE,method = 'loess',span=0.4)

  ppp<-  ggplot(LAND[LAND$CLASS=="NDB",], aes(x = classifiedWinEnd, y = LAND[LAND$CLASS=="NDB",stat], color= MOD,fill=MOD,ordered = FALSE))  +
    labs(title =LETTERS[stat-4], x = "Genome position (Mb)", y = ggplot2:::parse_safe(colnames(LAND)[stat])) +
    scale_y_continuous(labels = function(x) format(x, scientific = FALSE))  +
    scale_x_continuous(labels = function(x) x / 1000000) +
    theme_bw(base_size = 12) +
    theme(axis.text.x = element_text(angle=0,vjust=0.5, hjust = 0.5), plot.title = element_text(hjust = 0, face = "bold"),strip.background = element_blank(),panel.grid.major = element_blank(), panel.grid.minor = element_blank(), strip.text = element_text(size = rel(0.9)), strip.text.y = element_text(angle =0)) +
    theme(legend.position="top",legend.background=element_blank()) +
    facet_grid(mut ~ self , space = "free") +
    geom_vline(xintercept = c(1000000,2000000),colour = "#bfbfbf",size = 0.75,linetype = "dashed") +
    stat_summary(fun = mean, geom="line",size=0.75) +
    theme(plot.margin = unit(c(0.6, 0.46, 0.5, 0.05), "cm"))+
    scale_color_manual(values = colorsBEN,name = "") + scale_fill_manual(values = colorsBEN,name = "")

  tiff(filename=paste0(colnames(LAND)[stat],".SLiM_BEN_only.2.tiff"), width = 7.2,height = 4.77,res=300,units="in")
  plot(ppp)
  dev.off()
}


### Neutral + Deleterious
colorsDEL<-c("#E6B243", "#b38807", "#fac116")


for (stat in 5:16) {
  #stat_summary(fun = mean, geom="line",size=1.05)
  #geom_smooth(alpha=1,size=1.1,se =FALSE,method = 'loess',span=0.4)
  ppp<-  ggplot(LAND[LAND$CLASS=="ND",], aes(x = classifiedWinEnd, y = LAND[LAND$CLASS=="ND",stat], color= MOD,fill=MOD,ordered = FALSE))  +
    labs(title =LETTERS[stat-4], x = "Genome position (Mb)", y = ggplot2:::parse_safe(colnames(LAND)[stat])) +
    scale_y_continuous(labels = function(x) format(x, scientific = FALSE))  +
    scale_x_continuous(labels = function(x) x / 1000000) +
    theme_bw(base_size = 12) +
    theme(axis.text.x = element_text(angle=0,vjust=0.5, hjust = 0.5), plot.title = element_text(hjust = 0, face = "bold"),strip.background = element_blank(),panel.grid.major = element_blank(), panel.grid.minor = element_blank(), strip.text = element_text(size = rel(0.9)), strip.text.y = element_text(angle =0)) +
    theme(legend.position="top",legend.background=element_blank()) +
    facet_grid(mut ~ self , space = "free") +
    geom_vline(xintercept = c(1000000,2000000),colour = "#bfbfbf",size = 0.75,linetype = "dashed") +
    stat_summary(fun = mean, geom="line",size=0.75) +
    theme(plot.margin = unit(c(0.6, 0.46, 0.5, 0.05), "cm"))+
    scale_color_manual(values = colorsBEN,name = "") + scale_fill_manual(values = colorsBEN,name = "")

  tiff(filename=paste0(colnames(LAND)[stat],".SLiM_DEL_only.2.tiff"), width = 7.2,height = 4.77,res=300,units="in")
  plot(ppp)
  dev.off()

}



##################################################
########### Supplementary plots ###################
###################################################



for (stat in 5:16){
  #stat_summary(fun = mean, geom="line",size=1.05)
  #geom_smooth(alpha=1,size=1.1,se =FALSE,method = 'loess',span=0.4)
  pl<-ggplot(LAND2, aes(x = classifiedWinEnd, y = LAND2[, stat], color= class,fill=class,ordered = FALSE))  +
    labs(title =LETTERS[stat-4], x = "Genome position (Mb)", y = ggplot2:::parse_safe(colnames(LAND2)[stat])) +
    scale_y_continuous(labels = function(x) format(x, scientific = FALSE))  +
    scale_x_continuous(labels = function(x) x / 1000000) +
    theme_bw(base_size = 12) +
    theme(axis.text.x = element_text(angle=0,vjust=0.5, hjust = 0.5), plot.title = element_text(hjust = 0, face = "bold"),strip.background = element_blank(),panel.grid.major = element_blank(), panel.grid.minor = element_blank(), strip.text = element_text(size = rel(0.9)), strip.text.y = element_text(angle =0)) +
    theme(legend.position="top",legend.background=element_blank()) +
    theme(legend.position = c(0.62,1.22)) +
    facet_grid(mut ~ self , space = "free") +
    geom_vline(xintercept = c(1000000,2000000),colour = "#bfbfbf",size = 0.75,linetype = "dashed") +
    stat_summary(fun = mean, geom="line",size=0.75) +
    theme(plot.margin = unit(c(0.6, 0.46, 0.5, 0.05), "cm")) +
    scale_color_manual(values = colors,name = "",labels = c("Neutral","Neutral & Deleterious","Neutral & Deleterious & Balancing","Neutral & Deleterious & Beneficial")) +
    scale_fill_manual(values = colors,name = "",labels = c("Neutral","Neutral & Deleterious","Neutral & Deleterious & Balancing","Neutral & Deleterious & Beneficial")) +
    guides(col=guide_legend(ncol=2,byrow=FALSE),colour=guide_legend(override.aes=list(fill=NA)))


  tiff(filename=paste0(colnames(LAND2)[stat],".SLiM_subset.2.tiff"), width = 7.2,height = 4.77,res=300,units="in")
  plot(pl)
  dev.off()
}






###################################################
#############      Beta       #####################
###################################################

BETA<-read.csv("SLIM_SEL_BETA_TABLE_40_FIN.out",header=T,sep="\t")

BETA$scenario<-BETA$sample
BETA$scenario<-gsub("sim_[0-9]+_", "", perl=T, BETA$scenario)
BETA$scenario<-gsub("Ne_5000_", "", perl=T, BETA$scenario)
BETA$scenario<-gsub(".12stats_40kb", "", perl=T, BETA$scenario)
BETA$mut<-gsub(".*_Mut_([12.15]+)_FrD.*", "\\1", perl=T, BETA$scenario)
BETA$sel<-gsub(".*(FrD_.*_FrB.*_.*_SDA_.*_SDC_.*_SB.*A_.*_SB.*C_.*)","\\1",perl=T,BETA$scenario)
BETA$sel<-gsub(".bottleneck","",perl=T,BETA$sel)
BETA$class<-gsub("FrD_0_FrB_0_SDA_0_SDC_0_SBA_0_SBC_0","Neutral",perl=T,BETA$sel)
BETA$class<-gsub(".*FrBa_0.01_.*","Neutral_Del_Bal",perl=T,BETA$class)
BETA$class<-gsub(".*FrB_0.01_.*","Neutral_Del_Ben",perl=T,BETA$class)
BETA$class<-gsub(".*FrB.*","Neutral_Del",perl=T,BETA$class)

BETA$self<-gsub("Self_(.*)_Mut_.*","\\1",perl=T, BETA$scenario)
BETA$self<-as.numeric(as.character(BETA$self))
BETA$self<-BETA$self*100

BETA$self<-gsub("$","%\nselfing",perl=T,BETA$self)
BETA[grep("bottleneck",BETA$scenario,fixed=TRUE),]$self<-"outcrossing\ninbreeding"
BETA$self<-gsub("^0%\nselfing","outcrossing",perl=T,BETA$self)


BETA$self <- factor(BETA$self , levels = c("outcrossing","outcrossing\ninbreeding", "90%\nselfing", "98%\nselfing", "99.9%\nselfing", "100%\nselfing"))

BETA$mut <-gsub("1.15","1.15-1-1.15", perl=FALSE, BETA$mut)
BETA$mut <-gsub("1.5","1.5-1-1.5", perl=FALSE, BETA$mut)
BETA$mut <-gsub("^1$","relative\nmutation rate\n(arm-center-arm)\n1-1-1",perl=TRUE,BETA$mut)
BETA$mut <-gsub("2","2-1-2",BETA$mut)
BETA$mut <- factor(BETA$mut , levels = c("relative\nmutation rate\n(arm-center-arm)\n1-1-1","1.15-1-1.15","1.5-1-1.5","2-1-2"))

#BETA<-BETA[!(BETA$self=="0i" & BETA$sel=="FrD_0.1_FrBa_0.01_SDA_15_SDC_15_SBaA_15_SBaC_15"),]

BETA2<-BETA[BETA$sel %in% c("FrD_0_FrB_0_SDA_0_SDC_0_SBA_0_SBC_0","FrD_0.1_FrB_0_SDA_15_SDC_15_SBA_0_SBC_0","FrD_0.1_FrBa_0.01_SDA_15_SDC_15_SBaA_15_SBaC_15","FrD_0.1_FrB_0.01_SDA_15_SDC_15_SBA_15_SBC_15"),]

#stat_summary(fun = mean, geom="line",size=1.05)
#geom_smooth(alpha=1,size=1.1,se =FALSE,method = 'loess',span=0.4)
pl<-ggplot(BETA2, aes(x = end, y = BETA, color= class,fill=class,ordered = FALSE))  +
  labs(title =LETTERS[13], x = "Genome position (Mb)", y = ggplot2:::parse_safe("beta")) +
  scale_y_continuous(labels = function(x) format(x, scientific = FALSE))  +
  scale_x_continuous(labels = function(x) x / 1000000) +
  theme_bw(base_size = 12) +
  theme(axis.text.x = element_text(angle=0,vjust=0.5, hjust = 0.5), plot.title = element_text(hjust = 0, face = "bold"),strip.background = element_blank(),panel.grid.major = element_blank(), panel.grid.minor = element_blank(), strip.text = element_text(size = rel(0.9)), strip.text.y = element_text(angle =0)) +
  theme(legend.position="top",legend.background=element_blank()) +
  theme(legend.position = c(0.62,1.22)) +
  facet_grid(mut ~ self , space = "free") +
  geom_vline(xintercept = c(1000000,2000000),colour = "#bfbfbf",size = 0.75,linetype = "dashed") +
  stat_summary(fun = mean, geom="line",size=0.75) +
  theme(plot.margin = unit(c(0.6, 0.46, 0.5, 0.05), "cm")) +
  scale_color_manual(values = colors,name = "",labels = c("Neutral","Neutral & Deleterious","Neutral & Deleterious & Balancing","Neutral & Deleterious & Beneficial")) +
  scale_fill_manual(values = colors,name = "",labels = c("Neutral","Neutral & Deleterious","Neutral & Deleterious & Balancing","Neutral & Deleterious & Beneficial")) +
  guides(col=guide_legend(ncol=2,byrow=FALSE),colour=guide_legend(override.aes=list(fill=NA)))



tiff(filename="beta.SLiM_subset.2.tiff", width = 7.2,height = 4.77,res=300,units="in")
plot(pl)
dev.off()


####BEN and DEL


BETA$CLASS<-BETA$sel
BETA$CLASS<- gsub("FrD_0_FrB_0_SDA_0_SDC_0_SBA_0_SBC_0","N",BETA$CLASS)
BETA$CLASS<- gsub("FrD_0.1_FrB_0_SDA_15_SDC_15_SBA_0_SBC_0","ND",BETA$CLASS)
BETA$CLASS<- gsub("FrD_0.1_FrB_0_SDA_3_SDC_3_SBA_0_SBC_0","ND",BETA$CLASS)
BETA$CLASS<- gsub("FrD_0.1_FrB_0_SDA_7.5_SDC_15_SBA_0_SBC_0","ND",BETA$CLASS)
BETA$CLASS<- gsub("FrD_0.1_FrB_0.01_SDA_3_SDC_3_SBA_15_SBC_15","NDB",BETA$CLASS)
BETA$CLASS<- gsub("FrD_0.1_FrB_0.01_SDA_15_SDC_15_SBA_15_SBC_15","NDB",BETA$CLASS)
BETA$CLASS<- gsub("FrD_0.1_FrB_0.01_SDA_7.5_SDC_15_SBA_7.5_SBC_15","NDB",BETA$CLASS)
BETA$CLASS<- gsub("FrD_0.1_FrBa_0.01_SDA_15_SDC_15_SBaA_15_SBaC_15","NDB",BETA$CLASS)


BETA$MOD<-BETA$sel
BETA$MOD<- gsub("FrD_0_FrB_0_SDA_0_SDC_0_SBA_0_SBC_0","N",BETA$MOD)
BETA$MOD<- gsub("FrD_0.1_FrB_0_SDA_15_SDC_15_SBA_0_SBC_0","SD-SD",BETA$MOD)
BETA$MOD<- gsub("FrD_0.1_FrB_0_SDA_3_SDC_3_SBA_0_SBC_0","WD-WD",BETA$MOD)
BETA$MOD<- gsub("FrD_0.1_FrB_0_SDA_7.5_SDC_15_SBA_0_SBC_0","MD-SD",BETA$MOD)
BETA$MOD<- gsub("FrD_0.1_FrB_0.01_SDA_3_SDC_3_SBA_15_SBC_15","WD-WD-SB-SB",BETA$MOD)
BETA$MOD<- gsub("FrD_0.1_FrB_0.01_SDA_15_SDC_15_SBA_15_SBC_15","SD-SD-SB-SB",BETA$MOD)
BETA$MOD<- gsub("FrD_0.1_FrB_0.01_SDA_7.5_SDC_15_SBA_7.5_SBC_15","MD-SD-MB-SB",BETA$MOD)
BETA$MOD<- gsub("FrD_0.1_FrBa_0.01_SDA_15_SDC_15_SBaA_15_SBaC_15","SD-SD-SBa-SBa",BETA$MOD)


colorsBEN<-c("#5de3da","#73A8BD","#266e69","#343635")


ppp<-  ggplot(BETA[BETA$CLASS=="NDB",], aes(x = end, y = BETA, color= MOD,fill=MOD,ordered = FALSE))  +
  labs(title =LETTERS[13], x = "Genome position (Mb)", y = ggplot2:::parse_safe("beta")) +
  scale_y_continuous(labels = function(x) format(x, scientific = FALSE))  +
  scale_x_continuous(labels = function(x) x / 1000000) +
  theme_bw(base_size = 12) +
  theme(axis.text.x = element_text(angle=0,vjust=0.5, hjust = 0.5), plot.title = element_text(hjust = 0, face = "bold"),strip.background = element_blank(),panel.grid.major = element_blank(), panel.grid.minor = element_blank(), strip.text = element_text(size = rel(0.9)), strip.text.y = element_text(angle =0)) +
  theme(legend.position="top",legend.background=element_blank()) +
  facet_grid(mut ~ self , space = "free") +
  geom_vline(xintercept = c(1000000,2000000),colour = "#bfbfbf",size = 0.75,linetype = "dashed") +
  stat_summary(fun = mean, geom="line",size=0.75) +
  theme(plot.margin = unit(c(0.6, 0.46, 0.5, 0.05), "cm"))+
  scale_color_manual(values = colorsBEN,name = "") + scale_fill_manual(values = colorsBEN,name = "")

tiff(filename="beta.SLiM_BEN_only.2.tiff", width = 7.2,height = 4.77,res=300,units="in")
plot(ppp)
dev.off()

### Neutral + Deleterious
colorsDEL<-c("#E6B243", "#b38807", "#fac116")

ppp<-  ggplot(BETA[BETA$CLASS=="ND",], aes(x = end, y = BETA, color= MOD,fill=MOD,ordered = FALSE))  +
  labs(title =LETTERS[13], x = "Genome position (Mb)", y = ggplot2:::parse_safe("beta")) +
  scale_y_continuous(labels = function(x) format(x, scientific = FALSE))  +
  scale_x_continuous(labels = function(x) x / 1000000) +
  theme_bw(base_size = 12) +
  theme(axis.text.x = element_text(angle=0,vjust=0.5, hjust = 0.5), plot.title = element_text(hjust = 0, face = "bold"),strip.background = element_blank(),panel.grid.major = element_blank(), panel.grid.minor = element_blank(), strip.text = element_text(size = rel(0.9)), strip.text.y = element_text(angle =0)) +
  theme(legend.position="top",legend.background=element_blank()) +
  facet_grid(mut ~ self , space = "free") +
  geom_vline(xintercept = c(1000000,2000000),colour = "#bfbfbf",size = 0.75,linetype = "dashed") +
  stat_summary(fun = mean, geom="line",size=0.75) +
  theme(plot.margin = unit(c(0.6, 0.46, 0.5, 0.05), "cm"))+
  scale_color_manual(values = colorsBEN,name = "") + scale_fill_manual(values = colorsBEN,name = "")


#ggsave("beta.SLiM_DEL_only.2.png",ppp,width = 7.2,height = 4.77,dpi=300,scale=1)
#ggsave("beta.SLiM_DEL_only.2.tiff",ppp,width = 7.2,height = 4.77,dpi=300,scale=1)

tiff(filename="beta.SLiM_DEL_only.2.tiff", width = 7.2,height = 4.77,res=300,units="in")
plot(ppp)
dev.off()
###Then divergence and Fis


FIS<-read.csv("SLiM_SEL_REC_SELF_FIS40.txt",header=T,sep=",")


FIS$scenario<-FIS$sample
FIS$scenario<-gsub("sim_[0-9]+_", "", perl=T, FIS$scenario)
FIS$scenario<-gsub("Ne_5000_", "", perl=T, FIS$scenario)
FIS$scenario<-gsub(".12stats_40kb", "", perl=T, FIS$scenario)
FIS$scenario<-gsub(".MAF", "", perl=T, FIS$scenario)
FIS$mut<-gsub(".*_Mut_([12.15]+)_FrD.*", "\\1", perl=T, FIS$scenario)
FIS$sel<-gsub(".*(FrD_.*_FrB.*_.*_SDA_.*_SDC_.*_SB.*A_.*_SB.*C_.*)","\\1",perl=T,FIS$scenario)
FIS$sel<-gsub(".bottleneck","",perl=T,FIS$sel)
FIS$class<-gsub("FrD_0_FrB_0_SDA_0_SDC_0_SBA_0_SBC_0","Neutral",perl=T,FIS$sel)
FIS$class<-gsub(".*FrBa_0.01_.*","Neutral_Del_Bal",perl=T,FIS$class)
FIS$class<-gsub(".*FrB_0.01_.*","Neutral_Del_Ben",perl=T,FIS$class)
FIS$class<-gsub(".*FrB.*","Neutral_Del",perl=T,FIS$class)

FIS$self<-gsub("Self_(.*)_Mut_.*","\\1",perl=T, FIS$scenario)
FIS$self<-as.numeric(as.character(FIS$self))
FIS$self<-FIS$self*100

FIS$self<-gsub("$","%\nselfing",perl=T,FIS$self)
FIS[grep("bottleneck",FIS$scenario,fixed=TRUE),]$self<-"outcrossing\ninbreeding"
FIS$self<-gsub("^0%\nselfing","outcrossing",perl=T,FIS$self)


FIS$self <- factor(FIS$self , levels = c("outcrossing","outcrossing\ninbreeding", "90%\nselfing", "98%\nselfing", "99.9%\nselfing", "100%\nselfing"))

FIS$mut <-gsub("1.15","1.15-1-1.15", perl=FALSE, FIS$mut)
FIS$mut <-gsub("1.5","1.5-1-1.5", perl=FALSE, FIS$mut)
FIS$mut <-gsub("^1$","relative\nmutation rate\n(arm-center-arm)\n1-1-1",perl=TRUE,FIS$mut)
FIS$mut <-gsub("2","2-1-2",FIS$mut)
FIS$mut <- factor(FIS$mut , levels = c("relative\nmutation rate\n(arm-center-arm)\n1-1-1","1.15-1-1.15","1.5-1-1.5","2-1-2"))

FIS2<-FIS[FIS$sel %in% c("FrD_0_FrB_0_SDA_0_SDC_0_SBA_0_SBC_0","FrD_0.1_FrB_0_SDA_15_SDC_15_SBA_0_SBC_0","FrD_0.1_FrBa_0.01_SDA_15_SDC_15_SBaA_15_SBaC_15","FrD_0.1_FrB_0.01_SDA_15_SDC_15_SBA_15_SBC_15"),]

#stat_summary(fun = mean, geom="line",size=1.05)
#geom_smooth(alpha=1,size=1.1,se =FALSE,method = 'loess',span=0.4)
pl<-ggplot(FIS2, aes(x = end, y = FIS, color= class,fill=class,ordered = FALSE))  +
  labs(title =LETTERS[14], x = "Genome position (Mb)", y = ggplot2:::parse_safe("F[IS]")) +
  scale_y_continuous(labels = function(x) format(x, scientific = FALSE))  +
  scale_x_continuous(labels = function(x) x / 1000000) +
  theme_bw(base_size = 12) +
  theme(axis.text.x = element_text(angle=0,vjust=0.5, hjust = 0.5), plot.title = element_text(hjust = 0, face = "bold"),strip.background = element_blank(),panel.grid.major = element_blank(), panel.grid.minor = element_blank(), strip.text = element_text(size = rel(0.9)), strip.text.y = element_text(angle =0)) +
  theme(legend.position="top",legend.background=element_blank()) +
  theme(legend.position = c(0.62,1.22)) +
  facet_grid(mut ~ self , space = "free") +
  geom_vline(xintercept = c(1000000,2000000),colour = "#bfbfbf",size = 0.75,linetype = "dashed") +
  stat_summary(fun = mean, geom="line",size=0.75) +
  theme(plot.margin = unit(c(0.6, 0.46, 0.5, 0.05), "cm")) +
  scale_color_manual(values = colors,name = "",labels = c("Neutral","Neutral & Deleterious","Neutral & Deleterious & Balancing","Neutral & Deleterious & Beneficial")) +
  scale_fill_manual(values = colors,name = "",labels = c("Neutral","Neutral & Deleterious","Neutral & Deleterious & Balancing","Neutral & Deleterious & Beneficial")) +
  guides(col=guide_legend(ncol=2,byrow=FALSE),colour=guide_legend(override.aes=list(fill=NA)))


tiff(filename="FIS.SLiM_subset.2.tiff", width = 7.2,height = 4.77,res=300,units="in")
plot(pl)
dev.off()


####BEN and DEL


FIS$CLASS<-FIS$sel
FIS$CLASS<- gsub("FrD_0_FrB_0_SDA_0_SDC_0_SBA_0_SBC_0","N",FIS$CLASS)
FIS$CLASS<- gsub("FrD_0.1_FrB_0_SDA_15_SDC_15_SBA_0_SBC_0","ND",FIS$CLASS)
FIS$CLASS<- gsub("FrD_0.1_FrB_0_SDA_3_SDC_3_SBA_0_SBC_0","ND",FIS$CLASS)
FIS$CLASS<- gsub("FrD_0.1_FrB_0_SDA_7.5_SDC_15_SBA_0_SBC_0","ND",FIS$CLASS)
FIS$CLASS<- gsub("FrD_0.1_FrB_0.01_SDA_3_SDC_3_SBA_15_SBC_15","NDB",FIS$CLASS)
FIS$CLASS<- gsub("FrD_0.1_FrB_0.01_SDA_15_SDC_15_SBA_15_SBC_15","NDB",FIS$CLASS)
FIS$CLASS<- gsub("FrD_0.1_FrB_0.01_SDA_7.5_SDC_15_SBA_7.5_SBC_15","NDB",FIS$CLASS)
FIS$CLASS<- gsub("FrD_0.1_FrBa_0.01_SDA_15_SDC_15_SBaA_15_SBaC_15","NDB",FIS$CLASS)


FIS$MOD<-FIS$sel
FIS$MOD<- gsub("FrD_0_FrB_0_SDA_0_SDC_0_SBA_0_SBC_0","N",FIS$MOD)
FIS$MOD<- gsub("FrD_0.1_FrB_0_SDA_15_SDC_15_SBA_0_SBC_0","SD-SD",FIS$MOD)
FIS$MOD<- gsub("FrD_0.1_FrB_0_SDA_3_SDC_3_SBA_0_SBC_0","WD-WD",FIS$MOD)
FIS$MOD<- gsub("FrD_0.1_FrB_0_SDA_7.5_SDC_15_SBA_0_SBC_0","MD-SD",FIS$MOD)
FIS$MOD<- gsub("FrD_0.1_FrB_0.01_SDA_3_SDC_3_SBA_15_SBC_15","WD-WD-SB-SB",FIS$MOD)
FIS$MOD<- gsub("FrD_0.1_FrB_0.01_SDA_15_SDC_15_SBA_15_SBC_15","SD-SD-SB-SB",FIS$MOD)
FIS$MOD<- gsub("FrD_0.1_FrB_0.01_SDA_7.5_SDC_15_SBA_7.5_SBC_15","MD-SD-MB-SB",FIS$MOD)
FIS$MOD<- gsub("FrD_0.1_FrBa_0.01_SDA_15_SDC_15_SBaA_15_SBaC_15","SD-SD-SBa-SBa",FIS$MOD)


colorsBEN<-c("#5de3da","#73A8BD","#266e69","#343635")


ppp<-  ggplot(FIS[FIS$CLASS=="NDB",], aes(x = end, y = FIS, color= MOD,fill=MOD,ordered = FALSE))  +
  labs(title =LETTERS[14], x = "Genome position (Mb)", y = ggplot2:::parse_safe("F[IS]")) +
  scale_y_continuous(labels = function(x) format(x, scientific = FALSE))  +
  scale_x_continuous(labels = function(x) x / 1000000) +
  theme_bw(base_size = 12) +
  theme(axis.text.x = element_text(angle=0,vjust=0.5, hjust = 0.5), plot.title = element_text(hjust = 0, face = "bold"),strip.background = element_blank(),panel.grid.major = element_blank(), panel.grid.minor = element_blank(), strip.text = element_text(size = rel(0.9)), strip.text.y = element_text(angle =0)) +
  theme(legend.position="top",legend.background=element_blank()) +
  facet_grid(mut ~ self , space = "free") +
  geom_vline(xintercept = c(1000000,2000000),colour = "#bfbfbf",size = 0.75,linetype = "dashed") +
  stat_summary(fun = mean, geom="line",size=0.75) +
  theme(plot.margin = unit(c(0.6, 0.46, 0.5, 0.05), "cm"))+
  scale_color_manual(values = colorsBEN,name = "") + scale_fill_manual(values = colorsBEN,name = "")


tiff(filename="FIS.SLiM_BEN_only.2.tiff", width = 7.2,height = 4.77,res=300,units="in")
plot(ppp)
dev.off()

### Neutral + Deleterious
colorsDEL<-c("#E6B243", "#b38807", "#fac116")

ppp<-  ggplot(FIS[FIS$CLASS=="ND",], aes(x = end, y = FIS, color= MOD,fill=MOD,ordered = FALSE))  +
  labs(title =LETTERS[14], x = "Genome position (Mb)", y = ggplot2:::parse_safe("F[IS]")) +
  scale_y_continuous(labels = function(x) format(x, scientific = FALSE))  +
  scale_x_continuous(labels = function(x) x / 1000000) +
  theme_bw(base_size = 12) +
  theme(axis.text.x = element_text(angle=0,vjust=0.5, hjust = 0.5), plot.title = element_text(hjust = 0, face = "bold"),strip.background = element_blank(),panel.grid.major = element_blank(), panel.grid.minor = element_blank(), strip.text = element_text(size = rel(0.9)), strip.text.y = element_text(angle =0)) +
  theme(legend.position="top",legend.background=element_blank()) +
  facet_grid(mut ~ self , space = "free") +
  geom_vline(xintercept = c(1000000,2000000),colour = "#bfbfbf",size = 0.75,linetype = "dashed") +
  stat_summary(fun = mean, geom="line",size=0.75) +
  theme(plot.margin = unit(c(0.6, 0.46, 0.5, 0.05), "cm"))+
  scale_color_manual(values = colorsBEN,name = "") + scale_fill_manual(values = colorsBEN,name = "")



tiff(filename="FIS.SLiM_DEL_only.2.tiff", width = 7.2,height = 4.77,res=300,units="in")
plot(ppp)
dev.off()


####### Divergence and divergence over diversity


#SLiM_SEL_REC_SELF_DIVG40.txt
DIVG<-read.csv("SLiM_SEL_REC_SELF_DIVERG_40.txt",header=T,sep=",")
DIVG$scenario<-DIVG$sample
DIVG$scenario<-gsub("sim_[0-9]+_", "", perl=T, DIVG$scenario)
DIVG$scenario<-gsub("Ne_5000_", "", perl=T, DIVG$scenario)
DIVG$scenario<-gsub(".12stats_40kb", "", perl=T, DIVG$scenario)
DIVG$mut<-gsub(".*_Mut_([12.15]+)_FrD.*", "\\1", perl=T, DIVG$scenario)
DIVG$sel<-gsub(".*(FrD_.*_FrB.*_.*_SDA_.*_SDC_.*_SB.*A_.*_SB.*C_.*)","\\1",perl=T,DIVG$scenario)
DIVG$sel<-gsub(".bottleneck","",perl=T,DIVG$sel)
DIVG$class<-gsub("FrD_0_FrB_0_SDA_0_SDC_0_SBA_0_SBC_0","Neutral",perl=T,DIVG$sel)
DIVG$class<-gsub(".*FrBa_0.01_.*","Neutral_Del_Bal",perl=T,DIVG$class)
DIVG$class<-gsub(".*FrB_0.01_.*","Neutral_Del_Ben",perl=T,DIVG$class)
DIVG$class<-gsub(".*FrB.*","Neutral_Del",perl=T,DIVG$class)

DIVG$self<-gsub("Self_(.*)_Mut_.*","\\1",perl=T, DIVG$scenario)
DIVG$self<-as.numeric(as.character(DIVG$self))
DIVG$self<-DIVG$self*100

DIVG$self<-gsub("$","%\nselfing",perl=T,DIVG$self)
DIVG[grep("bottleneck",DIVG$scenario,fixed=TRUE),]$self<-"outcrossing\ninbreeding"
DIVG$self<-gsub("^0%\nselfing","outcrossing",perl=T,DIVG$self)
DIVG$self <- factor(DIVG$self , levels = c("outcrossing","outcrossing\ninbreeding", "90%\nselfing", "98%\nselfing", "99.9%\nselfing", "100%\nselfing"))


DIVG$mut <-gsub("1.15","1.15-1-1.15", perl=FALSE, DIVG$mut)
DIVG$mut <-gsub("1.5","1.5-1-1.5", perl=FALSE, DIVG$mut)
DIVG$mut <-gsub("^1$","relative\nmutation rate\n(arm-center-arm)\n1-1-1",perl=TRUE,DIVG$mut)
DIVG$mut <-gsub("2","2-1-2",DIVG$mut)
DIVG$mut <- factor(DIVG$mut , levels = c("relative\nmutation rate\n(arm-center-arm)\n1-1-1","1.15-1-1.15","1.5-1-1.5","2-1-2"))


DIVG$Divergence<-DIVG$Divergence/40000

DIVG2<-DIVG[DIVG$sel %in% c("FrD_0_FrB_0_SDA_0_SDC_0_SBA_0_SBC_0","FrD_0.1_FrB_0_SDA_15_SDC_15_SBA_0_SBC_0","FrD_0.1_FrBa_0.01_SDA_15_SDC_15_SBaA_15_SBaC_15","FrD_0.1_FrB_0.01_SDA_15_SDC_15_SBA_15_SBC_15"),]





pl<-ggplot(DIVG2, aes(x = end, y = Divergence, color= class,fill=class,ordered = FALSE))  +
  labs(title =LETTERS[15], x = "Genome position (Mb)", y = "Divergence") +
  scale_y_continuous(labels = function(x) format(x, scientific = FALSE))  +
  scale_x_continuous(labels = function(x) x / 1000000) +
  theme_bw(base_size = 12) +
  theme(axis.text.x = element_text(angle=0,vjust=0.5, hjust = 0.5), plot.title = element_text(hjust = 0, face = "bold"),strip.background = element_blank(),panel.grid.major = element_blank(), panel.grid.minor = element_blank(), strip.text = element_text(size = rel(0.9)), strip.text.y = element_text(angle =0)) +
  theme(legend.position="top",legend.background=element_blank()) +
  theme(legend.position = c(0.62,1.22)) +
  facet_grid(mut ~ self , space = "free") +
  geom_vline(xintercept = c(1000000,2000000),colour = "#bfbfbf",size = 0.75,linetype = "dashed") +
  stat_summary(fun = mean, geom="line",size=0.75) +
  theme(plot.margin = unit(c(0.6, 0.46, 0.5, 0.05), "cm")) +
  scale_color_manual(values = colors,name = "",labels = c("Neutral","Neutral & Deleterious","Neutral & Deleterious & Balancing","Neutral & Deleterious & Beneficial")) +
  scale_fill_manual(values = colors,name = "",labels = c("Neutral","Neutral & Deleterious","Neutral & Deleterious & Balancing","Neutral & Deleterious & Beneficial")) +
  guides(col=guide_legend(ncol=2,byrow=FALSE),colour=guide_legend(override.aes=list(fill=NA)))



tiff(filename="DIVG.SLiM_subset.2.tiff", width = 7.2,height = 4.77,res=300,units="in")
plot(pl)
dev.off()


DIVG2$sample<-gsub("$",".12stats_40kb",perl=T,DIVG2$sample)

colnames(DIVG2)[3]<-"classifiedWinEnd"
DIVG2$classifiedWinEnd <- DIVG2$classifiedWinEnd-20000
LAND2$sample<-as.character(LAND2$sample)

DIVG3<-merge(DIVG2,LAND2[,c("sample","classifiedWinEnd","pi")], by=c("sample","classifiedWinEnd"), all.x=T)

pl<-ggplot(DIVG3, aes(x = classifiedWinEnd, y = pi/Divergence, color= class,fill=class,ordered = FALSE))  +
  labs(title =LETTERS[16], x = "Genome position (Mb)", y = "Diversity/Divergence") +
  scale_y_continuous(labels = function(x) format(x, scientific = FALSE))  +
  scale_x_continuous(labels = function(x) x / 1000000) +
  theme_bw(base_size = 12) +
  theme(axis.text.x = element_text(angle=0,vjust=0.5, hjust = 0.5), plot.title = element_text(hjust = 0, face = "bold"),strip.background = element_blank(),panel.grid.major = element_blank(), panel.grid.minor = element_blank(), strip.text = element_text(size = rel(0.9)), strip.text.y = element_text(angle =0)) +
  theme(legend.position="top",legend.background=element_blank()) +
  theme(legend.position = c(0.62,1.22)) +
  facet_grid(mut ~ self , space = "free") +
  geom_vline(xintercept = c(1000000,2000000),colour = "#bfbfbf",size = 0.75,linetype = "dashed") +
  stat_summary(fun = mean, geom="line",size=0.75) +
  theme(plot.margin = unit(c(0.6, 0.46, 0.5, 0.05), "cm")) +
  scale_color_manual(values = colors,name = "",labels = c("Neutral","Neutral & Deleterious","Neutral & Deleterious & Balancing","Neutral & Deleterious & Beneficial")) +
  scale_fill_manual(values = colors,name = "",labels = c("Neutral","Neutral & Deleterious","Neutral & Deleterious & Balancing","Neutral & Deleterious & Beneficial")) +
  guides(col=guide_legend(ncol=2,byrow=FALSE),colour=guide_legend(override.aes=list(fill=NA)))



tiff(filename="DIVG-PI.SLiM_subset.2.tiff", width = 7.2,height = 4.77,res=300,units="in")
plot(pl)
dev.off()




#########################################################
##################TMRCA #################################
###############Only for a subset of simulations #########
#########################################################


TM<-read.table("/Users/anastasia/Documents/Phillips_lab/drafts/CR_popgen/TMRCA/SLIM_SEL_TABLE_40KB_TMRCA.out",header=FALSE,sep="\t")

TM$POS<-gsub("(.*) .*$","\\1", perl=T, TM$V1)
TM$TMRCA<-gsub("(.*) (.*)$","\\2", perl=T, TM$V1)
TM$scenario<-TM$V2
TM$scenario<-gsub("sim_[0-9]+_", "", perl=T, TM$scenario)
TM$scenario<-gsub("Ne_5000_", "", perl=T, TM$scenario)
TM$mut<-gsub(".*_Mut_([12.15]+)_FrD.*", "\\1", perl=T, TM$scenario)
TM$sel<-gsub(".*(FrD_.*_FrB.*_.*_SDA_.*_SDC_.*_SB.*A_.*_SB.*C_.*)","\\1",perl=T,TM$scenario)
TM$sel<-gsub(".bottleneck","",perl=T,TM$sel)

TM$class<-gsub("FrD_0_FrB_0_SDA_0_SDC_0_SBA_0_SBC_0","Neutral",perl=T,TM$sel)
TM$class<-gsub(".*FrBa_0.01_.*","Neutral_Del_Bal",perl=T,TM$class)
TM$class<-gsub(".*FrB_0.01_.*","Neutral_Del_Ben",perl=T,TM$class)
TM$class<-gsub(".*FrB.*","Neutral_Del",perl=T,TM$class)

TM$self<-gsub("Self_(.*)_Mut_.*","\\1",perl=T, TM$scenario)
TM$self<-as.numeric(as.character(TM$self))
TM$self<-TM$self*100

TM$self<-gsub("$","%\nselfing",perl=T,TM$self)
TM[grep("bottleneck",TM$scenario,fixed=TRUE),]$self<-"outcrossing\ninbreeding"
TM$self<-gsub("^0%\nselfing","outcrossing",perl=T,TM$self)


TM$self <- factor(TM$self , levels = c("outcrossing","outcrossing\ninbreeding", "90%\nselfing", "98%\nselfing", "99.9%\nselfing", "100%\nselfing"))

TM$mut <-gsub("1.15","1.15-1-1.15", perl=FALSE, TM$mut)
TM$mut <-gsub("1.5","1.5-1-1.5", perl=FALSE, TM$mut)
TM$mut <-gsub("^1$","relative\nmutation rate\n(arm-center-arm)\n1-1-1",perl=TRUE,TM$mut)
TM$mut <-gsub("2","2-1-2",TM$mut)
TM$mut <- factor(TM$mut , levels = c("relative\nmutation rate\n(arm-center-arm)\n1-1-1","1.15-1-1.15","1.5-1-1.5","2-1-2"))

TM$POS<-as.numeric(as.character(TM$POS))
TM$TMRCA<-as.numeric(as.character(TM$TMRCA))
TM2<-TM[TM$sel %in% c("FrD_0_FrB_0_SDA_0_SDC_0_SBA_0_SBC_0","FrD_0.1_FrB_0_SDA_15_SDC_15_SBA_0_SBC_0","FrD_0.1_FrBa_0.01_SDA_15_SDC_15_SBaA_15_SBaC_15","FrD_0.1_FrB_0.01_SDA_15_SDC_15_SBA_15_SBC_15"),]

pltm<- ggplot(TM2, aes(x = POS, y = TMRCA, color= class,fill=class,ordered = FALSE))  +
  labs(title =LETTERS[17], x = "Genome position (Mb)", y = "TMRCA") +
  scale_y_continuous(labels = function(x) format(x, scientific = FALSE))  +
  scale_x_continuous(labels = function(x) x / 1000000) +
  theme_bw(base_size = 12) +
  theme(axis.text.x = element_text(angle=0,vjust=0.5, hjust = 0.5), plot.title = element_text(hjust = 0, face = "bold"),strip.background = element_blank(),panel.grid.major = element_blank(), panel.grid.minor = element_blank(), strip.text = element_text(size = rel(0.9)), strip.text.y = element_text(angle =0)) +
  theme(legend.position="top",legend.background=element_blank()) +
  theme(legend.position = c(0.62,1.22)) +
  facet_grid(mut ~ self , space = "free") +
  geom_vline(xintercept = c(1000000,2000000),colour = "#bfbfbf",size = 0.75,linetype = "dashed") +
  stat_summary(fun = mean, geom="line",size=0.75) +
  theme(plot.margin = unit(c(0.6, 0.46, 0.5, 0.05), "cm")) +
  scale_color_manual(values = colors,name = "",labels = c("Neutral","Neutral & Deleterious","Neutral & Deleterious & Balancing","Neutral & Deleterious & Beneficial")) +
  scale_fill_manual(values = colors,name = "",labels = c("Neutral","Neutral & Deleterious","Neutral & Deleterious & Balancing","Neutral & Deleterious & Beneficial")) +
  guides(col=guide_legend(ncol=2,byrow=FALSE),colour=guide_legend(override.aes=list(fill=NA)))



tiff(filename="TMRCA.SLiM_subset.2.tiff", width = 7.2,height = 4.77,res=300,units="in")
plot(pltm)
dev.off()


#png(filename="TMRCA.SLiM_subset.2.png", width = 7.2,height = 4.77,res=300,units="in")
#plot(pltm)
#dev.off()




##########################################################
########## Decay of the ancestral polymorphism ############
##########################################################


DECAY<-read.csv("SLIM_DECAY_TABLE.out",header=TRUE,sep="\t")

DECAY$scenario<-DECAY$sample
DECAY$scenario<-gsub("sim_[0-9]+_", "", perl=T, DECAY$scenario)
DECAY$scenario<-gsub(".decay.12stats", "", perl=T, DECAY$scenario)
DECAY$scenario<-gsub("Ne_5000_", "", perl=T, DECAY$scenario)
colnames(DECAY)<-gsub("_win0","",colnames(DECAY))

DECAY$generation<-gsub(".*_[0-9]+_([0-9]+)$", "\\1", perl=TRUE, DECAY$scenario)

DECAY$mut<-gsub(".*_Mut_([12.15]+)_FrD.*", "\\1", perl=T, DECAY$scenario)
DECAY$sel<-gsub(".*(FrD_.*_FrB.*_.*_SDA_.*_SDC_.*_SB.*A_.*_SB.*C_.*)_[0-9]+$","\\1",perl=T,DECAY$scenario)
DECAY$class<-gsub("FrD_0_FrB_0_SDA_0_SDC_0_SBA_0_SBC_0","Neutral",perl=T,DECAY$sel)
DECAY$class<-gsub(".*FrBa_0.01_.*","Neutral_Del_Bal",perl=T,DECAY$class)
DECAY$class<-gsub(".*FrB_0.01_.*","Neutral_Del_Ben",perl=T,DECAY$class)
DECAY$class<-gsub(".*FrB.*","Neutral_Del",perl=T,DECAY$class)

DECAY$self<-gsub("Self_(.*)_Mut_.*","\\1",perl=T, DECAY$scenario)
DECAY$self<-as.numeric(as.character(DECAY$self))
DECAY$self<-DECAY$self*100

DECAY$self<-gsub("$","%\nselfing",perl=T,DECAY$self)
#DECAY[grep("bottleneck",DECAY$scenario,fixed=TRUE),]$self<-"outcrossing\ninbreeding"
#DECAY$self<-gsub("^0%\nselfing","outcrossing",perl=T,DECAY$self)


DECAY$self <- factor(DECAY$self , levels = c("98%\nselfing", "100%\nselfing"))

DECAY$mut <-gsub("1.15","1.15-1-1.15", perl=FALSE, DECAY$mut)
DECAY$mut <-gsub("1.5","1.5-1-1.5", perl=FALSE, DECAY$mut)
DECAY$mut <-gsub("^1$","relative\nmutation rate\n(arm-center-arm)\n1-1-1",perl=TRUE,DECAY$mut)
DECAY$mut <-gsub("2","2-1-2",DECAY$mut)
DECAY$mut <- factor(DECAY$mut , levels = c("relative\nmutation rate\n(arm-center-arm)\n1-1-1","1.15-1-1.15","1.5-1-1.5","2-1-2"))

DECAYSH<-c();
for (i in 1:length(levels(as.factor(DECAY$sample)))){
  print(i);
  A<-DECAY[DECAY$sample==levels(as.factor(DECAY$sample))[i],];

  DECAYSH<-rbind(DECAYSH,data.frame(sample=A$sample[1],scenario=A$scenario[1],self=A$self[1], generation=A$generation[1], mut=A$mut[1],class=A$class[1], ARMPI= mean(A[A$classifiedWinStart<1000000,]$pi,na.rm=TRUE), CENTPI=mean(A[A$classifiedWinStart>1000000 & A$classifiedWinStart<2000000 ,]$pi,na.rm=TRUE), stringsAsFactors = FALSE))


}

#expected values
LANDSH<-c();
for (i in 1:length(levels(as.factor(LAND2$sample)))){
  print(i);
  A<-LAND2[LAND2$sample==levels(as.factor(LAND2$sample))[i],];

  LANDSH<-rbind(LANDSH,data.frame(sample=A$sample[1],scenario=A$scenario[1],self=A$self[1], mut=A$mut[1],class=A$class[1], ARMPI= mean(A[A$classifiedWinStart<1000000,]$pi,na.rm=TRUE), CENTPI=mean(A[A$classifiedWinStart>1000000 & A$classifiedWinStart<2000000 ,]$pi,na.rm=TRUE), stringsAsFactors = FALSE))


}

TMP1<-LANDSH[,c(1:6)]
TMP2<-LANDSH[,c(1:5,7)]
colnames(TMP1)<-c(colnames(LANDSH)[c(1:5)], "PI")
colnames(TMP2)<-c(colnames(LANDSH)[c(1:5)], "PI")
TMP1$Type<-"Arm"
TMP2$Type<-"Center"
LANDSH2<-rbind(TMP1,TMP2)

LAND3<-LANDSH2[LANDSH2$self %in% c("98%\nselfing","100%\nselfing") & LANDSH2$Type=="Arm",]

LAND3 %>%
  group_by(class, self, mut, Type) %>%
  summarise(pi = mean(PI)) ->LAND4;

LAND4<-as.data.frame(LAND4)



DECAYSH$generation<-as.numeric(as.character(DECAYSH$generation))
DECAYSH$mut<-as.character(DECAYSH$mut)

colorsDEC<-c("#6e6b6b","#db9d18","#2381a6","#07434a")
DECAYSH$self<-factor(DECAYSH$self,levels=c("98%\nselfing", "100%\nselfing"))
DECAYSH$mut<-factor(DECAYSH$mut,levels=c("relative\nmutation rate\n(arm-center-arm)\n1-1-1","1.15-1-1.15","1.5-1-1.5","2-1-2"))


### here is the 95% CI for the median, for the CI of the mean "mean_cl_boot" instead of "median_hilow", for SD "mean_sdl"

dec<-ggplot(DECAYSH[DECAYSH$mut!="1.15-1-1.15",], aes(x = generation, y = ARMPI, color= class,fill=class,ordered = FALSE))  +
  labs(title ="", x = "Generation", y = bquote("Nucleotide diversity ("~pi~")")) +
  scale_y_continuous(labels = function(x) format(x, scientific = FALSE),breaks=c(0,0.0003,0.0006,0.0009))  +
  scale_x_continuous(labels = function(x) paste0(x / 5000 - 10, "N")) +
  theme_bw(base_size = 12) +
  theme(axis.text.x = element_text(angle=0,vjust=0.5, hjust = 0.5), plot.title = element_text(hjust = 0, face = "bold"),strip.background = element_blank(),panel.grid.major = element_blank(), panel.grid.minor = element_blank(), strip.text = element_text(size = rel(0.9)), strip.text.y = element_text(angle =0)) +
  theme(legend.position="bottom") +
  guides(col=guide_legend(ncol=2,byrow=FALSE)) +
  facet_grid( mut ~ self , space = "free") +
  stat_summary(fun.y =mean, geom = "line", size = 1.1,alpha = 0.5) +
  stat_summary(fun.data = "median_hilow",
               geom = 'ribbon',
               alpha = 0.25,
               colour = "NA"
  ) +
  geom_hline(data = LAND4[LAND4$mut!="1.15-1-1.15",],
             aes(yintercept = pi, color = class),
             linetype = "dashed",
             size = 0.50) +
  theme(plot.margin = unit(c(0.05, 0.05, 0.05, 0.05), "cm")) +
  scale_color_manual(values = colorsDEC,name = "",labels = c("Neutral","Neutral & Deleterious","Neutral & Deleterious & Balancing","Neutral & Deleterious & Beneficial")) +
  scale_fill_manual(values = colorsDEC,name = "",labels = c("Neutral","Neutral & Deleterious","Neutral & Deleterious & Balancing","Neutral & Deleterious & Beneficial"))




#ggsave(filename = "DECAY_time.png",dec,width = 7.4,height = 5.8,dpi=300,scale=1)
#ggsave(filename = "DECAY_time.tiff",dec,width = 7.4,height = 5.8,dpi=300,scale=1)

tiff(filename="Fig7.tiff", width = 7.4,height = 5.8,res=300,units="in")
plot(dec)
dev.off()



######################################################
######## Exponential growth ##########################
#####################################################

EXP<-read.csv("SLIM_EXP_TABLE.out",header=TRUE,sep="\t")

EXP$scenario<-EXP$sample
EXP$scenario<-gsub("sim_[0-9]+_", "", perl=T, EXP$scenario)
EXP$scenario<-gsub(".1.03exp.12stats", "", perl=T, EXP$scenario)
EXP$scenario<-gsub("Ne_5000_", "", perl=T, EXP$scenario)
EXP$sample<-gsub(".12stats", "", perl=F, EXP$sample)
colnames(EXP)<-gsub("_win0","",colnames(EXP))

EXP$generation<-gsub(".*_[0-9]+_([0-9]+)$", "\\1", perl=TRUE, EXP$scenario)

EXP$mut<-gsub(".*_Mut_([12.15]+)_FrD.*", "\\1", perl=T, EXP$scenario)
EXP$sel<-gsub(".*(FrD_.*_FrB.*_.*_SDA_.*_SDC_.*_SB.*A_.*_SB.*C_.*)_[0-9]+$","\\1",perl=T,EXP$scenario)
EXP$class<-gsub("FrD_0_FrB_0_SDA_0_SDC_0_SBA_0_SBC_0","Neutral",perl=T,EXP$sel)



EXP$self<-gsub("Self_(.*)_Mut_.*","\\1",perl=T, EXP$scenario)
EXP$self<-as.numeric(as.character(EXP$self))
EXP$self<-EXP$self*100

EXP$mut <-gsub("1.15","1.15-1-1.15", perl=FALSE, EXP$mut)
EXP$mut <-gsub("1.5","1.5-1-1.5", perl=FALSE, EXP$mut)
EXP$mut <-gsub("^1$","1-1-1",perl=TRUE,EXP$mut)
EXP$mut <-gsub("2","2-1-2",EXP$mut)
EXP$mut <- factor(EXP$mut , levels = c("1-1-1","1.15-1-1.15","1.5-1-1.5","2-1-2"))



EXP$POSITION<-"NONE"
EXP[EXP$classifiedWinStart<1000000,]$POSITION<-"Arm"
EXP[EXP$classifiedWinStart>1000000 & EXP$classifiedWinStart<2000000,]$POSITION<-"Center"
EXP<-EXP[EXP$POSITION!="NONE",]


#add beta
#100 kb windows
EXPB<-read.csv("SLIM_EXP_BETA_TABLE.out",header=TRUE,sep="\t")
colnames(EXPB)<-gsub("start","classifiedWinStart",colnames(EXPB))
EXPB$classifiedWinStart<- EXPB$classifiedWinStart-50000+1
EXPB<-merge(EXP,EXPB[,c("classifiedWinStart","sample","BETA")], by=c("classifiedWinStart","sample"),all.x = TRUE,sort=FALSE)


EXPON<-EXPB

library(dplyr)
EXPON[complete.cases(EXPON) & EXPON$POSITION == "Arm", ] %>%
  group_by(generation, self, mut, POSITION) %>%
  summarise(
    pi_Arm = mean(pi),
    theta_Arm = mean(thetaW),
    tajD_Arm = mean(tajD),
    distVar_Arm = mean(distVar),
    distSkew_Arm=mean(distSkew),
    distKurt_Arm=mean(distKurt),
    nDiplos_Arm=mean(nDiplos),
    diplo_H1_Arm = mean(diplo_H1),
    diplo_H12_Arm = mean(diplo_H12),
    diplo_H2.H1_Arm = mean(diplo_H2.H1),
    diplo_ZnS_Arm = mean(diplo_ZnS),
    diplo_Omega_Arm = mean(diplo_Omega),
    beta_Arm = mean(BETA)
  ) -> EXPSUM_Arm
EXPON[complete.cases(EXPON) & EXPON$POSITION=="Center",] %>%
  group_by(generation, self, mut, POSITION) %>%

  summarise(
    pi_Center = mean(pi),
    theta_Center = mean(thetaW),
    tajD_Center = mean(tajD),
    distVar_Center = mean(distVar),
    distSkew_Center=mean(distSkew),
    distKurt_Center=mean(distKurt),
    nDiplos_Center=mean(nDiplos),
    diplo_H1_Center = mean(diplo_H1),
    diplo_H12_Center = mean(diplo_H12),
    diplo_H2.H1_Center = mean(diplo_H2.H1),
    diplo_ZnS_Center = mean(diplo_ZnS),
    diplo_Omega_Center = mean(diplo_Omega),
    beta_Center = mean(BETA)
  ) -> EXPSUM_Center


EXPSUM_Center<-as.data.frame(EXPSUM_Center)
EXPSUM_Arm<-as.data.frame(EXPSUM_Arm)


EXPSUM<-merge(EXPSUM_Arm,EXPSUM_Center,by=c("generation","mut","self"),sort=FALSE)
EXPSUM$generation<-gsub("50000","Constant",EXPSUM$generation)
EXPSUM$generation<-gsub("50100","Exponent",EXPSUM$generation)
EXPSUM$NAME<-paste0(EXPSUM$generation,"_",EXPSUM$self, "_", EXPSUM$mut,"_", EXPSUM$POSITION)
rownames(EXPSUM)<-gsub("_$","",perl=T,rownames(EXPSUM))
EXPSUM$NAME<-gsub("_$","",perl=T,EXPSUM$NAME)

EXPSUM$CLASS<-gsub("Exponent_","",perl=T,EXPSUM$NAME)
EXPSUM$CLASS<-gsub("Constant_","",perl=T,EXPSUM$CLASS)



EXPODIFF<-c();

for (cl in 1:length(levels(as.factor(EXPSUM$CLASS)))){
  A <-
    as.numeric(as.character((EXPSUM[EXPSUM$CLASS == levels(as.factor(EXPSUM$CLASS))[cl] &
                                      EXPSUM$generation == "Exponent", colnames(EXPSUM)[grepl('Arm', colnames(EXPSUM))]] - EXPSUM[EXPSUM$CLASS == levels(as.factor(EXPSUM$CLASS))[cl] & EXPSUM$generation == "Constant", colnames(EXPSUM)[grepl('Arm', colnames(EXPSUM))]]) / abs(EXPSUM[EXPSUM$CLASS == levels(as.factor(EXPSUM$CLASS))[cl] &  EXPSUM$generation == "Constant", colnames(EXPSUM)[grepl('Arm', colnames(EXPSUM))]])))

  B <- as.numeric(as.character((EXPSUM[EXPSUM$CLASS == levels(as.factor(EXPSUM$CLASS))[cl] &
                                         EXPSUM$generation == "Exponent", colnames(EXPSUM)[grepl('Center', colnames(EXPSUM))]] - EXPSUM[EXPSUM$CLASS == levels(as.factor(EXPSUM$CLASS))[cl] & EXPSUM$generation == "Constant", colnames(EXPSUM)[grepl('Center', colnames(EXPSUM))]]) / abs(EXPSUM[EXPSUM$CLASS == levels(as.factor(EXPSUM$CLASS))[cl] &  EXPSUM$generation == "Constant", colnames(EXPSUM)[grepl('Center', colnames(EXPSUM))]])))


  A <- c(A,B)
  EXPODIFF<-rbind(EXPODIFF,A);


}



rownames(EXPODIFF)<- as.character(levels(as.factor(EXPSUM$CLASS)))
colnames(EXPODIFF) <- c(colnames(EXPSUM)[grepl('Arm', colnames(EXPSUM))],colnames(EXPSUM)[grepl('Center', colnames(EXPSUM))])

#with 3 landscapes
#EXPODIFF<-EXPODIFF[c(1:3,7:15,4:6),c(1,14,2,15,3,16,4,17,5,18,6,19,7,20,8,21,9,22,10,23,11,24,12,25,13,26)]

EXPODIFF<-EXPODIFF[c(1:4,9:20,5:8),c(1,14,2,15,3,16,4,17,5,18,6,19,7,20,8,21,9,22,10,23,11,24,12,25,13,26)]
EXPODIFF<-as.matrix(EXPODIFF)


#change in stats
reshape2::melt(EXPODIFF) -> EXPORDER
head(EXPORDER[order(EXPORDER$value),],30)

#most affected statistics:
#tajD
#distSkew
#distKurt
#diplo_Omega





library(corrplot)
library(dichromat)
library(gplots)


my_palette <-colorRampPalette(c("#183843",  "#73A8BD", "white", "#f2c76f","#594114"))(n = 299);

col_breaks = c(seq(-5,-0.4,length=120),  # for blue
               seq(-0.2,0.2,length=60),           # white
               seq(0.3,5,length=120))

png("EXP_to_CONSTANT_Hetmap3.png",    # create PNG for the heat map
    width = 8*300,        # 5 x 300 pixels
    height = 6*300,
    res = 300,            # 300 pixels per inch
    pointsize = 11)



heatmap.2(EXPODIFF,main="B",
          notecol="black",      # change font color of cell labels to black
          density.info="none",  # turns off density plot inside color legend
          trace="none",         # turns off trace lines inside the heat map
          margins =c(12,9),     # widens margins around plot
          col=my_palette,  breaks=sort(col_breaks),     # use on color palette defined earlier
          dendrogram="none",     # only draw a row dendrogram
          Colv=FALSE,Rowv=FALSE, key.xlab="Fold change", key.title = "",srtCol=45)

dev.off()

tiff("EXP_to_CONSTANT_Hetmap3.tiff",    # create PNG for the heat map
     width = 8*300,        # 5 x 300 pixels
     height = 6*300,
     res = 300,            # 300 pixels per inch
     pointsize = 11)



heatmap.2(EXPODIFF, main="B. Exponential growth",
          notecol="black",      # change font color of cell labels to black
          density.info="none",  # turns off density plot inside color legend
          trace="none",         # turns off trace lines inside the heat map
          margins =c(12,9),     # widens margins around plot
          col=my_palette,  breaks=sort(col_breaks),     # use on color palette defined earlier
          dendrogram="none",     # only draw a row dendrogram
          Colv=FALSE,Rowv=FALSE, key.xlab="Fold change", key.title = "",srtCol=45)

dev.off()


######################################################
######## Fluctuations growth ##########################
#####################################################

FLUC<-read.csv("SLIM_FLUCT_TABLE.out",header=TRUE,sep="\t")

FLUC$scenario<-FLUC$sample
FLUC$scenario<-gsub("sim_[0-9]+_", "", perl=T, FLUC$scenario)
FLUC$scenario<-gsub(".5-15-fluct.12stats", "", perl=T, FLUC$scenario)
FLUC$scenario<-gsub("Ne_5000_", "", perl=T, FLUC$scenario)
FLUC$sample<- gsub(".12stats","",FLUC$sample)
colnames(FLUC)<-gsub("_win0","",colnames(FLUC))

FLUC$generation<-gsub(".*_[0-9]+_([0-9]+)$", "\\1", perl=TRUE, FLUC$scenario)

FLUC$mut<-gsub(".*_Mut_([12.15]+)_FrD.*", "\\1", perl=T, FLUC$scenario)
FLUC$sel<-gsub(".*(FrD_.*_FrB.*_.*_SDA_.*_SDC_.*_SB.*A_.*_SB.*C_.*)_[0-9]+$","\\1",perl=T,FLUC$scenario)
FLUC$class<-gsub("FrD_0_FrB_0_SDA_0_SDC_0_SBA_0_SBC_0","Neutral",perl=T,FLUC$sel)

FLUC$self<-gsub("Self_(.*)_Mut_.*","\\1",perl=T, FLUC$scenario)
FLUC$self<-as.numeric(as.character(FLUC$self))
FLUC$self<-FLUC$self*100

FLUC$mut <-gsub("1.15","1.15-1-1.15", perl=FALSE, FLUC$mut)
FLUC$mut <-gsub("1.5","1.5-1-1.5", perl=FALSE, FLUC$mut)
FLUC$mut <-gsub("^1$","1-1-1",perl=TRUE,FLUC$mut)
FLUC$mut <-gsub("2","2-1-2",FLUC$mut)
FLUC$mut <- factor(FLUC$mut , levels = c("1-1-1","1.15-1-1.15","1.5-1-1.5","2-1-2"))



FLUC$POSITION<-"NONE"
FLUC[FLUC$classifiedWinStart<1000000,]$POSITION<-"Arm"
FLUC[FLUC$classifiedWinStart>1000000 & FLUC$classifiedWinStart<2000000,]$POSITION<-"Center"
FLUC<-FLUC[FLUC$POSITION!="NONE",]

#add beta
#100 kb windows
FLUCB<-read.csv("SLIM_FLUCT_BETA_TABLE.out",header=TRUE,sep="\t")
colnames(FLUCB)<-gsub("start","classifiedWinStart",colnames(FLUCB))
FLUCB$classifiedWinStart<- FLUCB$classifiedWinStart-50000+1
FLUCB<-merge(FLUC,FLUCB[,c("classifiedWinStart","sample","BETA")], by=c("classifiedWinStart","sample"),all.x = TRUE,sort=FALSE)




FLUCON<-FLUCB

library(dplyr)
FLUCON[complete.cases(FLUCON) & FLUCON$POSITION=="Arm",] %>%
  group_by(generation, self, mut, POSITION) %>%
  summarise(
    pi_Arm = mean(pi),
    theta_Arm = mean(thetaW),
    tajD_Arm = mean(tajD),
    distVar_Arm = mean(distVar),
    distSkew_Arm=mean(distSkew),
    distKurt_Arm=mean(distKurt),
    nDiplos_Arm=mean(nDiplos),
    diplo_H1_Arm = mean(diplo_H1),
    diplo_H12_Arm = mean(diplo_H12),
    diplo_H2.H1_Arm = mean(diplo_H2.H1),
    diplo_ZnS_Arm = mean(diplo_ZnS),
    diplo_Omega_Arm = mean(diplo_Omega),
    beta_Arm = mean(BETA)
  ) -> FLUCSUM_Arm

FLUCON[complete.cases(FLUCON) & FLUCON$POSITION=="Center",] %>%
  group_by(generation, self, mut, POSITION) %>%
  summarise(
    pi_Center = mean(pi),
    theta_Center = mean(thetaW),
    tajD_Center = mean(tajD),
    distVar_Center = mean(distVar),
    distSkew_Center=mean(distSkew),
    distKurt_Center=mean(distKurt),
    nDiplos_Center=mean(nDiplos),
    diplo_H1_Center = mean(diplo_H1),
    diplo_H12_Center = mean(diplo_H12),
    diplo_H2.H1_Center = mean(diplo_H2.H1),
    diplo_ZnS_Center = mean(diplo_ZnS),
    diplo_Omega_Center = mean(diplo_Omega),
    beta_Center = mean(BETA)
  ) -> FLUCSUM_Center
FLUCSUM_Center<-as.data.frame(FLUCSUM_Center)

FLUCSUM_Arm<-as.data.frame(FLUCSUM_Arm)
FLUCSUM<-merge(FLUCSUM_Arm,FLUCSUM_Center,by=c("generation","mut","self"),order=FALSE)
FLUCSUM$generation<-gsub("50000","Constant",FLUCSUM$generation)
FLUCSUM$generation<-gsub("50100","Exponent",FLUCSUM$generation)
FLUCSUM$NAME<-paste0(FLUCSUM$generation,"_",FLUCSUM$self, "_", FLUCSUM$mut,"_", FLUCSUM$POSITION)
rownames(FLUCSUM)<-FLUCSUM$NAME
rownames(FLUCSUM)<-gsub("_$","",perl=T,rownames(FLUCSUM))
FLUCSUM$NAME<-gsub("_$","",perl=T,FLUCSUM$NAME)

FLUCSUM$CLASS<-gsub("Exponent_","",perl=T,FLUCSUM$NAME)
FLUCSUM$CLASS<-gsub("Constant_","",perl=T,FLUCSUM$CLASS)


FLUCODIFF<-c();

for (cl in 1:length(levels(as.factor(FLUCSUM$CLASS)))){

  A <-
    as.numeric(as.character((FLUCSUM[FLUCSUM$CLASS == levels(as.factor(FLUCSUM$CLASS))[cl] &
                                       FLUCSUM$generation == "Exponent", colnames(FLUCSUM)[grepl('Arm', colnames(FLUCSUM))]] - FLUCSUM[FLUCSUM$CLASS == levels(as.factor(FLUCSUM$CLASS))[cl] & FLUCSUM$generation == "Constant", colnames(FLUCSUM)[grepl('Arm', colnames(FLUCSUM))]]) / abs(FLUCSUM[FLUCSUM$CLASS == levels(as.factor(FLUCSUM$CLASS))[cl] &  FLUCSUM$generation == "Constant", colnames(FLUCSUM)[grepl('Arm', colnames(FLUCSUM))]])))

  B <- as.numeric(as.character((FLUCSUM[FLUCSUM$CLASS == levels(as.factor(FLUCSUM$CLASS))[cl] &
                                          FLUCSUM$generation == "Exponent", colnames(FLUCSUM)[grepl('Center', colnames(FLUCSUM))]] - FLUCSUM[FLUCSUM$CLASS == levels(as.factor(FLUCSUM$CLASS))[cl] & FLUCSUM$generation == "Constant", colnames(FLUCSUM)[grepl('Center', colnames(FLUCSUM))]]) / abs(FLUCSUM[FLUCSUM$CLASS == levels(as.factor(FLUCSUM$CLASS))[cl] &  FLUCSUM$generation == "Constant", colnames(FLUCSUM)[grepl('Center', colnames(FLUCSUM))]])))

  A <- c(A,B)
  FLUCODIFF<-rbind(FLUCODIFF,A);

}

rownames(FLUCODIFF)<- as.character(levels(as.factor(FLUCSUM$CLASS)))
colnames(FLUCODIFF) <- c(colnames(FLUCSUM)[grepl('Arm', colnames(FLUCSUM))],colnames(FLUCSUM)[grepl('Center', colnames(FLUCSUM))])

reshape2::melt(FLUCODIFF) -> FLUCORDER
tail(FLUCORDER[order(FLUCORDER$value),],30) #same stats affected


FLUCODIFF<-FLUCODIFF[c(1:4,9:20,5:8),c(1,14,2,15,3,16,4,17,5,18,6,19,7,20,8,21,9,22,10,23,11,24,12,25,13,26)]
FLUCODIFF<-as.matrix(FLUCODIFF)


my_palette <-colorRampPalette(c("#183843",  "#73A8BD", "white", "#f2c76f","#594114"))(n = 299);

col_breaks = c(seq(-5,-0.4,length=120),  # for blue
               seq(-0.2,0.2,length=60),           # white
               seq(0.3,5,length=120))


png("FLUC_to_CONSTANT_Hetmap3.png",    # create PNG for the heat map
    width = 8*300,        # 5 x 300 pixels
    height = 6*300,
    res = 300,            # 300 pixels per inch
    pointsize = 11)

heatmap.2(FLUCODIFF, main="A",
          notecol="black",      # change font color of cell labels to black
          density.info="none",  # turns off density plot inside color legend
          trace="none",         # turns off trace lines inside the heat map
          margins =c(12,9),     # widens margins around plot
          col=my_palette,  breaks=sort(col_breaks),     # use on color palette defined earlier
          dendrogram="none",     # only draw a row dendrogram
          Colv=FALSE,Rowv=FALSE, key.xlab="Fold change", key.title = "",srtCol=45)

dev.off()

tiff("FLUC_to_CONSTANT_Hetmap3.tiff",    # create PNG for the heat map
     width = 8*300,        # 5 x 300 pixels
     height = 6*300,
     res = 300,            # 300 pixels per inch
     pointsize = 11)

heatmap.2(FLUCODIFF,main="A. Fluctuations in population size",
          notecol="black",      # change font color of cell labels to black
          density.info="none",  # turns off density plot inside color legend
          trace="none",         # turns off trace lines inside the heat map
          margins =c(12,9),     # widens margins around plot
          col=my_palette,  breaks=sort(col_breaks),     # use on color palette defined earlier
          dendrogram="none",     # only draw a row dendrogram
          Colv=FALSE,Rowv=FALSE, key.xlab="Fold change", key.title = "",srtCol=45)

dev.off()





#library(magick)
#getwd()
#[1] "/Users/anastasia/Documents/Phillips_lab/drafts/CR_popgen/CR_CE_pi_draft/SLiM/100kb_data"
setwd("/Users/anastasia/Documents/Phillips_lab/drafts/CR_popgen/CR_CE_pi_draft/SLiM/100kb_data")

TIF <- file.info(list.files(pattern="*.SLiM_subset.2.tiff"))
TIF <- TIF[with(TIF, order(as.POSIXct(mtime))), ]
TIF<-rownames(TIF)
TIF<-gsub("tiff","tif",TIF)
setwd("/Users/anastasia/Downloads/PACE Corrected-1/")

combo<-image_append(do.call("c", lapply(TIF, function(h){image_read(h)})), stack = TRUE)
image_write(combo, path = "S10.tiff", format = "tiff", quality = 100, compression="LZW")




setwd("/Users/anastasia/Documents/Phillips_lab/drafts/CR_popgen/CR_CE_pi_draft/SLiM/100kb_data")
TIF <- file.info(list.files(pattern="*_BEN_only.2.tiff"))
TIF <- TIF[with(TIF, order(as.POSIXct(mtime))), ]
TIF<-rownames(TIF)
TIF<-gsub("tiff","tif",TIF)

setwd("/Users/anastasia/Downloads/PACE Corrected-2/")
combo<-image_append(do.call("c", lapply(TIF, function(h){image_read(h)})), stack = TRUE)
image_write(combo, path = "S12.tiff", format = "tiff", quality = 100, compression="LZW")



setwd("/Users/anastasia/Documents/Phillips_lab/drafts/CR_popgen/CR_CE_pi_draft/SLiM/100kb_data")
TIF <- file.info(list.files(pattern="*_DEL_only.2.tiff"))
TIF <- TIF[with(TIF, order(as.POSIXct(mtime))), ]
TIF<-rownames(TIF)
TIF<-gsub("tiff","tif",TIF)

setwd("/Users/anastasia/Downloads/PACE Corrected-2/")
combo<-image_append(do.call("c", lapply(TIF, function(h){image_read(h)})), stack = TRUE)
image_write(combo, path = "S11.tiff", format = "tiff", quality = 100, compression="LZW")



###################

setwd("/Users/anastasia/Documents/Phillips_lab/drafts/CR_popgen/CR_CE_pi_draft/SLiM/100kb_data")

combo<-image_append(c(image_read("FLUC_to_CONSTANT_Hetmap3.tiff"),image_read("EXP_to_CONSTANT_Hetmap3.tiff")), stack = TRUE)
image_write(combo, path = "S13.tiff", format = "tiff", quality = 100)
