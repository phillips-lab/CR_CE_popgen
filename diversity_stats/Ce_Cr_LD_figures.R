#######################################################
################### Plot LD ###########################
###########  C. elegans & C.remanei   #################
#######################################################
###########       Fig 4, S5       ####################

setwd("/Users/anastasia/Documents/Phillips_lab/drafts/CR_popgen/LD")


library(reshape2)
library(ggplot2)
library(ggpubr)
library(gridExtra)
library(magick)
library(dichromat)


LDCE<- read.table("CE_r2_LD_windNOPRUNNING_INTER_FIN.summary",header=FALSE, sep=" ")
colnames(LDCE)<-c("chr_A","pos_A","chr_B","pos_B","Mean")
LDCE$species<-"C. elegans"
LDCE<-LDCE[LDCE$chr_A!="CHR_A",]
LDCE<-LDCE[LDCE$chr_A!="CHR_B",]

LDCE$chr_A<-as.character(LDCE$chr_A)
LDCE$chr_B<-as.character(LDCE$chr_B)
LDCE$chr_A<-gsub("23","X",LDCE$chr_A)
LDCE$chr_B<-gsub("23","X",LDCE$chr_B)


LDCE$Distance<-(LDCE$pos_B-LDCE$pos_A)/10
LDCE[LDCE$chr_A!=LDCE$chr_B,]$Distance<-"NA" #remove intrachromosomal LD


LDCR<- read.table("CR_r2_LD_windNOPRUNNING_INTER_FIN.summary",header=FALSE, sep=" ")
colnames(LDCR)<-c("chr_A","pos_A","chr_B","pos_B","Mean")
LDCR$species<-"C. remanei"
LDCR<-LDCR[LDCR$chr_A!="CHR_A",]
LDCR<-LDCR[LDCR$chr_A!="CHR_B",]
LDCR$chr_A<-gsub("23","X",LDCR$chr_A)
LDCR$chr_B<-gsub("23","X",LDCR$chr_B)
LDCR$Distance<-(LDCR$pos_B-LDCR$pos_A)/10
LDCR[LDCR$chr_A!=LDCR$chr_B,]$Distance<-"NA" #remove intrachromosomal LD


LDCRK<- read.table("CR_r2_LD_wind_2KB_FIN.summary",header=FALSE, sep=" ")
colnames(LDCRK)<-c("chr_A","pos_A","chr_B","pos_B","Mean")
LDCRK$species<-"C. remanei"
LDCRK<-LDCRK[LDCRK$chr_A!="CHR_A",]
LDCRK<-LDCRK[LDCRK$chr_A!="CHR_B",]
LDCRK$chr_A<-gsub("23","X",LDCRK$chr_A)
LDCRK$chr_B<-gsub("23","X",LDCRK$chr_B)
LDCRK$Distance<-((LDCRK$pos_B*100)-(LDCRK$pos_A*100))
LDCRK$chr_A<-as.character(LDCRK$chr_A)
LDCRK$chr_B<-as.character(LDCRK$chr_B)


LDCEK<- read.table("CE_r2_LD_wind_2KB_FIN.summary",header=FALSE, sep=" ")
colnames(LDCEK)<-c("chr_A","pos_A","chr_B","pos_B","Mean")
LDCEK$species<-"C. elegans"
LDCEK<-LDCEK[LDCEK$chr_A!="CHR_A",]
LDCEK<-LDCEK[LDCEK$chr_A!="CHR_B",]
LDCEK$chr_A<-gsub("23","X",LDCEK$chr_A)
LDCEK$chr_B<-gsub("23","X",LDCEK$chr_B)
LDCEK$Distance<-((LDCEK$pos_B*100)-(LDCEK$pos_A*100))
LDCEK$chr_A<-as.character(LDCEK$chr_A)
LDCEK$chr_B<-as.character(LDCEK$chr_B)


colorsMUT<-c("#183843", "#196770", "#73A8BD", "#C0BDBC", "#f2c76f", "#916b20")

my_palette <-colorRampPalette(c("white", "#07434a"))(n = 10);
col_breaks = c(seq(0,0.5,length=10))


ggplot(LDCRK, aes(x = Distance, y = Mean, color= chr_A, fill=chr_A, ordered = FALSE))  +
  labs( x = "bp", y = "") +
  scale_y_continuous(breaks = seq(0, 1, by = 0.2)) +
  coord_cartesian(ylim=c(0.1,0.46),clip="off")+
  theme_bw(base_size = 12) +
  theme(axis.text.x = element_text(angle=0,vjust=0.5, hjust = 0.5), plot.title = element_text(hjust = 0, face = "bold"),strip.background = element_blank(),panel.grid.major = element_blank(), panel.grid.minor = element_blank(), strip.text = element_text(size = rel(1))) +
  theme(legend.position = "none") +
  theme(strip.text.y = element_text(face = "italic"))+
  theme(strip.background = element_blank(),panel.grid.major = element_blank(), panel.grid.minor = element_blank(), strip.text = element_text(size = rel(1))) +
  stat_summary(fun = median, geom="line",size=0.75)+
  theme(plot.margin = unit(c(0, 0, 0, 0), "pt")) +
  scale_color_manual(values = colorsMUT,name="") +
  scale_fill_manual(values = colorsMUT,name="") +
  facet_grid(species~., space="free") + xlim(0,1000) +theme(strip.text.y = element_blank())  + theme(panel.background = element_blank()) + theme(panel.background = element_rect(fill = "transparent", colour = NA),  plot.background = element_rect(fill = "transparent", colour = NA)) ->plKB


ggplot(LDCEK, aes(x = Distance, y = Mean, color= chr_A, fill=chr_A, ordered = FALSE))  +
  labs( x = "bp", y = "") +
  scale_y_continuous(breaks = seq(0, 1, by = 0.2)) +
  coord_cartesian(ylim=c(0.65,1.01),clip="off")+
  theme_bw(base_size = 12) +
  theme(axis.text.x = element_text(angle=0,vjust=0.5, hjust = 0.5), plot.title = element_text(hjust = 0, face = "bold"),strip.background = element_blank(),panel.grid.major = element_blank(), panel.grid.minor = element_blank(), strip.text = element_text(size = rel(1))) +
  theme(legend.position = "none") +
  theme(strip.text.y = element_text(face = "italic"))+
  theme(strip.background = element_blank(),panel.grid.major = element_blank(), panel.grid.minor = element_blank(), strip.text = element_text(size = rel(1))) +
  stat_summary(fun = median, geom="line",size=0.75)+
  theme(plot.margin = unit(c(0, 0, 0, 0), "pt")) +
  scale_color_manual(values = colorsMUT,name="") +
  scale_fill_manual(values = colorsMUT,name="") +
  facet_grid(species~., space="free") + xlim(0,1000) +theme(strip.text.y = element_blank())  + theme(panel.background = element_blank()) + theme(panel.background = element_rect(fill = "transparent", colour = NA),  plot.background = element_rect(fill = "transparent", colour = NA)) ->plKBCE

annotation_custom2 <- function (grob, xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf, data)
{
  layer(data = data, stat = StatIdentity, position = PositionIdentity,
        geom = ggplot2:::GeomCustomAnn,
        inherit.aes = FALSE, params = list(grob = grob,
                                           xmin = xmin, xmax = xmax,
                                           ymin = ymin, ymax = ymax))
}



COMBO<-rbind(LDCE,LDCR)
COMBO$Distance<-as.numeric(as.character(COMBO$Distance))

pl3 <- ggplot(COMBO[COMBO$chr_A==COMBO$chr_B,], aes(x = Distance, y = Mean, color= chr_A, fill=chr_A, ordered = FALSE))  +
  labs(title ="", x = "Distance (Mb)", y = expression(paste("Linkage Disequilibrium (",r^{2},")"))) +
  scale_y_continuous(breaks = seq(0, 1, by = 0.2),limits=c(0,1)) +
  theme_bw(base_size = 12) +
  theme(axis.text.x = element_text(angle=0,vjust=0.5, hjust = 0.5), plot.title = element_text(hjust = 0, face = "bold"),strip.background = element_blank(),panel.grid.major = element_blank(), panel.grid.minor = element_blank(), strip.text = element_text(size = rel(1))) +
  theme(legend.position = "right") +
  theme(strip.text.y = element_text(face = "italic"))+
  theme(strip.background = element_blank(),panel.grid.major = element_blank(), panel.grid.minor = element_blank(), strip.text = element_text(size = rel(1))) +
  stat_summary(fun = median, geom="line",size=0.75)+
  theme(plot.margin = unit(c(1, 3, 1.5, 3), "pt")) +
  scale_color_manual(values = colorsMUT,name="") +
  scale_fill_manual(values = colorsMUT,name="") +
  facet_grid(species~., space="free")


pl3+ annotation_custom2(grob=ggplotGrob(plKB), data=LDCRK[LDCRK$Distance<1000,],
                        xmin = 16.5, xmax=26.5, ymin=0.35, ymax=0.99)  + annotation_custom2(grob=ggplotGrob(plKBCE), data=LDCEK[LDCEK$Distance<1000,],
                                                                                           xmin = 16.5, xmax=26.5, ymin=0.35, ymax=0.9) ->pl4

tiff(filename = "Fig4.tiff",width = 7.2,height = 3.7,res=300,units="in")
plot(pl4)
dev.off()




##########################
###  Intrachromosomal  ###
##########################

LDCE$TYPE<-"Inter"
LDCE[LDCE$chr_A==LDCE$chr_B,]$TYPE<-"Intra"
LDCE$TYPE<-as.factor(LDCE$TYPE)
cohensD(Mean ~ TYPE, data=LDCE,method="corrected")
#[1] 0.5778754
oneway_test(Mean ~ TYPE, data=LDCE,alternative="less",distribution=approximate(nresample=10000))
#Approximative Two-Sample Fisher-Pitman Permutation Test
#
#data:  Mean by TYPE (Inter, Intra)
#Z = -151.08, p-value < 1e-04
#alternative hypothesis: true mu is less than 0


################outcrossing Rate


#C. elegans population size
#2*4*N*2.3e-09=0.0006267073

Nce=0.0006267073/(2*4*2.3e-09)
Nce
#[1] 34060.18



mean(sapply(LDCE[LDCE$TYPE=="Inter",]$Mean, function(r) { return(1-(1-r)/(r*(1+(8*Nce*0.5))-1))}))*100
#[1] 99.99765
sd(sapply(LDCE[LDCE$TYPE=="Inter",]$Mean, function(r) { return(1-(1-r)/(r*(1+(8*Nce*0.5))-1))}))*100
#[1] 0.00491034


LDCEI2<-LDCE[,c(3,4,1,2,5,6,7)]
colnames(LDCEI2)<-colnames(LDCE)
LDCEI2<-rbind(LDCE,LDCEI2)


LDCRI2<-LDCR[,c(3,4,1,2,5,6,7)]
colnames(LDCRI2)<-colnames(LDCR)
LDCRI2<-rbind(LDCR,LDCRI2)
LDCEI2$chr_B = factor(LDCEI2$chr_B, levels=c("I","II","III","IV","V","X"))
LDCEI2$chr_A = factor(LDCEI2$chr_A, levels=c("X","V","IV","III","II","I"))

pl4<-ggplot(LDCEI2, aes(x = pos_B/10, y = pos_A/10, ordered = FALSE))  +
  labs(title ="C. elegans", x = "Genome position (Mb)", y = "Genome position (Mb)") +
  theme_bw(base_size = 12) +
  theme(axis.text.x = element_text(angle=0,vjust=0.5, hjust = 0.5), plot.title = element_text(hjust = 0, face = "bold"),strip.background = element_blank(),panel.grid.major = element_blank(), panel.grid.minor = element_blank(), strip.text = element_text(size = rel(1))) +
  theme(plot.title = element_text(face="bold.italic")) +
  theme(strip.background = element_blank(),panel.grid.major = element_blank(), panel.grid.minor = element_blank(), strip.text = element_text(size = rel(1)))  +
  facet_grid(chr_A ~ chr_B,scale="free") +
  stat_summary_2d(fun=median,aes(z=Mean))  +theme(legend.position = "right") +
  scale_fill_gradientn(colours = my_palette, limits=c(0, 1), name=expression(paste("LD (",r^{2},")")))


#ggsave(filename = "CE_LD_heatmap.png",pl4,width = 7.4,height = 7,dpi=300,scale=1)

tiff(filename = "CE_LD_heatmap.tiff",width = 7.4,height = 7,res=300,units="in")
plot(pl4)
dev.off()

LDCRI2$chr_B = factor(LDCRI2$chr_B, levels=c("I","II","III","IV","V","X"))
LDCRI2$chr_A = factor(LDCRI2$chr_A, levels=c("X","V","IV","III","II","I"))

pl5<-ggplot(LDCRI2, aes(x = pos_B/10, y = pos_A/10, ordered = FALSE))  +
  labs(title ="C. remanei", x = "Genome position (Mb)", y = "Genome position (Mb)") +
  theme_bw(base_size = 12) +
  theme(axis.text.x = element_text(angle=0,vjust=0.5, hjust = 0.5), plot.title = element_text(hjust = 0, face = "bold"),strip.background = element_blank(),panel.grid.major = element_blank(), panel.grid.minor = element_blank(), strip.text = element_text(size = rel(1))) +
  theme(plot.title = element_text(face="bold.italic")) +
  theme(strip.background = element_blank(),panel.grid.major = element_blank(), panel.grid.minor = element_blank(), strip.text = element_text(size = rel(1)))  +
  facet_grid(chr_A ~ chr_B, scale="free", space = "free") +
  stat_summary_2d(fun=median,aes(z=Mean))  +theme(legend.position = "right") +
  scale_fill_gradientn(colours = my_palette, name=expression(paste("LD (",r^{2},")")))


#ggsave(filename = "CR_LD_heatmap.png",pl5,width = 7.4,height = 7,dpi=300,scale=1)

tiff(filename = "CR_LD_heatmap.tiff",width = 7.4,height = 7,res=300,units="in")
plot(pl5)
dev.off()


TIF <- file.info(list.files(pattern="*_LD_heatmap.tiff"))
TIF <- TIF[with(TIF, order(as.POSIXct(mtime))), ]
#rownames(TIF)
TIF<-rownames(TIF)

combo<-image_append(do.call("c", lapply(TIF, function(h){image_read(h)})), stack = TRUE)
image_write(combo, path = "S5.tiff", format = "tiff", quality = 100)
