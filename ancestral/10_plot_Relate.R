#############################################
########### Plot Relate results #############
#############################################


setwd("~/Documents/Phillips_lab/drafts/CR_popgen/RELATE")

library(scales)
library(ggplot2)
library(lsr)
library(coin)
library(magick)
library(harmonicmeanp)

my_files<-list.files(pattern="*.coal*")

my_names<-gsub("_COAL_2it.coal","", perl=TRUE, my_files)

POP<-c();for (i in 1:(length(my_files)) ) { B=t(read.csv(file=my_files[i], header=F, sep=" ")); POP<-rbind(POP,data.frame(B,sample=my_names[i]));}

POP<-POP[,c(2,3,4)]
POP<-POP[complete.cases(POP), ]
POP$X3<-as.numeric(as.character(POP$X3))
POP$X2<-as.numeric(as.character(POP$X2))
#to get the rate for diploids:
POP$Ne<-0.5/POP$X3

POP$CHR<-gsub("CR_Relate_chr_","",POP$sample)


colorsMUT<-c("#183843", "#196770", "#73A8BD", "#C0BDBC", "#f2c76f", "#916b20")

FIG7<-  ggplot(POP, aes(x = X2, y = Ne, color= CHR, fill=CHR, ordered = FALSE))  +
  labs(title ="A", x = "Generations", y = "Effective population size") +
  scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),labels = trans_format(trans="log10", math_format(10^.x))) +
  scale_y_log10( breaks = c(10000,100000,1000000,10000000),labels = trans_format(trans="log10", math_format(10^.x))) +
  theme_bw(base_size = 12) +
  theme(axis.text.x = element_text(angle=0,vjust=0.5, hjust = 0.5), plot.title = element_text(hjust = 0, face = "bold"),strip.background = element_blank(),panel.grid.major = element_blank(), panel.grid.minor = element_blank(), strip.text = element_text(size = rel(1))) +
  theme(legend.position = "right") +
  theme(strip.text.y = element_text(face = "italic"))+
  theme(strip.background = element_blank(),panel.grid.major = element_blank(), panel.grid.minor = element_blank(), strip.text = element_text(size = rel(1))) +
  geom_line(aes(group=sample),size=1.1, alpha=0.75) +
  theme(plot.margin = unit(c(1, 3, 1.5, 3), "pt")) +
  scale_color_manual(values = colorsMUT,name="") + scale_fill_manual(values = colorsMUT,name="")

#ggsave(filename = "CR_Relate_demography.tiff",FIG7,width = 7.4,height = 2.7,dpi=300,scale=1)
#ggsave(filename = "CR_Relate_demography.png",FIG7,width = 7.4,height = 2.7,dpi=300,scale=1)

tiff(filename = "CR_Relate_demography2.tiff",width = 7.4,height = 2.7,res=300,units="in")
plot(FIG7)
dev.off()

##################################################
####### Positive selection #######################
#################################################


#from Relate
my_files2<-list.files(pattern="*_POSSEL_SEL.sele*")

my_names2<-gsub("_POSSEL_SEL.sele","", perl=TRUE, my_files2)



SEL<-c();
for (i in 1:(length(my_files2)) ) {
  A<-read.csv(file=my_files2[i], header=TRUE, sep=" ");
  #A$when_mutation_has_freq2<-NULL;
  A$rs_id<-NULL
  H<-gsub("X","",colnames(A))
  H<-H[-1]
  H<-H[-length(H)]
  H<-as.integer(as.numeric(as.character(H)))
  for(j in 1:length(H)){
    B<-A[,c(1,j+1,ncol(A))]
    colnames(B)<-c("pos","pval","totalpval")

    SEL<-rbind(SEL,data.frame(B,sample=my_names2[i],generation=H[j]));
  }
}
SEL$CHR<-gsub("CR_Relate_chr_","",SEL$sample)
SEL<-SEL[complete.cases(SEL),]



TOTALSEL<-c();

for (i in 1:(length(my_files2)) ) {
  A<-read.csv(file=my_files2[i], header=TRUE, sep=" ");
  B<-A[,c(1,ncol(A))];
  colnames(B)<-c("pos","totalpval");
  TOTALSEL<-rbind(TOTALSEL,data.frame(B,sample=my_names2[i]));
}


TOTALSEL$CHR<-gsub("CR_Relate_chr_","",TOTALSEL$sample)
domainsCR<-data.frame(pos=c(5803000,13007000,6828000, 14292000,5849000,11888000,8639000,17108000,6767000,13852000,9040000,16070000), chrom=c("I","I","II","II","III","III","IV","IV","V","V","X","X"), Species=rep("C. remanei",12))

domainsCR$POS<-as.integer(domainsCR$pos/1000000)

domainsCR2<-domainsCR
colnames(domainsCR2)[2]<-"CHR"

#https://cran.r-project.org/web/packages/harmonicmeanp/vignettes/harmonicmeanp.html
#library(harmonicmeanp)
# Specify the false positive rate
L=nrow(TOTALSEL)
alpha = 0.05
TOTALSEL$w <- 1/L
TOTALSEL$p <- 10^(TOTALSEL$totalpval)


chrom <-c("I","II","III","IV","V","X")
###chromosome sizes
endsCR<-c(17247545,19935723,17877849,25790997,22502457,21501900)

HMP<-c()

for (i in 1:6) {
#for (i in 1) {
  print(chrom[i])
  #bins
  for (r in 1:round(endsCR[i] / 1000000, 0)) {
    print((r) * 1000000 + 1000000);
    print((r - 1) * 1000000 - 1000000);

    A <-
      TOTALSEL[TOTALSEL$CHR == chrom[i] &
                 TOTALSEL$pos < (r) * 1000000 + 1000000 & TOTALSEL$pos >= (r - 1) * 1000000 - 1000000, ]
                 PV<-p.hmp(A$p,A$w,L)
    # Calculate sums of weights for each combined test
    W<-   sum(A$w)
    HMP<-rbind(HMP,data.frame(CHR=chrom[i],POS=r * 1000000, HMPval=PV/W))
  }
}



pl<-ggplot(HMP, aes(x = POS/1000000, y = -log10(HMPval), ordered = FALSE))  +
  labs(title ="C", x = "Genome position (Mb)", y = "Selection (HMP)") +
  scale_y_continuous(labels = function(x) format(x, scientific = FALSE)) +
  theme_bw(base_size = 12) +
  theme(axis.text.x = element_text(angle=0,vjust=0.5, hjust = 0.5), plot.title = element_text(hjust = 0, face = "bold"),strip.background = element_blank(),panel.grid.major = element_blank(), panel.grid.minor = element_blank(), strip.text = element_text(size = rel(1))) +
  theme(legend.position = "none") +
  theme(strip.text.y = element_text(face = "italic"))+
  theme(strip.background = element_blank(),panel.grid.major = element_blank(), panel.grid.minor = element_blank(), strip.text = element_text(size = rel(1))) +
  facet_grid(Species ~ CHR, scales = "free" ) +
  geom_vline(data = domainsCR2,aes(xintercept = pos/1000000),colour = "#bfbfbf",size = 0.75,linetype = "dashed") +
  geom_hline(yintercept = 1.5,colour = "#262626",size = 0.75,linetype = "dotdash") +
 geom_point(alpha=0.6,size=0.6, col="#196770") +
  theme(plot.margin = unit(c(1, 3, 1.5, 3), "pt"))

pl
#ggsave("Relate_CR_FDR.png",pl,width = 7.2,height = 2.77,dpi=300,scale=1)
#ggsave("Relate_CR_CE_FDR",pl,width = 7.2,height = 2.77,dpi=300,scale=1)


tiff(filename="Relate_CR_HMP.tiff",width = 7.2,height = 2.77,res=300,units="in")
plot(pl)
dev.off()






TOTALSEL$Type<-"Arm"
for (chr in 1:length(levels(as.factor(TOTALSEL$CHR)))){
  TOTALSEL[TOTALSEL$CHR==levels(as.factor(TOTALSEL$CHR))[chr] & TOTALSEL$pos>domainsCR2$pos[chr*2-1] & TOTALSEL$pos<domainsCR2$pos[chr*2],]$Type<-"Center"
}



unipval<-runif(n = length(TOTALSEL$totalpval))
unipvalA<-runif(n = length(TOTALSEL[TOTALSEL$Type=="Arm",]$totalpval))
unipvalC<-runif(n = length(TOTALSEL[TOTALSEL$Type=="Center",]$totalpval))


QQdata <-data.frame(Exp = c(sort(-log10(unipvalA)), sort(-log10(unipvalC))),Obs = c(sort(-TOTALSEL[TOTALSEL$Type == "Arm", ]$totalpval),sort(-TOTALSEL[TOTALSEL$Type == "Center", ]$totalpval)),Type = c(rep("Arm", as.numeric(table(TOTALSEL$Type)[1])), rep("Center", as.numeric(table(TOTALSEL$Type)[2]))), stringsAsFactors = FALSE)

##########################
###### stats #############
##########################


cohensD(totalpval ~Type, data=TOTALSEL, method="corrected")
#[1] 0.0839197


TOTALSEL$Type<-as.factor(TOTALSEL$Type)
independence_test(totalpval ~Type, data=TOTALSEL, alternative="less",distribution=approximate(nresample=10000))

#	Approximative General Independence Test
#
#data:  totalpval by Type (Arm, Center)
#Z = -23.473, p-value < 1e-04
#alternative hypothesis: less



###QQplot

ggplot(QQdata, aes(x = Exp, y = Obs, color=Type,ordered = FALSE))  +
  labs(title ="D", x = "Expected -log10(p-value)", y = "Observed -log10(p-value)") +
  scale_y_continuous(breaks = seq(0, 6, by = 2)) +
  scale_x_continuous(breaks = seq(0, 6, by = 2)) +
  theme_bw(base_size = 12) +
  theme(axis.text.x = element_text(angle=0,vjust=0.5, hjust = 0.5), plot.title = element_text(hjust = 0, face = "bold"),strip.background = element_blank(),panel.grid.major = element_blank(), panel.grid.minor = element_blank(), strip.text = element_text(size = rel(1))) +
  theme(legend.position = c(0.75,1.1)) +
  theme(legend.key =element_blank(),legend.background=element_blank()) +
  guides(colour = guide_legend(nrow = 1)) +
  theme(strip.background = element_blank(),panel.grid.major = element_blank(), panel.grid.minor = element_blank(), strip.text = element_text(size = rel(1))) +
  geom_abline(intercept = 0,colour = "#bfbfbf",size = 0.75,linetype = "dashed") +
  geom_point(alpha=0.45,size=0.75) +
  scale_color_manual(values=c("#f2c76f","#183843"), name="") +
  scale_fill_manual(values=c("#f2c76f","#183843"), name="") +
  theme(plot.margin = unit(c(1, 3, 1.5, 3), "pt")) ->qq

#ggsave(filename = "QQplot_Relate.tiff",qq,width = 3.5,height = 3.5,dpi=300,scale=1)
#ggsave(filename = "QQplot_Relate.png",qq,width = 3.5,height = 3.5,dpi=300,scale=1)

tiff(filename = "QQplot_Relate.tiff",width = 3.5,height = 3.5,res=300,units="in")
plot(qq)
dev.off()


############ trees from Relate
library(phangorn)
library(phytools)

setwd("~/Documents/Phillips_lab/drafts/CR_popgen/RELATE/TREES/")

treelist <- list.files(pattern = "tree.newick")

TR <- c()
COMBOTR <- c()

for (i in 1:(length(treelist))) {
  #for (i in 1) {
  name <- treelist[i]
  print(name)
  TR <- read.tree(treelist[i])
  TMRCA <- c()
  RTH <- c()
  TERM <- c()
  #  for (tr in 1) {
  for (tr in 1:length(TR)) {
    print(tr)
    #add TMRCA to the list
    TMRCA <- c(TMRCA, max(nodeHeights(TR[[tr]])))
    #add RTH to the list
    #2*n-1= 2*28 -1 = 55 - the number of nodes
    #n + n/2 - 1 # number of nodes when a half of tips coalesced if to order from tips, 28+28/2 -1 =41
    #or it's n/2 of to count from the root, #14
    #as they are diploids, the number of haplotypes is always odd


    #but I sort in the reverse order, so it should be 55-41=14
    RTH <-
      c(RTH, (max(node.depth.edgelength(TR[[tr]])) - rev(sort(node.depth.edgelength(TR[[tr]])))[41]) / max(node.depth.edgelength(TR[[tr]])))

    #add all terminal tips to the list
    tips <- TR[[tr]]$tip.label
    nodes <-
      sapply(tips, function(x, y)
        which(y == x), y = TR[[tr]]$tip.label)
    ## then get the edge lengths for those nodes
    TERM <-
      c(TERM, TR[[tr]]$edge.length[sapply(nodes, function(x, y)
        which(y == x), y = TR[[tr]]$edge[, 2])])

  }

  #estimate median for everything
  COMBOTR <-
    rbind(COMBOTR, data.frame(
      Name = name,
      Stat = "TMRCA",
      Value = median(TMRCA)
    ))
  COMBOTR <-
    rbind(COMBOTR,
          data.frame(
            Name = name,
            Stat = "Num.trees",
            Value = length(TR)
          ))
  COMBOTR <-
    rbind(COMBOTR, data.frame(
      Name = name,
      Stat = "RTH",
      Value = median(RTH)
    ))
  COMBOTR <-
    rbind(COMBOTR,
          data.frame(
            Name = name,
            Stat = "Terminal",
            Value = median(TERM)
          ))

}




COMBOTR$chrom<-gsub("CR_Relate_chr_(.*)_(.*)_(.*)_tree.newick","\\1",COMBOTR$Name)
COMBOTR$pos<-gsub("CR_Relate_chr_(.*)_(.*)_(.*)_tree.newick","\\2",COMBOTR$Name)
COMBOTR$pos<-as.numeric(as.character(COMBOTR$pos))
COMBOTR$Stat<-gsub("Terminal","Terminal\nbranches", COMBOTR$Stat)
COMBOTR$Stat<-gsub("Terminal\nbranches", "Terminal\nbranches (gen)", COMBOTR$Stat)
COMBOTR$Stat<-gsub("TMRCA", "TMRCA (gen)", COMBOTR$Stat)
COMBOTR$Stat<-factor(COMBOTR$Stat,levels=c("TMRCA (gen)","RTH","Terminal\nbranches (gen)","Num.trees"))
domains<-data.frame(pos=c(5803000,13007000,6828000, 14292000,5849000,11888000,8639000,17108000,6767000,13852000,9040000,16070000), chrom=c("I","I","II","II","III","III","IV","IV","V","V","X","X"))




treepl<-ggplot(COMBOTR[COMBOTR$Stat!="Num.trees",], aes(x = pos/1000000, y = Value, ordered = FALSE))  +
  labs(title ="B", x = "Genome position (Mb)", y = "Value") +
  scale_y_continuous(labels = function(x) format(x, scientific = FALSE)) +
  theme_bw(base_size = 12) +
  theme(axis.text.x = element_text(angle=0,vjust=0.5, hjust = 0.5), plot.title = element_text(hjust = 0, face = "bold"),strip.background = element_blank(),panel.grid.major = element_blank(), panel.grid.minor = element_blank(), strip.text = element_text(size = rel(1))) +
  theme(legend.position = "none") +
  theme(strip.background = element_blank(),panel.grid.major = element_blank(), panel.grid.minor = element_blank(), strip.text = element_text(size = rel(1))) +
  facet_grid(Stat ~ chrom,scales = "free" ) +
  geom_vline(data = domains,aes(xintercept = pos/1000000),colour = "#bfbfbf",size = 0.75,linetype = "dashed") +
  geom_point(alpha=0.2,size=0.6, col="#196770") +
  geom_smooth(col="#196770",alpha=1,size=1.1,se =FALSE,method = 'loess',span=0.4) + theme(plot.margin = unit(c(1, 3, 1.5, 3), "pt"))


setwd("~/Documents/Phillips_lab/drafts/CR_popgen/RELATE")

tiff(filename = "tree_stats_Relate.tiff",width = 7.2,height = 4.5,res=300,units="in")
plot(treepl)
dev.off()

#TIF<-c("CR_Relate_demography.tiff","tree_stats_Relate2.tiff", "Relate_CR_FDR.tiff","QQplot_Relate.tiff")
TIF<-c("CR_Relate_demography.tiff","tree_stats_Relate.tiff", "Relate_CR_HMP.tiff","QQplot_Relate.tiff")

combo<-image_append(do.call("c", lapply(TIF, function(h){image_read(h)})), stack = TRUE)
image_write(combo, path = "S8.tiff", format = "tiff", quality = 100)
