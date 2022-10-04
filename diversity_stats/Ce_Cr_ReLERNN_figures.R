################################
###### ReLERNN results #########
################################

setwd("/Users/anastasia/Documents/Phillips_lab/drafts/CR_popgen/ReLERNN")

library(ggplot2)
library(lsr)
library(coin)

my_files<-list.files(pattern="*BSCORRECTED")
my_names<-gsub("phased_combined.PREDICT.BSCORRECTED.","", perl=T, my_files)

LD<-c();
for (i in 1:(length(my_files)) ) {
  B=read.csv(file=my_files[i], header=T, sep="\t");
  LD<-rbind(LD,data.frame(B,sample=my_names[i]));
}




LD$Species<-gsub("^([CER]+)_.*","\\1",perl=T,LD$sample)
LD$type<-gsub("^[CER]+_(.*).txt","\\1",perl=T,LD$sample)
LD$Species<-gsub("CE","C. elegans", LD$Species)
LD$Species<-gsub("CR","C. remanei", LD$Species)


chrom <-c("I","II","III","IV","V","X")
###chromosome sizes
endsCR<-c(17247545,19935723,17877849,25790997,22502457,21501900)
endsCE<-c(15331301,15525148,14108536,17759200,21243235,18110855)

LD$NORMPOS<-888

for (i in 1:6){

  LD[LD$Species=="C. elegans" & LD$chrom==chrom[i],]$NORMPOS <- LD[LD$Species=="C. elegans" & LD$chrom==chrom[i],]$start/endsCE[i]
  LD[LD$Species=="C. remanei" & LD$chrom==chrom[i],]$NORMPOS <- LD[LD$Species=="C. remanei" & LD$chrom==chrom[i],]$start/endsCR[i]

}

# the central domains
domains<-data.frame(start=c(3858000,11040000,4879000,12020000,3722000,10340000,3896000,12970000,5897000,16550000,6137000,12480000,5803000,13007000,6828000, 14292000,5849000,11888000,8639000,17108000,6767000,13852000,9040000,16070000), chrom=rep(c("I","I","II","II","III","III","IV","IV","V","V","X","X"),2), Species=c(rep("C. elegans",12),rep("C. remanei",12)))


domains$NORMPOS<-888

for (i in 1:6){
  domains[domains$Species=="C. elegans" & domains$chrom==chrom[i],]$NORMPOS <- domains[domains$Species=="C. elegans" & domains$chrom==chrom[i],]$start/endsCE[i]
  domains[domains$Species=="C. remanei" & domains$chrom==chrom[i],]$NORMPOS <- domains[domains$Species=="C. remanei" & domains$chrom==chrom[i],]$start/endsCR[i]

}

ENDS<-data.frame(ends=c(endsCE,endsCR),Species=c(rep("C. elegans",6),rep("C. remanei", 6)))


#re-estimate if 100-kb window sizes
NEWLD <- c()
for (sp in c("C. elegans", "C. remanei")) {
  print(sp)
  for (i in 1:6) {
    print(i)
    #bins
    for (r in 1:round(ENDS[ENDS$Species == sp, "ends"][i] / 100000, 0)) {
      A <-
        LD[LD$Species == sp &
             LD$chrom == chrom[i] &
             LD$start < (r) * 100000 & LD$end >= (r - 1) * 100000, ]

      if (nrow(A) > 0) {
        newrate = 0
        SIZE=0
        for (wind in 1:nrow(A)) {
          if (A$end[wind] < r * 100000) {

            newrate <<-
              newrate + (A$end[wind] - A$start[wind]) * A$recombRate[wind]
            SIZE <<- SIZE + (A$end[wind] - A$start[wind])
          }
          else if (A$end[wind] >= r * 100000) {

            newrate <<-
              newrate + (r * 100000 - A$start[wind]) * A$recombRate[wind]
            SIZE <<- SIZE + (r * 100000 - A$start[wind])
          }
        }
        newrate <- newrate/SIZE
      }
      else{
        newrate = "NA"
      }
      NEWLD <-
        rbind(
          NEWLD,
          data.frame(
            Species = sp,
            chrom = chrom[i],
            start = as.integer((r - 1) * 100000),
            end = as.integer(r * 100000 - 1),
            recomb = newrate
          )
        )
    }
  }
}


NEWLD$NORMPOS<-888

for (i in 1:6){

  NEWLD[NEWLD$Species=="C. elegans" & NEWLD$chrom==chrom[i],]$NORMPOS <- NEWLD[NEWLD$Species=="C. elegans" & NEWLD$chrom==chrom[i],]$start/endsCE[i]
  NEWLD[NEWLD$Species=="C. remanei" & NEWLD$chrom==chrom[i],]$NORMPOS <- NEWLD[NEWLD$Species=="C. remanei" & NEWLD$chrom==chrom[i],]$start/endsCR[i]

}

NEWLD<-NEWLD[complete.cases(NEWLD),]
NEWLD$recomb<- as.numeric(as.character(NEWLD$recomb))

#filterout the windows with bad data!!






setwd("~/Documents/Phillips_lab/drafts/CR_popgen/CR_CE_pi_draft/NEW_HARD_MASK_CE_CR/")




#diploSHIC statistics
my_files<-list.files(pattern=".100K.0.1.stats")
my_names<-gsub(".100K.0.1.stats","", perl=T, my_files)

STATS<-c();for (i in 1:(length(my_files)) ) { B=read.csv(file=my_files[i], header=T, sep="\t"); STATS<-rbind(STATS,data.frame(B,sample=my_names[i]));}

STATS$Species<-gsub("^([CER]+)_.*","\\1",perl=T,STATS$sample)
STATS$Species<-gsub("CE","C. elegans", STATS$Species)
STATS$Species<-gsub("CR","C. remanei", STATS$Species)


STATS$classifiedWinStart<-STATS$classifiedWinEnd - 50000

colnames(STATS)[2]<-"start"


MDATA<-merge(STATS[,c("chrom","start","Species")],NEWLD,by=c("chrom","start","Species"),all.x = TRUE,sort = F)


### ReleRNN outputs r. not rho, so I should multiplu it by 4Ne
### Ne for C. remanei is 1.4x10^6
### Ne for C. elegans is 3.4x10^4

MDATA$rho<-0
MDATA[MDATA$Species=="C. elegans",]$rho <-MDATA[MDATA$Species=="C. elegans",]$recomb*4*3.4*10000
MDATA[MDATA$Species=="C. remanei",]$rho <-MDATA[MDATA$Species=="C. remanei",]$recomb*4*1.4*1000000



 pl<- ggplot(MDATA, aes(x = NORMPOS, y = rho, ordered = FALSE))  +
  labs(title ="", x = "Normalized genome position", y =  bquote("Recombination rate ("~ rho == 4*N[e]*r ~")")) +
  scale_y_continuous(labels = function(x) format(x, scientific = TRUE)) +
  scale_x_continuous(breaks=c(0,0.5,1),labels=c(0,0.5,1)) +
  theme_bw(base_size = 12) +
  theme(axis.text.x = element_text(angle=0,vjust=0.5, hjust = 0.5), plot.title = element_text(hjust = 0, face = "bold"),strip.background = element_blank(),panel.grid.major = element_blank(), panel.grid.minor = element_blank(), strip.text = element_text(size = rel(1))) +
  theme(legend.position = "none") +
  theme(strip.text.y = element_text(face = "italic"))+
  theme(strip.background = element_blank(),panel.grid.major = element_blank(), panel.grid.minor = element_blank(), strip.text = element_text(size = rel(1))) +
  facet_grid(Species ~ chrom,scales = "free" ) +
  geom_vline(data = domains,aes(xintercept = NORMPOS),colour = "#bfbfbf",size = 0.75,linetype = "dashed") +
  geom_point(alpha=0.2,size=0.6, col="#196770") +
  geom_smooth(col="#196770",alpha=1,size=1.1,se =FALSE,method = 'loess',span=0.4) + theme(plot.margin = unit(c(1, 3, 1.5, 3), "pt"))

pl

#ggsave("ReLERNN_CE_0.5M_2T_CR_M_T_100Kb_norm.png",pl,width = 7.2,height = 3.77,dpi=300,scale=1)
 setwd("/Users/anastasia/Documents/Phillips_lab/drafts/CR_popgen/ReLERNN")

tiff(filename ="S6.tiff",width = 7.2,height = 3.77,res=300,units="in" )
plot(pl)
dev.off()

NEWLD2<-NEWLD

NEWLD<-MDATA


###########################################
##### the difference between arms and centers
############################################


NEWLD$CHRTYPE<-"Arm"
for (chr in 1:length(levels(as.factor(as.character(NEWLD$chrom))))){
  NEWLD[NEWLD$Species=="C. elegans" & NEWLD$chrom==levels(as.factor(as.character(NEWLD$chrom)))[chr] & NEWLD$start>domains$start[chr*2-1] & NEWLD$start<domains$start[chr*2],]$CHRTYPE<-"Center"
  NEWLD[NEWLD$Species=="C. remanei" & NEWLD$chrom==levels(as.factor(as.character(NEWLD$chrom)))[chr] & NEWLD$start>domains$start[12+chr*2-1] & NEWLD$start<domains$start[12+chr*2],]$CHRTYPE<-"Center"

}

NEWLD<-NEWLD[complete.cases(NEWLD),]
NEWLD$CHRTYPEN<-as.factor(NEWLD$CHRTYPE)

mean(NEWLD[NEWLD$Species=="C. remanei",]$rho)/mean(NEWLD[NEWLD$Species=="C. elegans",]$rho)
#[1] 8.366912

mean(NEWLD[NEWLD$Species=="C. elegans",]$rho)
#[1] 0.0003844049
sd(NEWLD[NEWLD$Species=="C. elegans",]$rho)
#[1] 0.0002693414

mean(NEWLD[NEWLD$Species=="C. remanei",]$rho)
#[1] 0.003216282
sd(NEWLD[NEWLD$Species=="C. remanei",]$rho)
#[1] 0.007726469


cohensD(rho ~CHRTYPE, data=NEWLD[NEWLD$Species=="C. elegans",], method="corrected")
#[1] 0.6833575


cohensD(rho ~CHRTYPE, data=NEWLD[NEWLD$Species=="C. remanei",], method="corrected")
#[1] 0.1490841


oneway_test(rho ~ CHRTYPEN, data=NEWLD[NEWLD$Species=="C. elegans",],alternative="two.sided",distribution=approximate(nresample=10000))

#Approximative Two-Sample Fisher-Pitman Permutation Test
#
#data:  recomb by CHRTYPEN (Arm, Center)
#Z = -9.5447, p-value < 1e-04
#alternative hypothesis: true mu is not equal to 0



oneway_test(rho ~ CHRTYPEN, data=NEWLD[NEWLD$Species=="C. remanei",],alternative="two.sided",distribution=approximate(nresample=10000))

#Approximative Two-Sample Fisher-Pitman Permutation Test
#
#data:  recomb by CHRTYPEN (Arm, Center)
#Z = 2.1787, p-value = 0.0262
#alternative hypothesis: true mu is not equal to 0
