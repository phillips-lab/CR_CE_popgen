#########################################
##### C. remanei genetic map ############
##### Figures and Tables ################
#########################################



setwd("~/Documents/Phillips_lab/drafts/CR_popgen/MAP/NEWMAP_4FAM")

library(ggplot2)
library(dplyr)
library(segmented)

set.seed(789999)


### import the data
###chromosomes
chrom <-c("I","II","III","IV","V","X")
###chromosome sizes
ends<-c(17247545,19935723,17877849,25790997,22502457,21501900)

###files from LepMap3
POS<-c();for (i in chrom ) {B=read.csv(file=paste(i,".POSITIONS2.txt",sep=""), header=T, sep="\t", skip=6); B$NUM<-c(1:nrow(B)); POS<-rbind(POS,B);}
CM<-c();for (i in chrom ) {B=read.csv(file=paste(i,".REC.MAP.SHORT2.txt",sep=""), header=T, sep="\t", skip=2); B$CHR<-i; CM<-rbind(CM,B);}


colnames(CM)[1]<-"NUM"
MAP<-merge(CM,POS,by=c("NUM","CHR"),all.x=TRUE,sort=FALSE)


#only maternal map (sex-averaged map for autosomes, and female only for the sex chromosome)
MAP<-MAP[,c(2,5,4)]
colnames(MAP)<-c("CHR","POS","CM")
MAP$Type<-"Maternal"


MAP$POS<-as.numeric(as.character(MAP$POS))
MAP$CHR<-as.character(MAP$CHR)
MAP<-MAP[complete.cases(MAP),]

#########################################################
############### rate off recombinaation per domian ######
################### for Table 1 #########################
#########################################################

dompos<-data.frame(chrom=c("I","I","I","I",
                           "II","II","II","II",
                           "III","III","III","III",
                           "IV","IV", "IV","IV",
                           "V","V","V","V",
                           "X","X","X","X"),
                    pos=c(1,5803000,13007000,17247545,
                         1,6828000, 14292000,19935723,
                         1,5849000,11888000,17877849,
                         1,8639000,17108000,25790997,
                         1,6767000,13852000,22502457,
                         1,9040000,16070000,21501900
                         ),
                   stringsAsFactors = FALSE)

dompos$cM<-0
dompos$rate<-0

for (i in 1:nrow(dompos)){
  print(i);
  chr<-dompos[i,1]
  print(chr);
  SUB<-MAP[MAP$CHR==chr,]
  if(dompos$pos[i]==1){
    dompos$cM[i]<-0;
    dompos$rate[i]<-0;
  }
  else{
  dompos$cM[i]<-SUB$CM[(which.min(abs(SUB$POS - dompos[i,2])))];
  dompos$rate[i]<-((dompos$cM[i]-dompos$cM[i-1])*1000000)/(dompos$pos[i]-dompos$pos[i-1])
  }

}



####################################################
#dompos
#chrom      pos     cM      rate
#1      I        1  0.000 0.0000000
#2      I  5803000 30.159 5.1971403
#3      I 13007000 35.008 0.6730983
#4      I 17247545 49.973 3.5290275
#5     II        1  0.000 0.0000000
#6     II  6828000 26.516 3.8834218
#7     II 14292000 32.544 0.8076099
#8     II 19935723 49.126 2.9381314
#9    III        1  0.000 0.0000000
#10   III  5849000 28.024 4.7912472
#11   III 11888000 29.933 0.3161119
#12   III 17877849 51.502 3.6009255
#13    IV        1  0.000 0.0000000
#14    IV  8639000 12.504 1.4473899
#15    IV 17108000 16.911 0.5203684
#16    IV 25790997 48.301 3.6151112
#17     V        1  0.000 0.0000000
#18     V  6767000 25.015 3.6966165
#19     V 13852000 26.925 0.2695836
#20     V 22502457 48.817 2.5307333
#21     X        1  0.000 0.0000000
#22     X  9040000 27.545 3.0470136
#23     X 16070000 33.722 0.8786629
#24     X 21501900 40.998 1.3394945






##########################################################
#############  rate of recombination per Mb ##############
#############        for Relate      #####################
##########################################################


for (i in 1:length(chrom)) {
  CMw<-c();
  print(nrow(MAP[MAP$CHR == chrom[i], ]))
  MAP[MAP$CHR == chrom[i], ] %>% arrange(POS) %>% group_by(CM) %>% slice(n()) ->B
  B<-data.frame(B)

  B$POS<-as.integer(B$POS/1000000)

  B %>% arrange(CM) %>% group_by(POS) %>% slice(n()) ->B2
  print(B)
  B<-data.frame(B2)
  B$POS<-B$POS*1000000
  POS<- c(0,B[B$CHR == chrom[i], ]$POS);
  print(POS)
  for (j in 1:(length(POS)-1)) {
    CMw <-
      rbind(CMw,
            data.frame(
              pos = POS[j],
              COMBINED_rate = (B[B$CHR == chrom[i] & B$POS==POS[j+1],3] - max(B[B$CHR == chrom[i] & B$POS==POS[j],3],0)) /(POS[j+1]-POS[j]) * 1000000,
              Genetic_Map = B[B$CHR == chrom[i] & B$POS==POS[j+1],3]
            ))
    print(j)
  }
  CMw$COMBINED_rate<-as.numeric(as.character(CMw$COMBINED_rate))
  CMw$pos<-as.integer(as.character(CMw$pos))
  CMw<-CMw[complete.cases(CMw),]
  CMw$COMBINED_rate<-round(CMw$COMBINED_rate, digits = 4)
  write.table(CMw, file = paste0(chrom[i],"_genetic_map.txt"), append = FALSE, quote = FALSE, sep = " ", eol = "\n", na = "NA", dec = ".", row.names = FALSE,col.names = TRUE)
}


##########################################################
################## detect break points ###################
############### domains of low recombination #############
##########################################################


breakstable<-c()


for (j in 1:6){

  #remove the tips
  test<-MAP[ MAP$CHR==chrom[j] &  MAP$POS< 13*ends[j]/15 & MAP$POS> 2*ends[j]/15,]

  #suggest the region where to look
  breakleft<-as.integer(ends[j]*(3/7))
  breakright<-as.integer(ends[j]*(5/7))

  A<-test$CM
  x<-test$POS
  fit0 <- lm(A ~ x)
  fit1 <- segmented.glm(fit0, seg.Z = ~ x, psi = list(x = c(breakleft, breakright)))

  #get the CIs
  BR<-as.data.frame(confint(fit1))
  BR$CHR<-chrom[j]
  breakstable<-rbind(breakstable,BR)

}


colnames(breakstable)[1]<-"POS"

###position of the central domains with 95% CI
write.table(breakstable,"Domain_positions_genetic_map.txt", append = FALSE, quote = FALSE, sep = "\t", eol = "\n", na = "NA", dec = ".", row.names = FALSE,col.names = TRUE)


breakstable<-read.csv("Domain_positions_genetic_map.txt",header=T,sep="\t")

breakstable
#POS CI.95...low CI.95...up CHR
#1   5802890     5752740    5853050   I
#2  13006600    12942300   13070900   I
#3   6828110     6778710    6877500  II
#4  14292000    14060400   14523700  II
#5   5849070     5827940    5870200 III
#6  11888400    11862800   11914000 III
#7   8638790     8396690    8880890  IV
#8  17108200    17071500   17144800  IV
#9   6766730     6716890    6816580   V
#10 13851600    13737000   13966200   V
#11  9039880     8959800    9119950   X
#12 16070100    15497800   16642400   X


round(breakstable[,c(1,2,3)]/1000,digits=0) #Table 1
#     POS CI.95...low CI.95...up
#1   5803        5753       5853
#2  13007       12942      13071
#3   6828        6779       6878
#4  14292       14060      14524
#5   5849        5828       5870
#6  11888       11863      11914
#7   8639        8397       8881
#8  17108       17072      17145
#9   6767        6717       6817
#10 13852       13737      13966
#11  9040        8960       9120
#12 16070       15498      16642

##########################################################
################ Plot the genetic map ####################
##########################################################


plmapM<-ggplot(MAP, aes(x = POS/1000000, y = CM, color= Type, fill=Type, ordered = FALSE))  +
  facet_grid(~CHR,scales="free", space="free") +
  theme_bw(base_size = 12) +
  labs(title ="", x = "Physical position (Mb)", y = "Genetic position (cM)") +
  theme(axis.text.x = element_text(angle=0,vjust=0.5, hjust = 0.5), plot.title = element_text(hjust = 0.5), strip.text = element_text(size = rel(1))) +
  theme(legend.position = "none") +
  scale_colour_manual(values = c("#196770","black"),name="") +
  theme(strip.background = element_blank(),panel.grid.major = element_blank(), panel.grid.minor = element_blank())  +
  geom_vline(data = breakstable,aes(xintercept = POS/1000000),colour = "#bfbfbf",size = 0.75,linetype = "dashed") +
  geom_rug(data=MAP, sides = "t", alpha=0.2,size=0.15,colour = "black") +
  geom_line(size=1.1)  +
  xlim(0,NA)

ggsave(filename = "CR_map.png",plmapM,width = 7.2,height = 2.7,dpi=300,scale=1)

#ggsave should be rescaled and works fine only with pdf, so tiff() it is
tiff(filename = "Fig1.tiff",width = 7.2,height = 2.7,units="in",res=300)
plot(plmapM)
dev.off()


#################################
##### Number of recombinations ##
#################################

##LepMap3 log only lines that starts with "Individual"
REC<-read.csv("Final_NUMBER_of_Recombination_SHORT.txt",sep="\t",header=F)
REC<- REC[,-c(1,4,6)]

colnames(REC)<-c("Family","Individual","Value")

table(REC$Family)
#A1  A2  B1  B2
#260 295 321 354


REC$CHR<-"N"
REC$CHR[1:222]<-"I"
REC$CHR[223:438]<-"II"
REC$CHR[439:664]<-"III"
REC$CHR[665:883]<-"IV"
REC$CHR[884:1101]<-"V"
REC$CHR[1102:1230]<-"X"


table(REC$CHR)
#I  II III  IV   V   X
#222 216 226 219 218 129

colorsMUT<-c("#183843", "#196770", "#73A8BD", "#C0BDBC", "#e8b474", "#a55f22")

pl<-ggplot(REC, aes(x = Value, color= Family, fill=Family, ordered = FALSE))  +
  facet_grid(Family~CHR,scales="free", space="free") +
  theme_bw(base_size = 12) +
  labs(title ="Recombination in crosses\n♀PX506 x ♂PX530 (A1, A2) and ♀PX530 x ♂PX506 (B1, B2)", y = "Count", x = "Number of Crossovers") +
  theme(axis.text.x = element_text(angle=0,vjust=0.5, hjust = 0.5), plot.title = element_text(hjust = 0.5), strip.text = element_text(size = rel(1))) +
  theme(legend.position = "none") +
  scale_colour_manual(values = colorsMUT,name="") + scale_fill_manual(values = colorsMUT,name="") +
  theme(strip.background = element_blank(),panel.grid.major = element_blank(), panel.grid.minor = element_blank())  + geom_bar()

#ggsave("Number_of_recombinations.tiff",pl,width = 7.2,height = 5.77,dpi=300,scale=1)
#ggsave("Number_of_recombinations.png",pl,width = 7.2,height = 5.77,dpi=300,scale=1)

tiff(filename ="S1.tiff",width = 7.2,height = 5.77,res=300, units="in" )
plot(pl)
dev.off()
