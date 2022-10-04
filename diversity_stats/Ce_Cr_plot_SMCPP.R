#######################################################
######### Plot Demographic histories ##################
#########  C. elegans & C.remanei    ##################
#######################################################


setwd("/Users/anastasia/Documents/Phillips_lab/drafts/CR_popgen/SMCPP")


library(ggplot2)
library(scales)
library(gridExtra)
library(boot)


my_files<-list.files(pattern="*.csv")
my_names<-gsub(".csv","", perl=T, my_files)

SMC<-c();for (i in 1:(length(my_files)) ) { B=read.csv(file=my_files[i], header=TRUE, sep=","); SMC<-rbind(SMC,data.frame(B,sample=my_names[i]));}
SMC$Species<-gsub("plot_(C..)*_aut_.*","\\1",perl=TRUE,SMC$sample)
SMC$Species<-gsub("CE8","C. elegans",SMC$Species)
SMC$Species<-gsub("CR8","C. remanei",SMC$Species)

#!!!!!read the discussion here https://github.com/popgenmethods/smcpp/issues/82
SMC[SMC$Species=="C. elegans",]$x<- SMC[SMC$Species=="C. elegans",]$x*2


mean(SMC[SMC$Species=="C. elegans" & SMC$x ==0, ]$y)
#[1] 6517.983

sd(SMC[SMC$Species=="C. elegans" & SMC$x ==0, ]$y)
#[1] 10354.84

mean(SMC[SMC$Species=="C. remanei" & SMC$x ==0, ]$y)
#[1] 247514.1

sd(SMC[SMC$Species=="C. remanei" & SMC$x ==0, ]$y)
#[1] 140341.7



# function to obtain the mean
Bmean <- function(data, indices) {
  d <- data[indices] # allows boot to select sample
  return(mean(d))
}

# bootstrapping with 1000 replications
#library(boot)
resultscr <- boot(data=SMC[SMC$Species=="C. remanei" & SMC$x ==0, ]$y, statistic=Bmean, R=1000)

# view results
resultscr

#ORDINARY NONPARAMETRIC BOOTSTRAP
#
#
#Call:
#boot(data = SMC[SMC$Species == "C. remanei" & SMC$x == 0, ]$y,
#    statistic = Bmean, R = 1000)
#
#
#Bootstrap Statistics :
#    original    bias    std. error
#t1* 247514.1 -320.1215    13645.54
plot(resultscr)

# get 95% confidence interval
bcicr<-boot.ci(resultscr, type="basic")
bcicr
#BOOTSTRAP CONFIDENCE INTERVAL CALCULATIONS
#Based on 1000 bootstrap replicates
#
#CALL :
#boot.ci(boot.out = resultscr, type = "basic")
#
#Intervals :
#Level      Basic
#95%   (219716, 275372 )
#Calculations and Intervals on Original Scale



resultsce <- boot(data=SMC[SMC$Species=="C. elegans" & SMC$x ==0, ]$y, statistic=Bmean, R=1000)

# view results
resultsce

#ORDINARY NONPARAMETRIC BOOTSTRAP
#
#
#Call:
#boot(data = SMC[SMC$Species == "C. elegans" & SMC$x == 0, ]$y,
#    statistic = Bmean, R = 1000)
#
#
#Bootstrap Statistics :
#    original   bias    std. error
#t1* 6517.983 19.05775    1051.815
plot(resultsce)

# get 95% confidence interval
bcice<-boot.ci(resultsce, type="basic")
bcice
#BOOTSTRAP CONFIDENCE INTERVAL CALCULATIONS
#Based on 1000 bootstrap replicates
#
#CALL :
#boot.ci(boot.out = resultsce, type = "basic")
#
#Intervals :
#Level      Basic
#95%   (4232, 8287 )
#Calculations and Intervals on Original Scale





FIG6B<-ggplot(SMC, aes(x = x, y = y, color= Species, fill=Species, ordered = FALSE))  +
  labs(title ="B", x = "Generations", y = "Effective population size") +
  scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),labels = trans_format(trans="log10", math_format(10^.x))) +
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),labels = trans_format(trans="log10", math_format(10^.x)),limits = c(1,0.5e7)) +
  theme_bw(base_size = 12) +
  theme(axis.text.x = element_text(angle=0,vjust=0.5, hjust = 0.5), plot.title = element_text(hjust = 0, face = "bold"),strip.background = element_blank(),panel.grid.major = element_blank(), panel.grid.minor = element_blank(), strip.text = element_text(size = rel(1))) +
  theme(legend.position = "none") +
  theme(strip.text.y = element_text(face = "italic"))+
  theme(strip.background = element_blank(),panel.grid.major = element_blank(), panel.grid.minor = element_blank(), strip.text = element_text(size = rel(1))) +
  facet_grid(Species~ .) +
  geom_rect(aes(xmin=0, xmax=1000, ymin=0, ymax=Inf), col="#e3e3e3",fill="#e3e3e3") +
  geom_line(aes(group=sample),size=0.5, alpha=0.4) +
  theme(plot.margin = unit(c(1, 3, 1.5, 3), "pt")) +
  scale_color_manual(values = c("#196770","#196770"),name="") + scale_fill_manual(values = c("#225D5B","#225D5B"),name="")

#ggsave(filename = "CR_CE_SMCPP.png",FIG6,width = 7.2,height = 4.2,dpi=300,scale=1)
#tiff(filename = "Fig5.tiff",width = 7.2,height = 4.2,res=300,units="in")
plot(FIG6B)
#dev.off()


###GONE results
my_files2<-list.files(pattern="Output_Ne_")
my_names2<-gsub("Output_Ne_","", perl=T, my_files2)
my_names2<-gsub("_data_rec","", perl=T, my_names2)

GONE <-c();
for (i in 1:(length(my_files2))) {
  B = read.csv(file = my_files2[i],
               header = TRUE,
               sep = "\t", skip = 1)
  GONE <- rbind(GONE, data.frame(B, sample = my_names2[i]))
}

GONE$Species<-gsub("CE","C. elegans",perl=TRUE,GONE$sample)
GONE$Species<-gsub("CR","C. remanei",perl=TRUE,GONE$Species)


#!!!!!read the discussion here https://github.com/popgenmethods/smcpp/issues/82
#GONE[GONE$Species=="C. elegans",]$Generation[1] <-0
GONE[GONE$Species=="C. elegans",]$Generation<- GONE[GONE$Species=="C. elegans",]$Generation*2
GONE[GONE$Species=="C. elegans",]$Generation[1] <-1
GONE[GONE$Species=="C. elegans",]$Geometric_mean <-GONE[GONE$Species=="C. elegans",]$Geometric_mean/2

FIG6A<-ggplot(GONE, aes(x = Generation, y = Geometric_mean, color= Species, fill=Species, ordered = FALSE))  +
  labs(title ="A", x = "Generations", y = "Effective population size") +
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),labels = trans_format(trans="log10", math_format(10^.x)),limits = c(1,0.5e7)) +
  theme_bw(base_size = 12) +
  theme(axis.text.x = element_text(angle=0,vjust=0.5, hjust = 0.5), plot.title = element_text(hjust = 0, face = "bold"),strip.background = element_blank(),panel.grid.major = element_blank(), panel.grid.minor = element_blank(), strip.text = element_text(size = rel(1))) +
  theme(legend.position = "none") +
  theme(strip.text.y = element_text(face = "italic"))+
  theme(strip.background = element_blank(),panel.grid.major = element_blank(), panel.grid.minor = element_blank(), strip.text = element_text(size = rel(1))) +
  facet_grid(Species~ .) +
  geom_line(aes(group=sample),size=1.1, alpha=1) + xlim(c(0,501)) +
  theme(plot.margin = unit(c(1, 3, 1.5, 3), "pt")) +
  scale_color_manual(values = c("#196770","#196770"),name="") + scale_fill_manual(values = c("#225D5B","#225D5B"),name="")

grid.arrange(FIG6A,FIG6B, nrow=1,widths=c(1,2)) ->F5
tiff(filename = "Fig5.tiff",width = 7.2,height = 3.8,res=300,units="in")
plot(F5)
dev.off()


####estimetion off the "recent Ne" by averagingthe last 100 generations






resultscrGONE <- boot(data=GONE[GONE$Species=="C. remanei" & GONE$Generation<=100, ]$Geometric_mean, statistic=Bmean, R=1000)

# view results
resultscrGONE

#ORDINARY NONPARAMETRIC BOOTSTRAP
#
#
#Call:
#  boot(data = GONE[GONE$Species == "C. remanei" & GONE$Generation <=
#                     100, ]$Geometric_mean, statistic = Bmean, R = 1000)
#
#
#Bootstrap Statistics :
#  original    bias    std. error
#t1* 116133.9 -238.9152    9047.852


resultsceGONE <- boot(data=GONE[GONE$Species=="C. elegans" & GONE$Generation<=100, ]$Geometric_mean, statistic=Bmean, R=1000)

# view results
resultsceGONE

#ORDINARY NONPARAMETRIC BOOTSTRAP
#
#
#Call:
#  boot(data = GONE[GONE$Species == "C. elegans" & GONE$Generation <=
#                     100, ]$Geometric_mean, statistic = Bmean, R = 1000)
#
#
#Bootstrap Statistics :
#  original  bias    std. error
#t1* 911025.3 -864.3987    51324.62
