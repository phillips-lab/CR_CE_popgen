#######################################################
######### Plot Demographic histories ##################
#########  C. elegans & C.remanei    ##################
#########       per chromosome       ##################
#######################################################


setwd("/Users/anastasia/Documents/Phillips_lab/drafts/CR_popgen/SMCPP/CHR")

library(ggplot2)
library(scales)


my_files<-list.files(pattern="*.csv")
my_names<-gsub(".csv","", perl=T, my_files)

SMC<-c();for (i in 1:(length(my_files)) ) { B=read.csv(file=my_files[i], header=TRUE, sep=","); SMC<-rbind(SMC,data.frame(B,sample=my_names[i]));}
SMC$Species<-gsub("plot_(C..)*_aut_.*","\\1",perl=TRUE,SMC$sample)
SMC$Species<-gsub("CE8","C. elegans",SMC$Species)
SMC$Species<-gsub("CR8","C. remanei",SMC$Species)
SMC$CHR<-gsub(".*_([IVX]+)$", "\\1",perl=TRUE,SMC$sample)

#!!!!!read the discussion here https://github.com/popgenmethods/smcpp/issues/82
SMC[SMC$Species=="C. elegans",]$x<- SMC[SMC$Species=="C. elegans",]$x*2


mean(SMC[SMC$Species=="C. elegans" & SMC$x ==0, ]$y)
#[1] 223070.8
mean(SMC[SMC$Species=="C. remanei" & SMC$x ==0, ]$y)
#[1] 913306.9
sd(SMC[SMC$Species=="C. elegans" & SMC$x ==0, ]$y)
#[1] 2030904
sd(SMC[SMC$Species=="C. remanei" & SMC$x ==0, ]$y)
#[1] 2808573

colorsMUT<-c("#183843", "#196770", "#73A8BD", "#C0BDBC", "#f2c76f", "#916b20")

FIGS6<-
ggplot(SMC, aes(x = x, y = y, color= CHR, fill=CHR, ordered = FALSE))  +
  labs(title ="", x = "Generations", y = "Effective population size") +
  scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),labels = trans_format(trans="log10", math_format(10^.x))) +
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),labels = trans_format(trans="log10", math_format(10^.x))) +
  theme_bw(base_size = 12) +
  theme(axis.text.x = element_text(angle=0,vjust=0.5, hjust = 0.5), plot.title = element_text(hjust = 0, face = "bold"),strip.background = element_blank(),panel.grid.major = element_blank(), panel.grid.minor = element_blank(), strip.text = element_text(size = rel(1))) +
  theme(legend.position = "right") +
  theme(strip.text.y = element_text(face = "italic"))+
  theme(strip.background = element_blank(),panel.grid.major = element_blank(), panel.grid.minor = element_blank(), strip.text = element_text(size = rel(1))) +
  facet_grid(Species~ .) +
  geom_rect(aes(xmin=0, xmax=1000, ymin=0, ymax=Inf), col="#e3e3e3",fill="#e3e3e3") +
  geom_line(aes(group=sample),size=0.5, alpha=0.6) +
  theme(plot.margin = unit(c(1, 3, 1.5, 3), "pt")) +
  scale_color_manual(values = colorsMUT,name="") + scale_fill_manual(values = colorsMUT,name="")

tiff(filename = "S7.tiff",width = 7.2,height = 4.2,res=300,units="in")
plot(FIGS6)
dev.off()
