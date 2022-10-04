################################
###### popVAE figure #########
################################


setwd("/Users/anastasia/Documents/Phillips_lab/drafts/CR_popgen/POPVAE")


library(ggplot2)
library(ggrepel)

CR<-read.csv("CR_popvae_FIN_latent_coords.txt", header=TRUE, sep="\t")
CE<-read.csv("CE_popvae_FIN_latent_coords.txt", header=TRUE, sep="\t")

CR$species<-"C. remanei"
CR$sample<-gsub("HI.*Index_[0-9]+.(.*)$","\\1", perl=TRUE, CR$sampleID)
CE$species <-"C. elegans"
CE$sample<-gsub("(.*)-.*$","\\1", perl=TRUE, CE$sampleID)

COMBO<-rbind(CR,CE)


pl<-ggplot(COMBO, aes(x = mean1, y = mean2, ordered = FALSE))  +
  labs(x = "LD1", y = "LD2") +
  theme_bw(base_size = 12) +
  theme(axis.text.x = element_text(angle=0,vjust=0.5, hjust = 0.5), plot.title = element_text(hjust = 0, face = "bold"),strip.background = element_blank(),panel.grid.major = element_blank(), panel.grid.minor = element_blank(), strip.text = element_text(size = rel(1))) +
  theme(legend.position = "none") +
  theme(strip.text.x = element_text(face = "italic"))+
  theme(strip.background = element_blank(),panel.grid.major = element_blank(), panel.grid.minor = element_blank(), strip.text = element_text(size = rel(1))) +
  facet_wrap(.~species,scales = "free") +
  theme(plot.margin = unit(c(1, 3, 1.5, 3), "pt")) +
  geom_text_repel(aes(label = sample), box.padding = unit(0.57, "lines"),colour = "#8a8a8a",size=3) +
  geom_point(alpha=0.55,size=2, col="#196770")

#ggsave("POPVAE_CR_CE.tiff",pl,width = 7.2,height = 5.7,dpi=300,scale=1)

tiff(filename= "S4.tiff", width = 7.2,height = 5.7,res=300,units = "in")
plot(pl)
dev.off()
