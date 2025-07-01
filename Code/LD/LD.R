# Filter VCF -------------------------------------------------------------------
accession_info <- read.csv('/Users/kehan/Library/CloudStorage/OneDrive-UniversityofCalifornia,Davis/Monroe_Lab/Cassava_Center/Tables/S1_Cassava_accessions_meta_information.csv')

Wild_noHuKistler<-accession_info$ID[which((accession_info$Biological.status=="Wild" | accession_info$Biological.status=="Other Manihot") & accession_info$Source!="Hu et al. 2021" & accession_info$Source!="Kistler et al. 2025")]
South_america_breeding_noHuKistler<-accession_info$ID[which(accession_info$Continent=="South America" & (accession_info$Biological.status=="Breeding line" | accession_info$Biological.status=="Domesticated") & accession_info$Species=="esculenta" & accession_info$Source!="This study" & accession_info$Source!="Hu et al. 2021" & accession_info$Source!="Kistler et al. 2025")]
Colombia_landraces<-accession_info$ID[which(accession_info$Source=="This study" & accession_info$Species=="esculenta" & accession_info$Biological.status=="Landrace")]
Africa_breeding<-accession_info$ID[which(accession_info$Continent=="Africa" & (accession_info$Biological.status=="Breeding line" | accession_info$Biological.status=="Domesticated") & accession_info$Species=="esculenta" & accession_info$Source!="This study")]

write(Wild_noHuKistler, '/Users/kehan/Library/CloudStorage/OneDrive-UniversityofCalifornia,Davis/Monroe_Lab/Cassava_Center/Accession_list_filtering/Wild_noHuKistler.txt')
write(South_america_breeding_noHuKistler, '/Users/kehan/Library/CloudStorage/OneDrive-UniversityofCalifornia,Davis/Monroe_Lab/Cassava_Center/Accession_list_filtering/South_america_breeding_noHuKistler.txt')
write(Colombia_landraces, '/Users/kehan/Library/CloudStorage/OneDrive-UniversityofCalifornia,Davis/Monroe_Lab/Cassava_Center/Accession_list_filtering/Colombia_landraces.txt')
write(Africa_breeding, '/Users/kehan/Library/CloudStorage/OneDrive-UniversityofCalifornia,Davis/Monroe_Lab/Cassava_Center/Accession_list_filtering/Africa_breeding.txt')


# Run VCF filtering and LD calculation using Plink on HPC 

# LD Decay ---------------------------------------------------------------------
Colombia_landraces_LD<-read.table("/Users/kehan/Library/CloudStorage/OneDrive-UniversityofCalifornia,Davis/Monroe_Lab/Cassava_Center/Code/LD/LD_Decay/Colombia_landraces_accesions_filtered_window100kb_thincount50000_may2025.ld",header = TRUE)
South_America_breeding_noHuKistler_LD<-read.table("/Users/kehan/Library/CloudStorage/OneDrive-UniversityofCalifornia,Davis/Monroe_Lab/Cassava_Center/Code/LD/LD_Decay/South_america_breeding_noHuKistler_accesions_filtered_window100kb_thincount50000_may2025.ld",header = TRUE)
Africa_breeding_LD<-read.table("/Users/kehan/Library/CloudStorage/OneDrive-UniversityofCalifornia,Davis/Monroe_Lab/Cassava_Center/Code/LD/LD_Decay/Africa_breeding_accesions_filtered_window100kb_thincount50000_may2025.ld",header = TRUE)
Wild_noHuKistler_LD<-read.table("/Users/kehan/Library/CloudStorage/OneDrive-UniversityofCalifornia,Davis/Monroe_Lab/Cassava_Center/Code/LD/LD_Decay/Wild_noHuKistler_accesions_filtered_window100kb_thincount50000_may2025.ld",header = TRUE)

Colombia_landraces_LD$Distance <- abs(as.numeric(Colombia_landraces_LD$BP_B) - as.numeric(Colombia_landraces_LD$BP_A))
Colombia_landraces_LD$R2<-as.numeric(Colombia_landraces_LD$R2)

South_America_breeding_noHuKistler_LD$Distance <- abs(as.numeric(South_America_breeding_noHuKistler_LD$BP_B) - as.numeric(South_America_breeding_noHuKistler_LD$BP_A))
South_America_breeding_noHuKistler_LD$R2<-as.numeric(South_America_breeding_noHuKistler_LD$R2)

Africa_breeding_LD$Distance <- abs(as.numeric(Africa_breeding_LD$BP_B) - as.numeric(Africa_breeding_LD$BP_A))
Africa_breeding_LD$R2<-as.numeric(Africa_breeding_LD$R2)

Wild_noHuKistler_LD$Distance <- abs(as.numeric(Wild_noHuKistler_LD$BP_B) - as.numeric(Wild_noHuKistler_LD$BP_A))
Wild_noHuKistler_LD$R2<-as.numeric(Wild_noHuKistler_LD$R2)

# plot raw dots without cut windows 
library(ggplot2)
library(data.table)
# Combined plot with smoothed curves
Colombia_landraces_LD$src<-"Colombian landraces"
Africa_breeding_LD$src<-"African breeding lines"
South_America_breeding_noHuKistler_LD$src<-"South American breeding lines (no Hu & Kistler)"
Wild_noHuKistler_LD$src<-"Wild_noHuKistler"

all_LD<-rbindlist(list(Colombia_landraces_LD, Africa_breeding_LD, South_America_breeding_noHuKistler_LD, Wild_noHuKistler_LD))
all_LD$src <- factor(all_LD$src, 
                     levels = c("Colombian landraces", 
                                "African breeding lines", 
                                "South American breeding lines (no Hu & Kistler)",
                                "Wild_noHuKistler"))

pdf("/Users/kehan/Library/CloudStorage/OneDrive-UniversityofCalifornia,Davis/Monroe_Lab/Cassava_Center/Figures/Figures_for_publication/Figure1c_LD_decay_100kb_thincount50000_may2025_v2.pdf", width=2.4, height=2)
ggplot(all_LD, aes(x = Distance, y = R2, col = src)) +
  geom_smooth(method = "auto") +
  labs(x = "Distance (bp)", y = "R²") +
  scale_color_manual(name = "Category", 
                     values = c("Colombian landraces" = "pink", 
                                "African breeding lines" = "#984ea3", 
                                "South American breeding lines (no Hu & Kistler)" = "#ffb000",
                                "Wild_noHuKistler" = "#319aa1"),
                     labels = c("Colombian landraces", "African breeding lines", "South American breeding lines", "Wild relatives")) +
  theme_classic(base_size = 6) +
  theme(axis.line = element_line(size = 0.5, colour = "black"),
        legend.key.size = unit(0.2, "cm"),
        legend.spacing.y = unit(0.2, "cm"),
        legend.position = c(0.7, 0.8),
        plot.background = element_blank(),
        panel.background = element_blank())
dev.off()

