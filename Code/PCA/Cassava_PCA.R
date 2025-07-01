# Load meta information
accession_info<-read.csv('/Users/kehan/Library/CloudStorage/OneDrive-UniversityofCalifornia,Davis/Monroe_Lab/Cassava_Center/Tables/TableS1_Cassava_accessions_meta_information.csv')

# Run PCA using Plink on HPC
# Load results
library(data.table)
PCA_result <- fread('/Users/kehan/Library/CloudStorage/OneDrive-UniversityofCalifornia,Davis/Monroe_Lab/Cassava_Center/Code/PCA/Data/All_accessions_withKistler_noHerbarium_filtered_Qual30_maf5pct_xgboost.eigenvec')
PCA_result<-as.data.frame(PCA_result)

renaming_map<-read.csv('/Users/kehan/Downloads/Cassava_Master_Sheet_2_2024 - Sample renaming map.csv', header = F)
library(dplyr)
head(PCA_result)
# Join to get new names
PCA_result <- PCA_result %>%
  left_join(renaming_map, by = c("V1" = "V1")) %>%
  mutate(V1 = ifelse(!is.na(V2.y), V2.y, V1)) %>%
  select(-V2.y) %>%
  rename(V2 = V1)
colnames(PCA_result)[1:2] <- c("V1", "V2")
head(PCA_result)

eigenValues <-read.delim('/Users/kehan/Library/CloudStorage/OneDrive-UniversityofCalifornia,Davis/Monroe_Lab/Cassava_Center/Code/PCA/Data/All_accessions_withKistler_noHerbarium_filtered_Qual30_maf5pct_xgboost.eigenval',sep = " ",header = F)

PCA_result$Country<-accession_info$Country[match(PCA_result$V1,accession_info$ID)]
PCA_result$Continent<-accession_info$Continent[match(PCA_result$V1,accession_info$ID)]
PCA_result$Source<-accession_info$Source[match(PCA_result$V1,accession_info$ID)]
PCA_result$Status<-accession_info$Biological.status[match(PCA_result$V1,accession_info$ID)]
PCA_result$Species<-accession_info$Species[match(PCA_result$V1,accession_info$ID)]
PCA_result$Collection<-accession_info$Collection.status[match(PCA_result$V1,accession_info$ID)]

eigen_percent<-round((eigenValues/(sum(eigenValues))*100),2)

# Plot PCA
library(ggplot2)
# Set color
species_colors <- c("alutacea" = "#33a02c", 
                    "caerulescens" = "#fb9a99",
                    "carthaginensis" = "#b15928",
                    "carthaginensis subsp. glaziovii" = "#ff7f00", 
                    "esculenta" = "grey", 
                    "esculenta x glaziovii" = "#ffff99",
                    "esculenta subsp. flabellifolia" = "#b2df8a",
                    "esculenta x flabellifolia" = "green",
                    "glaziovii" = "blue",
                    "longepetiolata" = "#e31a1c",
                    "orbicularis" = "#6a3d9a",
                    "peruviana" = "#cab2d6",
                    "pseudoglaziovii" = "#fdbf6f",
                    "triphylla" = "#1f78b4",
                    "tristis" = "#a6cee3")

pdf("/Users/kehan/Library/CloudStorage/OneDrive-UniversityofCalifornia,Davis/Monroe_Lab/Cassava_Center/Figures/Figures_for_publication/Figure1a_PCA_All_PC1&2_recolored_withKistler_noHerbarium.pdf", width=2.3, height=2)
ggplot(data = PCA_result) +
  # Base points for all samples
  geom_point(mapping = aes(x = V3, y = V4, colour = Species), size = 0.3, show.legend = TRUE) +
  
  labs(x = paste0("PC 1 (", eigen_percent[1, 1], " %)"),
       y = paste0("PC 2 (", eigen_percent[2, 1], " %)")) +
  
  scale_colour_manual(values = species_colors) +

  theme_classic(base_size = 6) +
  theme(axis.line = element_line(size = 0.5, colour = "black"),
        legend.key.size = unit(0.2, "cm"),
        legend.spacing.y = unit(0.2, "cm"),
        legend.position = c(0.73, 0.68),
        plot.background = element_blank(),
        panel.background = element_blank())
dev.off()

# Plot only esculenta
library(data.table)
PCA_result <- fread('/Users/kehan/Library/CloudStorage/OneDrive-UniversityofCalifornia,Davis/Monroe_Lab/Cassava_Center/Code/PCA/Data/Esculenta_accesions_withKistler_noHerbarium_filtered.eigenvec')
PCA_result<-as.data.frame(PCA_result)

# Join to get new names
PCA_result <- PCA_result %>%
  left_join(renaming_map, by = c("V1" = "V1")) %>%
  mutate(V1 = ifelse(!is.na(V2.y), V2.y, V1)) %>%
  select(-V2.y) %>%
  rename(V2 = V1)
colnames(PCA_result)[1:2] <- c("V1", "V2")
head(PCA_result)

eigenValues <-read.delim('/Users/kehan/Library/CloudStorage/OneDrive-UniversityofCalifornia,Davis/Monroe_Lab/Cassava_Center/Code/PCA/Data/Esculenta_accesions_withKistler_noHerbarium_filtered.eigenval',sep = " ",header = F)

PCA_result$Country<-accession_info$Country[match(PCA_result$V1,accession_info$ID)]
PCA_result$Continent<-accession_info$Continent[match(PCA_result$V1,accession_info$ID)]
PCA_result$Source<-accession_info$Source[match(PCA_result$V1,accession_info$ID)]
PCA_result$Status<-accession_info$Biological.status[match(PCA_result$V1,accession_info$ID)]
PCA_result$Species<-accession_info$Species[match(PCA_result$V1,accession_info$ID)]
PCA_result$Collection<-accession_info$Collection.status[match(PCA_result$V1,accession_info$ID)]

eigen_percent<-round((eigenValues/(sum(eigenValues))*100),2)

library(ggplot2)
Continent_colors <- c("Africa" = "#984ea3", 
                      "Asia" = "#bd0026",
                      "Europe" = "#4daf4a",
                      "Oceania" = "#377eb8", 
                      "South America" = "#ffb000")

pdf("/Users/kehan/Library/CloudStorage/OneDrive-UniversityofCalifornia,Davis/Monroe_Lab/Cassava_Center/Figures/Figures_for_publication/Figure1b_PCA_Esculenta_PC1&2_recolored_withKistler_noHerbarium_v2.pdf", width=2.3, height=2)
ggplot(data = PCA_result) +
  geom_point(mapping = aes(x = V3, y = V4, colour = Continent), size = 0.3, show.legend = TRUE) +
  
  # Highlight points from "This study" with black circle outline
  geom_point(data = PCA_result[PCA_result$Source == "This study", ],
             mapping = aes(x = V3, y = V4),
             shape = 21, fill = NA, colour = "black", size = 0.9, stroke = 0.1) +
  
  labs(x = paste0("PC 1 (", eigen_percent[1, 1], " %)"),
       y = paste0("PC 2 (", eigen_percent[2, 1], " %)")) +
  
  scale_x_continuous(expand = expansion(mult = c(0.02, 0.02))) +  # Add space on both sides of x-axis
  scale_y_continuous(expand = expansion(mult = c(0.02, 0.02))) +  # Optional: do the same for y-axis
  
  scale_colour_manual(values = Continent_colors) +
  
  theme_classic(base_size = 6) +
  theme(
    axis.line = element_line(size = 0.5, colour = "black"),
    legend.key.size = unit(0.2, "cm"),
    legend.spacing.y = unit(0.2, "cm"),
    legend.position = c(0.18, 0.18),
    plot.background = element_blank(),
    panel.background = element_blank()
  )
dev.off()


