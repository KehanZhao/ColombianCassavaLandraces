library(ggplot2)
library(dplyr)
library(tidyr)
library(ggsci)
library(readr)
library(data.table)

# Load accession data
accession_info <- read.csv('/Users/kehan/Library/CloudStorage/OneDrive-UniversityofCalifornia,Davis/Monroe_Lab/Cassava_Center/Tables/S1_Cassava_accessions_meta_information.csv')

# Filter for 'esculenta' species and remove accessions from Hu et al. 2021
Esculenta_noHu<-accession_info$ID[which(accession_info$Species=="esculenta" & accession_info$Source!="Hu et al. 2021")]

# Function to process Q data and add location information
process_q_data <- function(file_path, num_clusters) {
  q_data <- read.table(file_path, header = FALSE)
  q_data <- cbind(Esculenta_noHu, q_data)
  colnames(q_data) <- c("sampleID", as.character(1:num_clusters))
  
  # Reshape data from wide to long format
  K <- gather(q_data, key = "popGroup", value = "prob", -sampleID)
  
  # Add continent information based on accession_info
  K$loc <- accession_info$Continent[match(K$sampleID, accession_info$ID)]
  return(K)
}

# Process Q files for K2 to K5
K2 <- process_q_data('/Users/kehan/Library/CloudStorage/OneDrive-UniversityofCalifornia,Davis/Monroe_Lab/Cassava_Center/Code/Admixture/Cassava_final_1152accessions_apr2025_Qual30_maf5pct_xgboost_esculenta_excludingHu_geno10pct.2.Q', 2)
K3 <- process_q_data('/Users/kehan/Library/CloudStorage/OneDrive-UniversityofCalifornia,Davis/Monroe_Lab/Cassava_Center/Code/Admixture/Cassava_final_1152accessions_apr2025_Qual30_maf5pct_xgboost_esculenta_excludingHu_geno10pct.3.Q', 3)
K4 <- process_q_data('/Users/kehan/Library/CloudStorage/OneDrive-UniversityofCalifornia,Davis/Monroe_Lab/Cassava_Center/Code/Admixture/Cassava_final_1152accessions_apr2025_Qual30_maf5pct_xgboost_esculenta_excludingHu_geno10pct.4.Q', 4)
K5 <- process_q_data('/Users/kehan/Library/CloudStorage/OneDrive-UniversityofCalifornia,Davis/Monroe_Lab/Cassava_Center/Code/Admixture/Cassava_final_1152accessions_apr2025_Qual30_maf5pct_xgboost_esculenta_excludingHu_geno10pct.5.Q', 5)

# Load classification data for different accession groups
South_america_breeding_noHu<-accession_info$ID[which(accession_info$Continent=="South America" & (accession_info$Biological.status=="Breeding line" | accession_info$Biological.status=="Domesticated") & accession_info$Species=="esculenta" & accession_info$Source!="This study" & accession_info$Source!="Hu et al. 2021")]
Colombia_landraces<-accession_info$ID[which(accession_info$Source=="This study" & accession_info$Species=="esculenta" & accession_info$Biological.status=="Landrace")]
Africa_breeding<-accession_info$ID[which(accession_info$Continent=="Africa" & (accession_info$Biological.status=="Breeding line" | accession_info$Biological.status=="Domesticated") & accession_info$Species=="esculenta" & accession_info$Source!="This study")]

# Function to assign group labels based on classification lists
assign_group_labels <- function(K) {
  K$grp <- NA
  K$grp[K$sampleID %in% Africa_breeding] <- "African breeding lines"
  K$grp[K$sampleID %in% South_america_breeding_noHu] <- "South American breeding lines"
  K$grp[K$sampleID %in% Colombia_landraces] <- "Colombian landraces"
  K <- K[!is.na(K$grp), ]  # Remove rows with NA group labels
  return(K)
}

# Apply group label assignments
K2 <- assign_group_labels(K2)
K3 <- assign_group_labels(K3)
K4 <- assign_group_labels(K4)
K5 <- assign_group_labels(K5)

# Reorder by PC
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

K2$PC1<-PCA_result$V3[match(K2$sampleID,PCA_result$V1)]
K2$PC2<-PCA_result$V4[match(K2$sampleID,PCA_result$V1)]

K3$PC1<-PCA_result$V3[match(K3$sampleID,PCA_result$V1)]
K3$PC2<-PCA_result$V4[match(K3$sampleID,PCA_result$V1)]

K4$PC1<-PCA_result$V3[match(K4$sampleID,PCA_result$V1)]
K4$PC2<-PCA_result$V4[match(K4$sampleID,PCA_result$V1)]

K5$PC1<-PCA_result$V3[match(K5$sampleID,PCA_result$V1)]
K5$PC2<-PCA_result$V4[match(K5$sampleID,PCA_result$V1)]

# Plotting ADMIXTURE
library(ggplot2)
library(forcats)
library(ggthemes)
library(patchwork)

# Define custom colors for popGroups
popGroup_colors2 <- c("1" = "#377EB8", 
                      "2" = "#4DAF4A") 

K2 <- K2 %>%
  mutate(grp = factor(grp, levels = c("Colombian landraces", 
                                      "South American breeding lines", 
                                      "African breeding lines"))) %>%
  arrange(grp, PC1) %>%
  mutate(sampleID_factor = factor(sampleID, levels = unique(sampleID)))

K2plot <- 
  ggplot(K2, aes(x = sampleID_factor, y = prob, fill = factor(popGroup))) +
  geom_col(color = "gray", size = 0.1) +
  facet_grid(~fct_inorder(grp), scales = "free", switch = "x", space = "free") +
  theme_minimal(base_size = 6) + 
  labs(x = NULL, title = "K=2", y = "Ancestry") +
  scale_y_continuous(expand = c(0, 0)) +
  scale_x_discrete(expand = expansion(add = 1)) +
  theme(
    panel.spacing.x = unit(0.1, "lines"),
    axis.text.x = element_blank(),
    panel.grid = element_blank(),
    plot.background = element_blank(),
    panel.background = element_blank(),
    strip.text = element_blank()
  ) +
  scale_fill_manual(values = popGroup_colors2, guide = FALSE)

popGroup_colors3 <- c("1" = "#984EA3", 
                      "2" = "#377EB8", 
                      "3" = "#4DAF4A")

K3 <- K3 %>%
  mutate(grp = factor(grp, levels = c("Colombian landraces", 
                                      "South American breeding lines", 
                                      "African breeding lines"))) %>%
  arrange(grp, PC1) %>%
  mutate(sampleID_factor = factor(sampleID, levels = unique(sampleID)))

K3plot <- 
  ggplot(K3, aes(x = sampleID_factor, y = prob, fill = factor(popGroup))) +
  geom_col(color = "gray", size = 0.1) +
  facet_grid(~fct_inorder(grp), scales = "free", switch = "x", space = "free") +
  theme_minimal(base_size = 6) + 
  labs(x = NULL, title = "K=3", y = "Ancestry") +
  scale_y_continuous(expand = c(0, 0)) +
  scale_x_discrete(expand = expansion(add = 1)) +
  theme(
    panel.spacing.x = unit(0.1, "lines"),
    axis.text.x = element_blank(),
    panel.grid = element_blank(),
    plot.background = element_blank(),
    panel.background = element_blank(),
    strip.text = element_blank()
  ) +
  scale_fill_manual(values = popGroup_colors3, guide = FALSE)

popGroup_colors4 <- c("1" = "#4DAF4A", 
                      "2" = "#E41A1C", 
                      "3" = "#984EA3",  
                      "4" = "#377EB8") 

K4 <- K4 %>%
  mutate(grp = factor(grp, levels = c("Colombian landraces", 
                                      "South American breeding lines", 
                                      "African breeding lines"))) %>%
  arrange(grp, PC1) %>%
  mutate(sampleID_factor = factor(sampleID, levels = unique(sampleID))) %>%
  mutate(popGroup = factor(popGroup, levels = c(4, 3, 2, 1)))

K4plot <- 
  ggplot(K4, aes(x = sampleID_factor, y = prob, fill = factor(popGroup))) +
  geom_col(color = "gray", size = 0.1) +
  facet_grid(~fct_inorder(grp), scales = "free", switch = "x", space = "free") +
  theme_minimal(base_size = 6) + 
  labs(x = NULL, title = "K=4", y = "Ancestry") +
  scale_y_continuous(expand = c(0, 0)) +
  scale_x_discrete(expand = expansion(add = 1)) +
  theme(
    panel.spacing.x = unit(0.1, "lines"),
    axis.text.x = element_blank(),
    panel.grid = element_blank(),
    plot.background = element_blank(),
    panel.background = element_blank(),
    strip.text = element_blank()
  ) +
  scale_fill_manual(values = popGroup_colors4, guide = FALSE)

popGroup_colors5 <- c("1" = "#984EA3", 
                      "2" = "#377EB8", 
                      "3" = "#E41A1C",  
                      "4" = "#FF7F00",  
                      "5" = "#4DAF4A") 

K5 <- K5 %>%
  mutate(grp = factor(grp, levels = c("Colombian landraces", 
                                      "South American breeding lines", 
                                      "African breeding lines"))) %>%
  arrange(grp, PC1) %>%
  mutate(sampleID_factor = factor(sampleID, levels = unique(sampleID)))

K5plot <- 
  ggplot(K5, aes(x = sampleID_factor, y = prob, fill = factor(popGroup))) +
  geom_col(color = "gray", size = 0.1) +
  facet_grid(~fct_inorder(grp), scales = "free", switch = "x", space = "free") +
  theme_minimal(base_size = 6) + 
  labs(x = "Individuals", title = "K=5", y = "Ancestry") +
  scale_y_continuous(expand = c(0, 0)) +
  scale_x_discrete(expand = expansion(add = 1)) +
  theme(
    panel.spacing.x = unit(0.1, "lines"),
    axis.text.x = element_blank(),
    panel.grid = element_blank(),
    plot.background = element_blank(),
    panel.background = element_blank()
  ) +
  scale_fill_manual(values = popGroup_colors5, guide = FALSE)

pdf("/Users/kehan/Library/CloudStorage/OneDrive-UniversityofCalifornia,Davis/Monroe_Lab/Cassava_Center/Figures/Figures_for_publication/Admixture_Esculenta_noHu_orderedbyPC1_may2025.pdf", width=7, height=3)
K2plot + K3plot + K4plot + K5plot + plot_layout(ncol = 1) + theme(
  plot.background = element_blank(),
  panel.background = element_blank()
)
dev.off()


