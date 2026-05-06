library(ggplot2)
library(dplyr)
library(tidyr)
library(ggsci)
library(readr)
library(data.table)

# Load accession data
accession_info <- read.csv("TableS1_Cassava_accessions_meta_information.csv")

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
K2 <- process_q_data("Cassava_final_1152accessions_apr2025_Qual30_maf5pct_xgboost_esculenta_excludingHu_geno10pct.2.Q", 2)
K3 <- process_q_data("Cassava_final_1152accessions_apr2025_Qual30_maf5pct_xgboost_esculenta_excludingHu_geno10pct.3.Q", 3)
K4 <- process_q_data("Cassava_final_1152accessions_apr2025_Qual30_maf5pct_xgboost_esculenta_excludingHu_geno10pct.4.Q", 4)
K5 <- process_q_data("Cassava_final_1152accessions_apr2025_Qual30_maf5pct_xgboost_esculenta_excludingHu_geno10pct.5.Q", 5)

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
PCA_result <- fread("Esculenta_accesions_withKistler_noHerbarium_filtered.eigenvec")
PCA_result<-as.data.frame(PCA_result)
renaming_map<-read.csv("Cassava_Master_Sheet_2_2024 - Sample renaming map.csv", header = F)

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

pdf("Admixture_Esculenta_noHu_orderedbyPC1_may2025.pdf", width=7, height=3)
K2plot + K3plot + K4plot + K5plot + plot_layout(ncol = 1) + theme(
  plot.background = element_blank(),
  panel.background = element_blank()
)
dev.off()

# CV of K choice

library(tidyverse)

txt <- "
CV error (K=1): 0.52154
CV error (K=10): 0.39407
CV error (K=11): 0.38778
CV error (K=12): 0.38404
CV error (K=13): 0.37985
CV error (K=14): 0.37750
CV error (K=15): 0.37392
CV error (K=16): 0.37201
CV error (K=17): 0.36927
CV error (K=19): 0.36373
CV error (K=2): 0.47995
CV error (K=20): 0.36271
CV error (K=21): 0.36088
CV error (K=3): 0.45766
CV error (K=4): 0.43793
CV error (K=5): 0.42567
CV error (K=6): 0.41574
CV error (K=7): 0.40754
CV error (K=8): 0.40174
CV error (K=9): 0.39582
"

# Convert to a data frame
df <- read_lines(txt) %>%
  tibble(line = .) %>%
  extract(
    line,
    into = c("K", "CV_error"),
    regex = "K=([0-9]+)\\): ([0-9.]+)"
  ) %>%
  mutate(
    K = as.numeric(K),
    CV_error = as.numeric(CV_error)
  ) %>%
  arrange(K)

df

# Plot
library(ggplot2)

pdf("FigureSx_ADMIXTURE_Cross_Validation_Error.pdf", width=3, height=2.5)
ggplot(df[1:15,], aes(x = K, y = CV_error)) +
  geom_line() +
  geom_point(size = 1) +
  theme_classic(base_size = 4) +
  labs(
    x = "K",
    y = "Cross-validation error"
  )+
  theme(
    plot.background = element_blank(),
    panel.background = element_blank()
  )
dev.off()


# K 6-10

# Process Q files for K6 to K10
K6 <- process_q_data("Cassava_final_1152accessions_apr2025_Qual30_maf5pct_xgboost_esculenta_excludingHu_geno10pct.6.Q", 6)
K7 <- process_q_data("Cassava_final_1152accessions_apr2025_Qual30_maf5pct_xgboost_esculenta_excludingHu_geno10pct.7.Q", 7)
K8 <- process_q_data("Cassava_final_1152accessions_apr2025_Qual30_maf5pct_xgboost_esculenta_excludingHu_geno10pct.8.Q", 8)
K9 <- process_q_data("Cassava_final_1152accessions_apr2025_Qual30_maf5pct_xgboost_esculenta_excludingHu_geno10pct.9.Q", 9)
K10 <- process_q_data("Cassava_final_1152accessions_apr2025_Qual30_maf5pct_xgboost_esculenta_excludingHu_geno10pct.10.Q", 10)

# Apply group label assignments
K6 <- assign_group_labels(K6)
K7 <- assign_group_labels(K7)
K8 <- assign_group_labels(K8)
K9 <- assign_group_labels(K9)
K10 <- assign_group_labels(K10)

K6$PC1<-PCA_result$V3[match(K6$sampleID,PCA_result$V1)]
K6$PC2<-PCA_result$V4[match(K6$sampleID,PCA_result$V1)]

K7$PC1<-PCA_result$V3[match(K7$sampleID,PCA_result$V1)]
K7$PC2<-PCA_result$V4[match(K7$sampleID,PCA_result$V1)]

K8$PC1<-PCA_result$V3[match(K8$sampleID,PCA_result$V1)]
K8$PC2<-PCA_result$V4[match(K8$sampleID,PCA_result$V1)]

K9$PC1<-PCA_result$V3[match(K9$sampleID,PCA_result$V1)]
K9$PC2<-PCA_result$V4[match(K9$sampleID,PCA_result$V1)]

K10$PC1<-PCA_result$V3[match(K10$sampleID,PCA_result$V1)]
K10$PC2<-PCA_result$V4[match(K10$sampleID,PCA_result$V1)]

# Plotting ADMIXTURE
library(ggplot2)
library(forcats)
library(ggthemes)
library(patchwork)

i <- 6

popGroup_colors6 <- c("1" = "#FF7F00",
                      "2" = "#E41A1C",
                      "3" = "#fb9a99",
                      "4" = "#377EB8",
                      "5" = "#984EA3",
                      "6" = "#4DAF4A")

K6 <- K6 %>%
  mutate(grp = factor(grp, levels = c("Colombian landraces",
                                      "South American breeding lines",
                                      "African breeding lines"))) %>%
  arrange(grp, PC1) %>%
  mutate(sampleID_factor = factor(sampleID, levels = unique(sampleID)))%>%
  mutate(popGroup = factor(popGroup, levels = c(5, 4, 2, 1, 3, 6)))

K6plot <-
  ggplot(K6, aes(x = sampleID_factor, y = prob, fill = factor(popGroup))) +
  geom_col(color = "gray", size = 0.1) +
  facet_grid(~fct_inorder(grp), scales = "free", switch = "x", space = "free") +
  theme_minimal(base_size = i) +
  labs(x = NULL, title = "K=6", y = NULL) +
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
  scale_fill_manual(values = popGroup_colors6, guide = FALSE)
K6plot

popGroup_colors7 <- c("1" = "#fdbf6f",
                      "2" = "#377EB8",
                      "3" = "#FF7F00",
                      "4" = "#E41A1C",
                      "5" = "#4DAF4A",
                      "6" = "#984EA3",
                      "7" = "#fb9a99")

K7 <- K7 %>%
  mutate(grp = factor(grp, levels = c("Colombian landraces",
                                      "South American breeding lines",
                                      "African breeding lines"))) %>%
  arrange(grp, PC1) %>%
  mutate(sampleID_factor = factor(sampleID, levels = unique(sampleID)))%>%
  mutate(popGroup = factor(popGroup, levels = c(1, 6, 2, 4, 3, 7, 5)))

K7plot <-
  ggplot(K7, aes(x = sampleID_factor, y = prob, fill = factor(popGroup))) +
  geom_col(color = "gray", size = 0.1) +
  facet_grid(~fct_inorder(grp), scales = "free", switch = "x", space = "free") +
  theme_minimal(base_size = i) +
  labs(x = NULL, title = "K=7", y = NULL) +
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
  scale_fill_manual(values = popGroup_colors7, guide = FALSE)
K7plot

popGroup_colors8 <- c("1" = "#fb9a99", # pink
                      "2" = "#FF7F00", # orange
                      "3" = "#377EB8", # blue
                      "4" = "#fdbf6f", # yellow
                      "5" = "#984EA3", # purple
                      "6" = "#cab2d6", # lavender
                      "7" = "#E41A1C", # red
                      "8" = "#4DAF4A") # green

K8 <- K8 %>%
  mutate(grp = factor(grp, levels = c("Colombian landraces",
                                      "South American breeding lines",
                                      "African breeding lines"))) %>%
  arrange(grp, PC1) %>%
  mutate(sampleID_factor = factor(sampleID, levels = unique(sampleID)))%>%
  mutate(popGroup = factor(popGroup, levels = c(4, 5, 3, 7, 6, 2, 1, 8)))

K8plot <-
  ggplot(K8, aes(x = sampleID_factor, y = prob, fill = factor(popGroup))) +
  geom_col(color = "gray", size = 0.1) +
  facet_grid(~fct_inorder(grp), scales = "free", switch = "x", space = "free") +
  theme_minimal(base_size = i) +
  labs(x = NULL, title = "K=8", y = NULL) +
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
  scale_fill_manual(values = popGroup_colors8, guide = FALSE)
K8plot

popGroup_colors9 <- c("1" = "#a6cee3",  # light blue
                      "2" = "#377EB8",  # blue
                      "3" = "#984EA3",  # purple
                      "4" = "#fb9a99",  # pink
                      "5" = "#E41A1C",  # red
                      "6" = "#fdbf6f",  # yellow
                      "7" = "#4DAF4A",  # green
                      "8" = "#FF7F00",  # orange
                      "9" = "#cab2d6")  # lavender

K9 <- K9 %>%
  mutate(grp = factor(grp, levels = c("Colombian landraces",
                                      "South American breeding lines",
                                      "African breeding lines"))) %>%
  arrange(grp, PC1) %>%
  mutate(sampleID_factor = factor(sampleID, levels = unique(sampleID)))%>%
  mutate(popGroup = factor(popGroup, levels = c(6, 3, 2, 5, 9, 8, 4, 1, 7)))

K9plot <-
  ggplot(K9, aes(x = sampleID_factor, y = prob, fill = factor(popGroup))) +
  geom_col(color = "gray", size = 0.1) +
  facet_grid(~fct_inorder(grp), scales = "free", switch = "x", space = "free") +
  theme_minimal(base_size = i) +
  labs(x = NULL, title = "K=9", y = NULL) +
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
  scale_fill_manual(values = popGroup_colors9, guide = FALSE)
K9plot

popGroup_colors10 <- c("1" = "#984EA3", # purple
                       "2" = "#377EB8", # blue
                       "3" = "#fdbf6f", # yellow
                       "4" = "#4DAF4A", # green
                       "5" = "#cab2d6", # lavender
                       "6" = "#fb9a99", # pink
                       "7" = "#a6cee3",
                       "8" = "#b2df8a",
                       "9" = "#E41A1C", # red
                       "10" = "#FF7F00")

K10 <- K10 %>%
  mutate(grp = factor(grp, levels = c("Colombian landraces",
                                      "South American breeding lines",
                                      "African breeding lines"))) %>%
  arrange(grp, PC1) %>%
  mutate(sampleID_factor = factor(sampleID, levels = unique(sampleID)))%>%
  mutate(popGroup = factor(popGroup, levels = c(8, 3, 1, 2, 9, 5, 10, 6, 7, 4)))

K10plot <-
  ggplot(K10, aes(x = sampleID_factor, y = prob, fill = factor(popGroup))) +
  geom_col(color = "gray", size = 0.1) +
  facet_grid(~fct_inorder(grp), scales = "free", switch = "x", space = "free") +
  theme_minimal(base_size = i) +
  labs(x = NULL, title = "K=10", y = NULL) +
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
  scale_fill_manual(values = popGroup_colors10, guide = FALSE)
K10plot

pdf("FigureSx_Admixture_Esculenta_noHu_orderedbyPC1_K6_10.pdf", width=12, height=8)
K6plot + K7plot + K8plot + K9plot + K10plot+ plot_layout(ncol = 1) + theme(
  plot.background = element_blank(),
  panel.background = element_blank()
)
dev.off()
