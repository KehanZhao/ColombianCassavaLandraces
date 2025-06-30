library(corrplot)
library(dplyr)
library(data.table)
setwd("D:/WorkSpace/Rstudio_Workspace/Traits")
print("Reading in Environmental Traits and scaling")
Env <- read.table("Colombia_Batch1234_Phenos.txt",header=T)
Env_id <- Env[,1]
Env <- Env[,2:ncol(Env)]
Env <- scale(Env, center=TRUE, scale=TRUE) # center=TRUE, scale=TRUE are the defaults for scale()
## Recovering scaling coefficients
scale_env <- attr(Env, 'scaled:scale')
center_env <- attr(Env, 'scaled:center')
## Climatic table
Env <- as.data.frame(Env)
row.names(Env) <- c(Env_id)
PCs <- read.table("Colombia_Batch1234_PC.txt",header=T)

Variables <- data.frame(PCs, Env)
LOFPCs <- read.table("CassavaV8_GeneLoF_Disordered_tails_PCA.txt",header=T)
colnames(LOFPCs) <- c("Taxa","LOF_PC1","LOF_PC2","LOF_PC3","LOF_PC4","LOF_PC5")
Variables <- data.frame(LOFPCs, Variables)

Variables <- Variables %>% select(-c("Taxa.1","Taxa",Plant.height..cm.,Dry.matter.content....,Total.amylose,Total.sugars,Total.HCN.content))

data.frame(colnames(Variables))


library(corrplot)
df <- na.omit(Variables)
y_vars <- df[, 1:10]

# Select the remaining columns as the x-axis variables
x_vars <- df[, 11:ncol(df)]

# Calculate the correlation matrix
corr_matrix <- cor(cbind(y_vars, x_vars), use = "pairwise.complete.obs")

# Create the correlation plot
corrplot(corr_matrix, 
         method = "color", 
         type = "upper", 
         order = "hclust", 
         addCoef.col = "black", 
         tl.pos = "lt", 
         tl.col = "black", 
         tl.cex = 0.8, 
         title = "Correlation Plot")




subcor <- corr_matrix[c(6:8,1:3),11:ncol(corr_matrix)]


corrplot(as.matrix(subcor), is.corr=TRUE, 
         tl.pos = "lt", 
         tl.col = "black", 
         tl.cex = 0.8, 
         tl.srt = 65)  # Rotate x-axis titles by 45 degrees)

subcor_noLOF <- subcor[1:3,]

colnames(subcor_noLOF) <- gsub("bio2","Mean diurnal range",colnames(subcor_noLOF))
colnames(subcor_noLOF) <- gsub("bio3","Isothermality",colnames(subcor_noLOF))
colnames(subcor_noLOF) <- gsub("bio4","Temperature seasonality",colnames(subcor_noLOF))
colnames(subcor_noLOF) <- gsub("bio5","Max temperature of warmest month",colnames(subcor_noLOF))
colnames(subcor_noLOF) <- gsub("bio6","Min temperature of coldest month",colnames(subcor_noLOF))
colnames(subcor_noLOF) <- gsub("bio7","Temperature annual range",colnames(subcor_noLOF))
colnames(subcor_noLOF) <- gsub("bio8","Mean temperature of wettest quarter",colnames(subcor_noLOF))
colnames(subcor_noLOF) <- gsub("bio9","Mean temperature of driest quarter",colnames(subcor_noLOF))
colnames(subcor_noLOF) <- gsub("bio10","Mean temperature of warmest quarter",colnames(subcor_noLOF))
colnames(subcor_noLOF) <- gsub("bio11","Mean temperature of coldest quarter",colnames(subcor_noLOF))
colnames(subcor_noLOF) <- gsub("bio12","Annual precipitation",colnames(subcor_noLOF))
colnames(subcor_noLOF) <- gsub("bio13","Precipitation of wettest month",colnames(subcor_noLOF))
colnames(subcor_noLOF) <- gsub("bio14","Precipitation of driest month",colnames(subcor_noLOF))
colnames(subcor_noLOF) <- gsub("bio15","Precipitation seasonality",colnames(subcor_noLOF))
colnames(subcor_noLOF) <- gsub("bio16","Precipitation of wettest quarter",colnames(subcor_noLOF))
colnames(subcor_noLOF) <- gsub("bio17","Precipitation of driest quarter",colnames(subcor_noLOF))
colnames(subcor_noLOF) <- gsub("bio18","Precipitation of warmest quarter",colnames(subcor_noLOF))
colnames(subcor_noLOF) <- gsub("bio19","Precipitation of coldest quarter",colnames(subcor_noLOF))
colnames(subcor_noLOF) <- gsub("bio1","Annual mean temperature",colnames(subcor_noLOF))
colnames(subcor_noLOF) <- gsub("Elevation_AWS","Elevation",colnames(subcor_noLOF))
subcor_noLOF <- subcor_noLOF[,2:ncol(subcor_noLOF)]


corrplot(as.matrix(subcor_noLOF), is.corr=TRUE, 
         tl.pos = "lt", 
         tl.col = "black", 
         tl.cex = 0.8, 
         tl.srt = 65)  # Rotate x-axis titles by 45 degrees)





pdf("Corrplot_rectangular.pdf",width=6,height = 1.5)
corrplot(as.matrix(subcor_noLOF), is.corr=TRUE, 
         tl.pos = "lt", 
         tl.col = "black", 
         tl.cex = 0.3, 
         tl.srt = 65,
         cl.cex = 0.4)  # Rotate x-axis titles by 45 degrees)
dev.off()


