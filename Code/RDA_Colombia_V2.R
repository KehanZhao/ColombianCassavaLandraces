library(pegas)
library(ggplot2)
library(raster)
library(rgdal)
library(LEA)
library(rnaturalearth)
library(rnaturalearthdata)
library(RColorBrewer)
library(ggpubr)
library(vegan)
library(qvalue)
library(robust)
library(WMDB)
#library(ggVennDiagram)
library(cowplot)
library(corrplot)
library(rgeos)
library(qqman)
library(corrplot)
library(rgeos)
library(ggrepel)
library(dplyr)
library(data.table)
args <- commandArgs(trailingOnly = T)
Genotypes <- fread(args[2])
Geno_id <- Genotypes[,1]
Genotypes <- Genotypes[,2:ncol(Genotypes)]
output_folder <- args[1]
rownames(Genotypes) <- as.data.frame(Geno_id)[,1]
Genotype_Positions <- read.table(paste0(args[2],".pos"))
Genotype_Positions <- Genotype_Positions[Genotype_Positions[,1] %in% 1:18,]  
colnames(Genotype_Positions) <- c("Chr","Pos","SNPID")
Genotypes <- Genotypes[,1:nrow(Genotype_Positions)]
#print(nrow(Genotype_Positions))
#print(ncol(Genotypes))
#Genotypes
rownames(Genotypes) <- as.data.frame(Geno_id)[,1]



#if(args[4]=="T"){
#	print("Performing Median imputation")
#	for (i in 1:ncol(Genotypes))
#		{
# 		Genotypes[which(is.na(Genotypes[,i])),i] <- median(Genotypes[-which(is.na(Genotypes[,i])),i], na.rm=TRUE)
#		}
#}
## Filtering on MAF
#print("Filtering to 5% MAF")
#freq_mean <- colMeans(Genotypes)
#if(max(freq_mean)>1)
#	{
#	freq_mean <- freq_mean/2
#	}	
#Genotypes <- Genotypes[,-which(freq_mean>=0.95 | freq_mean<=0.05)]
#Genotype_Positions <- Genotype_Positions[-which(freq_mean>=0.95 | freq_mean<=0.05),]
## Ordering loci based on their scaffold
#Genotypes <- Genotypes[,order(colnames(Genotypes))]



print("Reading in Environmental Traits and scaling")
Env <- read.table("Cassava_Data/Colombia_Batch1234_Phenos.txt",header=T)
Env_id <- Env[,1]
Env <- Env[,2:ncol(Env)]
Env <- scale(Env, center=TRUE, scale=TRUE) # center=TRUE, scale=TRUE are the defaults for scale()
## Recovering scaling coefficients
scale_env <- attr(Env, 'scaled:scale')
center_env <- attr(Env, 'scaled:center')
## Climatic table
Env <- as.data.frame(Env)
row.names(Env) <- c(Env_id)
PCs <- read.table("Cassava_Data/Colombia_Batch1234_PC.txt",header=T)

Variables <- data.frame(PCs, Env)
#print(Variables[])
## Null model
#RDA0 <- rda(Genotypes ~ 1, Variables)


Variables <- Variables %>%  arrange(as.character(rownames(.)))
Genotypes <- Genotypes %>%  arrange(as.character(rownames(.)))

#print(cbind.data.frame(rownames(Variables),rownames(Genotypes)))

## Full model
Genotypes <- Genotypes[!is.na(Variables$bio1),]
Variables <- Variables %>% filter(!is.na(bio1))
#RDAfull <- rda(Genotypes ~ bio1 + bio2 + bio3 + bio4 + bio5 + bio6 + bio7 + bio8 + bio9 + bio10 + bio11 + bio12 + bio13 + bio14 + bio15 + bio16 + bio17 + bio18 + bio19 + Altitude2 + EnvPC1 + EnvPC2 + EnvPC3 + EnvPC4 + EnvPC5, Variables)
system(paste0("mkdir -p ",output_folder))
setwd(output_folder)

#mod <- ordiR2step(RDA0, RDAfull, Pin = 0.01, R2permutations = 1000, R2scope = T)
#print(mod)
print("Now that we have selected variables, let's run some RDA")
print("First, piecewise RDA estimate variance capture by pop structure, geography, etc")
print(colnames(Variables))
pRDAfull <- rda(Genotypes ~ PC1 + PC2 + PC3 + PC4 + PC5 + Longitude + Latitude + EnvPC1 + bio1 + bio4 + bio5 + bio6 + bio10 + bio14 + bio15 + bio18,  Variables)
pRDAfull_rq <- RsquareAdj(pRDAfull)
pRDAfull_anova <- anova(pRDAfull)
print(pRDAfull_rq)

pRDAClimGeo <- rda(Genotypes ~  Longitude + Latitude + EnvPC1 + bio1 + bio4 + bio5 + bio6 + bio10 + bio14 + bio15 + bio18 + Condition(PC1 + PC2 + PC3 + PC4 + PC5),  Variables)
pRDAClimGeo_rq <- RsquareAdj(pRDAClimGeo)
pRDAClimGeo_anova <- anova(pRDAClimGeo)

pRDAclim <- rda(Genotypes ~ EnvPC1 + bio1 + bio4 + bio5 + bio6 + bio10 + bio14 + bio15 + bio18 + Condition(Longitude + Latitude + PC1 + PC2 + PC3 + PC4 + PC5),  Variables)
pRDAclim_rq <- RsquareAdj(pRDAclim)
pRDAclim_anova <- anova(pRDAclim)


## Pure neutral population structure model
pRDAstruct <- rda(Genotypes ~ PC1 + PC2 + PC3 + PC4 + PC5 + Condition(Longitude + Latitude + EnvPC1 + bio1 + bio4 + bio5 + bio6 + bio10 + bio14 + bio15 + bio18),  Variables)
pRDAstruct_rq <- RsquareAdj(pRDAstruct)
pRDAstruct_anova <- anova(pRDAstruct)

pRDAstructandGeo <- rda(Genotypes ~ PC1 + PC2 + PC3 + PC4 + PC5 + Longitude + Latitude + Condition(EnvPC1 + bio1 + bio4 + bio5 + bio6 + bio10 + bio14 + bio15 + bio18),  Variables)
pRDAstructandGeo_rq <- RsquareAdj(pRDAstructandGeo)
pRDAstructadnGeo_anova <- anova(pRDAstructandGeo)

pRDAstructandClim <- rda(Genotypes ~ PC1 + PC2 + PC3 + PC4 + PC5 + EnvPC1 + bio1 + bio4 + bio5 + bio6 + bio10 + bio14 + bio15 + bio18 + Condition(Longitude + Latitude),  Variables)
pRDAstructandClim_rq <- RsquareAdj(pRDAstructandClim)
pRDAstructandClim_anova <- anova(pRDAstructandClim)

##Pure geography model 
pRDAgeog <- rda(Genotypes ~ Longitude + Latitude + Condition(EnvPC1 + bio1 + bio4 + bio5 + bio6 + bio10 + bio14 + bio15 + bio18 + PC1 + PC2 + PC3 + PC4 + PC5),  Variables)
pRDAgeog_rq <- RsquareAdj(pRDAgeog)
pRDAgeog_anova <- anova(pRDAgeog)

#print(Variables[,c("EnvPC1","bio1","bio4","bio5","bio6","bio10","bio14","bio15","bio18","PC1","PC2","PC3","Longitude","Latitude", "Altitude2")])


#Get Some raw RDAs
RDA_raw_env <-  rda(Genotypes ~ EnvPC1 + EnvPC2 + EnvPC3 + EnvPC4 + EnvPC5,  Variables)
saveRDS(RDA_raw_env,"RDA_9V_raw_Env.rds")
RDA_env_condGeo <-  rda(Genotypes ~ EnvPC1 + EnvPC2 + EnvPC3 + EnvPC4 + EnvPC5 + Condition(Longitude + Latitude),  Variables)
saveRDS(RDA_env_condGeo,"RDA_9V_env_condGeo.rds")



##Write out partial models results
Partial_RDAtable <- cbind.data.frame(
Model=c("Full","Structure_And_Geography","Structure_And_Climate","Climate_And_Geography","Climate","Structure","Geography"),
Variance=c(pRDAfull_anova$Variance[1],
	   pRDAstructadnGeo_anova$Variance[1],
	   pRDAstructandClim_anova$Variance[1],
	   pRDAClimGeo_anova$Variance[1],
	   pRDAclim_anova$Variance[1],
	   pRDAstruct_anova$Variance[1],
	   pRDAgeog_anova$Variance[1]),
Rsquared=c(pRDAfull_rq$r.squared,
	   pRDAstructandGeo_rq$r,
	   pRDAstructandClim_rq$r,
	   pRDAClimGeo_rq$r,
	   pRDAclim_rq$r.squared,
	   pRDAstruct_rq$r.squared,
	   pRDAgeog_rq$r.squared),
Adjsquared=c(pRDAfull_rq$adj.r.squared,
	     pRDAstructandGeo_rq$adj.r.squared,
	     pRDAstructandClim_rq$adj.r.squared,
	     pRDAClimGeo_rq$adj.r.squared,
	     pRDAclim_rq$adj.r.squared,
	     pRDAstruct_rq$adj.r.squared,
	     pRDAgeog_rq$adj.r.squared),
pvalue=c(pRDAfull_anova[[4]][1],
	 pRDAstructadnGeo_anova[[4]][1],
	 pRDAstructandClim_anova[[4]][1],
	 pRDAClimGeo_anova[[4]][1],
	 pRDAclim_anova[[4]][1],
	 pRDAstruct_anova[[4]][1],
	 pRDAgeog_anova[[4]][1]))
write.table(Partial_RDAtable,"Partial_RDAtable_9V.txt",col.names=T,row.names=F,sep="\t",quote=F)


## Make correlation Plots
  png(file="Corrplot_selected.png",width=800, height=800)
  corrplot(cor(Variables[,c("EnvPC1","bio1","bio4","bio5","bio6","bio10","bio14","bio15","bio18","PC1","PC2","PC3","Longitude","Latitude", "Elevation_AWS")]), type="upper")
  dev.off()
  png(file="Corrplot_full.png",width=800, height=800)
  corrplot(cor(Variables[,2:ncol(Variables)]), type="upper")
  dev.off()
RDA_env <- rda(Genotypes ~ EnvPC1 + bio1 + bio4 + bio5 + bio6 + bio10 + bio14 + bio15 + bio18 + Condition(PC1 + PC2 + PC3 + PC4 + PC5 + Longitude + Latitude),  Variables)
saveRDS(RDA_env,"RDA_9V.rds")
png(file="Eigenvalues.png",width=800, height=800)
screeplot(RDA_env, main="Eigenvalues of constrained axes")
dev.off()

print("let's test different numbers of RDA dimensions to include")
source("/group/gmonroegrp2/evlong/Cassava/RDA/RDA-landscape-genomics/src/rdadapt.R")
for(RDAK in 3:5){
print(RDAK)
###First 5 RDA axes look good
## Function rdadapt

## Running the function with K = RDAK
rdadapt_env<-rdadapt(RDA_env, RDAK)
write.table(rdadapt_env,paste0("RDA_outputk",RDAK,".txt"),col.names=T,row.names=T,quote=F,sep="\t")

## P-values threshold after Bonferroni correction
thres_env <- 0.01/length(rdadapt_env$p.values)
## Identifying the loci that are below the p-value threshold
outliers <- data.frame(Loci = colnames(Genotypes)[which(rdadapt_env$p.values<thres_env)], p.value = rdadapt_env$p.values[which(rdadapt_env$p.values<thres_env)], contig = unlist(lapply(strsplit(colnames(Genotypes)[which(rdadapt_env$p.values<thres_env)], split = "_"), function(x) x[1])))

## Top hit outlier per contig
print(outliers)
if(nrow(outliers >=1)){
outliers <- outliers[order(outliers$contig, outliers$p.value),]
## List of outlier names
outliers_rdadapt_env <- as.character(outliers$Loci[!duplicated(outliers$contig)])

## Formatting table for ggplot
locus_scores <- scores(RDA_env, choices=c(1:RDAK), display="species", scaling="none") # vegan references "species", here these are the loci
TAB_loci <- data.frame(names = row.names(locus_scores), locus_scores)
TAB_loci$type <- "Neutral"
TAB_loci$type[TAB_loci$names%in%outliers$Loci] <- "All outliers"
TAB_loci$type[TAB_loci$names%in%outliers_rdadapt_env] <- "Top outliers"
TAB_loci$type <- factor(TAB_loci$type, levels = c("Neutral", "All outliers", "Top outliers"))
TAB_loci <- cbind.data.frame(Genotype_Positions,TAB_loci,rdadapt_env$p.values)
colnames(TAB_loci)[ncol(TAB_loci)] <- "p.value"
write.table(TAB_loci,paste0("TAB_loci_k",RDAK,"_9V.txt"),col.names=T,row.names=F,quote=F,sep="\t")
TAB_loci <- TAB_loci[order(TAB_loci$type),]
TAB_var <- as.data.frame(scores(RDA_env, choices=c(1:RDAK), display="bp")) # pull the biplot scores
write.table(TAB_var,"TAB_var_9V.txt",col.names=T,row.names=F,quote=F,sep="\t")
## Biplot of RDA loci and variables scores
print(colnames(TAB_loci))
print(colnames(TAB_var))
print("Hey let's plot!")
p1 <- ggplot() +
  geom_hline(yintercept=0, linetype="dashed", color = gray(.80), size=0.6) +
  geom_vline(xintercept=0, linetype="dashed", color = gray(.80), size=0.6) +
  geom_point(data = TAB_loci, aes(x=RDA1*20, y=RDA2*20, colour = type), size = 1.4) +
  scale_color_manual(values = c("gray90", "#F9A242FF", "#6B4596FF")) +
  geom_segment(data = TAB_var, aes(xend=RDA1, yend=RDA2, x=0, y=0), colour="black", size=0.15, linetype=1, arrow=arrow(length = unit(0.02, "npc"))) +
  geom_text(data = TAB_var, aes(x=1.1*RDA1, y=1.1*RDA2, label = row.names(TAB_var)), size = 2.5, family = "Times") +
  xlab("RDA 1") + ylab("RDA 2") +
  facet_wrap(~"RDA space") +
  guides(color=guide_legend(title="Locus type")) +
  theme_bw(base_size = 11, base_family = "Times") +
  theme(panel.background = element_blank(), legend.background = element_blank(), panel.grid = element_blank(), plot.background = element_blank(), legend.text=element_text(size=rel(.8)), strip.text = element_text(size=11))+
  geom_label_repel(data=subset(TAB_loci, grepl("outlier",type)),aes(x=RDA1*20,y=RDA2*20,label = SNPID))

p2 <- ggplot() +
  geom_hline(yintercept=0, linetype="dashed", color = gray(.80), size=0.6) +
  geom_vline(xintercept=0, linetype="dashed", color = gray(.80), size=0.6) +
  geom_point(data = TAB_loci, aes(x=RDA2*20, y=RDA3*20, colour = type), size = 1.4) +
  scale_color_manual(values = c("gray90", "#F9A242FF", "#6B4596FF")) +
  geom_segment(data = TAB_var, aes(xend=RDA2, yend=RDA3, x=0, y=0), colour="black", size=0.15, linetype=1, arrow=arrow(length = unit(0.02, "npc"))) +
  geom_text(data = TAB_var, aes(x=1.1*RDA2, y=1.1*RDA3, label = row.names(TAB_var)), size = 2.5, family = "Times") +
  xlab("RDA 2") + ylab("RDA 3") +
  facet_wrap(~"RDA space") +
  guides(color=guide_legend(title="Locus type")) +
  theme_bw(base_size = 11, base_family = "Times") +
  theme(panel.background = element_blank(), legend.background = element_blank(), panel.grid = element_blank(), plot.background = element_blank(), legend.text=element_text(size=rel(.8)), strip.text = element_text(size=11))+
  geom_label_repel(data=subset(TAB_loci, grepl("outlier",type)),aes(x=RDA2*20,y=RDA3*20,label = SNPID))

print("plot")
png(file=paste0("RDA_biplot1v2_k",RDAK,"_9V.png"),width=800, height=800)
	print(p1)
	dev.off()
png(file=paste0("RDA_biplot2v3_k",RDAK,"_9V.png"),width=800, height=800)
	print(p2)
	dev.off()

  png(file=paste0("RDA_Manhattan",RDAK,"_EnvPC.png"),width=800, height=800)
	par(mfrow = c(2, 1))
  	manhattan(TAB_loci, chr="Chr",snp="SNPID", bp="Pos", p="p.value", ylim=c(0,10),col = c("blue4", "orange3"), suggestiveline = F, genomewideline = F, annotatePval = 0.05/length(TAB_loci$p.value))
  	FDR_input <- cbind.data.frame(TAB_loci$SNPID,TAB_loci$Chr,TAB_loci$Pos,-log10(TAB_loci$p.value))
  	FDR <- RAINBOWR::CalcThreshold(FDR_input, sig.level = 0.05, method = "BH")
  	abline(h=FDR, lty=2, col="orange")
  	abline(h=-log10(0.05/length(TAB_loci$p.value)), lty=2, col="red")
  	qqman::qq(TAB_loci$p.value)

  dev.off()

}
}


####
#### Now do it with regard to environmental PCs
####

RDA_env <- rda(Genotypes ~ EnvPC1 + EnvPC2 + EnvPC3 + EnvPC4 + EnvPC5 + Condition(PC1 + PC2 + PC3 + PC4 + PC5 + Longitude + Latitude),  Variables)
saveRDS(RDA_env,"RDA_EnvPC.rds")
print("Now that we have selected variables, let's run some RDA")
pRDAfull <- rda(Genotypes ~ PC1 + PC2 + PC3 + PC4 + PC5 + Longitude + Latitude + EnvPC1 + EnvPC2 + EnvPC3 + EnvPC4 + EnvPC5,  Variables)
pRDAfull_rq <- RsquareAdj(pRDAfull)
pRDAfull_anova <- anova(pRDAfull)

pRDAClimGeo <- rda(Genotypes ~  Longitude + Latitude + EnvPC1 + EnvPC2 + EnvPC3 + EnvPC4 + EnvPC5 + Condition(PC1 + PC2 + PC3 + PC4 + PC5),  Variables)
pRDAClimGeo_rq <- RsquareAdj(pRDAClimGeo)
pRDAClimGeo_anova <- anova(pRDAClimGeo)

pRDAclim <- rda(Genotypes ~ EnvPC1 + EnvPC2 + EnvPC3 + EnvPC4 + EnvPC5 + Condition(Longitude + Latitude + PC1 + PC2 + PC3 + PC4 + PC5),  Variables)
pRDAclim_rq <- RsquareAdj(pRDAclim)
pRDAclim_anova <- anova(pRDAclim)


## Pure neutral population structure model
pRDAstruct <- rda(Genotypes ~ PC1 + PC2 + PC3 + PC4 + PC5 + Condition(Longitude + Latitude + EnvPC1 + EnvPC2 + EnvPC3 + EnvPC4 + EnvPC5),  Variables)
pRDAstruct_rq <- RsquareAdj(pRDAstruct)
pRDAstruct_anova <- anova(pRDAstruct)

pRDAstructandGeo <- rda(Genotypes ~ PC1 + PC2 + PC3 + PC4 + PC5 + Longitude + Latitude + Condition(EnvPC1 + EnvPC2 + EnvPC3 + EnvPC4 + EnvPC5),  Variables)
pRDAstructandGeo_rq <- RsquareAdj(pRDAstructandGeo)
pRDAstructadnGeo_anova <- anova(pRDAstructandGeo)

pRDAstructandClim <- rda(Genotypes ~ PC1 + PC2 + PC3 + PC4 + PC5 + EnvPC1 + EnvPC2 + EnvPC3 + EnvPC4 + EnvPC5 + Condition(Longitude + Latitude),  Variables)
pRDAstructandClim_rq <- RsquareAdj(pRDAstructandClim)
pRDAstructandClim_anova <- anova(pRDAstructandClim)


##Pure geography model 
pRDAgeog <- rda(Genotypes ~ Longitude + Latitude + Condition(EnvPC1 + EnvPC2 + EnvPC3 + EnvPC4 + EnvPC5 + PC1 + PC2 + PC3 + PC4 + PC5),  Variables)
pRDAgeog_rq <- RsquareAdj(pRDAgeog)
pRDAgeog_anova <- anova(pRDAgeog)

#Get Some raw RDAs
RDA_raw_env <-  rda(Genotypes ~ EnvPC1 + EnvPC2 + EnvPC3 + EnvPC4 + EnvPC5,  Variables)
saveRDS(RDA_raw_env,"RDA_EnvPC_raw_Env.rds")
RDA_env_condGeo <-  rda(Genotypes ~ EnvPC1 + EnvPC2 + EnvPC3 + EnvPC4 + EnvPC5 + Condition(Longitude + Latitude),  Variables)
saveRDS(RDA_env_condGeo,"RDA_EnvPC_env_condGeo.rds")


##Write out partial models results
Partial_RDAtable <- cbind.data.frame(
Model=c("Full","Structure_And_Geography","Structure_And_Climate","Climate_And_Geography","Climate","Structure","Geography"),
Variance=c(pRDAfull_anova$Variance[1],
	   pRDAstructadnGeo_anova$Variance[1],
	   pRDAstructandClim_anova$Variance[1],
	   pRDAClimGeo_anova$Variance[1],
	   pRDAclim_anova$Variance[1],
	   pRDAstruct_anova$Variance[1],
	   pRDAgeog_anova$Variance[1]),
Rsquared=c(pRDAfull_rq$r.squared,
	   pRDAstructandGeo_rq$r,
	   pRDAstructandClim_rq$r,
	   pRDAClimGeo_rq$r,
	   pRDAclim_rq$r.squared,
	   pRDAstruct_rq$r.squared,
	   pRDAgeog_rq$r.squared),
Adjsquared=c(pRDAfull_rq$adj.r.squared,
	     pRDAstructandGeo_rq$adj.r.squared,
	     pRDAstructandClim_rq$adj.r.squared,
	     pRDAClimGeo_rq$adj.r.squared,
	     pRDAclim_rq$adj.r.squared,
	     pRDAstruct_rq$adj.r.squared,
	     pRDAgeog_rq$adj.r.squared),
pvalue=c(pRDAfull_anova[[4]][1],
	 pRDAstructadnGeo_anova[[4]][1],
	 pRDAstructandClim_anova[[4]][1],
	 pRDAClimGeo_anova[[4]][1],
	 pRDAclim_anova[[4]][1],
	 pRDAstruct_anova[[4]][1],
	 pRDAgeog_anova[[4]][1]))
write.table(Partial_RDAtable,"Partial_RDAtable_EnvPC.txt",col.names=T,row.names=F,sep="\t",quote=F)



for(RDAK in 3:5){
###First 5 RDA axes look good
## Function rdadapt
print(RDAK)
## Running the function with K = RDAK
rdadapt_env<-rdadapt(RDA_env, RDAK)
write.table(rdadapt_env,paste0("RDA_outputk",RDAK,"_EnvPC.txt"),col.names=T,row.names=T,quote=F,sep="\t")

## P-values threshold after Bonferroni correction
thres_env <- 0.01/length(rdadapt_env$p.values)
## Identifying the loci that are below the p-value threshold
#print(rdadapt_env)
outliers <- data.frame(Loci = colnames(Genotypes)[which(rdadapt_env$p.values<thres_env)], p.value = rdadapt_env$p.values[which(rdadapt_env$p.values<thres_env)], contig = unlist(lapply(strsplit(colnames(Genotypes)[which(rdadapt_env$p.values<thres_env)], split = "_"), function(x) x[1])))

print(outliers)
if(nrow(outliers >=1)){
## Top hit outlier per contig
outliers <- outliers[order(outliers$contig, outliers$p.value),]
## List of outlier names
outliers_rdadapt_env <- as.character(outliers$Loci[!duplicated(outliers$contig)])

## Formatting table for ggplot
locus_scores <- scores(RDA_env, choices=c(1:RDAK), display="species", scaling="none") # vegan references "species", here these are the loci
TAB_loci <- data.frame(names = row.names(locus_scores), locus_scores)
TAB_loci$type <- "Neutral"
TAB_loci$type[TAB_loci$names%in%outliers$Loci] <- "All outliers"
TAB_loci$type[TAB_loci$names%in%outliers_rdadapt_env] <- "Top outliers"
TAB_loci$type <- factor(TAB_loci$type, levels = c("Neutral", "All outliers", "Top outliers"))
TAB_loci <- cbind.data.frame(Genotype_Positions,TAB_loci,rdadapt_env$p.values)
colnames(TAB_loci)[ncol(TAB_loci)] <- "p.value"
write.table(TAB_loci,paste0("TAB_loci_k",RDAK,"_EnvPC.txt"),col.names=T,row.names=T,quote=F,sep="\t")
TAB_loci <- TAB_loci[order(TAB_loci$type),]
TAB_var <- as.data.frame(scores(RDA_env, choices=c(1:RDAK), display="bp")) # pull the biplot scores
write.table(TAB_var,"TAB_var_EnvPC.txt",col.names=T,row.names=T,quote=F,sep="\t")

## Biplot of RDA loci and variables scores
print(colnames(TAB_loci))
print(colnames(TAB_var))
print("Hey let's plot!")
p1 <- ggplot() +
  geom_hline(yintercept=0, linetype="dashed", color = gray(.80), size=0.6) +
  geom_vline(xintercept=0, linetype="dashed", color = gray(.80), size=0.6) +
  geom_point(data = TAB_loci, aes(x=RDA1*20, y=RDA2*20, colour = type), size = 1.4) +
  scale_color_manual(values = c("gray90", "#F9A242FF", "#6B4596FF")) +
  geom_segment(data = TAB_var, aes(xend=RDA1, yend=RDA2, x=0, y=0), colour="black", size=0.15, linetype=1, arrow=arrow(length = unit(0.02, "npc"))) +
  geom_text(data = TAB_var, aes(x=1.1*RDA1, y=1.1*RDA2, label = row.names(TAB_var)), size = 2.5, family = "Times") +
  xlab("RDA 1") + ylab("RDA 2") +
  facet_wrap(~"RDA space") +
  guides(color=guide_legend(title="Locus type")) +
  theme_bw(base_size = 11, base_family = "Times") +
  theme(panel.background = element_blank(), legend.background = element_blank(), panel.grid = element_blank(), plot.background = element_blank(), legend.text=element_text(size=rel(.8)), strip.text = element_text(size=11))+
  geom_label_repel(data=subset(TAB_loci, grepl("outlier",type)),aes(x=RDA1*20,y=RDA2*20,label = SNPID))

p2 <- ggplot() +
  geom_hline(yintercept=0, linetype="dashed", color = gray(.80), size=0.6) +
  geom_vline(xintercept=0, linetype="dashed", color = gray(.80), size=0.6) +
  geom_point(data = TAB_loci, aes(x=RDA2*20, y=RDA3*20, colour = type), size = 1.4) +
  scale_color_manual(values = c("gray90", "#F9A242FF", "#6B4596FF")) +
  geom_segment(data = TAB_var, aes(xend=RDA2, yend=RDA3, x=0, y=0), colour="black", size=0.15, linetype=1, arrow=arrow(length = unit(0.02, "npc"))) +
  geom_text(data = TAB_var, aes(x=1.1*RDA2, y=1.1*RDA3, label = row.names(TAB_var)), size = 2.5, family = "Times") +
  xlab("RDA 2") + ylab("RDA 3") +
  facet_wrap(~"RDA space") +
  guides(color=guide_legend(title="Locus type")) +
  theme_bw(base_size = 11, base_family = "Times") +
  theme(panel.background = element_blank(), legend.background = element_blank(), panel.grid = element_blank(), plot.background = element_blank(), legend.text=element_text(size=rel(.8)), strip.text = element_text(size=11))+
  geom_label_repel(data=subset(TAB_loci, grepl("outlier",type)),aes(x=RDA2*20,y=RDA3*20,label = SNPID))


png(file=paste0("RDA_biplot1v2_k",RDAK,"_EnvPC.png"),width=800, height=800)
        print(p1)
        dev.off()
png(file=paste0("RDA_biplot2v3_k",RDAK,"_EnvPC.png"),width=800, height=800)
        print(p2)
	dev.off()


  png(file=paste0("RDA_Manhattan",RDAK,"_EnvPC.png"),width=800, height=800)
	par(mfrow = c(2, 1))
  	manhattan(TAB_loci, chr="Chr",snp="SNPID", bp="Pos", p="p.value", ylim=c(0,10),col = c("blue4", "orange3"), suggestiveline = F, genomewideline = F, annotatePval = 0.05/length(TAB_loci$p.value))
  	FDR_input <- cbind.data.frame(TAB_loci$SNPID,TAB_loci$Chr,TAB_loci$Pos,-log10(TAB_loci$p.value))
  	FDR <- RAINBOWR::CalcThreshold(FDR_input, sig.level = 0.05, method = "BH")
  	abline(h=FDR, lty=2, col="orange")
  	abline(h=-log10(0.05/length(TAB_loci$p.value)), lty=2, col="red")
  	qqman::qq(TAB_loci$p.value)

  dev.off()

}
}



