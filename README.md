This repository accompanies the paper:
Zhao et al. <I>Unlocking Genetic Diversity in Colombian Cassava Landraces for Accelerated Breeding</I>.
It contains figures, supplemental data, and code used in the analyses described in the manuscript.

This repository includes scripts and pipelines for:
1. <B>Variant Calling Pipeline</B>: Scripts for processing raw sequencing data and calling variants in cassava genomes.
2. <B>Population Structure Analyses</B>, including Principal Component Analysis (PCA), ADMIXTURE analyses, Linkage Disequilibrium analyses.
3. <B>Imputation of Missing Variants</B> - BEAGLE.sh: Wrapper script for genotype imputation using BEAGLE.
4. <B>Environmental Variation Correlation Plot</B> - Corplot_LOF_Rectangle.R: Generates a correlation heatmap used in environmental analyses.
5. <B>Genome-Wide Association Studies (GWAS)</B> - GLM_GWAS_LOF.sh: Runs GWAS in TASSEL, incorporating the first five principal components as covariates.
6. <B>Loss-of-Function Variant Parsing</B> - LoF_Parser_phased_disordtails.R: Processes SnpEff-annotated VCF files, integrates custom annotations of disordered protein regions, produces gene-level VCFs for loss-of-function (LoF) mutations, and reports mutation classes, effects, and positions relative to protein structure.
7. <B>Interactive Mapping of Colombian Cassava Diversity</B> - Mapping_shinyapp.R: R Shiny app for visualizing multiple dimensions of genetic and environmental variation across Colombian cassava populations.
8. <B>Redundancy Analysis (RDA)</B> - RDA_Colombia_V2.R: Performs RDA and partial RDA analyses testing the effects of population structure, geography, and climate on genetic variation.
Code authorship: Kehan Zhao & Evan Long

<B>Data Availability</B>:
Raw sequencing data generated in this study are available in the NCBI Sequence Read Archive (SRA) under BioProject accession PRJNA1228154.
Sample metadata for both new and previously published datasets are provided in Table S1 of the manuscript.
For access to additional post-processed data formats or results, please contact the corresponding authors.
