library(data.table)
library(dplyr)
library(parallel)
args<- commandArgs(trailingOnly=T)
threads=100
LOFVCF <- fread(args[1],sep="\t",header=T)
#print(head(LOFVCF))
Genes <- fread("PrimaryTranscripts.txt", header = F)
#Genes <- fread("transcript_test.txt", header = F)
Genes <- unlist(as.vector(Genes[,1]))

Disorder_DF <- fread("DisorderedTails.txt.gz",header=T)

###This function is used to parse allele identity from the SNPEFF annotated string
ParseAllele <- function(SNPEFFstring){
  AlleleIdent <- gsub(".*ANN=","",SNPEFFstring)
  AlleleIdent <- strsplit(AlleleIdent,"\\|")[[1]][1]
  AlleleIdent <- unname(gsub("*ANN=","",AlleleIdent))
  return(AlleleIdent)
}

ParseConsequence <- function(SNPEFFstring){
  AlleleIdent <- gsub(".*ANN=","",SNPEFFstring)
  AlleleIdent <- strsplit(AlleleIdent,"\\|")[[1]][2]
  AlleleIdent <- unname(gsub("*ANN=","",AlleleIdent))
  return(AlleleIdent)
}


###This function will take Alt alleles identity, SNPEFF string, and transcript identity and identify which alleles are deleterious ("High") in this transcript
Get_LoF_Alleles <- function(genostrings,gene){
  Chr <-  unlist(c(genostrings[1]))
  Pos <-   unlist(c(genostrings[2]))
  AltAllString <- unlist(c(genostrings[5]))
  snpEFF_string <- unlist(c(genostrings[8]))
  LoF_Vec <- strsplit(snpEFF_string,split = ",")[[1]]
  LoF_Vec_gene <- LoF_Vec[grepl(gene,LoF_Vec)]
  #print(LoF_Vec_gene)
  AllIdent <- sapply(LoF_Vec_gene, function(x) ParseAllele(x))
  AllMatch <- match(strsplit(AltAllString,",")[[1]],AllIdent)
  LoF_del <- grepl("HIGH",LoF_Vec_gene)[AllMatch]
  #print(LoF_del)
  LoF_del[is.na(LoF_del)] <- F
  if(sum(LoF_del)>=1){
  	#Evaluate Consequences and if they happen in disordered tails
  	disorder_subdf <- Disorder_DF[Disorder_DF$transcript==gene,]
  	if(nrow(disorder_subdf)>0){
	##Got to add the stop codon positions to the DF to make sure it's good.  Adding it to beginning and end to account of iether orientation
#	lower_cushion <- cbind.data.frame(Chr=disorder_subdf[1,1],Pos=c(min(disorder_subdf$Pos,na.rm=T)-3,min(disorder_subdf$Pos,na.rm=T)-2,min(disorder_subdf$Pos,na.rm=T)-1),gene=disorder_subdf[1,3],AminoAcid="STOP",fivep_disord=0,threep_disord=1)
#	upper_cushion <- cbind.data.frame(Chr=disorder_subdf[1,1],Pos=c(max(disorder_subdf$Pos,na.rm=T)+1,max(disorder_subdf$Pos,na.rm=T)+2,max(disorder_subdf$Pos,na.rm=T)+3),gene=disorder_subdf[1,3],AminoAcid="STOP",fivep_disord=0,threep_disord=1) 
#	disorder_subdf <- rbind.data.frame(lower_cushion,
#					   disorder_subdf,
#					   upper_cushion
#	)
#	print(disorder_subdf)
  	fivep <- as.logical(disorder_subdf$fivep_disord[disorder_subdf$Chr==Chr & disorder_subdf$Pos==Pos])
  	threep <- as.logical(disorder_subdf$threep[disorder_subdf$Chr==Chr & disorder_subdf$Pos==Pos])
	##This line accounts for LoF mutations that fall outside the CDS such as splice lost
  	if(length(threep)<1){threep <- F}
	else if(length(threep)>1){threep <- threep[1]}
  	if(length(fivep)<1){fivep <- F}
	else if(length(fivep)>1){fivep <- fivep[1]}
	print(threep)
	print(fivep)
    	#Get the actual consequence of the mutations
    	All_Consequence <- sapply(LoF_Vec_gene, function(x) ParseConsequence(x))[AllMatch]
	
	##Need to simplify mutations with multiple effects so am instituting a hierarchy with 
	##Frameshift > stop gain > start lost > splice donor variant > splice acceptor variant > stop lost
        All_Consequence[grepl("frameshift_variant",All_Consequence)] <- "frameshift_variant"
	All_Consequence[grepl("stop_gained",All_Consequence)] <- "stop_gained"
	All_Consequence[grepl("start_lost",All_Consequence)] <- "start_lost"
	All_Consequence[grepl("splice_donor_variant",All_Consequence)] <- "splice_donor_variant"
	All_Consequence[grepl("splice_acceptor_variant",All_Consequence)] <- "splice_acceptor_variant"
	All_Consequence[grepl("stop_lost",All_Consequence)] <- "stop_lost"
	###Threep
	##Scrub all LoF mutations in threep disordered tail
  	if(threep){
		for(cons in 1:length(All_Consequence)){
                        if(LoF_del[cons]==T){
                        print(paste(gene,genostrings[1],genostrings[2],genostrings[5],"Three_prime",All_Consequence[cons],"Recovered",sep=" "))
  		LoF_del[cons] <- F
                        }
		}

  	}
	###Fivep
  	else if(fivep){
    		for(cons in 1:length(All_Consequence)){
      			if(!is.na(All_Consequence[cons])){
        			#If it's one of these 2 mutation types, we'll look for another start codon within the 
        			#5' disordered region.  But if it's a frameshift etc, I'm still calling it LoF
      				if(All_Consequence[cons]=="start_lost" | All_Consequence[cons]=="stop_gained"){
					#>3 because its for each base, not each amino acid
        				check_starts <- sum(disorder_subdf$AminoAcid=="MET" & disorder_subdf$fivep_disord==1) > 3
        				if(check_starts){
          					LoF_del[cons] <- F
						print(paste(gene,genostrings[1],genostrings[2],genostrings[5],"Five_prime",All_Consequence[cons], "Recovered",sep=" "))
        				}
					else{
						print(paste(gene,genostrings[1],genostrings[2],genostrings[5],"Five_prime", All_Consequence[cons], "Unrecovered",sep=" "))
					   }
      				}
      			else if(LoF_del[cons]==T){
				print(paste(gene,genostrings[1],genostrings[2],genostrings[5],"Five_prime", All_Consequence[cons],"Unrecovered",sep=" "))
				}
      			}
    		}	
  	}
	###Central
  	else{
		for(cons in 1:length(All_Consequence)){
			if(LoF_del[cons]==T){
			print(paste(gene,genostrings[1],genostrings[2],genostrings[5],"Central",All_Consequence[cons],"Unrecovered",sep=" "))	
			}
		}
	}
}
}
  return(LoF_del)
}



###Takes vector of Booleans corresponding to which Alternate alleles are deleterious and a SNP call ("0|1") then returns a new SNP call based on delteriousness
Allele_del <- function(LoFAllele,allele){
  allele_1 <- substring(allele,1,1)  
  allele_2 <- substring(allele,3,3)
  if(allele_1=="."){return(".")}
  if(allele_2=="."){return(".")}
  if(allele_1 %in% which(LoFAllele)){allele_1 <- "1"}else{allele_1 <- "0"}
  if(allele_2 %in% which(LoFAllele)){allele_2 <- "1"}else{allele_2 <- "0"}
  return(paste0(allele_1,"|",allele_2))
}

###Parses delterious level of phased SNP across all individuals for this gene/SNP combination
siteparse <- function(SNPSite, gene){
  LoF <- Get_LoF_Alleles(SNPSite, gene)
  SNPSite[10:length(SNPSite)] <- sapply(SNPSite[10:length(SNPSite)], function(SNP) Allele_del(LoF,substring(SNP,1,3)))
  return(SNPSite)
}
###Checks all SNPs in a haplotype and verifies if all are functional or if a deletetrious variant exists
CheckallFunctional <- function(funcvec){
  Allele_df <- data.frame(LOF=funcvec)
  Allele_df_split <- data.frame(do.call('rbind', strsplit(as.character(Allele_df$LOF),'|',fixed=TRUE)))
  Allele_vec_1 <- Allele_df_split[,1]
  Allele_vec_2 <- Allele_df_split[,2]
  #return(sum(1 %in% Allele_vec_1,1 %in% Allele_vec_2))
  if(sum(1 %in% Allele_vec_1,1 %in% Allele_vec_2)==2){return("1/1")}
  else if(sum(1 %in% Allele_vec_1,1 %in% Allele_vec_2)==1){return("0/1")}
  else{return("0/0")}
}
### subsets vcf and ensures LoF variants exist before parsing 
geneparse <- function(VCF,gene){
  genevcf <-    as.data.frame(VCF[grepl(gene,LOFVCF$INFO),])
  Chromosome <- substr(gene,7,8)
  if(substr(Chromosome,1,1)=="0"){Chromosome <- substr(Chromosome,2,2)}
  if(nrow(genevcf)>=1){
    genevcf <- as.data.frame(t(apply(genevcf, 1, function(genevcfrow) siteparse(genevcfrow, gene))))
    genefunc <- apply(as.data.frame(genevcf[,10:ncol(genevcf)]), 2, function(genevcfcol) CheckallFunctional(genevcfcol))
    genereport <- c(Chromosome,"1",gene,"A","T","100",".",".","GT",genefunc)
  }
  else{
    genereport <- c(Chromosome,"1",gene,"A","T","100",".",".","GT",rep("0/0",ncol(genevcf)-9))
  }
  genemat <- t(as.data.frame(unname(genereport)))
  rownames(genemat) <- gene

  return(t(genemat))
}
###Runs LOF Parsing on all genes 
LoF_GeneDF <- t(as.data.frame(mclapply(Genes,function(x) geneparse(LOFVCF,x)),mc.cores = threads))
#LoF_GeneDF <- t(as.data.frame(lapply(Genes,function(x) geneparse(LOFVCF,x))))
LoF_GeneDF <- LoF_GeneDF[order(LoF_GeneDF[,3]),]
#print(LoF_GeneDF)
LoF_GeneDF[,2] <- 1:nrow(LoF_GeneDF)
colnames(LoF_GeneDF) <- colnames(LOFVCF)
write.table(LoF_GeneDF,args[2],sep="\t",quote=F,col.names=T,row.names = F)

