#=============================================================================#
# Project Period, Liver cholestasis data analysis         												  
#	
# Version: 1.0   															  
# Date: 9-1-2020											             	  
# Author:  Jip de Kok, Stefan Meier, Ariadna Fosch & Ravin Schmidl 
# Note: adapted from "Introduction to anamiR - 2019-05-09"
# Ref: Wang T, Lu T, Lee C, Lai L, Tsai M, Chuang E (2015).
#     "anamiR-an integrated analysis package of miRNA and mRNA expression." Manuscript submitted for publication.
#     https://bioconductor.org/packages/release/bioc/vignettes/anamiR/inst/doc/IntroductionToanamiR.html#phenotype-data
#=============================================================================#

#===================#
## Data Formatting ##
#===================#

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("anamiR")

library(anamiR)
library(rstudioapi)

format <- function(){
  
  # Get data in correct format
  current_path <- getActiveDocumentContext()$path 
  setwd(dirname(current_path ))
  mrna = read.delim('../Data/GeneExpressionNormalized.txt', check.names = FALSE) #Load expression data
  mirna = read.delim('../Data/miRNAexpression.txt', check.names = FALSE) #Load miRNA expression data
  
  
  key <- mrna[,1:2] # Key for maintaining gene symbol and entrez gene ID
  
  
  # Column names are different for miRNA and mRNA datasets, but are conserved with sample names in SampleGroups.xlsx
  #IE FGS_01 is the name of sample 1; 
  #where miRNA sample is US10063773_254606410403_S01_miRNA_107_Sep09_1_1 - 2;
  #and mRNA sample is FGS_01_410978_1_1
  
  #=== Here we change the columns names of mrna (mRNA samples) to the sample names ===
  mrna = subset(mrna, select = -2) # Remove gene symbol column
  columns = colnames(mrna)
  columns[2:29] <- substr(unlist(columns[2:29]), 1, 6) #Only keep first 6 characters of column names
  colnames(mrna) <- columns # Assign this to the actual column names
  
  #=== Here we change the columns names of mirna (miRNA samples) to the sample names ===
  labels = read.delim("../Data/colNames.csv", sep = ',', header = FALSE, colClasses = 'character')
  
  
  columns = colnames(mirna)
  columns[[1]] <- 'miRNA'
  for (i in rep(2:length(colnames(mirna)))){
    for ( j in rep(1:nrow(labels))){
      if(colnames(mirna[i]) == labels[j,1]){
        columns[[i]] = labels[j,3] # Assign correct column name to array
        # print(labels[j,1])
        # print(colnames(mirna[i]))
        # print(labels[j, 3])
        # print("======")
      }
    }
  }
  colnames(mirna) <- columns
  
  #=== Here we initialise the sample grouping ===
  sampleGroups <- read.delim("Data/SampleGroups.csv", sep = ',', header = TRUE, colClasses = 'character')
  sampleGroups <- sampleGroups[,6:7]
  
  index_cholestatic <- sampleGroups == 'cholestatic'
  sampleGroups$id[index_cholestatic[,1]] <- 1
  index_drained <- sampleGroups == 'drained'
  sampleGroups$id[index_drained[,1]] <- 2
  index_control <- sampleGroups == 'control'
  sampleGroups$id[index_control[,1]] <- 3
  
  #=============================#
  # Create phenotype data       #
  #=============================#
  pheno.mrna = labels
  pheno.mrna[[1]] <- colnames(mrna)[2:29]
  colnames(pheno.mrna) <- c("", "Subtype", "ER")
  
  # Set controls
  case <- pheno.mrna[,'Subtype'] == "cholestatic"
  pheno.mrna[case,3] <- 'case'
  
  # Set drained
  case <- pheno.mrna[,'Subtype'] == "drained"
  pheno.mrna[case,3] <- 'case'
  
  # set controls
  case <- pheno.mrna[,'Subtype'] == "control"
  pheno.mrna[case,3] <- 'control'
  
  # Order cases and controls
  pheno.mrna <- pheno.mrna[order(as.character(pheno.mrna$ER)),]
  
  # Set pheno.mirna which is identical to pheno.mrna
  pheno.mirna = pheno.mrna
  
  
  #=============================#
  # Summarised experiment class #
  #=============================#
  
  # First we change the dataframes into matrices as the miRrna package works with matrices.
  pheno.mrna <- as.matrix(pheno.mrna)
  pheno.mirna <- as.matrix(pheno.mirna)
  mrna = data.matrix(mrna)
  mirna = data.matrix(mirna)
  
  # Set rownames for mrna
  rownames(mrna) <- mrna[,1]
  mrna <- mrna[,2:29]
  
  # Set rownames for mirna
  rownames(mirna) <- mirna[,1]
  mirna <- mirna[,2:29]
  
  # Fill NA's with mean for the time-being
  mrna[is.na(mrna)] <- mean(mrna, na.rm = TRUE)
  mirna[is.na(mirna)] <- mean(mirna, na.rm = TRUE)
  
  return(mrna, mirna, labels, pheno.mrna, pheno.mirna, key)
}


