#=============================================================================#
# Project Period, Liver cholestasis data analysis         												  
# Version: 1.0   															  
# Date: 9-1-2020											             	  
# Author:  Jip de Kok, Stefan Meier, Ariadna Fosch & Ravin Schmidl                                    
#=============================================================================#

#===================#
## Data Formatting ##
#===================#

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("anamiR")

library(anamiR)
library(rstudioapi)



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

# Set rownames for mrna
rownames(mrna) <- mrna[,1]
mrna <- mrna[,2:29]

# remove N

se <- normalization(data = mrna2, method = "quantile")

mrna_se <- SummarizedExperiment::SummarizedExperiment(
  assays = S4Vectors::SimpleList(counts=mrna),
  colData = pheno.mrna)

mirna_se <- SummarizedExperiment::SummarizedExperiment(
  assays = S4Vectors::SimpleList(counts=mirna),
  colData = pheno.mirna)



mrna_d <- differExp_discrete(se = mrna_se,
                             class = "ER", method = "t.test",
                             t_test.var = FALSE, log2 = TRUE,
                             p_value.cutoff = 0.05,  logratio = 0.5
)
