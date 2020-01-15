#=============================================================================#
# Project Period, Liver cholestasis data analysis         												  
#	
# Version: 1.0   															  
# Date: 13-1-2020											             	  
# Author:  Jip de Kok, Stefan Meier, Ariadna Fosch & Ravin Schmidl 
#=============================================================================#

#=== Install and load required packages ===#

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

if (!requireNamespace("miRBaseConverter", quietly = TRUE))
  BiocManager::install("miRBaseConverter")

library(miRBaseConverter)

#==== miRNA version conversion ===#
miRNA_names <- row.names(mirna) #                                                 Store miRNA names
version = checkMiRNAVersion(row.names(mirna), verbose = FALSE) #                   Check miRNA name version
miRNA_names <- miRNA_NameToAccession(miRNA_names, version = version) #            Include Accession numbers
miRNA_names_v22 <- miRNA_AccessionToName(miRNA_names[,2],targetVersion = "v22")#  Convert to version 22

#========================================================================================================#
# This function calculates pearson correlations between micro- and messenger RNA expression              #
# You will input the RNA expression values for your group of interest (e.g. all cases)                   #
# The output is a correlation matrix, showing the correlation between each possible pair of micro- and   #
# messenger RNA and its significance value (p-value)                                                     #
#========================================================================================================#
miRNA_correlate <- function(mirna, mrna){
  numberOfPairs <- nrow(mrna) * nrow(mirna)
  corMatrix <- matrix(nrow = numberOfPairs, ncol = 4)
  colnames(corMatrix) <- c("miRNA", "mRNA", "pearson-cor", "p-value")
  
  # Initiate progress bar
  progress <- txtProgressBar(min = 0, max = numberOfPairs, initial = 0)
  
  step = 0
  for (i in 1:nrow(mirna.case)){
    for (j in 1:nrow(mrna.case)){
      step = step + 1
      # Calculate pearson correlation
      miRNA <- as.numeric(mirna.case[i,])
      mRNA <- as.numeric(mrna.case[j,])
      correlation <- stats::cor.test(miRNA, mRNA, method = "pearson", use = "complete.obs")
      
      # Store correlation results
      miRNA_name <- row.names(mirna)[i]
      mRNA_name <- row.names(mrna)[j]
      corMatrix[step,] <- c(miRNA_name, mRNA_name, correlation$estimate[[1]], correlation$p.value)
      
      # Update progress bar
      setTxtProgressBar(progress, step)
    }
  }
  
  # Correct for multiple testing (using the FDR method of Benjamini Hochberg)
  corMatrix$"p.adjust" <- p.adjust(corMatrix$p.value, method = "BH")
  
  
  return(corMatrix)
}

miRNA_target_query <- function(mirna, mrna){
  
  prefix <- "https://app1.bioinformatics.mdanderson.org/tarhub/_design/basic/_view/"
  link_temaplate <- paste(prefix, 'by_matureMIRcount?startkey=[%s,%s]&endkey=[%s,{}]',  sep = "")
  
  name = dQuote("hsa-mir-212") 
  name <- "%22hsa-mir-29a%22"
  minSources = 3
  miRNA_link <- sprintf(link_temaplate, name, minSources, name)
  
  targets <- fromJSON(miRNA_link)
  mirna_targets <- matrix(nrow = nrow(mirna), ncol = 1)
  
  for (i in 1:nrow(mirna)){
    name = tolower(paste("%22", row.names(mirna)[i], "%22", sep = ""))
    miRNA_link <- sprintf(link_temaplate, name, minSources, name)
    mirna_targets[i] <- fromJSON(miRNA_link)
  }
}

