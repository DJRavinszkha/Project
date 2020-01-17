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
library(jsonlite)

#========================================================================================================#
# This function calculates pearson correlations between micro- and messenger RNA expression              #
# You will input the RNA expression values for your group of interest (e.g. all cases)                   #
# The output is a correlation matrix, showing the correlation between each possible pair of micro- and   #
# messenger RNA and its significance value (p-value)                                                     #
#========================================================================================================#
miRNA_correlate <- function(mirna, mrna){
  numberOfPairs <- nrow(mrna) * nrow(mirna)                           # Calculate number of miRNA-mRNA pairs
  corMatrix <- matrix(nrow = numberOfPairs, ncol = 4)                 # Initialise corMatix which will store all correlations
  colnames(corMatrix) <- c("miRNA", "mRNA", "pearson-cor", "p-value") # Set column names
  
  progress <- txtProgressBar(min = 0, max = numberOfPairs, initial = 0) # Initiate progress bar
  
  step = 0 # Initialise step which will count total number of correlations made
  for (i in 1:nrow(mirna.case)){
    for (j in 1:nrow(mrna.case)){
      step = step + 1
      
      # Calculate pearson correlation
      miRNA <- as.numeric(mirna.case[i,])                   # Load miRNA expression values
      mRNA <- as.numeric(mrna.case[j,])                     # Load mRNA expression values
      correlation <- stats::cor.test(miRNA,                 # Calculate miRNA-mRNA correlation
                                     mRNA,
                                     method = "pearson",
                                     use = "complete.obs")
      
      # Store correlation results
      miRNA_name <- row.names(mirna)[i]                     # Load miRNA name
      mRNA_name <- row.names(mrna)[j]                       # Load mRNA name
      corMatrix[step,] <- c(miRNA_name,                     # store results in corMatrix
                            mRNA_name,
                            correlation$estimate[[1]],
                            correlation$p.value)
      
      
      setTxtProgressBar(progress, step)                     # Update progress bar
    }
  }
  
  # Correct for multiple testing (using the FDR method of Benjamini Hochberg)
  corMatrix$"p.adjust" <- p.adjust(corMatrix$p.value, method = "BH")
  
  
  return(corMatrix)
}

#========================================================================================================#
# This function lookup miRNA targets from various databases using the targetHUB API                      #
#                                                                                                        #
# https://app1.bioinformatics.mdanderson.org/tarhub/_design/basic/index.html                             #
#                                                                                                        #
# miRBase: integrating microRNA annotation and deep-sequencing data.                                     #
# Kozomara A, Griffiths-Jones S.                                                                         #
# Nucleic Acids Res. 2011 39:D152-D157                                                                   #
#========================================================================================================#
miRNA_target_query <- function(mirna, mrna){
  #==== miRNA version conversion - Currently not used!!! ===#
  miRNA_names <- row.names(mirna)                                       # Store miRNA names
  version = checkMiRNAVersion(row.names(mirna), verbose = TRUE)         # Check miRNA name version
  miRNA_names <- miRNA_NameToAccession(miRNA_names, version = version)  # Include Accession numbers
  miRNA_names_v18 <- miRNA_AccessionToName(miRNA_names[,2],
                                           targetVersion = "v18")       # Convert to version 18
  
  #=== miRNA targer prediction - query ===#
  prefix <- "https://app1.bioinformatics.mdanderson.org/tarhub/_design/basic/_view/"
  link_temaplate <- paste(prefix, 'by_matureMIRcount?startkey=[%s,%s]&endkey=[%s,{}]',  sep = "")
  
  minSources = 1          # Set minimum number of sources
  mirna_targets <- list() # Initialise the mina_targets matrix
  
  progress <- txtProgressBar(min = 0, max = nrow(miRNA_names), initial = 0) # Initiate progress bar
  
  options(timeout = 100000) # Set timeout to increase time for retrying queries
  error_list = list()
  error_count = 0
  step = 0
  nTargets = 0
  for (i in 1:nrow(miRNA_names)){
    name = paste("%22", miRNA_names[i,1], "%22", sep = "")        # Set name in correct URL format
    miRNA_link <- sprintf(link_temaplate, name, minSources, name) # Create URL link to miRNA targets
    target <- try(fromJSON(miRNA_link))                           # Query the miRNA - try returns error for server errors
    
    if (typeof(target) == 'char'){
      error_list <- append(error_list, c(name, target))           # Capture error
      error_count <- error_count + 1                              # Keep track of number of errors
      next                                                        # Jump to next iteration
    }
    
    mirna_targets <- append(mirna_targets, list(target$rows))     # Add the miRNA targets to one long list of targets
    nTargets <- nTargets +  nrow(as.matrix(target$rows))          # Keep track of the total number of targets
    
    step = step + 1
    setTxtProgressBar(progress, step)                             # Update progress bar
  }
  print(sprintf("miRNA_target_query executed succesfully with %s errors", error_count))
  
  #=== The section below converts the mirna_targets list into one long data frame ===#
  progress <- txtProgressBar(min = 0, max = step, initial = 0)              # Initiate progress bar
  mirna_targets_df <- data.frame(matrix(nrow = nTargets, ncol = 3))         # Initialise mirna_targets_df
  colnames(mirna_targets_df) <- c("miRNA", "mRNA", "sources")               # set the column names
  index = 1                                                                 # Initialise index
  for (i in 1:step){
    targets <- data.frame(mirna_targets[[i]])           # load targets from mirna_targets
    if (length(targets) == 0){                          # If no targets, then skip to next iteration
      next
    }
    index_end <- index + nrow(targets) - 1              # Set end of index range
    mirna_targets_df[index:index_end,1:2] <- 
      t(data.frame(strsplit(targets[,1], ":")))         # Set miRNA and mRNA names
    mirna_targets_df[index:index_end,3] <- 
      as.numeric(unlist(targets[,2])[c(FALSE, TRUE)])   # Set number of sources
    index <- index_end + 1
    setTxtProgressBar(progress, i)                      # Update progress bar
  }
}

