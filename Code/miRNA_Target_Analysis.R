#=============================================================================#
# Project Period, Liver cholestasis data analysis         							      #					  
#	miRNA Target Analysis                                                       #
# Version: 1.0   															                                #
# Date: 9-1-2020											             	                          #
# Authors: Ariadna Fosch i Muntané, ID: I6215203, Maastricht University       #
#          Jip de Kok, ID: I6119367 , Maastricht University                   #
#          Ravin Schmidl, ID: I6125866, Maastricht University                 #
#          Stefan Meier, ID: I6114194 , Maastricht University                 #
#=============================================================================#
#=========================================#
##         Install libraries             ##
#=========================================#

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

if (!requireNamespace("miRBaseConverter", quietly = TRUE))
  BiocManager::install("miRBaseConverter")

library(miRBaseConverter)
library(jsonlite)

#=========================================#
##            Initialise Data            ##
#=========================================#
Data <- format()
mrna <- data.frame(Data[1])   #mRNA expression data (contains entrez ID as index)
mirna <- data.frame(Data[2])  #miRNA expression data (contains miRNA name as index)
labels <- data.frame(Data[3]) #batch and treatment id/labels for samples
key <- data.frame(Data[4])    #entrezID to genesymbol key

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
  for (i in 1:nrow(mirna)){ #There was a .case; what is this?
    for (j in 1:nrow(mrna)){
      step = step + 1
      # Calculate pearson correlation
      miRNA <- as.numeric(mirna[i,])
      mRNA <- as.numeric(mrna[j,])
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
  #==== miRNA version conversion ===#
  miRNA_names <- row.names(mirna)                                       # Store miRNA names
  version = checkMiRNAVersion(row.names(mirna), verbose = TRUE)         # Check miRNA name version
  miRNA_names <- miRNA_NameToAccession(miRNA_names, version = version)  # Include Accession numbers
  miRNA_names_v18 <- miRNA_AccessionToName(miRNA_names[,2],
                                           targetVersion = "v18")       # Convert to version 18
  
  #=== miRNA targer prediction - query ===#
  prefix <- "https://app1.bioinformatics.mdanderson.org/tarhub/_design/basic/_view/"
  link_temaplate <- paste(prefix, 'by_matureMIRcount?startkey=[%s,%s]&endkey=[%s,{}]',  sep = "")
  
  minSources = 1 # Set minimum number of sources
  mirna_targets <- list() # Initialise the mina_targets matrix
  
  list = rep( list(list(rep(0, 3))), nrow(mirna) ) 
  target = list()
  
  # Initiate progress bar
  progress <- txtProgressBar(min = 0, max = length(miRNA_names), initial = 0)
  
  step = 0
  nTargets = 0
  for (i in 1:length(miRNA_names)){
    step = step + 1
    name = paste("%22", miRNA_names[i], "%22", sep = "")          # Set name in correct URL format
    miRNA_link <- sprintf(link_temaplate, name, minSources, name) # Create URL link to miRNA targets
    target <- try(fromJSON(miRNA_link))                           # Query the miRNA - try returns error for server errors
    #names(target)[3] <- miRNA_names[i]                           # Set list item name to miRNA name                    
    mirna_targets <- append(mirna_targets, list(as.matrix(target$rows)))# Add the miRNA targets to one long list of targets
    nTargets <- nTargets +  nrow(as.matrix(target$rows))          # Keep track of the total number of targets
    
    # substr(unlist(columns[2:29]), 1, 6)
    
    # Update progress bar
    setTxtProgressBar(progress, step)
  }
}

exp.Correlation <- miRNA_correlate(mirna, mrna)


colnames(mirna_targets) <- c("miRNA", "mRNA", "sources")
strsplit(target$rows[1,1], ":")