#=============================================================================#
# Project Period, Liver cholestasis data analysis         												  
#	
# Version: 1.0   															  
# Date: 13-1-2020											             	  
# Author:  Jip de Kok, Stefan Meier, Ariadna Fosch & Ravin Schmidl 
#=============================================================================#


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
  return(corMatrix)
}