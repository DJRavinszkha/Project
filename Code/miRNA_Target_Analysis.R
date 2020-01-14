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
  
  # Correct for multiple testing
  trial <- corMatrix[1:500,]
  trial$"p-adjust" <- p.adjust(trial["p-value"], method = "BH")
  
  
  return(corMatrix)
}