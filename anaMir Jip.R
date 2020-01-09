#=============================================================================#
# Assignment_projectper3                                                      #
#                                                                             #                                                               
# Date:                                                                       #
# Authors: Ariadna Fosch i Muntané, ID: I6215203, Maastricht University       #
#          Jip de Kok, ID: I6119367 , Maastricht University                   #
#          Ravin Schmidl, ID: I6125866, Maastricht University                 #
#          Stefan Meier, ID: I6114194 , Maastricht University                 #
#                                                                             #
#                                                                             #
#                                                                             #                                                            
#                                                                             #
#=============================================================================#

# =================================================================================================================================
# 0. Load required packages:
# =================================================================================================================================
# Note: Please install any packages that you are missing.

setRepositories(FALSE, 1:9)

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("anamiR")

library(anamiR)

options(stringsAsFactors = F)

# =================================================================================================================================
# 1. Import the file into your R session and define the main directories
# =================================================================================================================================
#Change the directories to your assigned folder.
DATA.DIR <- "C:\\Users\\patar\\OneDrive\\Documenten\\GitHub\\Project\\Data"

RESULTS.DIR <- "C:\\Users\\patar\\OneDrive\\Documenten\\GitHub\\Project\\Data"

#-----------------------------------------------------------------------------#
# 1.1 Import the text file containing the data into a new object
#-----------------------------------------------------------------------------#
# Provide an overview with the minimum value, maximum value, first quantile, median, mean and third quantile 
# values per measurement point per time point by means of the summary function.

setwd(DATA.DIR)

df = read.delim('GeneExpressionNormalized.txt', check.names = FALSE) # Load mRNA expression data.

df2 = read.delim('miRNAexpression.txt', check.names = FALSE) # Load miRNA expression data.

df[is.na(df)] <- NA # Change all the NaNs to NA in the mRNA object.

summary(df)

summary(df2)

key <- df[,1:2] # Key for maintaining gene symbol and entrez gene ID

mrna = df[,2:30] # Generate dataframe with expression values; no entrezID

# Column names are different for miRNA and mRNA datasets, but are conserved with sample names in SampleGroups.xlsx
#IE FGS_01 is the name of sample 1; 
  #where miRNA sample is US10063773_254606410403_S01_miRNA_107_Sep09_1_1 - 2;
  #and mRNA sample is FGS_01_410978_1_1

#=== Here we change the columns names of df (mRNA samples) to the sample names ===
    #...write code...

#=== Here we change the columns names of df2 (miRNA samples) to the sample names ===
labels = read.delim("colNames.csv", sep = ',', header = FALSE, colClasses = 'character')

columns = colnames(df2)
columns[[1]] <- 'miRNA'
for (i in rep(2:length(colnames(df2)))){
  for (j in rep(1:nrow(labels))){
    if(colnames(df2[i]) == labels[j,1]){
      columns[[i]] = labels[j,3] # Assign correct column name to array
      # print(labels[j,1])
      # print(colnames(df2[i]))
      # print(labels[j, 3])
      # print("======")
    }
  }
}
colnames(df2) <- columns

summary(df2)

mirna_se <- SummarizedExperiment::SummarizedExperiment(
  assays = S4Vectors::SimpleList(counts=df2),
  colData = labels)


