if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("anamiR")

library(anamiR)


# Get data in correct format
df = read.delim('/Users/Ariadna/Documents/GitHub/Project/wetransfer-cb1c0d/PRO4002_Data/DataProject/GeneExpressionNormalized.txt', check.names = FALSE) #Load expression data
df2 = read.delim('/Users/Ariadna/Documents/GitHub/Project/wetransfer-cb1c0d/PRO4002_Data/DataProject/miRNAexpression.txt', check.names = FALSE) #Load miRNA expression data


key <- df[,1:2] # Key for maintaining gene symbol and entrez gene ID
mrna = df[,2:30] # Genrate dataframe with expression values; no entrezID

# Column names are different for miRNA and mRNA datasets, but are conserved with sample names in SampleGroups.xlsx
#IE FGS_01 is the name of sample 1; 
  #where miRNA sample is US10063773_254606410403_S01_miRNA_107_Sep09_1_1 - 2;
  #and mRNA sample is FGS_01_410978_1_1

#=== Here we change the columns names of df (mRNA samples) to the sample names ===
    #...write code...

#=== Here we change the columns names of df2 (miRNA samples) to the sample names ===
labels = read.delim("Data/colNames.csv", sep = ',', header = FALSE, colClasses = 'character')

columns = colnames(df2)
columns[[1]] <- 'miRNA'
for (i in rep(2:length(colnames(df2)))){
  for ( j in rep(1:nrow(labels))){
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


df3 = merge(df2, labels, by.x = colnames(df2[,2:29]), by.y = "V1")
