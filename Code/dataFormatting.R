#=============================================================================#
# Project Period, Liver cholestasis data analysis         												  
#	
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


# Get data in correct format
setwd('/Users/Ariadna/Documents/GitHub/Project')
mrna = read.delim('Data/GeneExpressionNormalized.txt', check.names = FALSE) #Load expression data
mirna = read.delim('Data/miRNAexpression.txt', check.names = FALSE) #Load miRNA expression data


key <- mrna[,1:2] # Key for maintaining gene symbol and entrez gene ID
mrna = mrna[,2:30] # Genrate dataframe with expression values; no entrezID

# Column names are different for miRNA and mRNA datasets, but are conserved with sample names in SampleGroups.xlsx
#IE FGS_01 is the name of sample 1; 
  #where miRNA sample is US10063773_254606410403_S01_miRNA_107_Sep09_1_1 - 2;
  #and mRNA sample is FGS_01_410978_1_1

#=== Here we change the columns names of mrna (mRNA samples) to the sample names ===
    # ...write code...
    # Pretinding to do stuff.....

#=== Here we change the columns names of mirna (miRNA samples) to the sample names ===
labels = read.delim("Data/colNames.csv", sep = ',', header = FALSE, colClasses = 'character')

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

## PCA plots 

mirna[is.na(mirna)] <- NA # change Nan for NA
pcaRes <- pca(mirna[,2:29],nPcs = 10)  # perform PCA
PCAtry<- data.frame(c(pcaRes@scores[,1]),
                    pcaRes@scores[,2],
                    pcaRes@scores[,3],
                    pcaRes@scores[,4],
                    pcaRes@scores[,5],
                    treat=labels$V2) 
colnames(PCAtry) = c("PCA1", "PCA2", "PCA3", "PCA4", "PCA5","treat")

ggplot(PCAtry, aes(x = PCA1, y = PCA2)) +
  geom_point(aes(colour = PCAtry$treat)) +
  scale_colour_manual(values = c("#04179b", "#da9e00", "#198c19"),
                      aesthetics = "fill") +
  theme_light()

