#=============================================================================#
# Assignment_projectper3                                                      #
#                                                                             #                                                               
# Date:                                                                       #
# Authors: Ariadna Fosch i Muntan?, ID: I6215203, Maastricht University       #
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

install.packages("dplyr")
install.packages("tidyverse")

library(anamiR)
library(rstudioapi)
library(limma)
library(qvalue)
library(tidyverse)
# =================================================================================================================================
# 1. Import the file into your R session and define the main directories
# =================================================================================================================================
#Change the directories to your assigned folder.
DATA.DIR <- "C:/Users/patar/OneDrive/Documenten/Project period 3/Project/Data"

RESULTS.DIR <- "C:/Users/patar/OneDrive/Documenten/Project period 3/Project/Data"

#-----------------------------------------------------------------------------#
# 1.1 Import the text file containing the data into a new object
#-----------------------------------------------------------------------------#
# Provide an overview with the minimum value, maximum value, first quantile, median, mean and third quantile 
# values per measurement point per time point by means of the summary function.


# Get data in correct format
current_path <- getActiveDocumentContext()$path 

setwd(dirname(current_path ))

mrna = read.delim('../Data/GeneExpressionNormalized.txt', check.names = FALSE) #Load expression data
mirna = read.delim('../Data/miRNAexpression.txt', check.names = FALSE)

mrna[is.na(mrna)] <- NA # Change all the NaNs to NA in the mRNA object.

key <- mrna[,1:2] # Key for maintaining gene symbol and entrez gene ID


# Column names are different for miRNA and mRNA datasets, but are conserved with sample names in SampleGroups.xlsx
#IE FGS_01 is the name of sample 1; 
  #where miRNA sample is US10063773_254606410403_S01_miRNA_107_Sep09_1_1 - 2;
  #and mRNA sample is FGS_01_410978_1_1

#=== Here we change the columns names of mrna (mRNA samples) to the sample names ===
    # ...write code...
    # Pretinding to do stuff.....

#=== Here we change the columns names of mirna (miRNA samples) to the sample names ===
labels = read.delim('../Data/colNames.csv', sep = ',', header = FALSE, colClasses = 'character')

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

#=== Here we change the columns names of mrna (mRNA samples) to the sample names ===
mrna = subset(mrna, select = -2) # Remove gene symbol column
columns = colnames(mrna)
columns[2:29] <- substr(unlist(columns[2:29]), 1, 6) #Only keep first 6 characters of column names
colnames(mrna) <- columns # Assign this to the actual column names

#=== Here we initialise the sample grouping ===
sampleGroups <- read.delim('../Data/SampleGroups.csv', sep = ',', header = TRUE, colClasses = 'character')
sampleGroups <- sampleGroups[,6:7]

index_cholestatic <- sampleGroups == 'cholestatic'
sampleGroups$id[index_cholestatic[,1]] <- 1
index_drained <- sampleGroups == 'drained'
sampleGroups$id[index_drained[,1]] <- 2
index_control <- sampleGroups == 'control'
sampleGroups$id[index_control[,1]] <- 3

#========#

design_matrix <- model.matrix(~ 0 + factor(sampleGroups$id))

colnames(design_matrix) <- c("cholestasis", "drained", "control")

cont_matrix <- makeContrasts (drained_v_control = drained - control,
                              cholestasis_v_control = cholestasis - control,
                              cholestasis_v_drained = cholestasis - drained,
                              levels = design_matrix)

fit <- lmFit(mrna[, 2:29], design_matrix)

fit_contrast <- contrasts.fit(fit, cont_matrix)

fit_contrast <- eBayes(fit_contrast)

results <- decideTests(fit_contrast)

summary(results)

top_genes <- topTable (fit_contrast, p.value = "0.05", number = nrow(mrna), adjust = "BH")

for (i in 1:ncol(fit_contrast)){
  volcanoplot(fit_contrast[,i], main= colnames(fit_contrast)[i], col=ifelse(fit_contrast[,i]$p.value > 0.05,"red","black"))
  abline(-log10(0.05),0)
  abline(v=log2(2))
  abline(v=-log2(2))
}
#=======================================================================#

design_matrix <- model.matrix(~ 0 + factor(sampleGroups$id))

colnames(design_matrix) <- c("cholestasis", "drained", "control")

cont_matrix <- makeContrasts (drained_v_control = drained - control,
                              cholestasis_v_control = cholestasis - control,
                              cholestasis_v_drained = cholestasis - drained,
                              levels = design_matrix)

#omitmirna <- na.omit(mirna) #Do we need impute or omit or mean? 

fit2 <- lmFit(mirna[, 2:29], design_matrix)

fit_contrast2 <- contrasts.fit(fit2, cont_matrix)

fit_contrast2 <- eBayes(fit_contrast2)

results2 <- decideTests(fit_contrast2)

summary(results2)

top_genes2 <- topTable (fit_contrast2, p.value = "0.05", number = nrow(mirna), adjust = "BH") 
#cut top genes list based on adj. sign. threshold. 


for (i in 1:ncol(fit_contrast2)){
  volcanoplot(fit_contrast2[,i], main= colnames(fit_contrast2)[i], col=ifelse(fit_contrast2[,i]$p.value > 0.05,"red","black"))
  abline(-log10(0.05),0)
  abline(v=log2(2))
  abline(v=-log2(2))
}


#=========================================#
# The qvalue function is used to calculate the adjusted pvalues of all the different comparisons#
#==========================================#
# this is pairwise testing (I believe)

# Adjust all p values using qvalue() for the mRNA set

mrna_drained_v_control_pvalue <- qvalue(fit_contrast$p.value[,1])
mrna_cholestasis_v_control_pvalue <- qvalue(fit_contrast$p.value[,2])
mrna_cholestasis_v_drained_pvalue <- qvalue(fit_contrast$p.value[,3])

# Adjust all p values using qvalue() for the miRNA set
mirna_drained_v_control_pvalue <- qvalue(fit_contrast2$p.value[,1])
mirna_cholestasis_v_control_pvalue <- qvalue(fit_contrast2$p.value[,2])
mirna_cholestasis_v_drained_pvalue <- qvalue(fit_contrast2$p.value[,3])

#======== CHECK IF ADJUSTED P VALUE IS EQUAL TO Q VALUE======

cont_matrix <- makeContrasts (drained_v_control = drained - control,
                              levels = design_matrix)

#omitmirna <- na.omit(mirna) #Do we need impute or omit or mean? 

fit2 <- lmFit(mirna[, 2:29], design_matrix)

fit_contrast2 <- contrasts.fit(fit2, cont_matrix)

fit_contrast2 <- eBayes(fit_contrast2)

results2 <- decideTests(fit_contrast2)

summary(results2)

top_genes2 <- topTable (fit_contrast2, number = nrow(mirna), adjust = "BH") 

#==============
#cc <- cor(mrna[,2:29])
#dend <- as.dendrogram(hclust(as.dist(1-cc)))

#======== still fucking with this so don't mind this


mrna_batch_410978 <- mrna[, c(2, 5, 3, 4, 7, 6, 8, 9), drop=F]
mrna_batch_410979 <- mrna[, c(12, 15, 17, 10, 13, 16, 11, 14), drop=F]
mrna_batch_410980 <- mrna[, c(20, 23, 18, 21, 24, 19, 22, 25), drop=F]
mrna_batch_412287 <- mrna[, c(26, 28, 27, 29), drop=F]
  


# THis is just me trying to figure out how to make this work, but we probably should first do PCA and ANOVA to check
# for actual batch effects
y2 <- removeBatchEffect(batch_410978, batch_410979)
par(mfrow=c(1,2))
boxplot(as.data.frame(batch_410978),main="Original")
boxplot(as.data.frame(y2),main="Batch corrected")



controls<- subset(sampleGroups, Treatment== "control")
cholestatic<- subset(sampleGroups, Treatment== "cholestatic")
drained<- subset(sampleGroups, Treatment== "drained")

drained_mrna<-mrna[match(drained$SampleName, colnames(mrna))]
cholestatic_mrna<-mrna[match(cholestatic$SampleName, colnames(mrna))]
controls_mrna<-mrna[match(controls$SampleName, colnames(mrna))]

control_vs_chole<- pairwise.t.test(controls_mrna,cholestatic_mrna), p.adjust.method = "BH") #ask the group about how to make the comparison between the different groups


