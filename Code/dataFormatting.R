#=============================================================================#
# Project Period, Liver cholestasis data analysis         							      #					  
#	Data Formatting                                                             #
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
library(limma)
library(rstudioapi)

#===================#
## Data Formatting ##
#===================#

format <- function(){
  # Get data in correct format
  current_path <- getActiveDocumentContext()$path 
  setwd(dirname(current_path ))
  mrna = read.delim('../Data/GeneExpressionNormalized.txt', check.names = FALSE) #Load expression data
  mirna = read.delim('../Data/miRNAexpression.txt', check.names = FALSE) #Load miRNA expression data
  colnames(mrna)[10] <- "FGS_10_410978_2_3_H" #Change the name of one variable due to error
  
  #Initialise a dataframe called labels, containing all:
      # file/batch names and sample names 
      # batch numbers and an int id attached to it
      # treatments and an int id attached to it
  labels = read.delim("../Data/colNames.csv", sep = ',', header = FALSE, colClasses = 'character')
  
  labels$V4 <- colnames(mrna[3:30])
  key <- mrna[,1:2] # Key for maintaining gene symbol and entrez gene ID
  
  labels$V1 <- colnames(mirna[2:29]) #maintaining names for mirna filenames
  
  #Here we initialise batch number id's of mRNA's
  labels$V5 <- substr(unlist(labels$V4), 8, 13) 
  labels$V6[labels$V5 == "410978"] <- 1 #set id 1 to mRNA batch
  labels$V6[labels$V5 == "410979"] <- 2 #set id 2
  labels$V6[labels$V5 == "410980"] <- 3 #set id 3
  labels$V6[labels$V5 == "412287"] <- 4 #set id 4
  
  #Here we initialise batch number id's for miRNA's
  labels$V7 <- substr(unlist(labels$V1), 12, 23)
  labels$V8[labels$V7 == "254606410403"] <- 1 #set id 1 to miRNA batch
  labels$V8[labels$V7 == "254606410404"] <- 2 #set id 2               
  labels$V8[labels$V7 == "254606411109"] <- 5 #set id 3                
  labels$V8[labels$V7 == "254606410405"] <- 3 #set id 4                
  labels$V8[labels$V7 == "254606410413"] <- 4 #set id 5                
  
  #Change the names of the columns
  colnames(labels) <- c("miRNA.file", "treatment", "sample.name", "mRNA.file", "mRNA.batch", "mRNA.batch.id", "miRNA.batch", "miRNA.batch.id")
  
  #=== Here we initialise the treatment numbers ===
  sampleGroups <- read.delim("../Data/SampleGroups.csv", sep = ',', header = TRUE, colClasses = 'character')
  sampleGroups <- sampleGroups[,6:7]
  
  labels$treatment.id[labels$treatment == "cholestatic"] <- 1 #set id 1 to cholestatic treatment
  labels$treatment.id[labels$treatment == "drained"] <- 2 #set id 2 to drained treatment
  labels$treatment.id[labels$treatment == "control"] <- 3 #set id 3 to control treatment
  
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

setRepositories(FALSE, 1:9)

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("anamiR")

library(anamiR)
library(rstudioapi)
library(limma)
library(qvalue)
library(tidyverse)

# =================================================================================================================================
# 1. Import the file into your R session and define the main directories
# =================================================================================================================================

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
  colnames(mirna) <- columns
  
  #mrna
  rownames(mrna) <- mrna[,1] # Set entrezid as index for mrna; can use key to find gene symbols based on index
  mrna <- mrna[,2:29] #Set mrna to only contain expression values
  
  #mirna
  rownames(mirna) <- mirna[,1] # Set mirna names as index for mirna
  mirna <- mirna[,2:29] #set mirna to only contain expression values
  
  #=============================#
  ##   Create phenotype data   ##     # Do we need this?
  #=============================#
  # pheno.mrna = labels
  # pheno.mrna[[1]] <- colnames(mrna)
  # colnames(pheno.mrna) <- c("", "Subtype", "ER")
  # 
  # # Set controls
  # case <- pheno.mrna[,'Subtype'] == "cholestatic"
  # pheno.mrna[case,3] <- 'case'
  # 
  # # Set drained
  # case <- pheno.mrna[,'Subtype'] == "drained"
  # pheno.mrna[case,3] <- 'case'
  # 
  # # set controls
  # case <- pheno.mrna[,'Subtype'] == "control"
  # pheno.mrna[case,3] <- 'control'
  # 
  # # Order cases and controls
  # pheno.mrna <- pheno.mrna[order(as.character(pheno.mrna$ER)),]
  # 
  # # Set pheno.mirna which is identical to pheno.mrna
  # pheno.mirna = pheno.mrna
  # 
  #===============================#
  ## Summarised experiment class ##   #Do we need this?
  #===============================#
  
  # First we change the dataframes into matrices as the miRrna package works with matrices.
  # pheno.mrna <- as.matrix(pheno.mrna)
  # pheno.mirna <- as.matrix(pheno.mirna)
  # mrna = data.matrix(mrna)
  # mirna = data.matrix(mirna)
  # 
  # # Fill NA's with mean for the time-being
  # mrna[is.na(mrna)] <- mean(mrna, na.rm = TRUE)
  # mirna[is.na(mirna)] <- mean(mirna, na.rm = TRUE)
  # 
  # # Remove drained cases from mrna
  # index = rep(TRUE, length(colnames(mrna)))
  # for(i in 1:length(colnames(mrna))){
  #   # if current column name is present in a list of all column names of cases:
  #   if(colnames(mrna)[i] %in% pheno.mrna[pheno.mrna[,2]=="drained",1]){
  #     index[i] <- FALSE
  #   }
  # }
  # mrna <- mrna[,index]
  # 
  # # Remove drained cases from mirna
  # index = rep(TRUE, length(colnames(mirna)))
  # for(i in 1:length(colnames(mirna))){
  #   # if current column name is present in a list of all column names of cases:
  #   if(colnames(mirna)[i] %in% pheno.mirna[pheno.mirna[,2]=="drained",1]){
  #     index[i] <- FALSE
  #   }
  # }
  # mirna <- mirna[,index]
  # 
  # # Remove drained cases from the phenotype data
  # pheno.mrna <- pheno.mrna[!(pheno.mrna[,2]=="drained"),]
  # pheno.mirna <- pheno.mirna[!(pheno.mirna[,2]=="drained"),]
  

  
  
  return(list(mrna, mirna, labels, key))
}
=======
colnames(mirna) <- columns

#=== Here we initialise the sample grouping ===
sampleGroups <- read.delim("../Data/SampleGroups.csv", sep = ',', header = TRUE, colClasses = 'character')
sampleGroups <- sampleGroups[,6:7]

index_cholestatic <- sampleGroups == 'cholestatic'
sampleGroups$id[index_cholestatic[,1]] <- 1
index_drained <- sampleGroups == 'drained'
sampleGroups$id[index_drained[,1]] <- 2
index_control <- sampleGroups == 'control'
sampleGroups$id[index_control[,1]] <- 3

# =================================================================================================================================
# 2. DIFFERENTIALLY EXPRESSED GENES WITH LMFIT
# =================================================================================================================================

#======== For mRNA ===================#
design_matrix_mrna <- model.matrix(~ 0 + factor(c(1,1,1,1,1,1,1,1,1,
                                                  2,2,2,2,2,2,2,2,2,2,
                                                  3,3,3,3,3,3,3,3,3)))

colnames(design_matrix_mrna) <- c("cholestasis", "drained", "control")

cont_matrix <- makeContrasts (drained_v_control = drained - control,
                              cholestasis_v_control = cholestasis - control,
                              cholestasis_v_drained = cholestasis - drained,
                              levels = design_matrix_mrna)

fit <- lmFit(mrna[, 2:29], design_matrix_mrna)

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
#======== For miRNA ===================#

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
# 3. Q-value test:The qvalue function is used to calculate the adjusted pvalues of all the different comparisons#
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

#======== Ari's CHECK IF ADJUSTED P VALUE IS EQUAL TO Q VALUE======

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

y2 <- removeBatchEffect(mrna[,2:29], mrna_batch_410978[,2:length(mrna_batch_410978)], mrna_batch_410980[,2:length(mrna_batch_410980)])

par(mfrow=c(1,2))
boxplot(as.data.frame(mrna[,3:29]),main="Original")
boxplot(as.data.frame(y2),main="Batch corrected")

#=========================================#
# 4. PAIRWISE T TEST: Separation on groups and controls. (Still work in progress)
#==========================================#

controls<- subset(sampleGroups, Treatment== "control")
cholestatic<- subset(sampleGroups, Treatment== "cholestatic")
drained<- subset(sampleGroups, Treatment== "drained")

drained_mrna<-mrna[match(drained$SampleName, colnames(mrna))]
cholestatic_mrna<-mrna[match(cholestatic$SampleName, colnames(mrna))]
controls_mrna<-mrna[match(controls$SampleName, colnames(mrna))]

labels_2<-labels[,2:3] # get treatment and mrna id
data_long<-stack(mrna[,2:29],colnames(mrna[,2:29])) # data in long format
merged_data<- merge(data_long, labels_2$V2) # merge dataset. Can take a few min. 
merged_data<- subset(merged_data, select=c("values","y")) # then select only the colums of the expression and the treatment. 
# Can take a few min. 
control_vs_chole<- pairwise.t.test(merged_data$values,merged_data$y, p.adjust.method = "BH") # This doesn't give good results
