#=============================================================================#
# Project Period, Liver cholestasis data analysis         												  
#	
# Version: 1.0   															  
# Date: 9-1-2020											             	  
# Authors: Ariadna Fosch i Muntané, ID: I6215203, Maastricht University       #
#          Jip de Kok, ID: I6119367 , Maastricht University                   #
#          Ravin Schmidl, ID: I6125866, Maastricht University                 #
#          Stefan Meier, ID: I6114194 , Maastricht University                 #
# Note: adapted from "Introduction to anamiR - 2019-05-09"
# Ref: Wang T, Lu T, Lee C, Lai L, Tsai M, Chuang E (2015).
#     "anamiR-an integrated analysis package of miRNA and mRNA expression." Manuscript submitted for publication.
#     https://bioconductor.org/packages/release/bioc/vignettes/anamiR/inst/doc/IntroductionToanamiR.html#phenotype-data
#=============================================================================#

#===================#
## Data Formatting ##
#===================#


library(limma)
library(rstudioapi)

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
  
  labels$V1 <- colnames(mirna[2:29])
  
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
  labels$V8[labels$V7 == "254606411109"] <- 3 #set id 3 
  labels$V8[labels$V7 == "254606410405"] <- 4 #set id 4 
  labels$V8[labels$V7 == "254606410413"] <- 5 #set id 5 
  
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
    }
  }
  colnames(mirna) <- columns
  
  #mrna
  rownames(mrna) <- mrna[,1] # Set entrezid as index for mrna; can use key to find gene symbols based on index
  mrna <- mrna[,2:29] #Set mrna to only contain expression values
  
  #mirna
  rownames(mirna) <- mirna[,1] # Set mirna names as index for mirna
  mirna <- mirna[,2:29] #set mirna to only contain expression values
  
  
  #=== Here we initialise the batch numbers ===
  
  #...Need to add code
  
  #=============================#
  # Create phenotype data       #     # Do we need this?
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
  #=============================#
  # Summarised experiment class #   #Do we need this?
  #=============================#
  
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
  

  
  
  return(mrna, mirna, labels, pheno.mrna, pheno.mirna, key)
}

