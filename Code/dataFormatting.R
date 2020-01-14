#=============================================================================#
# Project Period, Liver cholestasis data analysis         												  
#	
# Version: 1.0   															  
# Date: 9-1-2020											             	  
<<<<<<< Updated upstream
# Author:  Jip de Kok, Stefan Meier, Ariadna Fosch & Ravin Schmidl 
# Note: adapted from "Introduction to anamiR - 2019-05-09"
# Ref: Wang T, Lu T, Lee C, Lai L, Tsai M, Chuang E (2015).
#     "anamiR-an integrated analysis package of miRNA and mRNA expression." Manuscript submitted for publication.
#     https://bioconductor.org/packages/release/bioc/vignettes/anamiR/inst/doc/IntroductionToanamiR.html#phenotype-data
=======
# Authors: Ariadna Fosch i Muntané, ID: I6215203, Maastricht University       #
#          Jip de Kok, ID: I6119367 , Maastricht University                   #
#          Ravin Schmidl, ID: I6125866, Maastricht University                 #
#          Stefan Meier, ID: I6114194 , Maastricht University                 #                           
>>>>>>> Stashed changes
#=============================================================================#

#===================#
## Data Formatting ##
#===================#

<<<<<<< Updated upstream
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager", "limma")

BiocManager::install("anamiR")


library(anamiR)
=======
setRepositories(FALSE, 1:9)


library(rstudioapi)
>>>>>>> Stashed changes
library(limma)
library(rstudioapi)

format <- function(){
  # Get data in correct format
  current_path <- getActiveDocumentContext()$path 
  setwd(dirname(current_path ))
  mrna = read.delim('../Data/GeneExpressionNormalized.txt', check.names = FALSE) #Load expression data
  mirna = read.delim('../Data/miRNAexpression.txt', check.names = FALSE) #Load miRNA expression data
  
  
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
  }
  colnames(mirna) <- columns
  
  # Set rownames for mrna
  rownames(mrna) <- mrna[,1]
  mrna <- mrna[,2:29]
  
  # Set rownames for mirna
  rownames(mirna) <- mirna[,1]
  mirna <- mirna[,2:29]
  
  #=== Here we initialise the sample grouping ===
  sampleGroups <- read.delim("../Data/SampleGroups.csv", sep = ',', header = TRUE, colClasses = 'character')
  sampleGroups <- sampleGroups[,6:7]
  
  index_cholestatic <- sampleGroups == 'cholestatic'
  sampleGroups$id[index_cholestatic[,1]] <- 1
  index_drained <- sampleGroups == 'drained'
  sampleGroups$id[index_drained[,1]] <- 2
  index_control <- sampleGroups == 'control'
  sampleGroups$id[index_control[,1]] <- 3
  
  #=============================#
  # Different gene expression  #
  #=============================#
  
  # First we identify differentially expressed genes for the mRNA's
  design_matrix <- model.matrix(~ 0 + factor(c(1,1,1,1,1,1,1,1,1,
                                               2,2,2,2,2,2,2,2,2,2,
                                               3,3,3,3,3,3,3,3,3)))
  
  colnames(design_matrix) <- c("cholestasis", "drained", "control")
  
  cont_matrix <- makeContrasts (drained_v_control = drained - control,
                                cholestasis_v_control = cholestasis - control,
                                cholestasis_v_drained = cholestasis - drained,
                                levels = design_matrix)
  
  fit <- lmFit(mrna, design_matrix)
  
  fit_contrast <- contrasts.fit(fit, cont_matrix)
  
  fit_contrast <- eBayes(fit_contrast)
  
  results <- decideTests(fit_contrast)
  
  summary(results)
  
  top_genes <- topTable (fit_contrast, number = nrow(mrna), adjust = "BH")
  
  # Subsequently we do the same for the miRNA data
  design_matrix <- model.matrix(~ 0 + factor(c(1,1,1,1,1,1,1,1,1,
                                               2,2,2,2,2,2,2,2,2,2,
                                               3,3,3,3,3,3,3,3,3)))
  
  colnames(design_matrix) <- c("cholestasis", "drained", "control")
  
  cont_matrix <- makeContrasts (drained_v_control = drained - control,
                                cholestasis_v_control = cholestasis - control,
                                cholestasis_v_drained = cholestasis - drained,
                                levels = design_matrix)
  #Here we change NA's with rownmeans SHOULD BE CHANGED LATER!!!
  omitmirna <- mirna
  
  omitmirna[is.na(omitmirna)] <- mean(as.matrix(mirna), na.rm = TRUE)
  
  fit2 <- lmFit(omitmirna, design_matrix)
  
  fit_contrast2 <- contrasts.fit(fit2, cont_matrix)
  
  fit_contrast2 <- eBayes(fit_contrast2)
  
  results2 <- decideTests(fit_contrast2)
  
  summary(results2)
  
  top_genes2 <- topTable (fit_contrast2, number = nrow(mirna), adjust = "BH")
  
  #=============================#
  # Create phenotype data       #
  #=============================#
  pheno.mrna = labels
  pheno.mrna[[1]] <- colnames(mrna)
  colnames(pheno.mrna) <- c("", "Subtype", "ER")
  
  # Set controls
  case <- pheno.mrna[,'Subtype'] == "cholestatic"
  pheno.mrna[case,3] <- 'case'
  
  # Set drained
  case <- pheno.mrna[,'Subtype'] == "drained"
  pheno.mrna[case,3] <- 'case'
  
  # set controls
  case <- pheno.mrna[,'Subtype'] == "control"
  pheno.mrna[case,3] <- 'control'
  
  # Order cases and controls
  pheno.mrna <- pheno.mrna[order(as.character(pheno.mrna$ER)),]
  
  # Set pheno.mirna which is identical to pheno.mrna
  pheno.mirna = pheno.mrna
  
  #=============================#
  # Summarised experiment class #
  #=============================#
  
  # First we change the dataframes into matrices as the miRrna package works with matrices.
  pheno.mrna <- as.matrix(pheno.mrna)
  pheno.mirna <- as.matrix(pheno.mirna)
  mrna = data.matrix(mrna)
  mirna = data.matrix(mirna)
  
  # Fill NA's with mean for the time-being
  mrna[is.na(mrna)] <- mean(mrna, na.rm = TRUE)
  mirna[is.na(mirna)] <- mean(mirna, na.rm = TRUE)
  
  # Remove drained cases from mrna
  index = rep(TRUE, length(colnames(mrna)))
  for(i in 1:length(colnames(mrna))){
    # if current column name is present in a list of all column names of cases:
    if(colnames(mrna)[i] %in% pheno.mrna[pheno.mrna[,2]=="drained",1]){
      index[i] <- FALSE
    }
  }
  mrna <- mrna[,index]
  
  # Remove drained cases from mirna
  index = rep(TRUE, length(colnames(mirna)))
  for(i in 1:length(colnames(mirna))){
    # if current column name is present in a list of all column names of cases:
    if(colnames(mirna)[i] %in% pheno.mirna[pheno.mirna[,2]=="drained",1]){
      index[i] <- FALSE
    }
  }
  mirna <- mirna[,index]
  
  # Remove drained cases from the phenotype data
  pheno.mrna <- pheno.mrna[!(pheno.mrna[,2]=="drained"),]
  pheno.mirna <- pheno.mirna[!(pheno.mirna[,2]=="drained"),]
  
  
  
  #=============================#
  # Change DE gene table format #
  #=============================#
  # First we change the layout for the mRNA data
  mrna_d <- top_genes[,c(2,6,7)]
  mrna_d[["mean_case"]] <- 0
  mrna_d[["mean_control"]] <- 0
  
  # Seperate controls from cases
  index = rep(FALSE, length(colnames(mrna)))
  for(i in 1:length(colnames(mrna))){
    # if current column name is present in a list of all column names of cases:
    if(colnames(mrna)[i] %in% pheno.mrna[pheno.mrna[,3]=="case",1]){
      index[i] <- TRUE
    }
  }
  
  mrna.case <- mrna[,index]
  mrna.control <- mrna[,!index]
  
  # Calculate means
  mean_control <- rowMeans(mrna.control)
  mean_case <- rowMeans(mrna.case)
  
  # Include means in mrna_d
  mrna_d[[4]] <- mean_case
  mrna_d[[5]] <- mean_control
  
  colnames(mrna_d) <- c("log-ratio", "P-Value", "P-adjust", "mean_case", "mean_control")
  
  #=============================================#
  # Now we change the format for the miRNA data #
  #=============================================#
  mirna_d <- top_genes2[,c(2,6,7)]
  mirna_d[["mean_case"]] <- 0
  mirna_d[["mean_control"]] <- 0
  
  # Seperate controls from cases
  index = rep(FALSE, length(colnames(omitmirna)))
  for(i in 1:length(colnames(omitmirna))){
    # if current column name is present in a list of all column names of cases:
    if(colnames(omitmirna)[i] %in% pheno.mirna[pheno.mirna[,3]=="case",1]){
      index[i] <- TRUE
    }
  }
  
  mirna.case <- omitmirna[,index]
  mirna.control <- omitmirna[,!index]
  
  # Calculate means
  mean_control <- rowMeans(mirna.control[,2:length(mirna.control)])
  mean_case <- rowMeans(mirna.case[,2:length(mirna.case)])
  
  # Include means in omitmirna
  mirna_d[[4]] <- mean_case
  mirna_d[[5]] <- mean_control
  
  colnames(mirna_d) <- c("log-ratio", "P-Value", "P-adjust", "mean_case", "mean_control")
  
  # # Set first column of phenotype data as rowname
  # rownames(pheno.mrna) <- pheno.mrna[,1]
  # pheno.mrna <- pheno.mrna[,2:3]
  # pheno.mirna <- pheno.mrna # pheno.mirna is identical to pheno.mrna (for clarity they both exist)
  # 
  
  return(mrna, mirna, labels, pheno.mrna, pheno.mirna, key)
}
<<<<<<< Updated upstream
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


>>>>>>> Stashed changes


