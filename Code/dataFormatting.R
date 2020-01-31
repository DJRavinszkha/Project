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

#===================#
## Data Formatting ##
#===================#

# Get data in correct format
current_path <- getActiveDocumentContext()$path 
setwd(dirname(current_path ))
mrna = read.delim('../Data/GeneExpressionNormalized.txt', check.names = FALSE) #Load expression data
mirna = read.delim('../Data/miRNAexpression.txt', check.names = FALSE)        #Load miRNA expression data

#Initialise a dataframe called labels, containing all:
# file/batch names and sample names 
# batch numbers and an integer id attached to it
# treatments and an integer id attached to it
labels = read.delim("../Data/colNames.csv", sep = ',', header = FALSE, colClasses = 'character')

labels$V4 <- colnames(mrna[3:30])
key <- mrna[,1:2]                   # Key for maintaining gene symbol and entrez gene ID

labels$V1 <- colnames(mirna[2:29])  #maintaining names for mirna filenames

#Here we initialise batch number id's of mRNA's
labels$V5 <- substr(unlist(labels$V4), 8, 13) 
labels$V6[labels$V5 == "410978"] <- 1 #set id 1 to mRNA batch
labels$V6[labels$V5 == "410979"] <- 2 #set id 2
labels$V6[labels$V5 == "410980"] <- 3 #set id 3
labels$V6[labels$V5 == "412287"] <- 4 #set id 4
labels$V6[labels$V5 == "422569"] <- 5 #set id 5

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
labels$treatment.id[labels$treatment == "drained"] <- 2     #set id 2 to drained treatment
labels$treatment.id[labels$treatment == "control"] <- 3     #set id 3 to control treatment

# Column names are different for miRNA and mRNA datasets, but are conserved with sample names in SampleGroups.xlsx
#IE FGS_01 is the name of sample 1; 
#where miRNA sample is US10063773_254606410403_S01_miRNA_107_Sep09_1_1 - 2;
#and mRNA sample is FGS_01_410978_1_1

#=== Here we change the columns names of mrna (mRNA samples) to the sample names ===
mrna = subset(mrna, select = -2)                     # Remove gene symbol column
columns = colnames(mrna)
columns[2:29] <- substr(unlist(columns[2:29]), 1, 6) # Only keep first 6 characters of column names
colnames(mrna) <- columns                            # Assign this to the actual column names

#mrna
rownames(mrna) <- mrna[,1]    # Set entrezid as index for mrna; can use key to find gene symbols based on index
mrna <- mrna[,2:29]           # Set mrna to only contain expression values

#mirna
rownames(mirna) <- mirna[,1]  # Set mirna names as index for mirna
mirna <- mirna[,2:29]         # Set mirna to only contain expression values
colnames(mirna) <- labels$sample.name[order(labels$treatment.id)] ##REQUIRED!!!! in order to preserve order of miRNA sample naming

#= Sample Groups =#
sampleGroups <- data.frame(treatment = labels$treatment, treatment.id = labels$treatment.id, sampleName = labels$sample.name)

# Get batch order mrna (copied from pca)
mrna.batches <- data.frame(batch = labels$mRNA.batch, batch.id = labels$mRNA.batch.id, file = labels$mRNA.file)

# Get treatment order mrna
mrna.treatmentOrder <- matrix(nrow = (ncol(mrna)), ncol = 1)
colnames(mrna.treatmentOrder) <- "sampleName"
mrna.treatmentOrder[,1] <- colnames(mrna)
mrna.treatmentOrder <- merge(mrna.treatmentOrder, sampleGroups, by = "sampleName", sort = FALSE)

# Get treatment order mirna
mirna.treatmentOrder <- matrix(nrow = (ncol(mirna)), ncol = 1)
colnames(mirna.treatmentOrder) <- "sampleName"
mirna.treatmentOrder[,1] <- colnames(mirna)
mirna.treatmentOrder <- merge(mirna.treatmentOrder, sampleGroups, by = "sampleName", sort = FALSE)
mirna.treatmentOrder <- mirna.treatmentOrder[order(mirna.treatmentOrder$treatment.id),] #order based on treatment ID

# Get batch order mrna
mrna.batches <- data.frame(batch = labels$mRNA.batch, batch.id = labels$mRNA.batch.id, file = labels$mRNA.file)

# Get batch order mirna
mirna.batches <- data.frame(batch = labels$miRNA.batch, batch.id = labels$miRNA.batch.id, file = labels$mRNA.file, sort= FALSE)