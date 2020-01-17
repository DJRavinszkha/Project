#=============================================================================#
# Project Period, Liver cholestasis data analysis         							      #					  
#	Expression Analysis                                                         #
# Version: 1.0   															                                #
# Date: 9-1-2020											             	                          #
# Authors: Ariadna Fosch i Muntan??, ID: I6215203, Maastricht University       #
#          Jip de Kok, ID: I6119367 , Maastricht University                   #
#          Ravin Schmidl, ID: I6125866, Maastricht University                 #
#          Stefan Meier, ID: I6114194 , Maastricht University                 #
#=============================================================================#

#=========================================#
##         Install libraries             ##
#=========================================#
library(limma)
library(rstudioapi)
BiocManager::install("qvalue")

library(qvalue)
library(limma)
#====================================#
## Differential Expression Analysis ##
#====================================#

Data <- format() #Remember to run dataFormatting.R first before initialising data
mrna <- data.frame(Data[1]) #mRNA expression data (contains entrez ID as index)
mirna <- data.frame(Data[2]) #miRNA expression data (contains miRNA name as index)
labels <- data.frame(Data[3]) #batch and treatment id/labels for samples
key <- data.frame(Data[4]) #entrezID to genesymbol key


##Link from the PCA.R to here to identify whether we need batch corrections
## Exploratroy PCA plots

## Initilaise Batch effect correction before model
## Set criteria for the q-values


#=======================================================================#
############################################################
##############################################################
##############################################################
#=======================================================================#
# STEFAN'S CODE : FOR TESTING PURPOSES

sampleGroups <- data.frame(treatment = labels$treatment, treatment.id = labels$treatment.id, sampleName = labels$sample.name)

# Get batch order mrna (copied from pca)
mrna.batches <- data.frame(batch = labels$mRNA.batch, batch.id = labels$mRNA.batch.id, file = labels$mRNA.file)

# Get batch order mirna (copied from pca)
mirna.batches <- data.frame(batch = labels$miRNA.batch, batch.id = labels$miRNA.batch.id, file = labels$mRNA.file)

# get mirna treatment order
mirna.treatmentOrder <- matrix(nrow = (ncol(mirna)), ncol = 1)
colnames(mirna.treatmentOrder) <- "sampleName"
mirna.treatmentOrder[,1] <- colnames(mirna)
mirna.treatmentOrder <- merge(mirna.treatmentOrder, sampleGroups, by = "sampleName", sort = FALSE)
mirna.treatmentOrder <- mirna.treatmentOrder[order(mirna.treatmentOrder$treatment.id),] #order based on treatment ID
#========#

#mrna.corrected <- removeBatchEffect(mrna, factor(labels$mRNA.batch.id), design=design1) Move to PCA; only useful for visualisation
#When we creqte this function, design be one of the outputs to be called in PCA.R

#Initialise metadata (age, sex and bilirubin lvels). Move to dataFormatting.R. 
#Refer to whatsapp files to download correct file, 
mrna.meta <- read.delim("../Data/metadata_mrna.csv", sep=',',header = TRUE, colClasses = 'character')
conf <- data.frame(sample.name = mrna.meta$Sample[1:28], sex = mrna.meta$gender[1:28], age = mrna.meta$age[1:28], bili.tot = mrna.meta$Bili..tot.[1:28])

# Sex adjustment to get factor for desing matrix
conf$sex.id[conf$sex == "M"] <- 1 #change male to 1
conf$sex.id[conf$sex == "F"] <- 2 #change female to 2
sex<- factor(conf$sex.id) #Create a factor of sex, as confounder

# Age for confounder. 
#age<- as.numeric(conf$age) #DidnÂ´t create factor due to number of unique ages; therefore kept it a vector
# WARNING: DIFFERENT ORDER THAN age mirna
age<- c(73.5, 76.4, 67.7, 63, 55.1, 84.5, 24.3, 72.1, 71.7, 70.3, 73.4, 51.5, 52.5, 49.5, 59.2, 57.8, 56.7, 48.1, 66.1, 67.2, 67.6,
        41.5, 74.7, 36.3, 75.4, 74.2, 77.8, 59.9) # change to get the data from conf$sex.id but having trouble due to integer/double type. 

design2 <- model.matrix(~ 0 + factor(labels$treatment.id) + factor(labels$mRNA.batch.id) + sex + age, 
                        contrast.arg=list(state=contrasts(state, contrasts=TRUE), batch=contrasts(batch, contrasts = TRUE)))

colnames(design2) <- c("cholestasis", "drained", "control",
                       "batch2","batch3","batch4","batch5",
                       "sex","age")

# Contrast matrices
conmat_DvC = makeContrasts(drained_v_control = drained - control,
                           levels = design2)
conmat_CHvC = makeContrasts(cholestasis_v_control = cholestasis - control,
                            levels = design2)
conmat_CHvD = makeContrasts(cholestasis_v_drained = cholestasis - drained,
                            levels = design2)

#-----------------------------------------------------------#
# linear model fitting
fit <- lmFit(mrna, design2)

fit_contrast_DvC = eBayes(contrasts.fit(fit, conmat_DvC))

fit_contrast_CHvC = eBayes(contrasts.fit(fit, conmat_CHvC))

fit_contrast_CHvD= eBayes(contrasts.fit(fit, conmat_CHvD))

#--------------------------------------------------------#
# Results for lmfit
result_DvC <- decideTests(fit_contrast_DvC)
summary(result_DvC)

result_CHvC <- decideTests(fit_contrast_CHvC)
summary(result_CHvC)

result_CHvD <- decideTests(fit_contrast_CHvD)
summary(result_CHvD)

#----------------------------------------------------------#
# Find significant DEGs with log fold change higher than (lfc).
top_genes_DvC = topTable(fit_contrast_DvC, lfc = 1.5, p.value=0.05, number = nrow(mrna), adjust = "BH")
top_genes_CHvC = topTable(fit_contrast_CHvC,lfc = 1.5, p.value=0.05, number = nrow(mrna), adjust = "BH")
top_genes_CHvD = topTable(fit_contrast_CHvD,lfc = 1.5, p.value=0.05, number = nrow(mrna), adjust = "BH")

#----------------------------------------------------------#

# Save Entrez ID of DEGs 
names_top_DvC <-rownames(top_genes_DvC)
names_top_CHvC<-rownames(top_genes_CHvC)
names_top_CHvD<-rownames(top_genes_CHvD)

# Find gene symbols 
symb.DvC<- key$Genesymbol[match(names_top_DvC,key$EntrezID)]
symb.CHvC<- key$Genesymbol[match(names_top_CHvC,key$EntrezID)]
symb.CHvD<- key$Genesymbol[match(names_top_CHvD,key$EntrezID)]

# Create DEG dataframe 
DEG_DvC<- data.frame(Gene.Symb=symb.DvC, logFC= top_genes_DvC$logFC,adj.p= top_genes_DvC$adj.P.Val)
#rownames(DEG_DvC)<- names_top_DvC  Uncomment if top_genes_DvC has length different than 0. 

DEG_CHvC<- data.frame(Gene.Symb=symb.CHvC,logFC= top_genes_CHvC$logFC,adj.p= top_genes_CHvC$adj.P.Val)
rownames(DEG_CHvC)<- names_top_CHvC

DEG_CHvD<- data.frame(Gene.Symb=symb.CHvD,logFC= top_genes_CHvD$logFC,adj.p= top_genes_CHvD$adj.P.Val)
rownames(DEG_CHvD)<- names_top_CHvD


#---------------------------------------------------------------------------------------#
##                                            MICRO RNA 
#---------------------------------------------------------------------------------------#

#mrna.corrected <- removeBatchEffect(mrna, factor(labels$mRNA.batch.id), design=design1) Move to PCA; only useful for visualisation
#When we creqte this function, design be one of the outputs to be called in PCA.R

#Initialise metadata (age, sex and bilirubin lvels). Move to dataFormatting.R. 
#Refer to whatsapp files to download correct file, 

mirna.meta <- read.delim("../Data/metadata_mirna.csv", sep=',',header = TRUE, colClasses = 'character')
conf.mirna <- data.frame(sample.name = mirna.meta$Sample[1:28], sex = mirna.meta$gender[1:28], age = mirna.meta$age[1:28], bili.tot = mirna.meta$Bili..tot.[1:28])

# Sex as confounder
conf.mirna$sex.id[conf$sex == "M"] <- 1 #change male to 1
conf.mirna$sex.id[conf$sex == "F"] <- 2 #change female to 2
sex.mirna<- factor(conf.mirna$sex.id) #Create a factor of sex, as confounder

# Age as confounder
#age<- as.numeric(conf$age) #DidnÂ´t create factor due to number of unique ages; therefore kept it a vector
age.mirna<- c(73.5, 63.0, 72.1, 73.4, 49.5, 57.8, 66.1, 41.5, 75.4, 76.4, 67.7, 84.5, 71.7, 51.5, 59.2, 56.7,
        67.2, 74.7, 77.8, 55.1, 24.3, 70.3, 52.5, 48.1, 67.6, 36.3, 74.2, 59.9)
# Warning: Don't use on design2 factor(labels$treatment.id) it is sorted on the mrna order not miRNA. 

# Create desing matrix mirna
design2 <- model.matrix(~ 0 + factor(mirna.treatmentOrder$treatment.id) + factor(labels$miRNA.batch.id) + sex.mirna + age.mirna, 
                        contrast.arg=list(state=contrasts(state, contrasts=TRUE), batch=contrasts(batch, contrasts = TRUE)))

colnames(design2) <- c("cholestasis", "drained", "control",
                       "batch2","batch3","batch4","batch5",
                       "sex","age")

# Create contrast matrices mirna 
mirna.conmat_DvC = makeContrasts(drained_v_control = drained - control,
                           levels = design2)
mirna.conmat_CHvC = makeContrasts(cholestasis_v_control = cholestasis - control,
                            levels = design2)
mirna.conmat_CHvD = makeContrasts(cholestasis_v_drained = cholestasis - drained,
                            levels = design2)

#----------------------------------------------
# Model fitting mirna 
fit.mirna <- lmFit(mirna, design2)

fit_contrast_DvC = eBayes(contrasts.fit(fit.mirna, mirna.conmat_DvC))

fit_contrast_CHvC = eBayes(contrasts.fit(fit.mirna, mirna.conmat_CHvC))

fit_contrast_CHvD= eBayes(contrasts.fit(fit.mirna, mirna.conmat_CHvD))

#----------------------------------------------
# Results lmFit mirna

result_DvC = decideTests(fit_contrast_DvC)
summary(result_DvC)

result_CHvC = decideTests(fit_contrast_CHvC)
summary(result_CHvC)

result_CHvD = decideTests(fit_contrast_CHvD)
summary(result_CHvD)
#----------------------------------------------
# Find DEGs on miRNA. 
mirna.top_genes_DvC = topTable(fit_contrast_DvC, lfc=1.2, p.value=0.05, number = nrow(mirna), adjust = "BH")
mirna.top_genes_CHvC = topTable(fit_contrast_CHvC,lfc=1.2, p.value=0.05, number = nrow(mirna), adjust = "BH")
mirna.top_genes_CHvD = topTable(fit_contrast_CHvD,lfc=1.2, p.value=0.05, number = nrow(mirna), adjust = "BH")

#----------------------------------------------

# Save Entrez ID of DEGs 
mirna.names_top_DvC <-rownames(mirna.top_genes_DvC)
mirna.names_top_CHvC<-rownames(mirna.top_genes_CHvC)
mirna.names_top_CHvD<-rownames(mirna.top_genes_CHvD)

# Find gene symbols 
mirna.symb.DvC<- key$Genesymbol[match(mirna.names_top_DvC, key$EntrezID)]
mirna.symb.CHvC<- key$Genesymbol[match(mirna.names_top_CHvC, key$EntrezID)]
mirna.symb.CHvD<- key$Genesymbol[match(mirna.names_top_CHvD, key$EntrezID)]

# Create DEG dataframe 
mirna.DEG_DvC<- data.frame(Gene.Symb=mirna.symb.DvC, logFC= mirna.top_genes_DvC$logFC,adj.p= mirna.top_genes_DvC$adj.P.Val)
#rownames(DEG_DvC)<- names_top_DvC  Uncomment if top_genes_DvC has length different than 0. 

mirna.DEG_CHvC<- data.frame(Gene.Symb=mirna.symb.CHvC,logFC= mirna.top_genes_CHvC$logFC, adj.p= mirna.top_genes_CHvC$adj.P.Val)
rownames(DEG_CHvC)<- mirna.names_top_CHvC

mirna.DEG_CHvD<- data.frame(Gene.Symb=mirna.symb.CHvD,logFC= mirna.top_genes_CHvD$logFC, adj.p= mirna.top_genes_CHvD$adj.P.Val)
rownames(DEG_CHvD)<- mirna.names_top_CHvD


#---------------------------------------------------------------------------------------#
##                                       GENE ONTOLOGY
#---------------------------------------------------------------------------------------#

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("topGO")
l
