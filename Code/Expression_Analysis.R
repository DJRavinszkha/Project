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


#====================================#
##   Function to initialise model   ##    #Note: stop using limma and convert to edgeR for batch correction
#====================================#

fit.eb <- function(data, design, cont){
  fit <- lmFit(data, design) #            Note: glmFit
  fit.con <- contrasts.fit(fit, cont) #         ignore
  fit.eb <- eBayes(fit.con) #                   ignore
  result <- decideTests(fit.eb) #               ignore
  print("Summary")
  print(summary(result))
  return(fit.eb)
}

#=========================================================#
##Function to identify top differentially expressed genes##
#=========================================================#
top.expressed <- function(fit.eb, p.value, number, adjust){
  top_expressed <- topTable(fit.eb, p.value, number, adjust) 
  return(top_expressed)
}

#==========================#
##Function to plot volcano##  #Note: Can keep the same
#==========================# 
volcano <- function(fit.eb){  
  for (i in 1:ncol(fit.eb)){
    volcanoplot(fit.eb[,i], main= colnames(fit.eb)[i], col=ifelse(fit.eb[,i]$p.value > 0.05,"red","black"))
    abline(-log10(0.05),0)
    abline(v=log2(2))
    abline(v=-log2(2))
  }
}

#================================#
##Function to calculate q values##
#================================#

q.value <- function(fit.eb){
  drained_v_control <- qvalue(fit.eb$p.value[,1])
  cholestasis_v_control <- qvalue(fit.eb$p.value[,2])
  cholestasis_v_drained <- qvalue(fit.eb$p.value[,3])
  
  #Init dataframe with all of the qvalues called q.values
  q.values <- data.frame(cbind(drained_v_control$qvalues, cholestasis_v_control$qvalues, cholestasis_v_drained$qvalues))
  colnames(q.values) <- c("drained_v_control_qvalue", "cholestasis_v_control_qvalue", "cholestasis_v_drained_qvalue")
  return(q.values)
}    

#============================================================#
## Initialise design matrices, colnames and contrast matrix ##
#============================================================#

#======     mRNA      ======#       Need to check design matrices, can use dataframe LABELS for initialisation
print("mRNA Design Matrix")
design_matrix_mrna <- model.matrix(~ 0 + factor(c("NEED", "TO", "CREATE")))

#======     miRNA      ======#
print("miRNA Design Matrix")
design_matrix_mirna <- model.matrix(~ 0 + factor(c("NEED", "TO", "CREATE")))

#======     other      ======#
colnames(design_matrix_mrna) <- c("cholestasis", "drained", "control")

cont_matrix <- makeContrasts (drained_v_control = drained - control,
                              cholestasis_v_control = cholestasis - control,
                              cholestasis_v_drained = cholestasis - drained,
                              levels = design_matrix_mrna)

#============================================#
## mRNA topgenes, volcano plot and q-values ##
#============================================#
print("Linear model of mRNA")
mrna.fit <- fit.eb(mrna, design_matrix_mrna, cont_matrix)
print("Top DE of mRNA")
mrna.topExpressed <- top.expressed(mrna.fit, 0.05, nrow(mrna), "BH")
print("mRNA volcano plot")
mrna.volcano <- volcano(mrna.fit)
print("mRNA q-values")
mrna.q <- q.value(mrna.fit)


#=============================================#
## miRNA topgenes, volcano plot and q-values ##
#=============================================#
print("Linear model of miRNA")
mirna.fit <- fit.eb(mirna, design_matrix_mirna, cont_matrix)
print("Top DE of miRNA")
mirna.topExpressed <- top.expressed(mirna.fit, 0.05, nrow(mirna), "BH")
print("miRNA volcano plot")
mirna.volcano <- volcano(mirna.fit)
print("miRNA q-values")
mirna.q <- q.value(mirna.fit)

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
#mrna_d = topgenes, calculated means for all genes.

#Error above in `[[<-.data.frame`(`*tmp*`, 4, value = c(`1` = 12.6283098666667,  : 
#replacement has 17581 rows, data has 1011
mrna_d[[5]] <- mean_control
#mrna_d = topgenes, calculated means for all genes.


#Error above in `[[<-.data.frame`(`*tmp*`, 5, value = c(`1` = 11.8289007905,  : 
#replacement has 17581 rows, data has 1011

colnames(mrna_d) <- c("log-ratio", "P-Value", "P-adjust", "mean_case", "mean_control")

#===============================================#
## Now we change the format for the miRNA data ##   #Need to be adapted, omitmirna
#===============================================#
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

#=========================================================#
## 4. PAIRWISE T TEST: Separation on groups and controls ##   (Still work in progress)
#=========================================================#

controls <- subset(sampleGroups, Treatment== "control")
cholestatic<- subset(sampleGroups, Treatment== "cholestatic")
drained<- subset(sampleGroups, Treatment== "drained")

drained_mrna<-mrna[match(drained$SampleName, colnames(mrna))]
cholestatic_mrna<-mrna[match(cholestatic$SampleName, colnames(mrna))]
controls_mrna<-mrna[match(controls$SampleName, colnames(mrna))]

control_vs_chole<- pairwise.t.test(controls_mrna,cholestatic_mrna, p.adjust.method = "BH") #ask the group about how to make the comparison between the different groups


#####################################################################

#========#
# design_matrix_mrna <- matrix(nrow = (ncol(mrna) -1), ncol = 1)
# colnames(design_matrix_mrna) <- "SampleName"
# design_matrix_mrna[,1] <- colnames(mrna[,2:length(mrna)])
# design_matrix_mrna <- merge(design_matrix_mrna, sampleGroups, by =  "SampleName", sort = FALSE)
# design_matrix_mrna <- model.matrix(~ 0 + factor(design_matrix_mrna[,3]) + factor(mrna.batches$batch.id))
# colnames(design_matrix_mrna) <- c("cholestasis", "drained", "control")
# 
# cont_matrix <- makeContrasts (drained_v_control = drained - control,
#                               cholestasis_v_control = cholestasis - control,
#                               cholestasis_v_drained = cholestasis - drained,
#                               levels = design_matrix_mrna)
# 
# fit <- lmFit(mrna[, 2:29], design_matrix_mrna)
# 
# fit_contrast <- contrasts.fit(fit, cont_matrix)
# 
# fit_contrast <- eBayes(fit_contrast)
# 
# results <- decideTests(fit_contrast)
# 
# summary(results)

top_genes <- topTable (fit_contrast, p.value = "0.05", number = nrow(mrna), adjust = "BH")

for (i in 1:ncol(fit_contrast)){
  volcanoplot(fit_contrast[,i], main= colnames(fit_contrast)[i], col=ifelse(fit_contrast[,i]$p.value > 0.05,"red","black"))
  abline(-log10(0.05),0)
  abline(v=log2(2))
  abline(v=-log2(2))
}
#=======================================================================#
# STEFAN'S CODE : FOR TESTING PURPOSES

sampleGroups <- data.frame(treatment = labels$treatment, treatment.id = labels$treatment.id, sampleName = labels$sample.name)

# Get batch order mrna (copied from pca)
mrna.batches <- data.frame(batch = labels$mRNA.batch, batch.id = labels$mRNA.batch.id, file = labels$mRNA.file)

# Get batch order mirna (copied from pca)
mirna.batches <- data.frame(batch = labels$miRNA.batch, batch.id = labels$miRNA.batch.id, file = labels$mRNA.file)
#========#

design <- matrix(nrow = (ncol(mrna)), ncol = 1)
colnames(design) <- "sampleName"
design[,1] <- colnames(mrna)
design <- merge(design, sampleGroups, by = "sampleName", sort = FALSE)

design1 <- model.matrix(~ 0 + factor(labels$treatment.id))
#mrna.corrected <- removeBatchEffect(mrna, factor(labels$mRNA.batch.id), design=design1) Move to PCA; only useful for visualisation
#When we creqte this function, design be one of the outputs to be called in PCA.R

#Initialise metadata (age, sex and bilirubin lvels). Move to dataFormatting.R. 
#Refer to whatsapp files to download correct file, 
mrna.meta <- read.delim("../Data/metadata_mrna.csv", sep=',',header = TRUE, colClasses = 'character')
conf <- data.frame(sample.name = mrna.meta$Sample[1:28], sex = mrna.meta$gender[1:28], age = mrna.meta$age[1:28], bili.tot = mrna.meta$Bili..tot.[1:28])

conf$sex.id[conf$sex == "M"] <- 1 #change male to 1
conf$sex.id[conf$sex == "F"] <- 2 #change male to 1

sex<- factor(conf$sex.id) #Create a factor of sex, as confounder
age<- factor(conf$age) #Didn´t create factor due to number of unique ages; therefore kept it a vector

design2 <- model.matrix(~ 0 + factor(labels$treatment.id) + factor(labels$mRNA.batch.id) + sex + age)

colnames(design2) <- c("cholestasis", "drained", "control", "batch2","batch3","batch4","sex","age")

acont_matrix <- makeContrasts (drained_v_control = drained - control,
                              cholestasis_v_control = cholestasis - control,
                              cholestasis_v_drained = cholestasis - drained,
                              levels = design2)

fit <- lmFit(mrna, design2)

fit_contrast <- contrasts.fit(fit, acont_matrix)

fit_contrast <- eBayes(fit_contrast)

results <- decideTests(fit_contrast)

summary(results)

top_genes <- topTable (fit_contrast, p.value=0.05, number = nrow(mrna), adjust = "BH")

for (i in 1:ncol(fit_contrast)){
  volcanoplot(fit_contrast[,i], main= colnames(fit_contrast)[i], col=ifelse(fit_contrast[,i]$p.value > 0.05,"red","black"))
  abline(-log10(0.05),0)
  abline(v=log2(2))
  abline(v=-log2(2))
}


#=======================================================================#
design_matrix_mirna <- matrix(nrow = (ncol(mirna)), ncol = 1)
colnames(design_matrix_mirna) <- "sampleName"
design_matrix_mirna[,1] <- colnames(mirna)
design_matrix_mirna <- merge(design_matrix_mirna, sampleGroups, by = "sampleName", sort = FALSE)
design_matrix_mirna <- model.matrix(~ 0 + factor(design_matrix_mirna[,3]))
colnames(design_matrix_mirna) <- c("cholestasis", "drained", "control")

cont_matrix <- makeContrasts (drained_v_control = drained - control,
                              cholestasis_v_control = cholestasis - control,
                              cholestasis_v_drained = cholestasis - drained,
                              levels = design_matrix_mirna)

#omitmirna <- na.omit(mirna) #Do we need impute or omit or mean? 

fit2 <- lmFit(mirna, design_matrix_mirna)

fit_contrast2 <- contrasts.fit(fit2, cont_matrix)

fit_contrast2 <- eBayes(fit_contrast2)

results2 <- decideTests(fit_contrast2)

summary(results2)

top_genes2 <- topTable (fit_contrast2, number = nrow(mirna), adjust = "BH") 
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

adj_q_mrna_drained_v_control <- data.frame(q_val= mrna_drained_v_control_pvalue$qvalues[which(mrna_drained_v_control_pvalue$qvalues <= 0.05)])
adj_q_mrna_cholestasis_v_control <- data.frame(q_val= mrna_cholestasis_v_control_pvalue$qvalues[which(mrna_cholestasis_v_control_pvalue$qvalues <= 0.05)])
adj_q_mrna_cholestasis_v_drained <- data.frame(q_val= mrna_cholestasis_v_drained_pvalue$qvalues[which(mrna_cholestasis_v_drained_pvalue$qvalues <= 0.05)])

adj_q_mirna_drained_v_control <- data.frame(q_val= mirna_drained_v_control_pvalue$qvalues[which(mirna_drained_v_control_pvalue$qvalues <= 0.05)])
adj_q_mirna_cholestasis_v_control <- data.frame(q_val= mirna_cholestasis_v_control_pvalue$qvalues[which(mirna_cholestasis_v_control_pvalue$qvalues <= 0.05)])
adj_q_mirna_cholestasis_v_drained <- data.frame(q_val= mirna_cholestasis_v_drained_pvalue$qvalues[which(mirna_cholestasis_v_drained_pvalue$qvalues <= 0.05)])

# Differentially expressed genes on top genes. 

deg_mirna_cholestasis_v_control <- match(rownames(adj_q_mirna_cholestasis_v_control),rownames(top_genes2))
deg_mirna_cholestasis_v_drained <- match(rownames(adj_q_mirna_cholestasis_v_drained),rownames(top_genes2))

deg_mrna_drained_v_control <- match(rownames(adj_q_mrna_drained_v_control),rownames(top_genes))
deg_mrna_cholestasis_v_control <- match(rownames(adj_q_mrna_cholestasis_v_control),rownames(top_genes))
deg_mrna_cholestasis_v_drained <- match(rownames(adj_q_mrna_cholestasis_v_drained),rownames(top_genes))

#miRNA diff exp genes
deg_mirna_drained_v_control<- top_genes2[match(rownames(adj_q_mirna_drained_v_control),rownames(top_genes2)),]
q_log_deg_mirna_dvco<- cbind(logFC= deg_mirna_drained_v_control$drained_v_control,adj_q_mirna_drained_v_control)

deg_mirna_cholestasis_v_control<- top_genes2[match(rownames(adj_q_mirna_cholestasis_v_control),rownames(top_genes2)),]
q_log_deg_mirna_chvco<- cbind(logFC= deg_mirna_cholestasis_v_control$cholestasis_v_control,adj_q_mirna_cholestasis_v_control)

deg_mirna_cholestasis_v_drained<- top_genes2[match(rownames(adj_q_mirna_cholestasis_v_drained),rownames(top_genes2)),]
q_log_deg_mirna_chvd<- cbind(logFC= deg_mirna_cholestasis_v_drained$cholestasis_v_drained,adj_q_mirna_cholestasis_v_drained)

#mRNA diff exp genes

deg_mrna_drained_v_control<- top_genes[match(rownames(adj_q_mrna_drained_v_control),rownames(top_genes)),]
q_log_deg_mrna_dvco<- cbind(logFC= deg_mrna_drained_v_control$drained_v_control,adj_q_mrna_drained_v_control)

deg_mrna_cholestasis_v_control<- top_genes[match(rownames(adj_q_mrna_cholestasis_v_control),rownames(top_genes)),]
q_log_deg_mrna_chvco<- cbind(logFC= deg_mrna_cholestasis_v_control$cholestasis_v_control,adj_q_mrna_cholestasis_v_control)

deg_mrna_cholestasis_v_drained<- top_genes[match(rownames(adj_q_mrna_cholestasis_v_drained),rownames(top_genes)),]
q_log_deg_mrna_chvd<- cbind(logFC= deg_mrna_cholestasis_v_drained$cholestasis_v_drained,adj_q_mrna_cholestasis_v_drained)

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
#========#

design <- matrix(nrow = (ncol(mrna)), ncol = 1)
colnames(design) <- "sampleName"
design[,1] <- colnames(mrna)
design <- merge(design, sampleGroups, by = "sampleName", sort = FALSE)

design1 <- model.matrix(~ 0 + factor(labels$treatment.id))
#mrna.corrected <- removeBatchEffect(mrna, factor(labels$mRNA.batch.id), design=design1) Move to PCA; only useful for visualisation
#When we creqte this function, design be one of the outputs to be called in PCA.R

#Initialise metadata (age, sex and bilirubin lvels). Move to dataFormatting.R. 
#Refer to whatsapp files to download correct file, 
mrna.meta <- read.delim("../Data/metadata_mrna.csv", sep=',',header = TRUE, colClasses = 'character')
conf <- data.frame(sample.name = mrna.meta$Sample[1:28], sex = mrna.meta$gender[1:28], age = mrna.meta$age[1:28], bili.tot = mrna.meta$Bili..tot.[1:28])

conf$sex.id[conf$sex == "M"] <- 1 #change male to 1
conf$sex.id[conf$sex == "F"] <- 2 #change female to 2

sex<- factor(conf$sex.id) #Create a factor of sex, as confounder
#age<- as.numeric(conf$age) #DidnÂ´t create factor due to number of unique ages; therefore kept it a vector
age<- c(73.5, 76.4, 67.7, 63, 55.1, 84.5, 24.3, 72.1, 71.7, 70.3, 73.4, 51.5, 52.5, 49.5, 59.2, 57.8, 56.7, 48.1, 66.1, 67.2, 67.6,
        41.5, 74.7, 36.3, 75.4, 74.2, 77.8, 59.9)
design2 <- model.matrix(~ 0 + factor(labels$treatment.id) + factor(labels$mRNA.batch.id) + sex + age, 
                        contrast.arg=list(state=contrasts(state, contrasts=TRUE),batch=contrasts(batch, contrasts = TRUE)))

colnames(design2) <- c("cholestasis", "drained", "control", "batch2","batch3","batch4","sex","age")

acont_matrix <- makeContrasts (drained_v_control = drained - control,
                               cholestasis_v_control = cholestasis - control,
                               cholestasis_v_drained = cholestasis - drained,
                               levels = design2)

conmat_DvC = makeContrasts(drained_v_control = drained - control,
                           levels = design2)
conmat_CHvC = makeContrasts(cholestasis_v_control = cholestasis - control,
                            levels = design2)
conmat_CHvD = makeContrasts(cholestasis_v_drained = cholestasis - drained,
                            levels = design2)


#-----------------------------------------------------------#


fit <- lmFit(mrna, design2)

fit_contrast_DvC = contrasts.fit(fit, conmat_DvC)
fit_contrast_DvC = eBayes(fit_contrast_DvC)

fit_contrast_CHvC = contrasts.fit(fit, conmat_CHvC)
fit_contrast_CHvC = eBayes(fit_contrast_CHvC)

fit_contrast_CHvD= contrasts.fit(fit, conmat_CHvD)
fit_contrast_CHvD = eBayes(fit_contrast_CHvD)

#--------------------------------------------------------#


result_DvC = decideTests(fit_contrast_DvC)
summary(result_DvC)

result_CHvC = decideTests(fit_contrast_CHvC)
summary(result_CHvC)

result_CHvD = decideTests(fit_contrast_CHvD)
summary(result_CHvD)

#-----------------------------------#
# Creating the pairwise comparisons in the contrast matrixs seperately.

conmat_DvC = makeContrasts(drained_v_control = drained - control,
                           levels = design2)
conmat_CHvC = makeContrasts(cholestasis_v_control = cholestasis - control,
                            levels = design2)
conmat_CHvD = makeContrasts(cholestasis_v_drained = cholestasis - drained,
                            levels = design2)


#-----------------------------------------------------------#


fit <- lmFit(mrna, design2)

fit_contrast_DvC = contrasts.fit(fit, conmat_DvC)
fit_contrast_DvC = eBayes(fit_contrast_DvC)

fit_contrast_CHvC = contrasts.fit(fit, conmat_CHvC)
fit_contrast_CHvC = eBayes(fit_contrast_CHvC)

fit_contrast_CHvD= contrasts.fit(fit, conmat_CHvD)
fit_contrast_CHvD = eBayes(fit_contrast_CHvD)

#--------------------------------------------------------#


result_DvC = decideTests(fit_contrast_DvC)
summary(result_DvC)

result_CHvC = decideTests(fit_contrast_CHvC)
summary(result_CHvC)

result_CHvD = decideTests(fit_contrast_CHvD)
summary(result_CHvD)


#----------------------------------------------------------#
top_genes_DvC = topTable(fit_contrast_DvC, p.value=0.05, number = nrow(mrna), adjust = "BH")
top_genes_CHvC = topTable(fit_contrast_CHvC, p.value=0.05, number = nrow(mrna), adjust = "BH")
top_genes_CHvD = topTable(fit_contrast_CHvD, p.value=0.05, number = nrow(mrna), adjust = "BH")

#---------------------------------------------------------------------------------------#