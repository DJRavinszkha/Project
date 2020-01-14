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

BiocManager::install("qvalue")

library(qvalue)
library(limma)
library(edgeR)
#====================================#
## Differential Expression Analysis ##
#====================================#

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

#=============================================#
# Now we change the format for the miRNA data #   #Need to be adapted, omitmirna
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

#=========================================#
# 4. PAIRWISE T TEST: Separation on groups and controls. (Still work in progress)
#==========================================#

controls <- subset(sampleGroups, Treatment== "control")
cholestatic<- subset(sampleGroups, Treatment== "cholestatic")
drained<- subset(sampleGroups, Treatment== "drained")

drained_mrna<-mrna[match(drained$SampleName, colnames(mrna))]
cholestatic_mrna<-mrna[match(cholestatic$SampleName, colnames(mrna))]
controls_mrna<-mrna[match(controls$SampleName, colnames(mrna))]

control_vs_chole<- pairwise.t.test(controls_mrna,cholestatic_mrna, p.adjust.method = "BH") #ask the group about how to make the comparison between the different groups










