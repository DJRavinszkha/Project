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

#====================================#
##            Functions             ##
#====================================#

#====================================#
##  Function to initialise design   ##
#====================================#
design <- function(object, contrasts.arg, batch, names){
  d <- model.matrix(object, 
                    contrast.arg=contrasts.arg, 
                    batch=batch)
  colnames(d) <- names
  
  return(d)
}

#====================================#
##   Function to initialise model   ##
#====================================#
fit.eb <- function(data, design, cont){
  fit <- lmFit(data, design)          
  fit.con <- contrasts.fit(fit, cont) 
  fit.eb <- eBayes(fit.con) 
  result <- decideTests(fit.eb)
  print("Summary")
  print(summary(result))
  return(fit.eb)
}

#=========================================================#
##Function to identify top differentially expressed genes##
#=========================================================#
top.expressed <- function(fit.eb, p.value, number, adjust){
  top_expressed <- topTable(fit = fit.eb, p.value = p.value, number = number, adjust = adjust) 
  return(top_expressed)
}

#==========================#
##Function to plot volcano##
#==========================# 
volcano <- function(fit.eb){  
  for (i in 1:ncol(fit.eb)){
    volcanoplot(fit.eb[,i], main= colnames(fit.eb)[i], col=ifelse(fit.eb[,i]$p.value > 0.05,"red","black"))
    abline(-log10(0.05),0)
    abline(v=log2(2))
    abline(v=-log2(2))
  }
}


#====================================#
## Differential Expression Analysis ##
#====================================#

#============================================#
##                   mRNA                   ##
#============================================#

init.mRNA <- function(){
  #mrna.corrected <- removeBatchEffect(mrna, factor(labels$mRNA.batch.id), design=design1) Moved to PCA; only used for visualisation
  
  #Initialise metadata (age, sex and bilirubin lvels). Move to dataFormatting.R.
  mrna.meta <- read.delim("../Data/metadata_mrna.csv", sep=',',header = TRUE, colClasses = 'character')
  conf.mrna <- data.frame(sample.name = mrna.meta$Sample[1:28], sex = mrna.meta$gender[1:28], age = mrna.meta$age[1:28], bili.tot = mrna.meta$Bili..tot.[1:28])
  
  #= Treatment Order =#
  mrna.treatmentOrder <- matrix(nrow = (ncol(mrna)), ncol = 1)
  colnames(mrna.treatmentOrder) <- "sampleName"
  mrna.treatmentOrder[,1] <- colnames(mrna)
  mrna.treatmentOrder <- merge(mrna.treatmentOrder, sampleGroups, by = "sampleName", sort = FALSE)
  
  ##==== Sex as Confounder ====##
  conf.mrna$sex.id[conf.mrna$sex == "M"] <- 1 #change male to 1
  conf.mrna$sex.id[conf.mrna$sex == "F"] <- 2 #change female to 2
  sex<- factor(conf.mrna$sex.id) #Create a factor of sex, as confounder
  
  ##==== Age as Confounder ====##
  age<- c(73.5, 76.4, 67.7, 63, 55.1, 84.5, 24.3, 72.1, 71.7, 70.3, 73.4, 51.5, 52.5, 49.5, 59.2, 57.8, 56.7, 48.1, 66.1, 67.2, 67.6,
          41.5, 74.7, 36.3, 75.4, 74.2, 77.8, 59.9) # change to get the data from conf$sex.id but having trouble due to integer/double type. 
  
  ##==== mRNA Design Matrix ====##
  colnames.mrnaDesign <- c("cholestasis", "drained", "control",
                           "batch2","batch3","batch4", "batch5",
                           "sex","age")
  mrna.design <- design(~ 0 + factor(labels$treatment.id) + factor(labels$mRNA.batch.id) + sex + age,
                        list(state=contrasts(state, contrasts=TRUE)),
                        contrasts(batch, contrasts = TRUE),
                        colnames.mrnaDesign)
  
  return(mrna.design)
}

##===== Drained Vs Control (DVC) =====##
mRNA_DVC <- function(save = TRUE){
  mrna.design <- init.mRNA()
  #= Contrasts =#
  contrasts.mrna.DVC = makeContrasts(drained_v_control = drained - control,
                                     levels = mrna.design)
  #= Linear Model =#
  fit.mrna.DVC <- fit.eb(mrna, mrna.design, contrasts.mrna.DVC)
  
  #= Top Expressed =#
  mrna.expressed.DVC = top.expressed(fit.mrna.DVC, 
                                     p.value=0.05, 
                                     number = nrow(mrna), 
                                     adjust = "BH")
  
  #= Identify Gene Symbols =#
  namesTop.mrna.DVC <-rownames(mrna.expressed.DVC)
  symb.mrna.DVC<- key$Genesymbol[match(namesTop.mrna.DVC,key$EntrezID)]
  
  #= Create Dataframe =#
  mrna.DEG.DVC<- data.frame(Gene.Symb=symb.mrna.DVC,
                            logFC= mrna.expressed.DVC$logFC,
                            adj.p= mrna.expressed.DVC$adj.P.Val)
  rownames(mrna.DEG.DVC)<- namesTop.mrna.DVC
  
  #= Volcano =#
  if(save){
    pdf("../Figures/Volcano_mRNA_DVC.pdf")
    volcano(fit.mrna.DVC)
    dev.off()
  } else {
    volcano(fit.mrna.DVC)
  }
  
  return(list(mrna.DEG.DVC, namesTop.mrna.DVC))
}

##===== Cholestatic Vs Control (CHVC) =====##

mRNA_CHVC <- function(save = TRUE){
  mrna.design <- init.mRNA()
  #= Contrasts =#
  contrasts.mrna.CHVC = makeContrasts(cholestasis_v_control = cholestasis - control,
                                      levels = mrna.design)
  
  #= Fit =#
  fit.mrna.CHVC <- fit.eb(mrna, mrna.design, contrasts.mrna.CHVC)
  
  #= Top Expressed =#
  mrna.expressed.CHVC = top.expressed(fit.mrna.CHVC, 
                                      p.value=0.05, 
                                      number = nrow(mrna), 
                                      adjust = "BH")
  
  #= Identify Gene Symbols =#
  namesTop.mrna.CHVC<-rownames(mrna.expressed.CHVC)
  symb.mrna.CHVC<- key$Genesymbol[match(namesTop.mrna.CHVC,key$EntrezID)]
  
  #= Create Dataframe =#
  mrna.DEG.CHVC <- data.frame(Gene.Symb=symb.mrna.CHVC, logFC= mrna.expressed.CHVC$logFC, adj.p= mrna.expressed.CHVC$adj.P.Val)
  rownames(mrna.DEG.CHVC)<- namesTop.mrna.CHVC
  
  #= Volcano =#
  if(save){
    pdf("../Figures/Volcano_mRNA_CHVC.pdf")
    volcano(fit.mrna.CHVC)
    dev.off()
  } else {
    volcano(fit.mrna.CHVC)
  }
  
  return(list(mrna.DEG.CHVC, namesTop.mrna.CHVC))
}
##===== Cholestatic Vs Drained (CHVD) =====##

mRNA_CHVD <- function(save = TRUE){
  mrna.design <- init.mRNA()
  #= Contrasts =#
  contrasts.mrna.CHVD = makeContrasts(cholestasis_v_drained = cholestasis - drained,
                                      levels = mrna.design)
  
  #= Fit =#
  fit.mrna.CHVD <- fit.eb(mrna, mrna.design, contrasts.mrna.CHVD)
  
  #= Top Expressed =#
  mrna.expressed.CHVD = top.expressed(fit.mrna.CHVD,
                                      p.value=0.05, 
                                      number = nrow(mrna), 
                                      adjust = "BH")
  
  #= Identify Gene Symbol =#
  namesTop.mrna.CHVD<-rownames(mrna.expressed.CHVD)
  symb.mrna.CHVD<- key$Genesymbol[match(namesTop.mrna.CHVD,key$EntrezID)]
  
  #= Create Dataframe =#
  mrna.DEG.CHVD<- data.frame(Gene.Symb=symb.mrna.CHVD,
                             logFC= mrna.expressed.CHVD$logFC,
                             adj.p= mrna.expressed.CHVD$adj.P.Val)
  rownames(mrna.DEG.CHVD)<- namesTop.mrna.CHVD
  
  #= Volcano =#
  if(save){
    pdf("../Figures/Volcano_mRNA_CHVD.pdf")
    volcano(fit.mrna.CHVD)
    dev.off()
  } else {
    volcano(fit.mrna.CHVD)
  }
  
  
  return(list(mrna.DEG.CHVD, namesTop.mrna.CHVD))
}

#============================================#
##                  miRNA                   ##
#============================================#
# Get batch order mirna (copied from pca)
init.miRNA <- function(){
  mirna.batches <- data.frame(batch = labels$miRNA.batch, batch.id = labels$miRNA.batch.id, file = labels$mRNA.file)
  
  #==== Treatment Order ====#
  mirna.treatmentOrder <- matrix(nrow = (ncol(mirna)), ncol = 1)
  colnames(mirna.treatmentOrder) <- "sampleName"
  mirna.treatmentOrder[,1] <- colnames(mirna)
  mirna.treatmentOrder <- merge(mirna.treatmentOrder, sampleGroups, by = "sampleName", sort = FALSE)
  mirna.treatmentOrder <- mirna.treatmentOrder[order(mirna.treatmentOrder$treatment.id),] #order based on treatment ID
  
  #==== Meta Data ====#
  mirna.meta <- read.delim("../Data/metadata_mrna.csv", sep=',',header = TRUE, colClasses = 'character')
  conf.mirna <- data.frame(sample.name = mirna.meta$Sample[1:28], sex = mirna.meta$gender[1:28], age = mirna.meta$age[1:28], bili.tot = mirna.meta$Bili..tot.[1:28])
  
  ##==== Sex as Confounder ====##
  conf.mirna$sex.id[conf.mirna$sex == "M"] <- 1 #change male to 1
  conf.mirna$sex.id[conf.mirna$sex == "F"] <- 2 #change female to 2
  sex.mirna<- factor(conf.mirna$sex.id) #Create a factor of sex, as confounder
  
  ##==== Age as Confounder ====##
  age.mirna<- c(73.5, 63.0, 72.1, 73.4, 49.5, 57.8, 66.1, 41.5, 75.4, 76.4, 67.7, 84.5, 71.7, 51.5, 59.2, 56.7,
                67.2, 74.7, 77.8, 55.1, 24.3, 70.3, 52.5, 48.1, 67.6, 36.3, 74.2, 59.9)
  # Warning: Don't use on design2 factor(labels$treatment.id) it is sorted on the mrna order not miRNA. 
  
  ##==== Design Matrix ====##
  colnames.mirnaDesign <- c("cholestasis", "drained", "control",
                            "batch2","batch3","batch4","batch5",
                            "sex","age")
  mirna.design <- design(~ 0 + factor(mirna.treatmentOrder$treatment.id) + factor(labels$miRNA.batch.id) + sex.mirna + age.mirna,
                         list(state=contrasts(state, contrasts=TRUE)),
                         contrasts(batch, contrasts = TRUE),
                         colnames.mirnaDesign)
  
  return(mirna.design)
}
##===== Drained Vs Control (DVC) =====##

miRNA_DVC <- function(save = TRUE){
  mirna.design <- init.miRNA()
  #= Contrasts =#
  contrasts.mirna.DVC = makeContrasts(drained_v_control = drained - control,
                                      levels = mirna.design)
  #= Fit =#
  fit.mirna.DVC <- fit.eb(mirna, mirna.design, contrasts.mirna.DVC)
  
  #= Top Expressed =#
  mirna.expressed.DVC = top.expressed(fit.mirna.DVC, 
                                      p.value=0.05, 
                                      number = nrow(mirna), 
                                      adjust = "BH")
  #= Identify Gene Symbol =#
  mirna.namesTop.DVC <-rownames(mirna.expressed.DVC)
  mirna.symb.DVC<- key$Genesymbol[match(mirna.namesTop.DVC, key$EntrezID)]
  
  #= Create Dataframe =#
  mirna.DEG.DVC<- data.frame(Gene.Symb=mirna.symb.DVC,
                             logFC= mirna.expressed.DVC$logFC,
                             adj.p= mirna.expressed.DVC$adj.P.Val)
  rownames(mirna.DEG.DVC)<- mirna.namesTop.DVC 
  
  #= Volcano =#
  if(save){
    pdf("../Figures/Volcano_miRNA_DVC.pdf")
    volcano(fit.mirna.DVC)
    dev.off()
  } else {
    volcano(fit.mirna.DVC)
  }
  
  return(list(mirna.DEG.DVC, mirna.namesTop.DVC))
}

##===== Cholestatic Vs Control (CHVC) =====##
miRNA_CHVC <- function(save = TRUE){
  mirna.design <- init.miRNA()
  #= Contrasts =#
  contrasts.mirna.CHVC = makeContrasts(cholestasis_v_control = cholestasis - control,
                                       levels = mirna.design)
  #= Fit =#
  fit.mirna.CHVC <- fit.eb(mirna, mirna.design, contrasts.mirna.CHVC)
  
  #= Top Expressed =#
  mirna.expressed.CHVC <- top.expressed(fit.mirna.CHVC, 
                                        p.value=0.05, 
                                        number = nrow(mirna), 
                                        adjust = "BH")
  
  #= Identify Gene Symbol =#
  mirna.namesTop.CHVC <- rownames(mirna.expressed.CHVC)
  
  #= Create Dataframe =#
  mirna.DEG.CHVC<- data.frame(mirna.Name = mirna.namesTop.CHVC,
                              logFC= mirna.expressed.CHVC$logFC, 
                              adj.p= mirna.expressed.CHVC$adj.P.Val)
  
  rownames(mirna.DEG.CHVC) <- mirna.namesTop.CHVC
  
  #= Volcano =#
  if(save){
    pdf("../Figures/Volcano_miRNA_CHVC.pdf")
    volcano(fit.mirna.CHVC)
    dev.off()
  } else {
    volcano(fit.mirna.CHVC)
  }
  
  
  return(list(mirna.DEG.CHVC, mirna.namesTop.CHVC))
}

##===== Cholestatic Vs Drained (CHVD) =====##
miRNA_CHVD <- function(save = TRUE){
  mirna.design <- init.miRNA()
  #= Contrasts =#
  contrasts.mirna.CHVD = makeContrasts(cholestasis_v_drained = cholestasis - drained,
                                       levels = mirna.design)
  #= Fit =#
  fit.mirna.CHVD <- fit.eb(mirna, mirna.design, contrasts.mirna.CHVD)
  
  #= Top Expressed =#
  mirna.expressed.CHVD <- top.expressed(fit.mirna.CHVD, 
                                        p.value=0.05, 
                                        number = nrow(mirna), 
                                        adjust = "BH")
  
  #= Identify Gene Symbol =#
  mirna.namesTop.CHVD<-rownames(mirna.expressed.CHVD)
  
  #= Create Dataframe =#
  mirna.DEG.CHVD<- data.frame(mirna.Name = mirna.namesTop.CHVD,
                              logFC= mirna.expressed.CHVD$logFC, 
                              adj.p= mirna.expressed.CHVD$adj.P.Val)
  rownames(mirna.DEG.CHVD)<- mirna.namesTop.CHVD
  
  #= Volcano =#
  if(save){
    pdf("../Figures/Volcano_miRNA_CHVD.pdf")
    volcano(fit.mirna.CHVD)
    dev.off()
  } else {
    volcano(fit.mirna.CHVD)
  }
  
  return(list(mirna.DEG.CHVD, mirna.namesTop.CHVD))
}





