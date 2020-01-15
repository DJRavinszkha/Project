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
library(pcaMethods)
library(gplots)
library(ggplot2)

#=========================================#
##           Initialise Data             ##
#=========================================#
Data <- format() #Remember to run dataFormatting.R first before initialising data
mrna <- data.frame(Data[1]) #mRNA expression data (contains entrez ID as index)
mirna <- data.frame(Data[2]) #miRNA expression data (contains miRNA name as index)
labels <- data.frame(Data[3]) #batch and treatment id/labels for samples
key <- data.frame(Data[4]) #entrezID to genesymbol key

#===========================================================================#
## 1. Get list of how treatments and batches are ordered on mRNA and miRNA ##
#===========================================================================#
        # Notes: 1. treat_order_mrna changed to mrna.treatmentOrder
        #        2. treat_order_mirna changed to mirna.treatmentOrder

#Initialise Sample Groups
sampleGroups <- data.frame(treatment = labels$treatment, treatment.id = labels$treatment.id, sampleName = labels$sample.name)

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

#=========================#
##    2. PCA Function    ## 
#=========================#


PCA <- function(data, 
                data.treatmentOrder, 
                treat, 
                data.batch, 
                title1, 
                title2){
  #Usage:
     #   PCA(mirna, 
     #       mirna.treatmentOrder$treatment, 
     #       0,
     #       mirna.batches$batch, 
     #       "miRNA pca coloured by Treatment", 
     #       "miRNA pca coloured by Batch"
     #   )
  pcaRes_data <- pca(t(data),nPcs = 10)  # perform PCA
  data.PCA <- data.frame(c(pcaRes_data@scores[,1]),
                         pcaRes_data@scores[,2],
                         pcaRes_data@scores[,3],
                         pcaRes_data@scores[,4],
                         pcaRes_data@scores[,5],
                         treat = data.treatmentOrder) 
  colnames(data.PCA) = c("PCA1", "PCA2", "PCA3","PCA4", "PCA5","treat")
  
  if(treat != 0){ 
    #need this if loop because we either use a value generated in function for x (as in PCA for mRNA and miRNA)
        # - please refer to lines 16 and 58 in PCA_longversin.R
    # or we prescribe it a value (as in PCA for corrected mRNA)
        # - pleasw refer to line 101 in PCA_longversin.R
    Treatment = treat
  } else{
    Treatment = data.PCA$treat
  }
  ##==== Colouring by treatment ====##
  plot1 <- ggplot(data.PCA, aes(x = PCA1, y = PCA2)) +
              geom_point(aes(colour = Treatment)) +
              scale_colour_manual(values = c("#04179b", "#da9e00", "#198c19"),
                        aesthetics = "fill") + theme_light() + 
              ggtitle(title1)
  
  ##==== Colouring by batch ====##
  plot2 <- ggplot(data.PCA, aes(x = PCA1, y = PCA2)) +
              geom_point(aes(colour = data.batch)) +
              scale_colour_manual(values = c("#04179b", "#da9e00", "#198c19"),
                        aesthetics = "fill") + theme_light() + 
              ggtitle(title2)
  
  ## !!!!!!! ## 
  #install.packages("cowplot") ## !!!!!! ## <<<<<------ need to install
  library(cowplot)
  plot_grid(plot1, plot2, labels = "AUTO")
  
  
}

#=========================#
## 3. PCA for mRNA DATA  ## 
#=========================#
# Notes: 1. Batch_number_mrna changed to mrna.batches
#        2. Batch_num changed to batch
#        3. PCA_mrna changed to mrna.PCA
#        4. Treatment changed to treatment
#MRNA for 28 samples
mrna[is.na(mrna)] <- NA # change Nan for NA
PCA(mrna, #data
    mrna.treatmentOrder$treatment, #treatment order
    0, #treatment color
    mrna.batches$batch, #batch color
    "mRNA pca coloured by Treatment", #title1 
    "mRNA pca coloured by Batch"#title2
    ) 

#=========================#
## 4. PCA for miRNA DATA ##
#=========================#
                        # Notes: 1. Batch_number_mirna changed to mirna.batches
                        #        2. Batch_num changed to batch
                        #        3. PCA_mirna changed to mirna.PCA
                        #        4. Treatment changed to treatment

PCA(mirna, #data
    mirna.treatmentOrder$treatment, #treatment order
    0, #treatment color
    mirna.batches$batch, #batch color
    "miRNA pca coloured by Treatment",#title1
    "miRNA pca coloured by Batch" #title2
    )

#================#
## 5. Heat maps ##
#================#
# Don't remove only silenced because it takes too much time to run. 

#heatmap.2(as.matrix(mrna), trace = "none", main="mRNA heatmap")
#heatmap.2(as.matrix(na.omit(mirna)), trace = "none", main="miRNA heatmap")


#========================================================================#
## 6. PCA showing batches and treatment simultaneously (symbol & color) ## 
#========================================================================#
        #Don't remove just in case we need it. 
        #Note: Remember to change variables as outlined in parts 2 and 3


# pcaRes2 <- pca(t(mrna),nPcs = 10)  # perform PCA
# PCA_28mrna<- data.frame(c(pcaRes2@scores[,1]),
#                         pcaRes2@scores[,2],
#                         pcaRes2@scores[,3],
#                         pcaRes2@scores[,4],
#                         pcaRes2@scores[,5],
#                         batch= batch_number_mrna) 
# colnames(PCA_28mrna) = c("PCA1", "PCA2", "PCA3","PCA4", "PCA5", "Batches")
# 
# ggplot(PCA_28mrna, aes(x = PCA1, y = PCA2)) +
#   geom_point(aes(colour = PCA_mrna$treat, shape=PCA_28mrna$Batches)) +
#   scale_colour_manual(values = c("#04179b", "#da9e00", "#198c19","#66049b"),
#                       aesthetics = "fill") +
#   theme_light()
# 
# mirna = read.delim('../Data/miRNAexpression.txt', check.names = FALSE)
# batch_number_mirna <- data.frame(Batch_num= substr(unlist(colnames(mirna[,2:length(mirna)])), 12, 23)) #Only keep the batch number 
# name_batches_mirna<- unique(batch_number_mirna)
# 
# pcaRes <- pca(t(mirna[,2:29]),nPcs = 10)  # perform PCA
# PCA_28mirna<- data.frame(c(pcaRes@scores[,1]),
#                         pcaRes@scores[,2],
#                         pcaRes@scores[,3],
#                         pcaRes@scores[,4],
#                         pcaRes@scores[,5],
#                         batch= batch_number_mirna) 
# colnames(PCA_28mirna) = c("PCA1", "PCA2", "PCA3","PCA4", "PCA5", "Batches")
# 
# ggplot(PCA_28mirna, aes(x = PCA1, y = PCA2)) +
#   geom_point(aes(colour = PCA_28mirna$Batches)) +
#   scale_colour_manual(values = c("#04179b", "#da9e00", "#198c19","#66049b"),
#                       aesthetics = "fill") +
#   theme_light()

#=============================#
## 7. ANOVA (mrna and mirna) ##
#=============================#

# Anova testing mRNA
mrna.means <- data.frame(means= rowMeans(t(mrna),na.rm = TRUE))
mrna.Anova <- aov(mrna.means[,1]~mrna.batches[,1])
TukeyHSD(mrna.Anova, ordered = FALSE, conf.level = 0.95)

# Anova Testing miRna
mirna.means <- data.frame(means= rowMeans(t(mirna),na.rm = TRUE))
mirna.Anova<- aov(mirna.means[,1]~mirna.batches[,1])
TukeyHSD(mirna.Anova, ordered = FALSE, conf.level = 0.95)
#The anova shows that the mirna doesn't need correction as there is no significant batch variance

#==============================#
## 8. Batch effect correction ##
#==============================#
# Perform correction
##====  mRNA corrected  ====##
mrna.corrected <- removeBatchEffect(mrna, batch = factor(mrna.batches$batch.id)) 

##====  miRNA corrected  ====##
mirna.corrected <- removeBatchEffect(mirna, batch = factor(mirna.batches$batch.id)) 

# Correction boxplot
boxplot(as.data.frame(mrna),main="Original")
boxplot(as.data.frame(mrna.corrected),main="Batch corrected")

# remove batch effects mirna: NOT NECESSARY
# mirna_corrected <- removeBatchEffect(mirna2[, 2:29], batch = factor(batch_number_mirna$id))
# boxplot(as.data.frame(mirna2[,2:29]),main="Original miRNA")
# boxplot(as.data.frame(mirna_corrected),main="Batches miRNA corrected")

#=============================#
## 9. PCA for corrected data ##
#=============================#

##====  mRNA corrected  ====##
PCA(mrna.corrected, #data
    mrna.batches$batch, #batch order
    mrna.treatmentOrder$treatment, #treatment color
    mrna.batches$batch, #batch color
    "mRNA corrected pca coloured by Treatment",  #title1
    "mRNA corrected pca coloured by Batch" #title2
)
##====  miRNA corrected  ====##
PCA(mirna.corrected, #data
    mirna.batches$batch, #batch order
    mirna.treatmentOrder$treatment, #treatment color
    mirna.batches$batch, #batch color
    "miRNA corrected pca coloured by Treatment",  #title1
    "miRNA corrected pca coloured by Batch" #title2
)

# Interbatch variability was very high before the coorection. This didn't allow us to see the 
# Intrabatch variability. After the correction all batches are more similar and we can see the
# intrabatch variability. 