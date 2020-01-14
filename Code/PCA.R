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
mirna.batches <- data.frame(batch = labels$miRNA.batch, batch.id = labels$miRNA.batch.id, file = labels$mRNA.file)

#=======================#
# 2. PCA for mRNA DATA  # 
#=======================#
                        # Notes: 1. Batch_number_mrna changed to mrna.batches
                        #        2. Batch_num changed to batch
                        #        3. PCA_mrna changed to mrna.PCA
                        #        4. Treatment changed to treatment
#MRNA for 28 samples
mrna[is.na(mrna)] <- NA # change Nan for NA
pcaRes_mrna <- pca(t(mrna),nPcs = 10)  # perform PCA
mrna.PCA <- data.frame(c(pcaRes_mrna@scores[,1]),
                       pcaRes_mrna@scores[,2],
                       pcaRes_mrna@scores[,3],
                       pcaRes_mrna@scores[,4],
                       pcaRes_mrna@scores[,5],
                       treat = mrna.treatmentOrder$treatment) 
colnames(mrna.PCA) = c("PCA1", "PCA2", "PCA3","PCA4", "PCA5","treat")

##==== Colouring by treatment ====##
ggplot(mrna.PCA, aes(x = PCA1, y = PCA2)) +
  geom_point(aes(colour = mrna.PCA$treat)) +
  scale_colour_manual(values = c("#04179b", "#da9e00", "#198c19"),
                      aesthetics = "fill") + theme_light() + 
  ggtitle("mRNA pca coloured by Treatment")

##==== Colouring by batch ====##
ggplot(mrna.PCA, aes(x = PCA1, y = PCA2)) +
  geom_point(aes(colour = mrna.batches$batch)) +
  scale_colour_manual(values = c("#04179b", "#da9e00", "#198c19"),
                      aesthetics = "fill") + theme_light() + 
  ggtitle("mRNA pca coloured by Batch")

#MRNA for all microRNA (17581)
# mrna[is.na(mrna)] <- NA # change Nan for NA
# pcaRes_mrna17581_treat <- pca(mrna[,2:29],nPcs = 10)  # perform PCA
# PCA_1758mrna<- data.frame(c(pcaRes_mrna17581_treat@scores[,1]),
#                         pcaRes_mrna17581_treat@scores[,2],
#                         pcaRes_mrna17581_treat@scores[,3],
#                         pcaRes_mrna17581_treat@scores[,4],
#                         pcaRes_mrna17581_treat@scores[,5]) 
# colnames(PCA_1758mrna) = c("PCA1", "PCA2", "PCA3","PCA4", "PCA5")
# 
# ggplot(PCA_1758mrna, aes(x = PCA1, y = PCA2)) +
#   geom_point() +
#   scale_colour_manual(values = c("#04179b", "#da9e00", "#198c19"),
#                       aesthetics = "fill") + theme_light() + 
#                       ggtitle("mRNA Treatment pca")

#=======================#
# 3. PCA for miRNA DATA #
#=======================#
                        # Notes: 1. Batch_number_mirna changed to mirna.batches
                        #        2. Batch_num changed to batch
                        #        3. PCA_mirna changed to mirna.PCA
                        #        4. Treatment changed to treatment
pcaRes_mirna <- pca(t(mirna[,1:28]),nPcs = 10)  # perform PCA
mirna.PCA <- data.frame(c(pcaRes_mirna@scores[,1]),
                        pcaRes_mirna@scores[,2],
                        pcaRes_mirna@scores[,3],
                        pcaRes_mirna@scores[,4],
                        pcaRes_mirna@scores[,5],
                        treat = mirna.treatmentOrder$treatment) 
colnames(mirna.PCA) = c("PCA1", "PCA2", "PCA3","PCA4", "PCA5", "treat")

##==== Colouring by treatment ====##
ggplot(mirna.PCA, aes(x = PCA1, y = PCA2)) +
  geom_point(aes(colour = mirna.PCA$treat)) +
  scale_colour_manual(values = c("#04179b", "#da9e00", "#198c19"),
                      aesthetics = "fill") +
  theme_light()+ ggtitle("miRNA pca coloured by treatment") # for the main title

##==== Colouring by batch ====## 
ggplot(mirna.PCA, aes(x = PCA1, y = PCA2)) +
  geom_point(aes(colour = mirna.batches$batch)) +
  scale_colour_manual(values = c("#04179b", "#da9e00", "#198c19"),
                      aesthetics = "fill") +
  theme_light() + ggtitle("miRNA pca coloured by batch") # for the main title

# MIRNA pca for all microRNA (1758)
# mirna[is.na(mirna)] <- NA # change Nan for NA
# pcaRes_mirna17581_treat <- pca(mirna[,2:29],nPcs = 10)  # perform PCA
# PCA_1758mirna<- data.frame(c(pcaRes_mirna17581_treat@scores[,1]),
#                     pcaRes_mirna17581_treat@scores[,2],
#                     pcaRes_mirna17581_treat@scores[,3],
#                     pcaRes_mirna17581_treat@scores[,4],
#                     pcaRes_mirna17581_treat@scores[,5]) 
# colnames(PCA_1758mirna) = c("PCA1", "PCA2", "PCA3","PCA4", "PCA5","treat")
# 
# ggplot(PCA_1758mirna, aes(x = PCA1, y = PCA2)) +
#   geom_point(aes()) +
#   scale_colour_manual(values = c("#04179b", "#da9e00", "#198c19"),
#                       aesthetics = "fill") +
#   theme_light()

##########

#==============#
# 4. Heat maps #
#==============#
# Don't remove only silenced because it takes too much time to run. 

#heatmap.2(as.matrix(mrna), trace = "none", main="mRNA heatmap")
#heatmap.2(as.matrix(na.omit(mirna)), trace = "none", main="miRNA heatmap")


#======================================================================#
# 5. PCA showing batches and treatment simultaneously (symbol & color) #
#======================================================================#
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

#===========================#
# 6. ANOVA (mrna and mirna) #
#===========================#

# Anova testing mRNA
means_mrna<- data.frame(means= rowMeans(t(mrna2[,2:29]),na.rm = TRUE))
batch_mrna<-batch_number_mrna
anova_Test_mrna<- aov(means_mrna[,1]~batch_mrna[,1])
TukeyHSD(anova_Test_mrna, ordered = FALSE, conf.level = 0.95)

# Anova Testing miRna
means_mirna<- data.frame(means= rowMeans(t(mirna2[,2:29]),na.rm = TRUE))
batch_mirna<-batch_number_mirna
anova_Test_mirna<- aov(means_mirna[,1]~batch_mirna[,1])
TukeyHSD(anova_Test_mirna, ordered = FALSE, conf.level = 0.95)
#The anova shows that the mirna doesn't need correction as there is no significant batch variance

#=============================#
# 7. Batch effect correction. #
#=============================#

# assign numbers 1-4 to batches
batch1 <- batch_number_mrna == "410978"
batch_number_mrna$id[batch1[,1]] <- 1
batch2 <- batch_number_mrna == "410979"
batch_number_mrna$id[batch2[,1]] <- 2
batch3 <- batch_number_mrna == "410980"
batch_number_mrna$id[batch3[,1]] <- 3
batch4 <- batch_number_mrna == "412287"
batch_number_mrna$id[batch4[,1]] <- 4

# Perform correction
mrna_corrected <- removeBatchEffect(mrna[, 2:29], batch = factor(batch_number_mrna$id)) 

# Correction boxplot
boxplot(as.data.frame(mrna2[,2:29]),main="Original")
boxplot(as.data.frame(mrna_corrected),main="Batch corrected")

# remove batch effects mirna: NOT NECESSARY
# mirna_corrected <- removeBatchEffect(mirna2[, 2:29], batch = factor(batch_number_mirna$id))
# boxplot(as.data.frame(mirna2[,2:29]),main="Original miRNA")
# boxplot(as.data.frame(mirna_corrected),main="Batches miRNA corrected")

#===========================#
# 7. PCA for corrected data #
#===========================#

pcaRes_mrna_corrected <- pca(t(mrna_corrected),nPcs = 10)  # perform PCA
PCA_mrna_corrected <- data.frame(c(pcaRes_mrna_corrected@scores[,1]),
                                 pcaRes_mrna_corrected@scores[,2],
                                 pcaRes_mrna_corrected@scores[,3],
                                 pcaRes_mrna_corrected@scores[,4],
                                 pcaRes_mrna_corrected@scores[,5],
                                 batch=batch_number_mrna$Batch_num) 
colnames(PCA_mrna_corrected) = c("PCA1", "PCA2", "PCA3","PCA4", "PCA5","Batches")

# 
ggplot(PCA_mrna_corrected, aes(x = PCA1, y = PCA2)) +
  geom_point(aes(colour = treat_order_mrna$Treatment)) +
  scale_colour_manual(values = c("#04179b", "#da9e00", "#198c19","#66049b"),
                      aesthetics = "fill") + theme_light() + 
  ggtitle("mRNA PCA coloured by treatment after correction")

ggplot(PCA_mrna_corrected, aes(x = PCA1, y = PCA2)) +
  geom_point(aes(colour = batch_number_mrna$Batch_num)) +
  scale_colour_manual(values = c("#04179b", "#da9e00", "#198c19","#66049b"),
                      aesthetics = "fill") + theme_light() + 
  ggtitle("mRNA PCA coloured by batch after correction")

# Interbatch variability was very high before the coorection. This didn't allow us to see the 
# Intrabatch variability. After the correction all batches are more similar and we can see the
# intrabatch variability. 