require(limma)
require(pcaMethods)
require(gplots)
require(ggplot2)
library(FactoMineR)
library(factoextra)

#=========================================#
# 1. Get list of how treatments and batches are ordered on mRNA and miRNA. 
#==========================================#

# Get treatment order mirna
treat_order_mirna <- matrix(nrow = (ncol(mirna) -1), ncol = 1)
colnames(treat_order_mirna) <- "SampleName"
treat_order_mirna[,1] <- colnames(mirna[,2:length(mirna)])
treat_order_mirna <- merge(treat_order_mirna, sampleGroups, by =  "SampleName", sort = FALSE)

# Get treatment order mrna
treat_order_mrna <- matrix(nrow = (ncol(mrna) -1), ncol = 1)
colnames(treat_order_mrna) <- "SampleName"
treat_order_mrna[,1] <- colnames(mrna[,2:length(mrna)])
treat_order_mrna <- merge(treat_order_mrna, sampleGroups, by =  "SampleName", sort = FALSE)


# Get batch order mrna
mrna2 = read.delim('../Data/GeneExpressionNormalized.txt', check.names = FALSE) #Load expression data
mrna2[is.na(mrna2)] <- NA
colnames(mrna2)[10]<-"FGS_10_410978_2_3_H" #Don't move this line, since this is the correction of the original data set
mrna2 = subset(mrna2, select = -2) # Remove gene symbol column
names_mrna2 = colnames(mrna2[,2:length(mrna2)]) #
batch_number_mrna2 <- data.frame(Batch_num= substr(unlist(names_mrna2), 8, 13)) #Only keep the batch number 
name_batches2<- unique(batch_number_mrna2)

# Get batch order mirna
mirna2 = read.delim('../Data/miRNAexpression.txt', check.names = FALSE)
batch_number_mirna2 <- data.frame(Batch_num= substr(unlist(colnames(mirna2[,2:length(mirna)])), 12, 23)) #Only keep the batch number 
name_batches_mirna2<- unique(batch_number_mirna2)

#=========================================#
# 2. PCA for mirna DATA
#==========================================#

pcaRes_mirna28_treat <- pca(t(mirna[,2:29]),nPcs = 10)  # perform PCA
PCA_28mirna_treat<- data.frame(c(pcaRes_mirna28_treat@scores[,1]),
                    pcaRes_mirna28_treat@scores[,2],
                    pcaRes_mirna28_treat@scores[,3],
                    pcaRes_mirna28_treat@scores[,4],
                    pcaRes_mirna28_treat@scores[,5],
                    treat=treat_order_mirna$Treatment) 
colnames(PCA_28mirna_treat) = c("PCA1", "PCA2", "PCA3","PCA4", "PCA5", "treat")

ggplot(PCA_28mirna_treat, aes(x = PCA1, y = PCA2)) +
  geom_point(aes(colour = PCA_28mirna_treat$treat)) +
  scale_colour_manual(values = c("#04179b", "#da9e00", "#198c19"),
                      aesthetics = "fill") +
  theme_light()+ ggtitle("miRNA pca coloured by treatment") # for the main title


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

#=========================================#
# 3. PCA for mRNA DATA coloured by treatment
#==========================================#

#MRNA for 28 samples
mrna[is.na(mrna)] <- NA # change Nan for NA
pcaRes_mrna28_treat <- pca(t(mrna[,2:29]),nPcs = 10)  # perform PCA
PCA_28mrna_treat<- data.frame(c(pcaRes_mrna28_treat@scores[,1]),
                    pcaRes_mrna28_treat@scores[,2],
                    pcaRes_mrna28_treat@scores[,3],
                    pcaRes_mrna28_treat@scores[,4],
                    pcaRes_mrna28_treat@scores[,5],
                    treat=treat_order_mrna$Treatment) 
colnames(PCA_28mrna_treat) = c("PCA1", "PCA2", "PCA3","PCA4", "PCA5","treat")

ggplot(PCA_28mrna_treat, aes(x = PCA1, y = PCA2)) +
  geom_point(aes(colour = PCA_28mrna_treat$treat)) +
  scale_colour_manual(values = c("#04179b", "#da9e00", "#198c19"),
                      aesthetics = "fill") + theme_light() + 
                      ggtitle("mRNA pca coloured by Treatment")

ggplot(PCA_28mrna_treat, aes(x = PCA1, y = PCA2)) +
  geom_point(aes(colour = batch_number_mrna2$Batch_num)) +
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

#=========================================#
# 4. Heat maps
#==========================================#
# Don't remove only silenced because it takes too much time to run. 

# heatmap.2(as.matrix(mrna[,2:29]), trace = "none", main="mRNA heatmap")
# heatmap.2(as.matrix(na.omit(mirna[,2:29])), trace = "none", main="miRNA heatmap")


#=========================================#
# 4. PCA showing batches and treatment simultaneously (symbol & color)
#==========================================#
#Don't remove just in case we need it. 

# mrna = read.delim('../Data/GeneExpressionNormalized.txt', check.names = FALSE) #Load expression data
# mrna = subset(mrna, select = -2) # Remove gene symbol column
# names_mrna = colnames(mrna[,2:length(mrna)]) #
# batch_number_mrna <- data.frame(Batch_num= substr(unlist(names_mrna), 8, 13)) #Only keep the batch number 
# name_batches<- unique(batch_number_mrna)
# 
# mrna[is.na(mrna)] <- NA # change Nan for NA
# pcaRes2 <- pca(t(mrna[,2:29]),nPcs = 10)  # perform PCA
# PCA_28mrna<- data.frame(c(pcaRes2@scores[,1]),
#                         pcaRes2@scores[,2],
#                         pcaRes2@scores[,3],
#                         pcaRes2@scores[,4],
#                         pcaRes2@scores[,5],
#                         batch= batch_number_mrna) 
# colnames(PCA_28mrna) = c("PCA1", "PCA2", "PCA3","PCA4", "PCA5", "Batches")
# 
# ggplot(PCA_28mrna, aes(x = PCA1, y = PCA2)) +
#   geom_point(aes(colour = PCA_28mrna_treat$treat, shape=PCA_28mrna$Batches)) +
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



#=========================================#
# 5. PCA on batches coloring by treatment BATCH EFFECT: NO CORRECTION
#==========================================#

pcaRes_mrna28_batch <- pca(t(mrna[,2:29]),nPcs = 10)  # perform PCA
PCA_28mrna_batch<- data.frame(c(pcaRes_mrna28_batch@scores[,1]),
                        pcaRes_mrna28_batch@scores[,2],
                        pcaRes_mrna28_batch@scores[,3],
                        pcaRes_mrna28_batch@scores[,4],
                        pcaRes_mrna28_batch@scores[,5],
                        batch= batch_number_mrna2$Batch_num)
colnames(PCA_28mrna_batch) = c("PCA1", "PCA2", "PCA3","PCA4", "PCA5","Batches")

ggplot(PCA_28mrna_batch, aes(x = PCA1, y = PCA2)) +
  geom_point(aes(colour = treat_order_mrna$Treatment)) +
  scale_colour_manual(values = c("#04179b", "#da9e00", "#198c19","#66049b"),
                      aesthetics = "fill") +
  theme_light() + ggtitle("mRNA PCA on batches coloured by treatment BEFORE correction")

#=========================================#
# 6. PCA on batches coloring by treatment BATCH EFFECT: CORRECTION
#==========================================#

# assign numbers 1-4 to batches
batch1 <- batch_number_mrna2 == "410978"
batch_number_mrna2$id[batch1[,1]] <- 1
batch2 <- batch_number_mrna2 == "410979"
batch_number_mrna2$id[batch2[,1]] <- 2
batch3 <- batch_number_mrna2 == "410980"
batch_number_mrna2$id[batch3[,1]] <- 3
batch4 <- batch_number_mrna2 == "412287"
batch_number_mrna2$id[batch4[,1]] <- 4

# Perform correction
y2 <- removeBatchEffect(mrna[, 2:29], batch = factor(batch_number_mrna2$id)) 

# Correction boxplot
#boxplot(as.data.frame(mrna2[,2:29]),main="Original")
#boxplot(as.data.frame(y2),main="Batch corrected")

# PCA corrected
pcaRes_mrna28_batch_corrected <- pca(t(y2),nPcs = 10)  # perform PCA
PCA_28mrna_batch_corrected<- data.frame(c(pcaRes_mrna28_batch_corrected@scores[,1]),
                              pcaRes_mrna28_batch_corrected@scores[,2],
                              pcaRes_mrna28_batch_corrected@scores[,3],
                              pcaRes_mrna28_batch_corrected@scores[,4],
                              pcaRes_mrna28_batch_corrected@scores[,5],
                              batch=batch_number_mrna2$Batch_num) 
colnames(PCA_28mrna_batch_corrected) = c("PCA1", "PCA2", "PCA3","PCA4", "PCA5","Batches")

ggplot(PCA_28mrna_batch_corrected, aes(x = PCA1, y = PCA2)) +
  geom_point(aes(colour = treat_order_mrna$Treatment)) +
  scale_colour_manual(values = c("#04179b", "#da9e00", "#198c19","#66049b"),
                      aesthetics = "fill") + theme_light() + 
  ggtitle("mRNA PCA on batches coloured by treatment AFTER correction")

ggplot(PCA_28mrna_batch_corrected, aes(x = PCA1, y = PCA2)) +
  geom_point(aes(colour = batch_number_mrna2$Batch_num)) +
  scale_colour_manual(values = c("#04179b", "#da9e00", "#198c19","#66049b"),
                      aesthetics = "fill") + theme_light() + 
  ggtitle("mRNA PCA coloured by batch AFTER correction")

# Interbatch variability was very high before the coorection. This didn't allow us to see the 
# Intrabatch variability. After the correction all batches are more 
#=========================================#
# 7. ANOVA 
#==========================================#

means<- data.frame(means= rowMeans(t(mrna2[,2:29]),na.rm = TRUE))
batch<-batch_number_mrna2
anova_Test<- aov(means[,1]~batch[,1])
summary(anova_Test)
TukeyHSD(anova_Test, ordered = FALSE, conf.level = 0.95)
###########################################################################
# anova in mirna
batch1 <- batch_number_mirna2 == "254606410403"
batch_number_mirna2$id[batch1[,1]] <- 1
batch2 <- batch_number_mirna2 == "254606410404"
batch_number_mirna2$id[batch2[,1]] <- 2
batch3 <- batch_number_mirna2 == "254606410405"
batch_number_mirna2$id[batch3[,1]] <- 3
batch4 <- batch_number_mirna2 == "254606410413"
batch_number_mirna2$id[batch4[,1]] <- 4
batch5 <- batch_number_mirna2 == "254606411109"
batch_number_mirna2$id[batch5[,1]] <- 5

remove_be_mirna <- removeBatchEffect(mirna2[, 2:29], batch = factor(batch_number_mirna2$id))
boxplot(as.data.frame(mirna2[,2:29]),main="Original miRNA")
boxplot(as.data.frame(remove_be_mirna),main="Batches miRNA corrected")

means_mirna<- data.frame(means= rowMeans(t(mirna2[,2:29]),na.rm = TRUE))
batch_mirna<-batch_number_mirna2
anova_Test_mirna<- aov(means_mirna[,1]~batch_mirna[,1])
summary(anova_Test_mirna)
TukeyHSD(anova_Test_mirna, ordered = FALSE, conf.level = 0.95)


# The anova shows that the mirna doesn't need correction as there is no significant batch variance
