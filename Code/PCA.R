require(limma)
require(pcaMethods)
require(gplots)
require(ggplot2)
library(FactoMineR)
library(factoextra)
# MIRNA pca for 28 samples
mirna[is.na(mirna)] <- NA # change Nan for NA
pcaRes <- pca(t(mirna[,2:29]),nPcs = 10)  # perform PCA
PCA_28mirna_treat<- data.frame(c(pcaRes@scores[,1]),
                    pcaRes@scores[,2],
                    pcaRes@scores[,3],
                    pcaRes@scores[,4],
                    pcaRes@scores[,5],
                    treat=labels$V2) 
colnames(PCA_28mirna_treat) = c("PCA1", "PCA2", "PCA3","PCA4", "PCA5", "treat")

ggplot(PCA_28mirna_treat, aes(x = PCA1, y = PCA2)) +
  geom_point(aes(colour = PCA_28mirna_treat$treat)) +
  scale_colour_manual(values = c("#04179b", "#da9e00", "#198c19"),
                      aesthetics = "fill") +
  theme_light()

# MIRNA pca for all microRNA (1758)
mirna[is.na(mirna)] <- NA # change Nan for NA
pcaRes <- pca(mirna[,2:29],nPcs = 10)  # perform PCA
PCA_1758mirna<- data.frame(c(pcaRes@scores[,1]),
                    pcaRes@scores[,2],
                    pcaRes@scores[,3],
                    pcaRes@scores[,4],
                    pcaRes@scores[,5],
                    treat=labels$V2) 
colnames(PCA_1758mirna) = c("PCA1", "PCA2", "PCA3","PCA4", "PCA5","treat")

ggplot(PCA_1758mirna, aes(x = PCA1, y = PCA2)) +
  geom_point(aes(colour = PCA_1758mirna$treat)) +
  scale_colour_manual(values = c("#04179b", "#da9e00", "#198c19"),
                      aesthetics = "fill") +
  theme_light()


#MRNA for 28 samples
mrna[is.na(mrna)] <- NA # change Nan for NA
pcaRes2 <- pca(t(mrna[,2:29]),nPcs = 10)  # perform PCA
PCA_28mrna_treat<- data.frame(c(pcaRes2@scores[,1]),
                    pcaRes2@scores[,2],
                    pcaRes2@scores[,3],
                    pcaRes2@scores[,4],
                    pcaRes2@scores[,5],
                    treat=labels$V2) 
colnames(PCA_28mrna_treat) = c("PCA1", "PCA2", "PCA3","PCA4", "PCA5","treat")

ggplot(PCA_28mrna_treat, aes(x = PCA1, y = PCA2)) +
  geom_point(aes(colour = PCA_28mrna_treat$treat)) +
  scale_colour_manual(values = c("#04179b", "#da9e00", "#198c19"),
                      aesthetics = "fill") +
  theme_light()

#MRNA for all microRNA (17581)
mrna[is.na(mrna)] <- NA # change Nan for NA
pcaRes2 <- pca(mrna[,2:29],nPcs = 10)  # perform PCA
PCA_1758mrna<- data.frame(c(pcaRes2@scores[,1]),
                        pcaRes2@scores[,2],
                        pcaRes2@scores[,3],
                        pcaRes2@scores[,4],
                        pcaRes2@scores[,5]) 
colnames(PCA_1758mrna) = c("PCA1", "PCA2", "PCA3","PCA4", "PCA5")

ggplot(PCA_1758mrna, aes(x = PCA1, y = PCA2)) +
  geom_point() +
  scale_colour_manual(values = c("#04179b", "#da9e00", "#198c19"),
                      aesthetics = "fill") +
  theme_light()

############# HEAT MAP ######################

heatmap.2(as.matrix(mrna[,2:29]), trace = "none", main="mRNA heatmap")
heatmap.2(as.matrix(na.omit(mirna[,2:29])), trace = "none", main="miRNA heatmap")

############# Use Factominer to complete PCA and HCPC ######################
res.pca <- PCA(mrna[,2:29], ncp = 3, graph = FALSE)
# Compute hierarchical clustering on principal components
res.hcpc <- HCPC(res.pca, graph = FALSE)

fviz_dend(res.hcpc, 
          cex = 0.7,                     # Label size
          palette = "jco",               # Color palette see ?ggpubr::ggpar
          rect = TRUE, rect_fill = TRUE, # Add rectangle around groups
          rect_border = "jco",           # Rectangle color
          labels_track_height = 0.8      # Augment the room for labels
)

#############  PCA FOR BATCHES  ######################

mrna = read.delim('../Data/GeneExpressionNormalized.txt', check.names = FALSE) #Load expression data
mrna = subset(mrna, select = -2) # Remove gene symbol column
names_mrna = colnames(mrna[,2:length(mrna)]) #
batch_number_mrna <- data.frame(Batch_num= substr(unlist(names_mrna), 8, 13)) #Only keep the batch number 
name_batches<- unique(batch_number_mrna)

mrna[is.na(mrna)] <- NA # change Nan for NA
pcaRes2 <- pca(t(mrna[,2:29]),nPcs = 10)  # perform PCA
PCA_28mrna<- data.frame(c(pcaRes2@scores[,1]),
                        pcaRes2@scores[,2],
                        pcaRes2@scores[,3],
                        pcaRes2@scores[,4],
                        pcaRes2@scores[,5],
                        batch= batch_number_mrna) 
colnames(PCA_28mrna) = c("PCA1", "PCA2", "PCA3","PCA4", "PCA5", "Batches")

ggplot(PCA_28mrna, aes(x = PCA1, y = PCA2)) +
  geom_point(aes(colour = PCA_28mrna_treat$treat, shape=PCA_28mrna$Batches)) +
  scale_colour_manual(values = c("#04179b", "#da9e00", "#198c19","#66049b"),
                      aesthetics = "fill") +
  theme_light()

mirna = read.delim('../Data/miRNAexpression.txt', check.names = FALSE)
batch_number_mirna <- data.frame(Batch_num= substr(unlist(colnames(mirna[,2:length(mirna)])), 12, 23)) #Only keep the batch number 
name_batches_mirna<- unique(batch_number_mirna)

pcaRes <- pca(t(mirna[,2:29]),nPcs = 10)  # perform PCA
PCA_28mirna<- data.frame(c(pcaRes@scores[,1]),
                        pcaRes@scores[,2],
                        pcaRes@scores[,3],
                        pcaRes@scores[,4],
                        pcaRes@scores[,5],
                        batch= batch_number_mirna) 
colnames(PCA_28mirna) = c("PCA1", "PCA2", "PCA3","PCA4", "PCA5", "Batches")

ggplot(PCA_28mirna, aes(x = PCA1, y = PCA2)) +
  geom_point(aes(colour = PCA_28mirna$Batches)) +
  scale_colour_manual(values = c("#04179b", "#da9e00", "#198c19","#66049b"),
                      aesthetics = "fill") +
  theme_light()

#### correction for batch that is separated###

#=========================================#
# 5. ANOVA 
#==========================================#
means<- data.frame(rowMeans(t(mrna[,2:29]),na.rm = TRUE))
batch<-batch_number_mrna
anova_Test<- aov(means[,1]~batch[,1])
summary(anova_Test)
TukeyHSD(anova_Test, ordered = FALSE, conf.level = 0.95)

