require(limma)
library(pcaMethods)
require(gplots)
require(ggplot2)


# MIRNA pca 
#Expression of variance of all samples (28 points)
#Interesting to see how samples are clustered based on treatment
mirna[is.na(mirna)] <- NA # change Nan for NA
pcaRes <- pca(t(mirna[,2:29]),nPcs = 10)  # perform PCA
PCA_28mirna<- data.frame(c(pcaRes@scores[,1]),
                    pcaRes@scores[,2],
                    pcaRes@scores[,3],
                    pcaRes@scores[,4],
                    pcaRes@scores[,5],
                    treat=labels$V2) 
colnames(PCA_28mirna) = c("PCA1", "PCA2", "PCA3","PCA4", "PCA5", "treat")

ggplot(PCA_28mirna, aes(x = PCA1, y = PCA2)) +
  geom_point(aes(colour = PCA_28mirna$treat)) +
  scale_colour_manual(values = c("#04179b", "#da9e00", "#198c19"),
                      aesthetics = "fill") +
  theme_light()

#Expression of variance of all miRNA (17581 points)
#These PCAS may not be that informative
#Warning: Need to remove coloring
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


#MRNA 
#Expression variance of samples (28 samples)
mrna[is.na(mrna)] <- NA # change Nan for NA
pcaRes2 <- pca(t(mrna[,2:29]),nPcs = 10)  # perform PCA
PCA_28mrna<- data.frame(c(pcaRes2@scores[,1]),
                    pcaRes2@scores[,2],
                    pcaRes2@scores[,3],
                    pcaRes2@scores[,4],
                    pcaRes2@scores[,5],
                    treat=labels$V2) 
colnames(PCA_28mrna) = c("PCA1", "PCA2", "PCA3","PCA4", "PCA5","treat")

ggplot(PCA_28mrna, aes(x = PCA1, y = PCA2)) +
  geom_point(aes(colour = PCA_28mrna$treat)) +
  scale_colour_manual(values = c("#04179b", "#da9e00", "#198c19"),
                      aesthetics = "fill") +
  theme_light()

#MRNA for all mRNA (17581 observations)
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

#############  PCA FOR BATCHES  ######################
#instead of doing PCA and coloring them based on treatment (drain, control, cholestatic); color
# them based on their batches.
  # The idea is to observe intrinsic variability between batches. 
  # PCA plot will highlight these changes....hopefully.

#Data formatting
#PCA
#saw some weird shit
#Differential Expression Analysis
#Tried to understand it 
#Realised there was something wrong
#Identified there was a batch effect
#Filtered for batch effect
#Conduct PCA including batch effect

#instead of labelling samples based on treatment, we want to identify 

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
  geom_point(aes(colour = PCA_28mrna$Batches)) +
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

