
# MIRNA pca for 28 samples
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

