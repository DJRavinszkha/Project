

mirna[is.na(mirna)] <- NA # change Nan for NA
pcaRes <- pca(mirna[,2:29],nPcs = 10)  # perform PCA
PCAtry<- data.frame(c(pcaRes@scores[,1]),
                    pcaRes@scores[,2],
                    pcaRes@scores[,3],
                    pcaRes@scores[,4],
                    pcaRes@scores[,5],
                    treat=labels$V2) 
colnames(PCAtry) = c("PCA1", "PCA2", "PCA3", "PCA4", "PCA5","treat")

ggplot(PCAtry, aes(x = PCA1, y = PCA2)) +
  geom_point(aes(colour = PCAtry$treat)) +
  scale_colour_manual(values = c("#04179b", "#da9e00", "#198c19"),
                      aesthetics = "fill") +
  theme_light()