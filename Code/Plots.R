#================================#
##       Adjacency Matrix       ##
#================================#
adj.Mat <- function(d){
  A <- matrix(data = 0,
              ncol = length(unique(d[,"mRNA"])),
              nrow = length(unique(d[,"mirna"])))
  colnames(A) <- unique(d[,"mRNA"])
  rownames(A) <- unique(d[,"mirna"])
  
  # Create adjacency matrix
  for (i in 1:nrow(d)){
    row <- as.character(d[i,"mirna"])
    col <- as.character(d[i,"mRNA"])
    A[row, col] <- A[row, col] + 1
  }
  return(A)
}



#================================#
##          Normalise           ##
#================================#
normalise <- function(x) {
  return ((x - min(x)) / (max(x) - min(x)))
}

#=========================================#
##              mRNA heatmap             ##
#=========================================#
##==== mRNA CHVC ====##
#= Heatmap =#
heatmap <- function(data, 
                    data.deg, 
                    treatment, 
                    title, 
                    save = TRUE, 
                    title_main){
  
  treatment <- data.frame(sample = labels[,3], treat = labels[,9])
  
  treatment <- treatment[order(treatment$treat),]
  
  data.deg$mRNA <- rownames(data.deg)
  data$mRNA <- rownames(mrna)
  data.heat <- merge(data, data.deg, by = "mRNA", sort = FALSE)
  dh2 <- data.heat[,2:29]
  
  dh3 <- dh2[, c(treatment$sample)]
  dh3[,29:31] <- data.heat[,30:32]
  dh3$mRNA <- data.heat$mRNA
  colnames(data.heat[,2:29]) <- c(treatment$sample)
  
  
  dh3$Mean.chol <- rowMeans(dh3[,1:9])
  dh3$Mean.drai <- rowMeans(dh3[,10:19])
  dh3$Mean.cont <- rowMeans(dh3[,20:28])
  
  if(save){
    pdf(title_main)
    heatmap.2(as.matrix(unique(dh3[,33:35])),
              trace = "none", 
              main = title,
              labRow = dh3$Gene.Symb,
              labCol = c("Cholestatic", "Control", "Drained"),
              srtCol = 90,
              cexCol = 0.75,
              col = bluered)
    dev.off()
  } else {
    heatmap(as.matrix(unique(dh3[,33:35])),
            trace = "none", 
            main = title,
            labRow = dh3$Gene.Symb,
            labCol = c("Cholestatic", "Control", "Drained"),
            srtCol = 90,
            col = bluered)
    dev.off()
  }
}


#=========================================#
##           miRNA Dendrogram            ##
#=========================================#

dendrogram <- function(miRNA, targets, save = TRUE, title){
  
  adj.matrix <- adj.Mat(targets)
  #Normalise and Scale
  # adj.matrix <- t(A)
  norm.adj.matrix <- normalise(adj.matrix)
  scale.norm.adj.matrix <- scale(norm.adj.matrix)
  summary(scale.norm.adj.matrix)
  
  #Dist matrix, cluster and dendrogram
  dist.mat <- dist(scale.norm.adj.matrix, method = 'euclidean')
  hc <- hclust(dist.mat, method = "ward.D2")
  hcd <- as.dendrogram(hc)
  #Colors
  colors <- c("#E69F00", "#56B4E9")
  
  miRNA.down <- miRNA[miRNA[,2]<0,] #red
  miRNA.up <- miRNA[miRNA[,2]>0,] #blue
  miRNA.down$groupCodes <- 1
  miRNA.up$groupCodes <- 2
  
  color_Codes <- rbind(miRNA.down, miRNA.up)
  hcd_labels <- data.frame(mirna.Name = hcd %>% labels)
  
  color_labels <- merge(hcd_labels, color_Codes, by = "mirna.Name", sort = FALSE)
  labels_colors(hcd) <- colors[color_labels$groupCodes]
  
  hcd <- set(hcd, "labels_cex", 1.2)
  
  par(mar = c(10,4,0.1,1))
  plot(hcd)
  mtext(side = 1, text = "miRNA", line = 8, cex = 1.2)
  mtext(side = 2, text = "Distance", line = 2, cex = 1.2)
  legend("topright",
         legend = c("Up Regulated" , "Down Regulated"), 
         col = c("#56B4E9", "#E69F00"), 
         pch = c(20,20,4,4,4), bty = "n",  pt.cex = 1.5, cex = 0.8 , 
         text.col = "black", horiz = FALSE, inset = c(0, 0.1))
  
  if(save){
    pdf(title)
    plot(hcd)
    mtext(side = 1, text = "miRNA", line = 3.5)
    mtext(side = 2, text = "Euclidean Distance", line = 2)
    legend(20,56,
           legend = c("Up Regulated" , "Down Regulated"), 
           col = c("#56B4E9", "#E69F00"), 
           pch = c(20,20,4,4,4), bty = "n",  pt.cex = 1.5, cex = 0.8 , 
           text.col = "black", horiz = FALSE, inset = c(0, 0.1))
    dev.off()
    
  }
}