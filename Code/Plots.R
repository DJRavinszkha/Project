#================================#
##       Adjacency Matrix       ##
#================================#
adj.Mat <- function(data){
  # Initiliase adjacency matrix
  adj.matrix <- matrix(data = 0,
                       ncol = length(unique(data[,"mirna"])),
                       nrow = length(unique(data[,"mRNA"])))
  colnames(adj.matrix) <- unique(data[,"mirna"])
  rownames(adj.matrix) <- unique(data[,"mRNA"])
  
  # Create adjacency matrix
  for (i in 1:nrow(data)){
    row <- data[i,"mRNA"]
    col <- data[i,"mirna"]
    adj.matrix[row, col] <- adj.matrix[row, col] + 1
  }
  return(adj.matrix)
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
heatmap <- function(data, data.deg, title){
  data.deg$mRNA <- rownames(data.deg)
  data$mRNA <- rownames(mrna)
  data.heat <- merge(data, data.deg, by = "mRNA", sort = FALSE)
  heatmap.2(as.matrix(unique(data.heat[,2:29])),
            trace = "none", 
            main = title,
            labRow = data.heat$Gene.Symb)
}


#=========================================#
##           miRNA Dendrogram            ##
#=========================================#

#= miRNA CHVC =#
mRNA.CHVC <- mRNA.CHVC[[1]]
miRNA.CHVC <- miRNA.CHVC[[1]]

adj.matrix <- adj.Mat(miRNA.targets(mRNA.CHVC, miRNA.CHVC))
#Normalise and Scale
adj.matrix <- t(adj.matrix)
norm.adj.matrix <- normalise(adj.matrix)
scale.norm.adj.matrix <- as.data.frame(scale(norm.adj.matrix))
summary(scale.norm.adj.matrix)

#Dist matrix, cluster and dendrogram
dist.mat <- dist(scale.norm.adj.matrix, method = 'euclidean')
hc <- hclust(dist.mat, method = "ward.D2")
hcd <- as.dendrogram(hc)

#Nodes
nodePar <- list(lab.cex = 0.6, 
                pch = c(NA, 19), 
                cex = 0.7, 
                col = "steelblue")

#Plot and add labels
plot(hcd, nodePar = nodePar)

#Colors
colors <- c("#E69F00", "#56B4E9")

miRNA.CHVC.down <- miRNA.CHVC[miRNA.CHVC[,2]<0,] #red
miRNA.CHVC.up <- miRNA.CHVC[miRNA.CHVC[,2]>0,] #blue
miRNA.CHVC.down$groupCodes <- 1
miRNA.CHVC.up$groupCodes <- 2

color_Codes <- rbind(miRNA.CHVC.down, miRNA.CHVC.up)
hcd_labels <- data.frame(mirna.Name = hcd %>% labels)

color_labels <- merge(hcd_labels, color_Codes, by = "mirna.Name", sort = FALSE)
labels_colors(hcd) <- colors[color_labels$groupCodes]

plot(hcd, nodePar = nodePar)
mtext(side = 1, text = "miRNA", line = 3.5)
mtext(side = 2, text = "Euclidean Distance", line = 2)
legend("topright", 
       legend = c("Up Regulated" , "Down Regulated"), 
       col = c("#56B4E9", "#E69F00"), 
       pch = c(20,20,4,4,4), bty = "n",  pt.cex = 1.5, cex = 0.8 , 
       text.col = "black", horiz = FALSE, inset = c(0, 0.1))

