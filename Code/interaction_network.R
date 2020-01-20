miRNA_mRNA_interaction_network <- function(miRNA.targets.CHVC){ # Still have to incorporate other  comparison groups
  data <- miRNA.targets.CHVC@data                                                                 # load data
  data <- data[!duplicated(data[,c("mature_mirna_id", "target_entrez", "pubmed_id")]),]  # Filter out duplicates
  data <- data[data[,"target_entrez"] != "",] # Remove rows with missing entrez ID's
  
  adj.matrix <- matrix(data = 0,
                       ncol = length(unique(data[,"mature_mirna_id"])),
                       nrow = length(unique(data[,"target_entrez"])))
  colnames(adj.matrix) <- unique(data[,"mature_mirna_id"])
  rownames(adj.matrix) <- unique(data[,"target_entrez"])
  
  for (i in 1:nrow(data)){
    row <- data[i,"target_entrez"]
    col <- data[i,"mature_mirna_id"]
    adj.matrix[row, col] <- adj.matrix[row, col] + 1
  }
  
  
  
}