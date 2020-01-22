#========================================================================================================#
# This function converts data into two datasets required for visualising a network                       #
#                                                                                                        #
# Input:
#   data: a matrix of which only the first two columns are used, it should contain interactions.
#   Columns 1 contain sources, and column two contains targets that the source interact with.
#   Column 3 can contain a value to weight the interaction --- CURRENTLY NON-FUNCTIONAL
#
# Output:
# network_node:
#   This data file contains the information on the different nodes in the  network, it contains the ID's #
#   starting from 0, their according name that they correspond to, this can for example be entrez gene ID#
#   or miRNA name. Also it contains a Group column where you can distinguish groups. By default this 
#========================================================================================================#

generate_network_data <- function(data){

  # Create Source ID's
  names <- unique(data[,1])
  source_ID <- matrix(ncol = 3, nrow = length(names))   # Initialise ID matrix
  source_ID[,2] <- names                                # Get all unqiue source names
  source_ID[,1] <- seq(0,length(names) - 1)             # Create the ID's
  source_ID[,3] <- 1                                    # Set group to Source
  colnames(source_ID) <- c("ID", "name", "Group")       # Create column names
  
  # Create Target ID's
  names <- unique(data[,2])
  target_ID <- matrix(ncol = 3, nrow = length(names))   # Initialise ID matrix
  target_ID[,2] <- names                                # Get all unqiue source names
  target_ID[,1] <- seq(0,length(names) - 1)             # Create the ID's
  target_ID[,3] <- 1                                    # Set group to Source
  colnames(target_ID) <- c("ID", "name", "Group")       # Create column names
  
  network_data <- merge(data, source_ID, by.x = 1, by.y = "name", sort = FALSE)         # Put source ID's in network_data
  network_data <- merge(network_data, target_ID, by.x = 2, by.y = "name", sort = FALSE) # Put target ID's in network_data
  
  # Create link data for network
  network_link <- network_data[,c("ID.x", "ID.y")]
  colnames(network_link) <- c("Source", "Target")
  
  # Create node data for network
  network_node <- data.frame(do.call(rbind, list(source_ID, target_ID))) # Create node data by stacking the two ID matrices on top of eachother
  
  return(list(network_link, network_node))
}



generate_miRNA_mRNA_interaction_network_network <- function(data){ # Still have to incorporate other  comparison groups
  data <- miRNA.targets.CHVC.DEG                                                                 # load data
  data <- data[!duplicated(data[,c("mature_mirna_id", "target_entrez", "pubmed_id")]),]  # Filter out duplicates
  data <- data[data[,"target_entrez"] != "",] # Remove rows with missing entrez ID's
  
  # Initiliase adjacency matrix
  adj.matrix <- matrix(data = 0,
                       ncol = length(unique(data[,"mature_mirna_id"])),
                       nrow = length(unique(data[,"target_entrez"])))
  colnames(adj.matrix) <- unique(data[,"mature_mirna_id"])
  rownames(adj.matrix) <- unique(data[,"target_entrez"])
  
  # Create adjacency matrix
  for (i in 1:nrow(data)){
    row <- data[i,"target_entrez"]
    col <- data[i,"mature_mirna_id"]
    adj.matrix[row, col] <- adj.matrix[row, col] + 1
  }
  
  
  #=== Create 3 columned matrix for network visualisation ===#
  data[,"Occurance"] <- 1
  network_data <- aggregate(Occurance~mature_mirna_id + target_entrez, data=data,FUN=sum)
  
  # Create miRNA ID's
  IDs <- unique(network_data[,"mature_mirna_id"])
  miRNA_ID <- matrix(ncol = 3, nrow = length(IDs))
  miRNA_ID[,2] <- IDs # Get all unqiue miRNA names
  miRNA_ID[,1] <- seq(0,length(IDs) - 1)
  miRNA_ID[,3] <- 1
  colnames(miRNA_ID) <- c("ID", "name", "Group")
  
  # Create mRNA ID's
  IDs <- unique(network_data[,"target_entrez"])
  mRNA_ID <- matrix(ncol = 3, nrow = length(IDs))
  mRNA_ID[,2] <- IDs # Get all unqiue miRNA names
  mRNA_ID[,1] <- seq(nrow(miRNA_ID),nrow(miRNA_ID) + length(IDs) - 1)
  mRNA_ID[,3] <- 2
  colnames(mRNA_ID) <- c("ID", "name", "Group")
  
  
  network_data <- merge(network_data, miRNA_ID, by.x = "mature_mirna_id", by.y = "name", sort = FALSE) # Put mirna ID's in network_data
  network_data <- merge(network_data, mRNA_ID, by.x = "target_entrez", by.y = "name", sort = FALSE) # Put mirna ID's in network_data
  
  # Create link data for network
  network_link <- network_data[,c("ID.x", "ID.y", "Occurance")]
  colnames(network_link) <- c("Source", "Target", "Value")
  network_link_filtered <- network_link[network_link[,"Value"] > 5,]
  
  # Create node data for network
  network_node <- data.frame(do.call(rbind, list(miRNA_ID, mRNA_ID))) # Create node data by stacking the two ID matrices on top of eachother
  
}

filter_interactions <- function(){
  source("Expression_Analysis.R")
  
  mRNA.CHVC <- mRNA.CHVC()
}

visualise_network <- function(network_link, network_node, silent = TRUE){
  if (!requireNamespace("networkD3", quietly = TRUE))
    install.packages("networkD3")
  library(networkD3)
  
  miRNA_interaction_network <- forceNetwork(Links   = network_link,
                                             Nodes   = network_node,
                                             Source  = "Source",
                                             Target  = "Target",
                                             Value   = "Value",
                                             NodeID  = "name",
                                             Group   = "Group",
                                             linkWidth = 1,
                                             opacity = 0.8,
                                             zoom    = TRUE,
                                             fontSize = 30
                                             )
  if (!silent){
    miRNA_interaction_network
  }
  
  saveNetwork(miRNA_interaction_network, "miRNA_interaction_network.html")
  
  forceNetwork(Links = MisLinks, Nodes = MisNodes, Source = "source",Target = "target", Value = "value", NodeID = "name",Group = "group", opacity = 0.4, zoom = TRUE)
  
  return(miRNA_interaction_network)
  }