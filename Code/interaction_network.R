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
##              Functions                ##
#=========================================#

#================================#
##        General Network       ##
#================================#
generate_network_data <- function(data){
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
  
  # Create Source ID's
  names <- unique(data[,1])
  source_ID <- matrix(ncol = 3, nrow = length(names))   # Initialise ID matrix
  source_ID[,2] <- as.character(names)                                # Get all unqiue source names
  source_ID[,1] <- seq(0,length(names) - 1)             # Create the ID's
  source_ID[,3] <- 1                                    # Set group to Source
  colnames(source_ID) <- c("ID", "name", "Group")       # Create column names
  
  # Create Target ID's
  names <- unique(data[,2])
  target_ID <- matrix(ncol = 3, nrow = length(names))   # Initialise ID matrix
  target_ID[,2] <- as.character(names)                                # Get all unqiue source names
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

#========================================#
##      miRNA Interaction Network       ##
#========================================#
miRNA_mRNA_interaction_network <- function(mRNA, 
                                           miRNA, 
                                           filter = FALSE){ # Still have to incorporate other  comparison groups
  
  #Inputs:
    # mRNA = mRNA data file
    # miRNA = miRNA data file
    # filter = (default is FALSE); value specifying minimum number of sources
  
  
  data <- miRNA.targets(mRNA, miRNA)                               # load data
  data <- data[!duplicated(data[,c("mirna", "mRNA", "source")]),]  # Filter out duplicates
  data <- data[data[,"mRNA"] != "",]                               # Remove rows with missing entrez ID's
  
  #=== Create 3 columned matrix for network visualisation ===#
  data[,"Occurance"] <- 1
  network_data <- aggregate(Occurance~mirna + mRNA, data=data,FUN=sum)
  
  # Create miRNA ID's
  IDs <- unique(network_data[,"mirna"])
  miRNA_ID <- matrix(ncol = 3, nrow = length(IDs))
  miRNA_ID[,2] <- as.character(IDs)                               # Get all unqiue miRNA names
  miRNA_ID[,1] <- seq(0,length(IDs) - 1)
  miRNA_ID[,3] <- 1
  colnames(miRNA_ID) <- c("ID", "name", "Group")
  
  # Create mRNA ID's
  IDs <- unique(network_data[,"mRNA"])
  mRNA_ID <- matrix(ncol = 3, nrow = length(IDs))
  mRNA_ID[,2] <- as.character(IDs)                               # Get all unqiue miRNA names
  mRNA_ID[,1] <- seq(nrow(miRNA_ID),nrow(miRNA_ID) + length(IDs) - 1)
  mRNA_ID[,3] <- 2
  colnames(mRNA_ID) <- c("ID", "name", "Group")
  
  
  network_data <- merge(network_data, miRNA_ID, by.x = "mirna", by.y = "name", sort = FALSE) # Put mirna ID's in network_data
  network_data <- merge(network_data, mRNA_ID, by.x = "mRNA", by.y = "name", sort = FALSE)   # Put mirna ID's in network_data
  
  # Create link data for network
  network_link <- network_data[,c("ID.x", "ID.y", "Occurance")]
  colnames(network_link) <- c("Source", "Target", "Value")
  
  #filter
  if(filter != FALSE){
    network_link_filtered <- network_link[network_link[,"Value"] > 5,]
  }
  
  # Create node data for network
  network_node <- data.frame(do.call(rbind, list(miRNA_ID, mRNA_ID))) # Create node data by stacking the two ID matrices on top of eachother
  
  return(list(network_link, network_node))
}


#========================================#
##          Visualise Network           ##
#========================================#
visualise_network <- function(network_link, 
                              network_node, 
                              silent = TRUE){
  
  miRNA_interaction_network <- forceNetwork(Links   = network_link,
                                            Nodes   = network_node,
                                            Source  = "Source",
                                            Target  = "Target",
                                            NodeID  = "name",
                                            Group   = "Group",
                                            linkWidth = 1,
                                            opacity = 0.8,
                                            zoom    = TRUE,
                                            fontSize = 30
                                            )
  saveNetwork(miRNA_interaction_network, "miRNA_interaction_network.html")
  
  if(!silent){
    return(miRNA_interaction_network)
  }
}