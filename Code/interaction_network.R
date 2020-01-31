#=============================================================================#
# Project Period, Liver cholestasis data analysis         							      #					  
#	Interaction Network Generation and visualisation                            #
# Version: 1.0   															                                #
# Date: 9-1-2020											             	                          #
# Authors: Ariadna Fosch i Muntané, ID: I6215203, Maastricht University      #
#          Jip de Kok, ID: I6119367 , Maastricht University                   #
#          Ravin Schmidl, ID: I6125866, Maastricht University                 #
#          Stefan Meier, ID: I6114194 , Maastricht University                 #
#=============================================================================#

#=========================================#
##              Functions                ##
#=========================================#

#====================================#
##        create  Network data      ##
#====================================#
generate_network_data <- function(data,
                                  contains_genes = FALSE,
                                  group1 = "mRNA",
                                  group2 = "miRNA"){
  #===========================================================================================================#
  # This function converts data into two datasets required for visualising a network                          #
  #                                                                                                           #
  # Input:                                                                                                    #
  #   data: a matrix of which only the first two columns are used, it should contain interactions.            #
  #   Column 1 contain sources, and column 2 contains targets that the source interacts with.                 #
  #   Column 3 can contain a value to weight the interaction --- CURRENTLY NON-FUNCTIONAL                     #
  #                                                                                                           #
  # Output:                                                                                                   #
  #   link_node:                                                                                              #
  #     This matrix contains all the interactions, column 1 has the source and column 2 has the target.       #
  #     Column three has a value with which the interactions can be weighted.                                 #
  #   network_node:                                                                                           #
  #     This data file contains the information on the different nodes in the  network, it contains the ID's  #
  #     starting from 0, their according name that they correspond to, this can for example be entrez gene ID #
  #     or miRNA name. Also it contains a Group column where you can distinguish groups.                      #
  #===========================================================================================================#
  
  # Create Source ID's
  if (contains_genes == TRUE){
    names <- unique(data[,1])
    source_ID <- matrix(ncol = 5, nrow = length(names))   # Initialise ID matrix
    source_ID[,2] <- as.character(names)                  # Get all unqiue source names
    source_ID[,1] <- seq(0,length(names) - 1)             # Create the ID's
    source_ID[,3] <- group1                               # Set group to Source
    colnames(source_ID) <- c("ID", "name", "Group", "Symbol", "Genename")   # Create column names
  } else {
    names <- unique(data[,1])
    source_ID <- matrix(ncol = 3, nrow = length(names))   # Initialise ID matrix
    source_ID[,2] <- as.character(names)                  # Get all unqiue source names
    source_ID[,1] <- seq(0,length(names) - 1)             # Create the ID's
    source_ID[,3] <- group1                               # Set group to Source
    colnames(source_ID) <- c("ID", "name", "Group")       # Create column names
  }
  
  # Create Target ID's
  names <- unique(data[,2])
  target_ID <- matrix(ncol = 3, nrow = length(names))                         # Initialise ID matrix
  target_ID[,2] <- as.character(names)                                        # Get all unqiue source names
  target_ID[,1] <- seq(nrow(source_ID), nrow(source_ID) + length(names) - 1)  # Create the ID's
  target_ID[,3] <- group2                                                     # Set group to Target
  colnames(target_ID) <- c("ID", "name", "Group")                             # Create column names
  
  # Replace entrezIDs with genenames in the Source_ID
  if (contains_genes == TRUE){
    cols <- c("SYMBOL", "GENENAME", "SYMBOL")
    name_mapping <- select(org.Hs.eg.db,
                           keys = source_ID[,"name"],
                           columns=cols, keytype="ENTREZID")
    source_ID[,4] <- name_mapping[,"SYMBOL"]
    source_ID[,5] <- name_mapping[,"GENENAME"]
  }
  
  # Combine source and target data
  network_data <- merge(data, source_ID, by.x = 1, by.y = "name", sort = FALSE)         # Put source ID's in network_data
  network_data <- merge(network_data, target_ID, by.x = 2, by.y = "name", sort = FALSE) # Put target ID's in network_data
  
  # Create link data for network
  network_link <- network_data[,c("ID.x", "ID.y")]
  colnames(network_link) <- c("Source", "Target")
  
  # Create node data for network
  if (contains_genes == TRUE){
    network_node <- data.frame(do.call(rbind,
                                       list(source_ID[,c("ID", "Symbol", "Group")],
                                            target_ID)))            # Create node data by stacking the two ID matrices on top of eachother
    network_node[,4] <- NA
    network_node[1:nrow(source_ID),4] <- source_ID[,"Genename"]     # Add gene names
    network_node[,5] <- NA
    network_node[1:nrow(source_ID),5] <- source_ID[,"name"]         # Add entrez gene IDs
    colnames(network_node) <- c("ID", "Name", "Group", "Genename", "EntrezID")  # Set correct column names
    network_node[,"Name"] <- as.character(network_node[,"Name"])    # Convert gene symbols to characters instead of factors
    network_node[is.na(network_node$Name), "Name"] <- network_node[is.na(network_node$Name), "EntrezID"]  # Replace missing gene symbols with entrez IDs
  } else {
    network_node <- data.frame(do.call(rbind,
                                       list(source_ID,
                                            target_ID)))            # Create node data by stacking the two ID matrices on top of eachother
    colnames(network_node) <- c("ID", "Name", "Group")              # Set correct column names
  }
  # Calculate node degrees
  network_node[,"nodeDegree"] <- 0
  for (i in 1:nrow(network_node)){
    id <- i - 1
    network_node[i,"nodeDegree"] <- length(network_link[network_link == id])
  }
  
  # Add value - this can be used to weight interactions
  network_link[,"Value"] <- 1
  
  return(list(network_link, network_node))
}

#========================================#
##      miRNA Interaction Network       ##
#========================================#
miRNA_mRNA_interaction_network <- function(mRNA, 
                                           miRNA, 
                                           filter = FALSE){
  
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
  miRNA_ID <- matrix(ncol = 5, nrow = length(IDs))
  miRNA_ID[,2] <- as.character(IDs)                               # Get all unqiue miRNA names
  miRNA_ID[,1] <- seq(0,length(IDs) - 1)
  miRNA_ID[,3] <- "miRNA"
  colnames(miRNA_ID) <- c("ID", "name", "Group", "Symbol", "Genename")
  
  # Create mRNA ID's
  IDs <- unique(network_data[,"mRNA"])
  mRNA_ID <- matrix(ncol = 3, nrow = length(IDs))
  mRNA_ID[,2] <- as.character(IDs)                               # Get all unqiue miRNA names
  mRNA_ID[,1] <- seq(nrow(miRNA_ID),nrow(miRNA_ID) + length(IDs) - 1)
  mRNA_ID[,3] <- "mRNA"
  colnames(mRNA_ID) <- c("ID", "name", "Group")
  
  
  network_data <- merge(network_data, miRNA_ID, by.x = "mirna", by.y = "name", sort = FALSE) # Put mirna ID's in network_data
  network_data <- merge(network_data, mRNA_ID, by.x = "mRNA", by.y = "name", sort = FALSE)   # Put mirna ID's in network_data
  
  # Create link data for network
  network_link <- network_data[,c("ID.x", "ID.y", "Occurance")]
  colnames(network_link) <- c("Source", "Target", "Value")
  
  #filter
  if(filter != FALSE){
    network_link_filtered <- network_link[network_link[,"Value"] > filter,]
  }
  
  # Create node data for network
  network_node <- data.frame(do.call(rbind, list(miRNA_ID, mRNA_ID))) # Create node data by stacking the two ID matrices on top of eachother
  colnames(network_node) <- c("ID", "Name", "Group")
  
  # Calculate node degrees
  network_node[,"nodeDegree"] <- 0
  for (i in 1:nrow(network_node)){
    id <- i - 1
    network_node[i,"nodeDegree"] <- length(network_link[network_link == id])
  }
  
  
  return(list(network_link, network_node))
}


#========================================#
##         Visualise the Network        ##
#========================================#
visualise_network <- function(network_link, 
                              network_node, 
                              silent = TRUE,
                              nodeSizeDiff = 0.5,
                              edgeWidth = 3,
                              nameVisiblity = 0,
                              width = NULL,
                              height = NULL,
                              bounded = FALSE,
                              force = -30,
                              linkDistance = 150,
                              saveName = "interaction_network.html"){
  
  # MyClickScript <- 'alert(Object.getOwnPropertyNames(d));' # use this to show all available node properties
  
  # Set JavaScript to display a message when clicking on a node
  MyClickScript <-   'alert("You clicked " + d.name + " which has a degree of: " + d.nodesize);'
  
  # Set how to calculate the node size
  nodeSize <- sprintf("Math.pow(d.nodesize, %f) + 4", nodeSizeDiff)
  
  # Set save directory
  saveName <- paste("../Figures/" ,saveName, sep = "")
  
  # Build the force directed network
  miRNA_interaction_network <- forceNetwork(Links   = network_link,
                                            Nodes   = network_node,
                                            Source  = "Source",
                                            Target  = "Target",
                                            NodeID  = "Name",
                                            Group   = "Group",
                                            Value   = "Value",
                                            Nodesize = "nodeDegree",
                                            radiusCalculation = nodeSize,
                                            opacity = 0.95,
                                            zoom    = TRUE,
                                            fontSize = 20,
                                            clickAction = MyClickScript,
                                            colourScale = JS("d3.scaleOrdinal(d3.schemeCategory10);"),
                                            legend = TRUE,
                                            linkWidth = edgeWidth,
                                            opacityNoHover = nameVisiblity,
                                            linkDistance = linkDistance,
                                            height = height,
                                            width = width,
                                            bounded = bounded,
                                            charge = force
                                            )
  # Save the network as html
  saveNetwork(miRNA_interaction_network, saveName)
  
  if(!silent){
    # Visualise the network in R
    print(miRNA_interaction_network)
    return(miRNA_interaction_network)
  }
  return(miRNA_interaction_network)
}