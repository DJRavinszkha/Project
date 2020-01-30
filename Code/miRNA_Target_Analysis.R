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

#=========================================#
##             Correlation               ##
#=========================================#
miRNA_correlate <- function(miRNA, mRNA,
                            only_DEG = FALSE,
                            miRNA_DEG = NULL,
                            mRNA_DEG = NULL,
                            mrna.treatmentOrder = NULL){
  #========================================================================================================#
  # This function calculates pearson correlations between micro- and messenger RNA expression              #
  # You will input the RNA expression values for your group of interest (e.g. all cases)                   #
  # The output is a correlation matrix, showing the correlation between each possible pair of micro- and   #
  # messenger RNA and its significance value (p-value)                                                     #
  #========================================================================================================#
  
  if (only_DEG == TRUE){ # Filter only differential expressed genes and miRNA's
    miRNA_DEG <- data.frame(miRNA_DEG)
    miRNA <-merge(miRNA, miRNA_DEG, by.x = "row.names", by.y = "miRNA_DEG")
    rownames(miRNA) <- miRNA[,1]
    miRNA <- miRNA[,2:ncol(miRNA)]
    
    mRNA_DEG <- data.frame(mRNA_DEG)
    mRNA <-merge(mRNA, mRNA_DEG, by.x = "row.names", by.y = "mRNA_DEG")
    rownames(mRNA) <- mRNA[,1]
    mRNA <- mRNA[,2:ncol(mRNA)]
    
  }
  # Filter only the cases
  cases <- mrna.treatmentOrder[mrna.treatmentOrder$treatment.id == 1,]
  miRNA_case_index <- colnames(miRNA) %in% cases$sampleName
  miRNA <- miRNA[,miRNA_case_index]
  mRNA_case_index <- colnames(mRNA) %in% cases$sampleName
  mRNA <- mRNA[,mRNA_case_index]
  
  # ===
  numberOfPairs <- nrow(mRNA) * nrow(miRNA)                           # Calculate number of miRNA-mRNA pairs
  corMatrix <- matrix(nrow = numberOfPairs, ncol = 7)                 # Initialise corMatix which will store all correlations
  colnames(corMatrix) <- c("miRNA", "mRNA", "pearson-cor", "p-value", "conf.int L", "conf.int R","p-adjust")  # Set column names
  
  progress <- txtProgressBar(min = 0, max = numberOfPairs, initial = 0) # Initiate progress bar
  
  print("Calculating miRNA-mRNA interactions, depeding on your data this can take a while...")
  step = 0 # Initialise step which will count total number of correlations made
  for (i in 1:nrow(miRNA)){
    for (j in 1:nrow(mRNA)){
      step = step + 1
      
      # Calculate pearson correlation
      miRNA_val <- as.numeric(miRNA[i,])                   # Load miRNA expression values
      mRNA_val <- as.numeric(mRNA[j,])                     # Load mRNA expression values
      correlation <- stats::cor.test(miRNA_val,            # Calculate miRNA-mRNA correlation
                                     mRNA_val,
                                     method = "pearson",
                                     use = "complete.obs")
      
      # Store correlation results
      miRNA_name <- row.names(miRNA)[i]                     # Load miRNA name
      mRNA_name <- row.names(mRNA)[j]                       # Load mRNA name
      corMatrix[step,] <- c(miRNA_name,                     # store results in corMatrix
                            mRNA_name,
                            correlation$estimate[[1]],
                            correlation$p.value,
                            correlation$conf.int[1:2],
                            NA)
      
      
      setTxtProgressBar(progress, step)                     # Update progress bar
    }
  }
  # Set column types as numeric
  corMatrix_df <- as.data.frame(corMatrix)
  corMatrix_df[,3:7] <- as.numeric(corMatrix[,3:7])
  
  # Correct for multiple testing (using the FDR method of Benjamini Hochberg)
  corMatrix_df[,"p-adjust"] <- p.adjust(corMatrix_df[,"p-value"], method = "BH")
  
  # Only take significnat correlation - without multiple testing
  corMatrix_df_sign <- corMatrix_df[corMatrix_df[,"p-value"] <0.05,] # p<0.05
  corMatrix_df_sign <- corMatrix_df_sign[corMatrix_df_sign[,"pearson-cor"] < -0.7 | corMatrix_df_sign[,"pearson-cor"] > 0.7,] 
  
  return(corMatrix_df_sign)
}

#=========================================#
##            Query Targets              ##
#=========================================#
miRNA_target_query_targetHUB <- function(mirna, mrna){
  #========================================================================================================#
  # This function lookup miRNA targets from various databases using the targetHUB API                      #
  #                                                                                                        #
  # https://app1.bioinformatics.mdanderson.org/tarhub/_design/basic/index.html                             #
  #                                                                                                        #
  # miRBase: integrating microRNA annotation and deep-sequencing data.                                     #
  # Kozomara A, Griffiths-Jones S.                                                                         #
  # Nucleic Acids Res. 2011 39:D152-D157                                                                   #
  #========================================================================================================#
  
  #==== miRNA version conversion - Currently not used!!! ===#
  miRNA_names <- row.names(mirna)                                       # Store miRNA names
  version = checkMiRNAVersion(row.names(mirna), verbose = TRUE)         # Check miRNA name version
  miRNA_names <- miRNA_NameToAccession(miRNA_names, version = version)  # Include Accession numbers
  miRNA_names_v18 <- miRNA_AccessionToName(miRNA_names[,2],
                                           targetVersion = "v18")       # Convert to version 18
  
  #=== miRNA targer prediction - query ===#
  prefix <- "https://app1.bioinformatics.mdanderson.org/tarhub/_design/basic/_view/"
  link_temaplate <- paste(prefix, 'by_matureMIRcount?startkey=[%s,%s]&endkey=[%s,{}]',  sep = "")
  
  minSources = 1          # Set minimum number of sources
  mirna_targets <- list() # Initialise the mina_targets matrix
  
  progress <- txtProgressBar(min = 0, max = nrow(miRNA_names), initial = 0) # Initiate progress bar
  
  options(timeout = 100000) # Set timeout to increase time for retrying queries
  error_list = list()
  error_count = 0
  step = 0
  nTargets = 0
  for (i in 1:nrow(miRNA_names)){
    name = paste("%22", miRNA_names[i,1], "%22", sep = "")        # Set name in correct URL format
    miRNA_link <- sprintf(link_temaplate, name, minSources, name) # Create URL link to miRNA targets
    target <- try(fromJSON(miRNA_link))                           # Query the miRNA - try returns error for server errors
    
    if (typeof(target) == 'char'){
      error_list <- append(error_list, c(name, target))           # Capture error
      error_count <- error_count + 1                              # Keep track of number of errors
      next                                                        # Jump to next iteration
    }
    
    mirna_targets <- append(mirna_targets, list(target$rows))     # Add the miRNA targets to one long list of targets
    nTargets <- nTargets +  nrow(as.matrix(target$rows))          # Keep track of the total number of targets
    
    step = step + 1
    setTxtProgressBar(progress, step)                             # Update progress bar
  }
  print(sprintf("miRNA_target_query executed succesfully with %s errors", error_count))
  
  #=== The section below converts the mirna_targets list into one long data frame ===#
  progress <- txtProgressBar(min = 0, max = step, initial = 0)              # Initiate progress bar
  mirna_targets_df <- data.frame(matrix(nrow = nTargets, ncol = 3))         # Initialise mirna_targets_df
  colnames(mirna_targets_df) <- c("miRNA", "mRNA", "sources")               # set the column names
  index = 1                                                                 # Initialise index
  for (i in 1:step){
    targets <- data.frame(mirna_targets[[i]])           # load targets from mirna_targets
    if (length(targets) == 0){                          # If no targets, then skip to next iteration
      next
    }
    index_end <- index + nrow(targets) - 1              # Set end of index range
    mirna_targets_df[index:index_end,1:2] <- 
      t(data.frame(strsplit(targets[,1], ":")))         # Set miRNA and mRNA names
    mirna_targets_df[index:index_end,3] <- 
      as.numeric(unlist(targets[,2])[c(FALSE, TRUE)])   # Set number of sources
    index <- index_end + 1
    setTxtProgressBar(progress, i)                      # Update progress bar
  }
}

#================================================#
##       Query miRNA Targets with multiMiR      ##
#================================================#
miRNA.targets <- function(mRNA, miRNA, removeDuplicates = TRUE){
  mirna.targets <- get_multimir(mirna = miRNA$mirna.Name, summary = TRUE)
  targets.names <- data.frame(mRNA = mirna.targets@data$target_entrez, 
                                   mirna = mirna.targets@data$mature_mirna_id,
                                   source = mirna.targets@data$pubmed_id)
  
  mRNA$mRNA <- rownames(mRNA)
  targets.DEG <- merge(targets.names, mRNA, by = "mRNA", sort = FALSE)
  targets.DEG <- targets.DEG[!duplicated(
    targets.DEG[,c("mRNA", "mirna", "source")]),] # Remove duplicates from same publication
  
  if (removeDuplicates == TRUE){
    targets.DEG[,"source"] <- 1
    targets.DEG <- aggregate(source ~ mRNA + mirna, data = targets.DEG, sum)
  }
  
  return(targets.DEG)
}

#=========================================#
##     find overlapping interactions     ##
#=========================================#
miRNA_target_overlap <- function(corMatrix, interactions){
  interactions_overlap <- merge(interactions, corMatrix[,c("mRNA", "miRNA")],
                                by.x = c("mRNA", "mirna"),
                                by.y = c("mRNA", "miRNA"))
  return(interactions_overlap)
}
