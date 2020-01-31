#=============================================================================#
# Project Period, Liver cholestasis data analysis         							      #					  
#	Gene Ontology                                                               #
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

#===================================================#
##    function to perform gene ontology analysis    #
#===================================================# 
GO <- function(p,
               names, 
               Onto){
  allGO2genes <-annFUN.org(whichOnto = Onto, 
                           feasibleGenes = NULL, 
                           mapping = "org.Hs.eg.db", 
                           ID = "entrez")
  all <- p
  names(all)<-names
  
  GO <- new("topGOdata", 
            ontology = Onto,
            allGenes = all,
            annot = annFUN.GO2genes,
            geneSel = function(x)x,
            GO2genes = allGO2genes,
            nodeSize = 10)
  
  print(GO)
  return(GO)
}

#==============================================#
##    function generate gene ontology plots    #
#==============================================# 
plotGO <- function(GO, 
                   show, #number of topNodes to be shown
                   title,
                   save = TRUE,
                   main_Title,
                   lab
){
  
  results.ks <- runTest(GO, algorithm = "classic", statistic = "ks")
  goEnrichment <- GenTable(GO, KS = results.ks, orderBy = "KS", topNodes = show)
  enrichment<- data.frame(GO.id=goEnrichment$GO.ID, BP=goEnrichment$Term)
  goEnrichment <- goEnrichment[goEnrichment$KS < 0.05,]
  goEnrichment <- goEnrichment[,c("GO.ID","Term","KS")]
  
  goEnrichment$Term <- gsub(" [a-z]*\\.\\.\\.$", "", goEnrichment$Term)
  goEnrichment$Term <- gsub("\\.\\.\\.$", "", goEnrichment$Term)
  goEnrichment$Term <- paste(goEnrichment$GO.ID, goEnrichment$Term, sep = ", ")
  goEnrichment$Term <- factor(goEnrichment$Term, levels = rev(goEnrichment$Term))
  goEnrichment$KS <- as.numeric(goEnrichment$KS)
  
  goplot <- ggplot(goEnrichment, aes(x = Term, y = -log10(KS))) +
    stat_summary(geom = "bar", fun.y = mean, position = "dodge") +
    xlab(lab) +
    ylab("Enrichment (-Log10(pValue))") +
    ggtitle(title) +
    scale_y_continuous(breaks = round(seq(0, max(-log10(goEnrichment$KS)), by = 2), 1)) +
    theme_bw(base_size = 24) +
    theme(
      legend.position='none',
      legend.background=element_rect(),
      plot.title=element_text(angle = 0, size = 36, face="bold", vjust = 1),
      axis.text.x=element_text(angle = 0, size = 32, face="bold", hjust = 1.10),
      axis.text.y=element_text(angle = 0, size = 32, face = "bold", vjust = 0.5),
      axis.title=element_text(size = 40, face = "bold"),
      legend.key=element_blank(),     #removes the border
      legend.key.size=unit(1, "cm"),      #Sets overall area/size of the legend
      legend.text=element_text(size = 32),  #Text size
      title=element_text(size = 32)) +
    guides(colour=guide_legend(override.aes=list(size = 2.5))) +
    coord_flip()
  
  if(save){
    ggsave(main_Title, plot = goplot, path = "../Figures/", 
           width = 640, height = 480, units = "mm")
  }
  print(goplot)
  return(list(results.ks, enrichment))
}

GO.adj.mat <-function(GO, DEG){
  
  ann.genes <- genesInTerm(GO, usedGO(GO)) ## get the annotations
  data <- rownames(DEG[[1]])
  
  # Initialize adj matrix
  adj.matrix<- data.frame(matrix(0, nrow = length(data), ncol = length(usedGO(GO))))
  colnames(adj.matrix)<- usedGO(GO)
  rownames(adj.matrix)<- data
  
  for (a in 1:length(usedGO(GO))){
    genes<-ann.genes[[a]]
    for (b in 1:length(genes)){
      current<-genes[b]
      adj.matrix[which(rownames(adj.matrix)== current),a]<- 1}
  }
  adj.matrix<- data.frame(ids=rownames(adj.matrix),adj.matrix) # Only for testing purposes
  return(adj.matrix)
}

#=========================================#
##              get_GO                   ##
#=========================================#

get_GO <- function(mRNA, 
                   mirna.targets, 
                   use_target = TRUE, 
                   prefix, 
                   title,
                   Onto,
                   save = TRUE,
                   main_Title,
                   lab){
  if(use_target){ 
    #miRNA
    uniqueMirna.targets <- mirna.targets[!duplicated(mirna.targets$mRNA),]
    
    GO <- GO(uniqueMirna.targets$adj.p, 
             uniqueMirna.targets$mRNA, Onto)
    
    A <- plotGO(GO, 
                20, 
                title,
                save,
                main_Title,
                lab
    )
    result.ks <- A[[1]]
    enrichment <- A[[2]]
    
  } else {
    #mRNA
    GO <- GO(mRNA[[1]]$adj.p, 
             mRNA[[2]], Onto)
    
    A <- plotGO(GO, 
                20, 
                title,
                save,
                main_Title,
                lab
    )
    result.ks <- A[[1]]
    enrichment <- A[[2]]
  }
  
  printGraph(GO, 
             result.ks, 
             firstSigNodes = 20, 
             fn.prefix = prefix,
             useInfo = "all", 
             pdfSW = TRUE
  )
  return(list(GO, result.ks, enrichment))
}

#=========================================#
##      GO Graph Initialisation          ##
#=========================================#

longFormat.Enrichment <- function(GO, 
                                  mRNA,
                                  enrichment
){
  adj.mat <- GO.adj.mat(GO, mRNA) 
  data_long <- melt(adj.mat) # change data to format long
  colnames(data_long)<- c("ids","GO.id","values")
  data.long <- data_long[data_long[,3] != 0,]
  data.long <- subset(data.long, select=c("ids","GO.id"))
  
  enrichment$GO.id<-gsub(":", ".", enrichment$GO.id)
  output <- merge(data.long,enrichment, by="GO.id")
  
  return(output)
}