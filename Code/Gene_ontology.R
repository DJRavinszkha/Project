#=============================================================================#
# Project Period, Liver cholestasis data analysis         							      #					  
#	Expression Analysis                                                         #
# Version: 1.0   															                                #
# Date: 9-1-2020											             	                          #
# Authors: Ariadna Fosch i Muntané, ID: I6215203, Maastricht University       #
#          Jip de Kok, ID: I6119367 , Maastricht University                   #
#          Ravin Schmidl, ID: I6125866, Maastricht University                 #
#          Stefan Meier, ID: I6114194 , Maastricht University                 #
#=============================================================================#

#=========================================#
##         Install libraries             ##
#=========================================#
BiocManager::install(c("multiMiR", "topGo", "org.Hs.eg.db", "Rgraphviz","edgeR"))
install.packages('padr')

library(qvalue)
library(limma)
library(multiMiR)
library(Rgraphviz)
library(topGO)
library(edgeR)
library(org.Hs.eg.db)
library(ggplot2)
library(rstudioapi)
library(reshape2)
library(plyr)

#=========================================#
##          Initialise Files             ##
#=========================================#

current_path <- getActiveDocumentContext()$path 
setwd(dirname(current_path))
source("dataFormatting.R")
source("Expression_Analysis.R")

mRNA.DVC <- mRNA.DVC()
mRNA.CHVC <- mRNA.CHVC()
mRNA.CHVD <- mRNA.CHVD()

miRNA.DVC <- miRNA.DVC()
miRNA.CHVC <- miRNA.CHVC()
miRNA.CHVD <- miRNA.CHVD()

#=========================================#
##              Functions                ##
#=========================================#

GO <- function(p, #mrna.DEG.CHVC$adj.P
               names){
  allGO2genes <-annFUN.org(whichOnto = "BP", 
                           feasibleGenes = NULL, 
                           mapping = "org.Hs.eg.db", 
                           ID = "entrez")
  
  all <- p
  names(all)<-names
  
    
  GO <- new("topGOdata", 
            ontology = "BP",
            allGenes = all,
            annot = annFUN.GO2genes,
            geneSel = function(x)x,
            GO2genes = allGO2genes,
            nodeSize = 10
            )

  print(GO)
  return(GO)
}

plotGO <- function(GO, 
                   show, #number of topNodes to be shown
                   title
                   ){
  
  results.ks <- runTest(GO, algorithm = "classic", statistic = "ks")
  goEnrichment <- GenTable(GO.mrna.CHVC, KS = results.ks, orderBy = "KS", topNodes = show)
  enrichment<- data.frame(GO.id=goEnrichment$GO.ID, BP=goEnrichment$Term)
  goEnrichment <- goEnrichment[goEnrichment$KS < 0.05,]
  goEnrichment <- goEnrichment[,c("GO.ID","Term","KS")]
 
  goEnrichment$Term <- gsub(" [a-z]*\\.\\.\\.$", "", goEnrichment$Term)
  goEnrichment$Term <- gsub("\\.\\.\\.$", "", goEnrichment$Term)
  goEnrichment$Term <- paste(goEnrichment$GO.ID, goEnrichment$Term, sep = ", ")
  goEnrichment$Term <- factor(goEnrichment$Term, levels = rev(goEnrichment$Term))
  goEnrichment$KS <- as.numeric(goEnrichment$KS)
 
  
  plot <- ggplot(goEnrichment, aes(x = Term, y = -log10(KS))) +
            stat_summary(geom = "bar", fun.y = mean, position = "dodge") +
            xlab("Biological process") +
            ylab("Enrichment") +
            ggtitle(title) +
            scale_y_continuous(breaks = round(seq(0, max(-log10(goEnrichment$KS)), by = 2), 1)) +
            theme_bw(base_size = 24) +
            theme(
              legend.position='none',
              legend.background=element_rect(),
              plot.title=element_text(angle = 0, size = 24, face="bold", vjust = 1),
              axis.text.x=element_text(angle = 0, size = 18, face="bold", hjust = 1.10),
              axis.text.y=element_text(angle = 0, size = 18, face = "bold", vjust = 0.5),
              axis.title=element_text(size = 24, face = "bold"),
              legend.key=element_blank(),     #removes the border
              legend.key.size=unit(1, "cm"),      #Sets overall area/size of the legend
              legend.text=element_text(size = 18),  #Text size
              title=element_text(size = 18)) +
            guides(colour=guide_legend(override.aes=list(size = 2.5))) +
            coord_flip()
  
  print(plot)
  return(list(results.ks, enrichment))
}

adj.mat <-function(GO, DEG, type){
 
  ann.genes <- genesInTerm(GO, usedGO(GO)) ## get the annotations
  data<- DEG[[1]]$mRNA
  
  # Initialize adj matrix
  adj.matrix<- data.frame(matrix(0, nrow = length(data), ncol = length(usedGO(GO.mrna.CHVC))))
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
##                mRNA                   ##
#=========================================#
##===== Cholestatic Vs Control (CHVC) =====##
GO.mrna.CHVC <- GO(mRNA.CHVC[[1]]$adj.p, mRNA.CHVC[[2]]) ## Error becauae havent initialised namesTop
result.ks.CHVC<-plotGO(GO.mrna.CHVC, 
                       20, 
                      "mRNA1: GO enrichment Cholestatic vs Control")
enrichment.CHVC <- result.ks.CHVC[[2]]

#showSigOfNodes(GOdata_CHvC, score(results.ks), firstSigNodes = 8, useInfo = 'def') #View network in panel
printGraph(GO.mrna.CHVC, result.ks.CHVC[[1]], firstSigNodes = 20, fn.prefix = "tGO", useInfo = "all", pdfSW = TRUE) # Save network in pdf

# Matrix for igraph 
adj.mat.GO_CHVC<-adj.mat(GO.mrna.CHVC,mRNA.CHVC) 
data_long<-melt(adj.mat.GO_CHVC) # change data to format long
colnames(data_long)<- c("ids","GO.id","values")
data.long_CHVC<- data_long[data_long[,3] != 0,]
data.long_CHVC<- subset(data.long_CHVC,select=c("ids","GO.id"))

enrichment.CHVC$GO.id<-gsub(":", ".", enrichment.CHVC$GO.id)
data.long_CHVC_enrich<- merge(data.long_CHVC,enrichment.CHVC, by="GO.id")

##===== Cholestatic Vs Drained (CHVD) =====##
GO.mrna.CHVD <- GO(mRNA.CHVD[[1]]$adj.p, mRNA.CHVD[[2]])
result.ks.CHVD<-plotGO(GO.mrna.CHVD, 
                        20, 
                        "mRNA1: GO enrichment Cholestatic vs Drained")
enrichment.CHVD <- result.ks.CHVD[[2]]

printGraph(GO.mrna.CHVD, result.ks.CHVD[[1]], firstSigNodes = 20, fn.prefix = "tGO", useInfo = "all", pdfSW = TRUE) # Save network in pdf

adj.mat.GO_CHVD<-adj.mat(GO.mrna.CHVD,mRNA.CHVD)
data_long<-melt(adj.mat.GO_CHVD) # change data to format long
colnames(data_long)<- c("ids","GO.id","values")
data.long_CHVD<- data_long[data_long[,3] != 0,]
data.long_CHVD<- subset(data.long_CHVD,select=c("ids","GO.id"))
enrichment.CHVD$GO.id<-gsub(":", ".", enrichment.CHVD$GO.id)
data.long_CHVD_enrich<- merge(data.long_CHVD,enrichment.CHVD, by="GO.id")

#=========================================#
##               miRNA                   ##
#=========================================#

#= Cholestatic V Control (CHVC) =#
# Get targets
mirna.targets.CHVC <- get_multimir(mirna = miRNA.CHVC[[1]]$mirna.Name, summary = TRUE)
targets.names.CHVC <- data.frame(mRNA = mirna.targets.CHVC@data$target_entrez, 
                                 mirna = mirna.targets.CHVC@data$mature_mirna_id,
                                 source = mirna.targets.CHVC@data$pubmed_id)

# Find which targets are present on our mRNA list
mirna.targets.DEG.CHVC <- merge(mRNA.CHVC[[1]], targets.names.CHVC, by = "mRNA", sort = FALSE)

#Find unique targets (Avoid double testing)
uniqueMirna.targets.DEG.CHVC <- mirna.targets.DEG.CHVC[!duplicated(mirna.targets.DEG.CHVC$mRNA),]

#GO enrichment for the mRNA targets
GO.mirna.CHVC <- GO(uniqueMirna.targets.DEG.CHVC$adj.p, uniqueMirna.targets.DEG.CHVC$mRNA)
result.mirna.ks.CHVC <-plotGO(GO.mirna.CHVC, 
                              20, 
                             "Targets miRNA: GO enrichment Cholestatic vs Control")
enrichment.mirna.CHVC <- result.mirna.ks.CHVC[[2]]

printGraph(GO.mirna.CHVC, result.mirna.ks.CHVC[[1]], firstSigNodes = 20, fn.prefix = "tGO", useInfo = "all", pdfSW = TRUE) # Save network in pdf

# Adjacency matrix
adj.mat.mirna.GO_CHVC<-adj.mat(GO.mirna.CHVC,mRNA.CHVC)
data_long<-melt(adj.mat.mirna.GO_CHVC) # change data to format long
colnames(data_long)<- c("ids","GO.id","values")
data.long.mirna_CHVC<- data_long[data_long[,3] != 0,]
data.long.mirna_CHVC<- subset(data.long.mirna_CHVC, select=c("ids","GO.id"))
enrichment.mirna.CHVC$GO.id<-gsub(":", ".", enrichment.mirna.CHVC$GO.id)
data.long.mirna_CHVC_enrich<- merge(data.long.mirna_CHVC,enrichment.mirna.CHVC, by="GO.id")

#= Cholestatic V Drained (CHVD) =#
# Get targets
mirna.targets.CHVD <- get_multimir(mirna = miRNA.CHVD[[1]]$mirna.Name, summary = TRUE)
targets.names.CHVD <- data.frame(mRNA = mirna.targets.CHVD@data$target_entrez, 
                                 mirna = mirna.targets.CHVD@data$mature_mirna_id,
                                 source = mirna.targets.CHVD@data$pubmed_id)

# Find which targets are present on our mRNA list
mirna.targets.DEG.CHVD <- merge(mRNA.CHVD[[1]], targets.names.CHVD, by = "mRNA", sort = FALSE)

#Find unique mRNAs
uniqueMirna.targets.DEG.CHVD <- mirna.targets.DEG.CHVD[!duplicated(mirna.targets.DEG.CHVD$mRNA),]

# GO enrichment for mRNA targets

GO.mirna.CHVD <- GO(uniqueMirna.targets.DEG.CHVD$adj.p, uniqueMirna.targets.DEG.CHVD$mRNA)

# Test not significant, therefrore we don't plot
#result.mirna.ks.CHVD <-plotGO(GO.mirna.CHVD, 
#                              20, 
#                            "Targets miRNA: GO enrichment Cholestatic vs Drained") 

# Adjacency matrix
adj.mat.mirna.GO_CHVD <- adj.mat(GO.mirna.CHVD,mRNA.CHVD)
data_long <- melt(adj.mat.mirna.GO_CHVD) # change data to format long
colnames(data_long)<- c("ids","GO.id","values")
data.long.mirna_CHVD <- data_long[data_long[,3] != 0,]
data.long.mirna_CHVD <- subset(data.long.mirna_CHVD, select=c("ids","variable"))
