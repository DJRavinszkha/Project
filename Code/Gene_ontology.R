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
  goEnrichment <- GenTable(GO, KS = results.ks, orderBy = "KS", topNodes = show)
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
  return(results.ks)
}

adj.mat <-function(GO, mRNA_DEG){
  
  ann.genes <- genesInTerm(GO, usedGO(GO)) ## get the annotations
  
  # Initialize adj matrix
  adj.matrix<- data.frame(matrix(0, nrow = length(mRNA_DEG[[1]]$mRNA), ncol = length(usedGO(GO.mrna.CHVC))))
  colnames(adj.matrix)<- usedGO(GO)
  rownames(adj.matrix)<- mRNA_DEG[[1]]$mRNA
  
  for (a in 1:length(usedGO(GO.mrna.CHVC))){
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

#showSigOfNodes(GOdata_CHvC, score(results.ks), firstSigNodes = 8, useInfo = 'def') #View network in panel
printGraph(GO.mrna.CHVC, result.ks.CHVC, firstSigNodes = 20, fn.prefix = "tGO", useInfo = "all", pdfSW = TRUE) # Save network in pdf

# Matrix for igraph 
adj.mat.GO_CHVC<-adj.mat(GO.mrna.CHVC,mRNA.CHVC)
data_long<-melt(adj.mat.GO_CHVC) # change data to format long
data_long_CHVC<- data_long[data_long[,3] != 0,]


##===== Cholestatic Vs Drained (CHVD) =====##
GO.mrna.CHVD <- GO(mRNA.CHVD[[1]]$adj.p, mRNA.CHVD[[2]])
result.ks.CHVD<-plotGO(GO.mrna.CHVD, 
                        20, 
                        "mRNA1: GO enrichment Cholestatic vs Drained")

printGraph(GO.mrna.CHVD, result.ks.CHVD, firstSigNodes = 20, fn.prefix = "tGO", useInfo = "all", pdfSW = TRUE) # Save network in pdf

adj.mat.GO_CHVD<-adj.mat(GO.mrna.CHVD,mRNA.CHVD)
data_long<-melt(adj.mat.GO_CHVD) # change data to format long
data_long_CHVD<- data_long[data_long[,3] != 0,]

#### 
# CODE FOR ADJ MATRIX WITHOUT MESSING THINGS UP 
# ann.genes <- genesInTerm(GO.mrna.CHVC, usedGO(GO.mrna.CHVC)) ## get the annotations
# 
# # Initialize adj matrix
# adj.mat.GO_CHVC<- data.frame(matrix(0, nrow = length(mRNA.CHVC[[1]]$mRNA), ncol = length(usedGO(GO.mrna.CHVC))))
# colnames(adj.mat.GO_CHVC)<- usedGO(GO.mrna.CHVC)
# rownames(adj.mat.GO_CHVC)<- mRNA.CHVC[[1]]$mRNA
# 
# for (a in 1:length(usedGO(GO.mrna.CHVC))){
#   genes<-ann.genes[[a]]
#   for (b in 1:length(genes)){
#     current<-genes[b]
#     adj.mat.GO_CHVC[which(rownames(adj.mat.GO_CHVC)== current),a]<- 1}
# }
# adj.mat.GO_CHVC<- data.frame(ids=rownames(adj.mat.GO_CHVC),adj.mat.GO_CHVC) # Only for testing purposes


#=========================================#
##               miRNA                   ##
#=========================================#

#= Cholestatic V Control (CHVC) =#
# Get targets
mirna.targets.CHVC <- get_multimir(mirna = miRNA.CHVC$mirna.Name, summary = TRUE)
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

# Adjacency matrix
adj.mat.mirna.GO_CHVC<-adj.mat(GO.mirna.CHVC,miRNA.CHVC)
data_long<-melt(adj.mat.GO_CHVD) # change data to format long
data_long_CHVD<- data_long[data_long[,3] != 0,]

#= Cholestatic V Drained (CHVD) =#
# Get targets
mirna.targets.CHVD <- get_multimir(mirna = miRNA.CHVD$mirna.Name, summary = TRUE)
targets.names.CHVD <- data.frame(mRNA = mirna.targets.CHVD@data$target_entrez, 
                                 mirna = mirna.targets.CHVD@data$mature_mirna_id,
                                 source = mirna.targets.CHVD@data$pubmed_id)

# Find which targets are present on our mRNA list
mirna.targets.DEG.CHVD <- merge(mRNA.CHVD[[1]], targets.names.CHVD, by = "mRNA", sort = FALSE)

#Find unique mRNAs
uniqueMirna.targets.DEG.CHVD <- mirna.targets.DEG.CHVD[!duplicated(mirna.targets.DEG.CHVD$mRNA),]

# GO enrichment for mRNA targets

GO.mirna.CHVD <- GO(uniqueMirna.targets.DEG.CHVD$adj.p, uniqueMirna.targets.DEG.CHVD$mRNA)
result.mirna.ks.CHVD <-plotGO(GO.mirna.CHVD, 
                              20, 
                            "Targets miRNA: GO enrichment Cholestatic vs Drained") 



####======####

# # Get genes corresponding to each GO term
# ann.genes <- genesInTerm(GO.mrna.CHVC, usedGO(GO.mrna.CHVC)) ## get the annotations
# 
# # Initialize adj matrix
# adj.mat.GO_CHVC<- data.frame(matrix(0, nrow = length(mRNA.CHVC[[1]]$mRNA), ncol = length(usedGO(GO.mrna.CHVC))))
# colnames(adj.mat.GO_CHVC)<- usedGO(GO.mrna.CHVC)
# rownames(adj.mat.GO_CHVC)<- mRNA.CHVC[[1]]$mRNA
# 
# for (a in 1:length(usedGO(GO.mrna.CHVC))){
#   genes<-ann.genes[[a]]
#   for (b in 1:length(genes)){
#   current<-genes[b]
#   adj.mat.GO_CHVC[which(rownames(adj.mat.GO_CHVC)== current),a]<- 1}
# }
# adj.mat.GO_CHVC<- data.frame(ids=rownames(adj.mat.GO_CHVC),adj.mat.GO_CHVC) # Only for testing purposes
# 
# 
# showSigOfNodes(GOdata_CHvC, score(results.ks), firstSigNodes = 8, useInfo = 'def')
# printGraph(GO.mrna.CHVC, result.ks.CHVC, firstSigNodes = 20, fn.prefix = "tGO", useInfo = "all", pdfSW = TRUE)

