#---------------------------------------------------------------------------------------#
##                                       GENE ONTOLOGY
#---------------------------------------------------------------------------------------#
# WARNING: Run dataFormatting.R and  Expression_Analysis.R before running the gene ontology file. 

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("topGO")
BiocManager::install("org.Hs.eg.db")
BiocManager::install("miRNAtap")

library ("topGO")
library("org.Hs.eg.db")
library(ggplot2)
library(miRNAtap)

# Info form the ontology used: BP=biological processes, ID= entrez gene ID, mapping= the ontology to use. 
# org.Hs.eg.db = Genome wide annotation for Human, primarily based on mapping using Entrez Gene identifiers.
allGO2genes <- annFUN.org(whichOnto="BP", feasibleGenes=NULL, mapping="org.Hs.eg.db", ID="entrez")

# Input for GO: p values and name for the genes in a vector (all.genes).
all.genes<- DEG_CHvC$adj.p
names(all.genes)<- names_top_CHvC

# Create a topGOdata object that stores the info on the ontology. 
# ontology= BP type process. 
# all genes= genes to search (all.genes)
# annot= Info from the ontology. 
# GO2genes= Conversion used
# nodeSize= you select genes that have only nodeSize match on the GO. 
# (Small nodesize gives more nodes).

GOdata_CHvC <- new("topGOdata",
                   ontology="BP",
                   allGenes= all.genes,
                   annot=annFUN.GO2genes,
                   geneSel= function(x)x,
                   GO2genes=allGO2genes,
                   nodeSize=10)
GOdata_CHvC

results.ks <- runTest(GOdata_CHvC, algorithm="classic", statistic="ks")
goEnrichment <- GenTable(GOdata_CHvC, KS=results.ks, orderBy="KS", topNodes=20)
goEnrichment <- goEnrichment[goEnrichment$KS<0.05,]
goEnrichment <- goEnrichment[,c("GO.ID","Term","KS")]
goEnrichment$Term <- gsub(" [a-z]*\\.\\.\\.$", "", goEnrichment$Term)
goEnrichment$Term <- gsub("\\.\\.\\.$", "", goEnrichment$Term)
goEnrichment$Term <- paste(goEnrichment$GO.ID, goEnrichment$Term, sep=", ")
goEnrichment$Term <- factor(goEnrichment$Term, levels=rev(goEnrichment$Term))
goEnrichment$KS <- as.numeric(goEnrichment$KS)

ggplot(goEnrichment, aes(x=Term, y=-log10(KS))) +
  stat_summary(geom = "bar", fun.y = mean, position = "dodge") +
  xlab("Biological process") +
  ylab("Enrichment") +
  ggtitle("GO enrichment Cholestatic vs control 1") +
  scale_y_continuous(breaks = round(seq(0, max(-log10(goEnrichment$KS)), by = 2), 1)) +
  theme_bw(base_size=24) +
  theme(
    legend.position='none',
    legend.background=element_rect(),
    plot.title=element_text(angle=0, size=24, face="bold", vjust=1),
    axis.text.x=element_text(angle=0, size=18, face="bold", hjust=1.10),
    axis.text.y=element_text(angle=0, size=18, face="bold", vjust=0.5),
    axis.title=element_text(size=24, face="bold"),
    legend.key=element_blank(),     #removes the border
    legend.key.size=unit(1, "cm"),      #Sets overall area/size of the legend
    legend.text=element_text(size=18),  #Text size
    title=element_text(size=18)) +
  guides(colour=guide_legend(override.aes=list(size=2.5))) +
  coord_flip()

#printGraph(GOdata_CHvC, resultWeight, firstSigNodes = 5, fn.prefix = "tGO",
#              pdfSW = TRUE) # Not working


### WARNING NEED MULTIPLE TESTING CORRECTION 



# Mirna
mir = 'miR-10b'
predictions <- getPredictedTargets(mir, species = 'hsa',
                                   method = 'geom', min_src = 1)

# option 1 
library(miRNAtap)
library(topGO)
library(org.Hs.eg.db)

mir = 'miR-10b'
predictions = getPredictedTargets(mir, species = 'hsa',
                                  method = 'geom', min_src = 2)
# option 2 
BiocManager::install("multiMiR")
library(multiMiR)

targets.mirna_CHvC <- get_multimir(mirna=mirna.names_top_CHvC, summary=TRUE)
genes.mirna_CHvC<- targets.mirna_CHvC@data$target_entrez



all.genes<- targets.mirna_CHvC@summary$validated.sum
names(all.genes)<- genes.mirna_CHvC
allGO2genes <- annFUN.org(whichOnto="BP", feasibleGenes=NULL, mapping="org.Hs.eg.db", ID="entrez")


GOdata_mirna.CHvC <- new("topGOdata",
                   ontology="BP",
                   allGenes= all.genes,
                   annot=annFUN.GO2genes,
                   geneSel= function(x)x,
                   GO2genes=allGO2genes,
                   nodeSize=5)
GOdata_mirna.CHvC

results.ks <- runTest(GOdata_mirna.CHvC, algorithm="classic", statistic="ks")
goEnrichment <- GenTable(GOdata_mirna.CHvC, KS=results.ks, orderBy="KS", topNodes=20)
#goEnrichment <- goEnrichment[goEnrichment$KS<0.05,]
goEnrichment <- goEnrichment[,c("GO.ID","Term","KS")]
goEnrichment$Term <- gsub(" [a-z]*\\.\\.\\.$", "", goEnrichment$Term)
goEnrichment$Term <- gsub("\\.\\.\\.$", "", goEnrichment$Term)
goEnrichment$Term <- paste(goEnrichment$GO.ID, goEnrichment$Term, sep=", ")
goEnrichment$Term <- factor(goEnrichment$Term, levels=rev(goEnrichment$Term))
goEnrichment$KS <- as.numeric(goEnrichment$KS)

ggplot(goEnrichment, aes(x=Term, y=-log10(KS))) +
  stat_summary(geom = "bar", fun.y = mean, position = "dodge") +
  xlab("Biological process") +
  ylab("Enrichment") +
  ggtitle("GO enrichment Cholestatic vs control 1") +
  scale_y_continuous(breaks = round(seq(0, max(-log10(goEnrichment$KS)), by = 2), 1)) +
  theme_bw(base_size=24) +
  theme(
    legend.position='none',
    legend.background=element_rect(),
    plot.title=element_text(angle=0, size=24, face="bold", vjust=1),
    axis.text.x=element_text(angle=0, size=18, face="bold", hjust=1.10),
    axis.text.y=element_text(angle=0, size=18, face="bold", vjust=0.5),
    axis.title=element_text(size=24, face="bold"),
    legend.key=element_blank(),     #removes the border
    legend.key.size=unit(1, "cm"),      #Sets overall area/size of the legend
    legend.text=element_text(size=18),  #Text size
    title=element_text(size=18)) +
  guides(colour=guide_legend(override.aes=list(size=2.5))) +
  coord_flip()

#######################################
BiocManager::install("Rgraphviz") 
library(Rgraphviz)

sel.terms <- sample(usedGO(GOdata_CHvC), length(GOdata_CHvC@graph@edgeL))
ann.genes <- genesInTerm(GOdata_CHvC, usedGO(GOdata_CHvC)) ## get the annotations
head(ann.genes)

# showSigOfNodes(GOdata_CHvC, score(results.ks), firstSigNodes = 8, useInfo = 'def')
printGraph(GOdata_CHvC, results.ks, firstSigNodes = 20, fn.prefix = "tGO", useInfo = "all", pdfSW = TRUE)
