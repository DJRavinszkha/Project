library(rstudioapi)
current_path <- getActiveDocumentContext()$path 
setwd(dirname(current_path))
source("dataFormatting.R")
source("Expression_Analysis.R")
source("PCA.R")
source("miRNA_Target_Analysis.R")
source("interaction_network.R")
source("Gene_ontology.R")
source("Plots.R")

if (!requireNamespace("rstudioapi", quietly = TRUE))
  install.packages("rstudioapi")
if (!requireNamespace("networkD3", quietly = TRUE))
  install.packages("networkD3")
if (!requireNamespace("multiMiR", quietly = TRUE))
  install.packages("multiMiR")
if (!requireNamespace("org.Hs.eg.db", quietly = TRUE))
  install.packages("org.Hs.eg.db")
if (!requireNamespace("Rgraphviz", quietly = TRUE))
  install.packages("Rgraphviz")
if (!requireNamespace("edgeR", quietly = TRUE))
  install.packages("edgeR")
if (!requireNamespace("limma", quietly = TRUE))
  install.packages("limma")
if (!requireNamespace("topGO", quietly = TRUE))
  install.packages("topGO")
if (!requireNamespace("ggplot2", quietly = TRUE))
  install.packages("ggplot2")
if (!requireNamespace("ggdendro", quietly = TRUE))
  install.packages("ggdendro")
if (!requireNamespace("plyr", quietly = TRUE))
  install.packages("plyr")
if (!requireNamespace("reshape2", quietly = TRUE))
  install.packages("reshape2")
if (!requireNamespace("gplots", quietly = TRUE))
  install.packages("gplots")
if (!requireNamespace("pcaMethods", quietly = TRUE))
  install.packages("pcaMethods")


library(networkD3)
library(limma)
library(multiMiR)
library(Rgraphviz)
library(topGO)
library(org.Hs.eg.db)
library(ggplot2)
library(rstudioapi)
library(ggdendro)
library(reshape2)
library(plyr)
library(gplots)
library(pcaMethods)

#============================================#
##                   mRNA                   ##
#============================================#
#= Drain V Control (DVC) =#
mRNA.DVC <- mRNA_DVC()


#= Cholestatic V Control (CHVC) =#
mRNA.CHVC <- mRNA_CHVC()
GO.mRNA.CHVC <- get_GO(mRNA.CHVC, 1, FALSE, "mRNA_CHVC")
mRNA.CHVC.GO.Network <- longFormat.Enrichment(GO.mRNA.CHVC[[1]], mRNA.CHVC, GO.mRNA.CHVC[[3]])

#= Cholestatic V Drain (CHVD) =#
mRNA.CHVD <- mRNA_CHVD()
GO.mRNA.CHVD <- get_GO(mRNA.CHVD, 1, FALSE, "mRNA_CHVC")
mRNA.CHVD.GO.Network <- longFormat.Enrichment(GO.mRNA.CHVD[[1]], mRNA.CHVD, GO.mRNA.CHVD[[3]])

#== Plots ==#
#= PCA =#
PCA.mrna()
PCA.mrnaCorrected()

#= Box Plots =#
mrna.boxPlot()

#= Heatmap =#
# CHVC #
heatmap(mrna, mRNA.CHVC[[1]], "mRNA heatmap (CHVC)")

# CHVD #
heatmap(mrna, mRNA.CHVD[[1]], "mRNA heatmap (CHVD)")

#= Anova =#
mRNA.Anova <- mrna.Anova()


#============================================#
##                  miRNA                   ##
#============================================#
#= Drain V Control (DVC) =#
miRNA.DVC <- miRNA_DVC()

#= Cholestatic V Control (CHVC) =#
miRNA.CHVC <- miRNA_CHVC()
mirna.CHVC.targets <- miRNA.targets(mRNA.CHVC[[1]], miRNA.CHVC[[1]])
GO.miRNA.CHVC <- get_GO(mRNA.CHVC, mirna.CHVC.targets, TRUE, "miRNA_CHVC")
miRNA.CHVC.GO.Network <- longFormat.Enrichment(GO.miRNA.CHVC[[1]], miRNA.CHVC, GO.miRNA.CHVC[[3]])
miRNA.GO.network.data <- generate_network_data(mRNA.CHVC.GO.Network)
link <- miRNA.GO.network.data[[1]]
node <- miRNA.GO.network.data[[2]]
visualise_network(link, node, silent = FALSE)

#= Cholestatic V Drain (CHVD) =#
miRNA.CHVD <- miRNA_CHVD()
mirna.CHVD.targets <- miRNA.targets(mRNA.CHVD[[1]], miRNA.CHVD[[1]])
GO.miRNA.CHVD <- get_GO(mRNA.CHVD, mirna.CHVD.targets, TRUE, "miRNA_CHVD")
miRNA.CHVD.GO.Network <- longFormat.Enrichment(GO.miRNA.CHVD[[1]], mRNA.CHVD, GO.miRNA.CHVD[[3]])

#== Plots ==#
#= PCA =#
PCA.mirna()
PCA.mirnaCorrected()

#= Box Plots =#
mirna.boxPlot()

#= Anova =#
miRNA.Anova <- mirna.Anova()

#============================================#
##           Interaction Network            ##
#============================================#

##==== miRNA-mRNA Interaction Network ====##

mirna_mrna.CHVC <- miRNA_mRNA_interaction_network(mRNA.CHVC[[1]], miRNA.CHVC[[1]])
visualise_network(mirna_mrna.CHVC[[1]], mirna_mrna.CHVC[[2]], silent = TRUE)

mirna_mrna.CHVD <- miRNA_mRNA_interaction_network(mRNA.CHVD[[1]], miRNA.CHVD[[1]])
visualise_network(mirna_mrna.CHVD[[1]], mirna_mrna.CHVD[[2]])

##==== GO-mRNA Interaction Network ====##






#============================================#
##       Interaction Network - Jip          ##
#============================================#

corMatrix <- miRNA_correlate(mirna, mrna,
                             only_DEG = TRUE,
                             miRNA_DEG = miRNA.CHVC[[2]],
                             mRNA_DEG = mRNA.CHVC[[2]],
                             mrna.treatmentOrder = mrna.treatmentOrder)
mirna.CHVC.targets_overlap <- miRNA_target_overlap(corMatrix, mirna.CHVC.targets)

target_overlap_network <- generate_network_data(mirna.CHVC.targets_overlap,
                                                contains_genes = TRUE)
target_overlap_network_visualisation <- visualise_network(target_overlap_network[[1]],
                                                          target_overlap_network[[2]],
                                                          silent = FALSE)

