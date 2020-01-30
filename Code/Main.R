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

#===========================================================#
##           load all required packages functions          ##
#===========================================================#
library(rstudioapi)
current_path <- getActiveDocumentContext()$path 
setwd(dirname(current_path))          # Set working directory
source("package_Manager.R")
source("dataFormatting.R")
source("Expression_Analysis.R")
source("PCA.R")
source("miRNA_Target_Analysis.R")
source("interaction_network.R")
source("Gene_ontology.R")
source("Plots.R")

dir.create("../Figures", showWarnings = F) # Initialise folder for figures

#============================================#
##                   mRNA                   ##
#============================================#
#== Drain V Control (DVC) ==#
mRNA.DVC <- mRNA_DVC(save = TRUE)


#== Cholestatic V Control (CHVC) ==#
mRNA.CHVC <- mRNA_CHVC(save = TRUE)

#= GO - Biological Process =#
GO.mRNA.CHVC.BP <- get_GO(mRNA.CHVC, 
                          1, 
                          FALSE, 
                          "mRNA_CHVC", 
                          "mRNA: Cholestatic V Control (CHVC) \n Enriched Biological Processes", 
                          Onto = "BP",
                          save = TRUE,
                          "GO_mRNA_CHVC_BP.pdf",
                          "Biological Processes")
mRNA.CHVC.GO.Network.BP <- longFormat.Enrichment(GO.mRNA.CHVC.BP[[1]], mRNA.CHVC, GO.mRNA.CHVC.BP[[3]])
mRNA.GO.network.data.CHVC.BP <- generate_network_data(mRNA.CHVC.GO.Network.BP,
                                                      group1 = "GO: Biological Process",
                                                      group2 = "mRNA (entrez ID")

visualise_network(mRNA.GO.network.data.CHVC.BP[[1]], 
                  mRNA.GO.network.data.CHVC.BP[[2]],
                  edgeWidth = 1,
                  saveName = "GO biological processess network.html",
                  silent = FALSE)

#= GO - Molecular Function =#
GO.mRNA.CHVC.MF <- get_GO(mRNA.CHVC, 
                          1, 
                          FALSE, 
                          "mRNA_CHVC", 
                          "mRNA: Cholestatic V Control (CHVC) \n Enriched Molecular Functions", 
                          Onto = "MF",
                          save = TRUE,
                          "GO_mRNA_CHVC_MF.pdf",
                          "Molecular Functions")

mRNA.CHVC.GO.Network.MF <- longFormat.Enrichment(GO.mRNA.CHVC.MF[[1]], mRNA.CHVC, GO.mRNA.CHVC.MF[[3]])

mRNA.GO.network.data.CHVC.MF <- generate_network_data(mRNA.CHVC.GO.Network.MF,
                                                      group1 = "GO: Molecular Function",
                                                      group2 = "mRNA (entrez ID")

visualise_network(mRNA.GO.network.data.CHVC.MF[[1]], 
                  mRNA.GO.network.data.CHVC.MF[[2]],
                  edgeWidth = 1,
                  saveName = "GO molecular function network.html",
                  silent = FALSE)

#= GO - Cellular Components =#
GO.mRNA.CHVC.CC <- get_GO(mRNA.CHVC, 
                          1, 
                          FALSE, 
                          "mRNA_CHVC", 
                          "mRNA: Cholestatic V Control (CHVC) \n Enriched Cellular Components", 
                          Onto = "CC",
                          save = TRUE,
                          "GO_mRNA_CHVC_CC.pdf",
                          "Cellular Components")

mRNA.CHVC.GO.Network.CC <- longFormat.Enrichment(GO.mRNA.CHVC.CC[[1]], mRNA.CHVC, GO.mRNA.CHVC.CC[[3]])

mRNA.GO.network.data.CHVC.CC <- generate_network_data(mRNA.CHVC.GO.Network.CC,
                                                      group1 = "GO: Cellular Component",
                                                      group2 = "mRNA (entrez ID")
visualise_network(mRNA.GO.network.data.CHVC.CC[[1]], 
                  mRNA.GO.network.data.CHVC.CC[[2]],
                  edgeWidth = 1,
                  saveName = "GO cellular component network.html",
                  silent = FALSE)

#= Heatmap =#
heatmap(mrna, mRNA.CHVC[[1]], labels,
        "mRNA heatmap (CHVC)", 
        save = TRUE,
        "../Figures/Heatmap_mRNA_CHVC.pdf")


#= Cholestatic V Drain (CHVD) =#
mRNA.CHVD <- mRNA_CHVD(save = TRUE)
#= Heatmap =#
heatmap(mrna, 
        mRNA.CHVD[[1]],
        labels,
        "mRNA heatmap (CHVD)", 
        save = TRUE,
        "../Figures/Heatmap_mRNA_CHVD.pdf")

#== Plots ==#
#= PCA =#
PCA.mrna("PCA_mRNA_treatment.pdf", save = TRUE)
PCA.mrnaCorrected("PCA_mRNA_Corrected.pdf", save = TRUE)

#= Box Plots =#
mrna.boxPlot(save = TRUE, "../Figures/Boxplot_mRNA.pdf")

#= Anova =#
mRNA_Anova <- mrna.Anova()

#============================================#
##                  miRNA                   ##
#============================================#
#============================================#
##                  miRNA                   ##
#============================================#
#= Drain V Control (DVC) =#
miRNA.DVC <- miRNA_DVC(save = TRUE)

#=== Cholestatic V Control (CHVC) ===#
miRNA.CHVC <- miRNA_CHVC(save = TRUE)
mirna.CHVC.targets <- miRNA.targets(mRNA.CHVC[[1]], miRNA.CHVC[[1]], removeDuplicates = FALSE)

#= miRNA dendrogram =#
dendrogram(miRNA.CHVC[[1]], mirna.CHVC.targets, TRUE, "../Figures/Dendrogram_CHVC.pdf")

#= GO - Biological Process =#
GO.miRNA.CHVC.BP <- get_GO(miRNA.CHVC[[1]], 
                           mirna.CHVC.targets, 
                           use_target = TRUE, 
                           "miRNA_CHVC", 
                           "miRNA: Cholestatic V Control (CHVC) \n Enriched Biological Processes", 
                           Onto = "BP",
                           save = TRUE,
                           "GO_miRNA_CHVC_BP.pdf",
                           "Biological Processes"
                            )

#= GO - Molecular Function =#
GO.miRNA.CHVC.MF <- get_GO(miRNA.CHVC[[1]], 
                           mirna.CHVC.targets, 
                           TRUE, 
                           "miRNA_CHVC", 
                           "miRNA: Cholestatic V Control (CHVC) \n Enriched Molecular Functions", 
                           Onto = "MF",
                           save = TRUE,
                           "GO_miRNA_CHVC_MF.pdf",
                           "Molecular Functions"
                            )

#= GO - Cellular Component =#
GO.miRNA.CHVC.CC <- get_GO(miRNA.CHVC[[1]], 
                           mirna.CHVC.targets, 
                           TRUE, 
                           "miRNA_CHVC", 
                           "miRNA: Cholestatic V Control (CHVC) \n Enriched Cellular Components", 
                           Onto = "CC",
                           save = TRUE,
                           "GO_miRNA_CHVC_CC.pdf",
                           "Cellular Components"
                            )

#=== Cholestatic V Drain (CHVD) ===#
miRNA.CHVD <- miRNA_CHVD(save = TRUE)
mirna.CHVD.targets <- miRNA.targets(mRNA.CHVD[[1]], miRNA.CHVD[[1]])

#= miRNA dend«rogram =#
dendrogram(miRNA.CHVD[[1]], mirna.CHVD.targets, TRUE, "../Figures/Dendrogram_CHVD.pdf")

#== Plots ==#
#= PCA =#
PCA.mirna("PCA_miRNA.pdf", save = TRUE)
PCA.mirnaCorrected("PCA_miRNA_Corrected.pdf", save = TRUE)

#= Box Plots =#
mirna.boxPlot(save = TRUE, "../Figures/Boxplot_miRNA.pdf")

#= Anova =#
miRNA_Anova <- mirna.Anova()

#============================================#
##           Interaction Network            ##
#============================================#
mirna.CHVC.targets.unique <- miRNA.targets(mRNA.CHVC[[1]], miRNA.CHVC[[1]], removeDuplicates = TRUE)


#=== multiMiR and correlation filter network ===#
#= Compute the miRNA-mRNA correlations  =#
corMatrix <- miRNA_correlate(mirna, mrna,
                             only_DEG = TRUE,
                             miRNA_DEG = miRNA.CHVC[[2]],
                             mRNA_DEG = mRNA.CHVC[[2]],
                             mrna.treatmentOrder = mrna.treatmentOrder)

#= only keep interactions also present from the multiMiR query =#
mirna.CHVC.targets.unique_overlap <- miRNA_target_overlap(corMatrix, mirna.CHVC.targets.unique)

#= Generate link and node data for the network visualisation =#
target_overlap_network <- generate_network_data(mirna.CHVC.targets.unique_overlap,
                                                contains_genes = TRUE)

#= Visualise the multiMiR with correlation filter network =#
target_overlap_network_visualisation <- visualise_network(target_overlap_network[[1]],
                                                          target_overlap_network[[2]],
                                                          silent = FALSE,
                                                          nodeSizeDiff = 0.95,
                                                          nameVisiblity  = 0.9,
                                                          force = -30,
                                                          linkDistance = 120,
                                                          bounded = TRUE,
                                                          width = 1490,
                                                          height = 720,
                                                          saveName = "multiMiR and correlation filter network.html"
)

#=== Correlation based network ===#
#= load the miRNA-mRNA interactions based on the pearson correlations =#
correlation_data <- corMatrix[,c("mRNA", "miRNA", "p-value")]

#= Generate link and node data for the network visualisation =#
target_overlap_network <- generate_network_data(correlation_data,
                                                contains_genes = TRUE)

#= Visualise the correlation based network =#
correlation_network_visualisation <- visualise_network(target_overlap_network[[1]],
                                                       target_overlap_network[[2]],
                                                       silent = FALSE,
                                                       nodeSizeDiff = 0.7,
                                                       nameVisiblity  = 0,
                                                       force = -30,
                                                       linkDistance = 150,
                                                       bounded = TRUE,
                                                       width = 1490,
                                                       height = 840,
                                                       edgeWidth = 1,
                                                       saveName = "Correlation network.html")

##==== multiMiR query miRNA-mRNA Interaction Network ====##
# === This part is obselete and unfunctional == #
mirna_mrna.CHVC <- miRNA_mRNA_interaction_network(mRNA.CHVC[[1]], miRNA.CHVC[[1]])
mirna_mrna.CHVC_network <- visualise_network(mirna_mrna.CHVC[[1]], mirna_mrna.CHVC[[2]],
                                             silent = FALSE,
                                             edgeWidth = 1,
                                             width = 1920,
                                             height = 2080,
                                             bounded = FALSE,
                                             saveName = "miRNA-mRNA interaction network CHvC.html")

mirna_mrna.CHVD <- miRNA_mRNA_interaction_network(mRNA.CHVD[[1]], miRNA.CHVD[[1]])
visualise_network(mirna_mrna.CHVD[[1]], mirna_mrna.CHVD[[2]])

##==== GO-mRNA Interaction Network ====##






#============================================#
##       Interaction Network - Jip          ##
#============================================#


