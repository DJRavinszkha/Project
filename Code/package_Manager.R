#=============================================================================#
# Project Period, Liver cholestasis data analysis         							      #					  
#	Package Manager                                                             #
# Version: 1.0   															                                #
# Date: 24-1-2020											             	                          #
# Authors: Ariadna Fosch i Muntané, ID: I6215203, Maastricht University       #
#          Jip de Kok, ID: I6119367 , Maastricht University                   #
#          Ravin Schmidl, ID: I6125866, Maastricht University                 #
#          Stefan Meier, ID: I6114194 , Maastricht University                 #
#=============================================================================#

#=========================================#
##       Install missing packages        ##
#=========================================#

# Principal Component Analysis
if (!requireNamespace("pcaMethods", quietly = TRUE))
  install.packages("pcaMethods")

# Expression Analysis
if (!requireNamespace("edgeR", quietly = TRUE))
  install.packages("edgeR")
if (!requireNamespace("limma", quietly = TRUE))
  install.packages("limma")

#miRNA Targets
if (!requireNamespace("multiMiR", quietly = TRUE))
  install.packages("multiMiR")

#Network
if (!requireNamespace("networkD3", quietly = TRUE))
  install.packages("networkD3")
if (!requireNamespace("Rgraphviz", quietly = TRUE))
  install.packages("Rgraphviz")

#Ontology
if (!requireNamespace("topGO", quietly = TRUE))
  install.packages("topGO")
if (!requireNamespace("org.Hs.eg.db", quietly = TRUE))
  install.packages("org.Hs.eg.db")

#Formatting
if (!requireNamespace("plyr", quietly = TRUE))
  install.packages("plyr")
if (!requireNamespace("reshape2", quietly = TRUE))
  install.packages("reshape2")
if (!requireNamespace("formattable", quietly = TRUE))
  install.packages("formattable")
if (!requireNamespace("dplyr", quietly = TRUE))
  install.packages("dplyr")

#Graphing
if (!requireNamespace("ggplot2", quietly = TRUE))
  install.packages("ggplot2")
if (!requireNamespace("ggdendro", quietly = TRUE))
  install.packages("ggdendro")
if (!requireNamespace("cowplot", quietly = TRUE))
  install.packages("cowplot")
if (!requireNamespace("gplots", quietly = TRUE))
    install.packages("gplots")
if (!requireNamespace("dendextend", quietly = TRUE))
  install.packages("dendextend")

#=========================================#
##      Load all required packages       ##
#=========================================#

# Principal Component Analysis
library("pcaMethods")

#Expression Analysis
library(limma)
library(edgeR)

#Formatting
library(formattable)
library(dplyr)
library(reshape2)
library(plyr)

#Graphing
library(gplots)
library(ggplot2)
library(ggdendro)
library(cowplot)
library(dendextend)

#Networks
library(Rgraphviz)
library(networkD3)

#miRNA targets
library(multiMiR)

#Ontology
library(topGO)
library(org.Hs.eg.db)