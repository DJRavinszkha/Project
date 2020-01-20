library(rstudioapi)
current_path <- getActiveDocumentContext()$path 
setwd(dirname(current_path))
source("dataFormatting.R")
source("Expression_Analysis.R")
source("PCA.R")

#============================================#
##                   mRNA                   ##
#============================================#
#= Drain V Control (DVC) =#
mRNA.DVC <- mRNA.DVC()

#= Cholestatic V Control (CHVC) =#
mRNA.CHVC <- mRNA.CHVC()

#= Cholestatic V Drain (CHVD) =#
mRNA.CHVD <- mRNA.CHVD()

#= PCA =#
PCA.mrna()
PCA.mrnaCorrected()

#= Box Plots =#
mrna.boxPlot()

#= Anova =#
mRNA.Anova <- mrna.Anova()

#============================================#
##                  miRNA                   ##
#============================================#
#= Drain V Control (DVC) =#
miRNA.DVC <- miRNA.DVC()

#= Cholestatic V Control (CHVC) =#
miRNA.CHVC <- miRNA.CHVC()

#= Cholestatic V Drain (CHVD) =#
miRNA.CHVD <- miRNA.CHVD()

#= PCA =#
PCA.mirna()
PCA.mirnaCorrected()

#= Box Plots =#
mirna.boxPlot()

#= Anova =#
miRNA.Anova <- mirna.Anova()
