#=============================================================================#
# Project Period, Liver cholestasis data analysis         												  
#	
# Version: 1.0   															  
# Date: 9-1-2020											             	  
# Author:  Jip de Kok, Stefan Meier, Ariadna Fosch & Ravin Schmidl 
# Note: adapted from "Introduction to anamiR - 2019-05-09"
# Ref: Wang T, Lu T, Lee C, Lai L, Tsai M, Chuang E (2015).
#     "anamiR-an integrated analysis package of miRNA and mRNA expression." Manuscript submitted for publication.
#     https://bioconductor.org/packages/release/bioc/vignettes/anamiR/inst/doc/IntroductionToanamiR.html#phenotype-data
#=============================================================================#

process <- function(mrna, mirna, pheno.mrna, pheno.mirna){
  
  #====================================#
  ## Differential Expression Analysis ##
  #====================================#
  
  mrna_se <- SummarizedExperiment::SummarizedExperiment(
    assays = S4Vectors::SimpleList(counts=mrna),
    colData = pheno.mrna)
  
  mirna_se <- SummarizedExperiment::SummarizedExperiment(
    assays = S4Vectors::SimpleList(counts=mirna),
    colData = pheno.mirna)
  
    mrna_d <- differExp_discrete(se = mrna_se,
                               class = "ER", method = "wilcox.test",
                               t_test.var = FALSE, log2 = FALSE,
                               p_value.cutoff = 0.05,  logratio = 0.5
    )
    
    
    mirna_d <- differExp_discrete(se = mirna_se,
                                  class = "ER", method = "t.test",
                                  t_test.var = FALSE, log2 = FALSE,
                                  p_value.cutoff = 0.05,  logratio = 0.5
    )
    
    
    # Perform manual t-test
    index = rep(FALSE, length(colnames(mrna)))
    for(i in 1:length(colnames(mrna))){
      # if current column name is present in a list of all column names of cases:
      if(colnames(mrna)[i] %in% pheno.mrna[pheno.mrna[,3]=="case",1]){
        index[i] <- TRUE
      }
    }
    
    mrna.case <- mrna[,index]
    mrna.control <- mrna[,!index]
    t.results = t.test(mrna.case, mrna.control)
    
    
    #====================================#
    ## Correlation Analysis             ##
    #====================================#
    cor <- negative_cor(mrna_data = as.matrix(mrna_d), mirna_data = as.matrix(mirna_d),
                        method = c("pearson"), cut.off = -0.5)
    
 
    mirna_21 <- miR_converter(data = mirna_d, remove_old = TRUE,
                              original_version = 17, latest_version = 21)
    
    
}






