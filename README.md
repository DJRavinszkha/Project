# CorQBase

CorQBase is a R pipeline for the identification of miRNA-mRNA interactions for complex diseases. 

## main.R
     * The main pipeline that runs all of the files
     * Initialises lists using dataFormatting.R and identified differentially expressed (DE) miRNA and mRNA using Expression_Analysis.R
     * Performs Gene Ontology enrichment using gene_ontology.R
     * Visualises plots using PCA.R and Plots.R
     * Creates network using interaction_network.R
     
## package_Manager.R
     * The main pipeline that runs all of the files
        * edgeR and limma
        * multiMiR
        * networkD3 and Rgraphviz
        * topGO and org.Hs.eg.db
        * plyr, reshape2, formattable and dplyr
        * ggplot2, gplots, cowplot, ggdendro and dendextend

## dataFormatting.R
     * Contains a function called format()
     * Collects Raw data and produces 4 objects:
        * mRNA is the expression profile of mRNA levels
        * miRNA is the expression profile of miRNA levels
           * labels contains all of the ids and names regarding samples, batches and treatments
           * A Key is generated connecting entrezID to the gene symbol
           
## PCA.R
     * dataFormatting.R objects are collected using format()
     * plots the PCA of mRNA and miRNA expression and color codes the plots based on batches and treatment
     * It was identified that batch correction needed to be implemented, and thus is implemented in PCA.R
     * plots for corrected expression levels were then created using PCA.R
     
## Expression_Analysis.R
     * Calls lists using format()
     * Utilises limma to create model on data. 
     * Calculates top DE mRNA and miRNA
     * Plots a volcano plot of mRNA and miRNA data expression 
     * Finally identifies significantly DE mRNA and miRNA
     
## miRNA_Target_Analysis.R
     * miRNA_correlate():
        * Calculates pearson correlation between all possible pairs of miRNA and mRNA that you input
     * miRNA_target_query_targetHUB()
        * Uses the targetHUB API to query miRNA targets
    * miRNA.targets()
        * uses multiMiR which accesses miRecords, miRTarBase and TarBase databases to query miRNA targets
     * miRNA_target_overlap:
        * Find overlapping interactions between two sets of interactions

## Gene_ontology.R
    * Inputs:
       * DE mRNA and miRNA from Expression_Analysis.R
       * miRNA targets from miRNA_Target_Analysis.R      
    * Outputs:
       * Top Gene Ontology terms related to DE mRNA's and mRNA targets
       * Plots Gene Ontology network
    
## interaction_network.R
     * Generates and builds interaction network
     * Possible inputs:
        * miRNA-mRNA interactions
        * Top DE miRNAs and mRNAs from Expression_Analysis.R
     * Outputs: 
        * miRNA-mRNA interaction network 
     
## Plots.R
     * Creates heatmaps of gene expression and hierarchically clusters miRNAs
     * Inputs:
        * mRNA/miRNA expression data from dataFormatting.R
        * Top DE miRNAs/mRNAs from Expression_Analysis.R
        * miRNA targets from miRNA_Target_Analysis.R
     * Outputs: 
        * Heatplot of mRNA expression
        * Hierarchical cluster dendrogram of miRNA similarity
        
## Demographic_Analysis.R
     * Creates table containing confounders and demographic data

## Citing
Adriane Fosch-Mutané, Jip de Kok, Stefain Meier, Ravin Schmidl (2019). CorQBase: miRNA-mRNA interactions from complex disease expression profiles in R. https://github.com/DJRavinszkha/Project/
