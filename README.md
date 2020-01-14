# Project_Period
 Cholestasis Systems BiologyProject

index:
1) dataFormatting.R
     - Is a function called format()
     - Collects Raw data and produces 4 objects:
      - mRNA is the expression profile of mRNA levels
      - miRNA is the expression profile of miRNA levels
       - labels contains all of the ids and names regarding samples, batches and treatments
      - A Key is generated connecting entrezID to the gene symbol
      
2) PCA.R
     - dataFormatting.R objects are collected using format()
     - plots the PCA of mRNA and miRNA expression and color codes the plots based on batches and treatment
     - It was identified that batch correctin needed to be implemented, and thus is implemented in PCA.R
     - plots for corrected expression levels were then created using PCA.R
     
3) Expression_Analysis.R
     - Calls objects using format()
     - Utilises limma to create model on data. 
     - Calculates top differentially expressed mRNA and miRNA
     - Plots a volcano plot of mRNA and miRNA data expression
     - Utilises qvalue to calculate q-values of expression levels. 
     - And finally identify significantly differentially expression mRNA and miRNA
     - Future:
      - Ideally should return significantly differentially expressed mRNA and miRNA to allow miRNA_Target_Analysis.R,     Network.R and Enrichment.R to use these values
       
4) miRNA_Target_Analysis.R
     - Calculates pearson correlations between mRNA and miRNA expression
     - Requires mRNA and/or miRNA expression values
     - Outputs correlation matrix of each possible mRNA-miRNA pairs and respective p-value

5) Enrichment.R
     - Notes: 
      a) https://www.biostars.org/p/222081/
      b) https://www.biostars.org/p/277303/
      c) https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-016-0971-3 (machine learning :O)
      d) https://bioconductor.org/packages/release/bioc/vignettes/topGO/inst/doc/topGO.pdf (topGO documentation)
      e) https://bioconductor.org/packages/devel/bioc/vignettes/GOexpress/inst/doc/GOexpress-UsersGuide.pdf (GOexpress documentation)
       
6) Network.R
     - hmmmmmmmm
      
