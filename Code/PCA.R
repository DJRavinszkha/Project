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

#=========================#
##    1. PCA Function    ## 
#=========================#

PCA <- function(data, 
                data.treatmentOrder, 
                treat, 
                data.batch, 
                title1, 
                title2, 
                main_Title, 
                save = TRUE
) {
  #= PCA =#
  pcaRes_data <- pca(t(data),nPcs = 10)  # perform PCA
  data.PCA <- data.frame(c(pcaRes_data@scores[,1]),
                         pcaRes_data@scores[,2],
                         pcaRes_data@scores[,3],
                         pcaRes_data@scores[,4],
                         pcaRes_data@scores[,5],
                         treat = data.treatmentOrder) 
  colnames(data.PCA) = c("PC1", "PC2", "PC3","PC4", "PC5","treat")
  
  if(treat != 0){ 
    #need this if loop because we either use a value generated in function for x (as in PCA for mRNA and miRNA)
    # - please refer to lines 16 and 58 in PCA_longversin.R
    # or we prescribe it a value (as in PCA for corrected mRNA)
    # - pleasw refer to line 101 in PCA_longversin.R
    Treatment = treat
  } else{
    Treatment = data.PCA$treat
  }
  ##==== Colouring by treatment ====##
  plot1 <- ggplot(data.PCA, aes(x = PC1, y = PC2)) +
    geom_point(aes(colour = Treatment), size = 5) +
    scale_colour_manual(values = c("#04179b", "#da9e00", "#198c19"),
                        aesthetics = "fill") + theme_light(base_size = 20) + 
    ggtitle("miRNA Corrected PCA Coloured by Treatment") 
  
  
  ##==== Colouring by batch ====##
  plot2 <- ggplot(data.PCA, aes(x = PC1, y = PC2)) +
    geom_point(aes(colour = data.batch), size = 5) +
    scale_colour_manual(values = c("#04179b", "#da9e00", "#198c19"),
                        aesthetics = "fill") + theme_light(base_size = 20) + 
    ggtitle("miRNA Corrected PCA Coloured by Batch")
  final_Plot <- plot_grid(plot1, plot2, labels = "AUTO")
  
  if(save){
    ggsave(main_Title, plot = final_Plot, path = "../Figures/", dpi = 300)
    
  }
}

#=========================#
## 2. PCA for mRNA DATA  ## 
#=========================#
PCA.mrna <- function(main_Title, save = TRUE){
  # Notes: 1. Batch_number_mrna changed to mrna.batches
  #        2. Batch_num changed to batch
  #        3. PCA_mrna changed to mrna.PCA
  #        4. Treatment changed to treatment
  #MRNA for 28 samples
  mrna[is.na(mrna)] <- NA # change Nan to NA
  plot <- PCA(mrna, #data
              mrna.treatmentOrder$treatment, #treatment order
              0, #treatment color
              mrna.batches$batch, #batch color
              "mRNA pca coloured by Treatment", #title1 
              "mRNA pca coloured by Batch",#title2
              main_Title, 
              save = save
  ) 
  
  return(plot)
}
#=========================#
## 3. PCA for miRNA DATA ##
#=========================#
# Notes: 1. Batch_number_mirna changed to mirna.batches
#        2. Batch_num changed to batch
#        3. PCA_mirna changed to mirna.PCA
#        4. Treatment changed to treatment

PCA.mirna <- function(main_Title, save = TRUE) {
  
  mirna[is.na(mirna)] <- NA # change Nan to NA
  plot <- PCA(mirna, #data
              mirna.treatmentOrder$treatment, #treatment order
              0, #treatment color
              mirna.batches$batch, #batch color
              "miRNA pca coloured by Treatment",#title1
              "miRNA pca coloured by Batch", #title2
              main_Title, 
              save = save
  )
  
  return(plot)
}

#=============================#
## 4. ANOVA (mrna and mirna) ##
#=============================#

# Anova testing mRNA batch variability
mrna.Anova <- function(){
  mrna.means <- data.frame(means= rowMeans(t(mrna),na.rm = TRUE))
  mrna.Anova <- aov(mrna.means[,1]~mrna.batches[,1])
  mrna.Tukey <- TukeyHSD(mrna.Anova, ordered = FALSE, conf.level = 0.95)
}

# Anova Testing miRna batch variability
mirna.Anova <- function(){
  mirna.means <- data.frame(means= rowMeans(t(mirna),na.rm = TRUE))
  mirna.Anova <- aov(mirna.means[,1]~mirna.batches[,1])
  mirna.Tukey <- TukeyHSD(mirna.Anova, ordered = FALSE, conf.level = 0.95)
}
#The anova shows that the mirna doesn't need correction as there is no significant batch variance

#==============================#
## 5. Batch effect correction ##
#==============================#

# Correction boxplot
mrna.boxPlot <- function(save = TRUE, title){
  mrna.corrected <- removeBatchEffect(mrna, batch = factor(mrna.batches$batch.id))
  boxplot(as.data.frame(mrna),main="Original")
  boxplot(as.data.frame(mrna.corrected),main="Batch corrected")
  
  if(save){
    pdf(title)
    boxplot(as.data.frame(mrna),main="Original")
    boxplot(as.data.frame(mrna.corrected),main="Batch corrected")
    dev.off()
  }
}
# remove batch effects mirna: NOT NECESSARY
mirna.boxPlot <- function(save = TRUE, title){
  mirna.corrected <- removeBatchEffect(mirna, batch = factor(mirna.batches$batch.id))
  boxplot(as.data.frame(mirna),main="Original miRNA")
  boxplot(as.data.frame(mirna.corrected),main="Batches miRNA corrected")
  
  if(save){
    pdf(title)
    boxplot(as.data.frame(mirna),main="Original miRNA")
    boxplot(as.data.frame(mirna.corrected),main="Batches miRNA corrected")
    dev.off()
  }
}
#=============================#
## 6. PCA for corrected data ##
#=============================#

##====  mRNA corrected  ====##
PCA.mrnaCorrected <- function(main_Title, save = TRUE){
  
  #=  mRNA corrected  =#
  mrna[is.na(mrna)] <- NA # change Nan to NA
  mrna.corrected <- removeBatchEffect(mrna, batch = factor(mrna.batches$batch.id)) 
  
  #= PCA =#
  plot <- PCA(mrna.corrected, #data
              mrna.batches$batch, #batch order
              mrna.treatmentOrder$treatment, #treatment color
              mrna.batches$batch, #batch color
              "mRNA corrected pca coloured by Treatment",  #title1
              "mRNA corrected pca coloured by Batch", #title2
              main_Title, 
              save = save
  )
  return(plot)
}
##====  miRNA corrected  ====##
PCA.mirnaCorrected <- function(main_Title, save = save){
  
  #=  miRNA corrected  =#
  mirna[is.na(mirna)] <- NA # change Nan to NA
  mirna.corrected <- removeBatchEffect(mirna, batch = factor(mirna.batches$batch.id)) 
  
  #= PCA =#
  plot <- PCA(mirna.corrected, #data
              mirna.batches$batch, #batch order
              mirna.treatmentOrder$treatment, #treatment color
              mirna.batches$batch, #batch color
              "miRNA corrected pca coloured by Treatment",  #title1
              "miRNA corrected pca coloured by Batch", #title2
              main_Title, 
              save = save
  )
  return(plot)
}