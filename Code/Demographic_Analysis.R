#=============================================================================#
# Project Period, Liver cholestasis data analysis         							      #					  
#	Demographic Analysis                                                        #
# Version: 1.0   															                                #
# Date: 9-1-2020											             	                          #
# Authors: Ariadna Fosch i Muntané, ID: I6215203, Maastricht University       #
#          Jip de Kok, ID: I6119367 , Maastricht University                   #
#          Ravin Schmidl, ID: I6125866, Maastricht University                 #
#          Stefan Meier, ID: I6114194 , Maastricht University                 #
#=============================================================================#

metadata <- read.delim('../Data/metadata_mrna_adapt.csv', sep =",", check.names = FALSE)

metadata_chol <- metadata[metadata[,"Group"] == "Cholestatic",]
malecount_chol <- nrow(metadata_chol[metadata_chol[, "Gender"] == "M",])
femalecount_chol <- nrow(metadata_chol[metadata_chol[, "Gender"] == "F",])
gender_ratio_chol <- paste(malecount_chol,"/",femalecount_chol)
mean_age_chol <- signif(mean(metadata_chol$Age), digits = 3)
std_dev_age_chol <- signif(sd(metadata_chol$Age), digits = 3)
mean_bili_tot_chol <- signif(mean(metadata_chol$Bilirubintot, na.rm = TRUE), digits = 3)
std_dev_bili_tot_chol <- signif(sd(metadata_chol$Bilirubintot, na.rm = TRUE), digits = 3)
mean_bili_dir_chol <- signif(mean(metadata_chol$Bilirubindir, na.rm = TRUE), digits = 3)
std_dev_bili_dir_chol <- signif(sd(metadata_chol$Bilirubindir, na.rm = TRUE), digits = 3)
mean_bileacid_chol <- signif(mean(metadata_chol$Bileacid, na.rm = TRUE), digits = 3)
std_dev_bileacid_chol <- signif(sd(metadata_chol$Bileacid, na.rm = TRUE), digits = 3)


metadata_drain <- metadata[metadata[,"Group"] == "Drained",]
malecount_drain <- nrow(metadata_drain[metadata_drain[, "Gender"] == "M",])
femalecount_drain <- nrow(metadata_drain[metadata_drain[, "Gender"] == "F",])
gender_ratio_drain <- paste(malecount_drain,"/",femalecount_drain)
mean_age_drain <- signif(mean(metadata_drain$Age), digits = 3)
std_dev_age_drain <- signif(sd(metadata_drain$Age), digits = 3)
mean_bili_tot_drain <- signif(mean(metadata_drain$Bilirubintot, na.rm = TRUE), digits = 3)
std_dev_bili_tot_drain <- signif(sd(metadata_drain$Bilirubintot, na.rm = TRUE), digits = 3)
mean_bili_dir_drain <- signif(mean(metadata_drain$Bilirubindir, na.rm = TRUE), digits = 3)
std_dev_bili_dir_drain <- signif(sd(metadata_drain$Bilirubindir, na.rm = TRUE), digits = 3)
mean_bileacid_drain <- signif(mean(metadata_drain$Bileacid, na.rm = TRUE), digits = 3)
std_dev_bileacid_drain <- signif(sd(metadata_drain$Bileacid, na.rm = TRUE), digits = 3)
# be aware that there were many NAs in the drainaige set with regards to bilirubin levels and bile acid levels


metadata_control <- metadata[metadata[,"Group"] == "Control",]
malecount_control <- nrow(metadata_control[metadata_control[,"Gender"] == "M",])
femalecount_control <- nrow(metadata_control[metadata_control[,"Gender"] == "F",])
gender_ratio_control <- paste(malecount_control,"/",femalecount_control)
mean_age_control <- signif(mean(metadata_control$Age), digits = 3)
std_dev_age_control <- signif(sd(metadata_control$Age), digits = 3)
mean_bili_tot_control <- signif(mean(metadata_control$Bilirubintot, na.rm = TRUE), digits = 3)
std_dev_bili_tot_control <- signif(sd(metadata_control$Bilirubintot, na.rm = TRUE), digits = 3)
mean_bili_dir_control <- signif(mean(metadata_control$Bilirubindir, na.rm = TRUE), digits = 3)
std_dev_bili_dir_control <- signif(sd(metadata_control$Bilirubindir, na.rm = TRUE), digits = 3)
mean_bileacid_control <- signif(mean(metadata_control$Bileacid, na.rm = TRUE), digits = 3)
std_dev_bileacid_control <- signif(sd(metadata_control$Bileacid, na.rm = TRUE), digits = 3)
# be aware that there were many NAs in the drainaige set with regards to bilirubin levels and bile acid level

metadata_analysis <- metadata
colnames(metadata_analysis) <- c("Sample", "Gender", "Age", "Bilitot", "Bilidir", "Bileacid", "Group")

aov_age <- aov(Age~Group, data = metadata_analysis)
aov_age_p <- summary(aov_age)[[1]][1, 4:5]
aov_age_correct <- TukeyHSD(aov_age, ordered = FALSE, conf.level = 0.95)
CHvC_aov_age <- aov_age_correct$Group[1,4]
CHvD_aov_age <- aov_age_correct$Group[2,4]
DvC_aov_age <- aov_age_correct$Group[3,4]


TAB <- table(metadata_analysis$Gender, metadata_analysis$Group)
TAB_CHvC <- TAB[2:3, 2:3]
TAB_CHvD <- TAB[2:3, c(2,4)]
TAB_DvC <- TAB[2:3, 3:4]
TAB_all <- TAB[2:3, 2:4]
Fisher_gender_CHvC <- fisher.test(TAB_CHvC, conf.int = T, conf.level = 0.95)
Fisher_gender_CHvD <- fisher.test(TAB_CHvD, conf.int = T, conf.level = 0.95)
Fisher_gender_DvC <- fisher.test(TAB_DvC, conf.int = T, conf.level = 0.95)
Fisher_gender <- fisher.test(TAB_all, conf.int = T, conf.level = 0.95)
Fisher_gender_CHvC
Fisher_gender_CHvD
Fisher_gender_DvC
Fisher_gender # The probability of the observed array of cell frequencies plus the sum of the 
# probabilities of all other cell-frequency arrays (such as would be consistent with the observed marginal totals) 
# that are equal to or smaller than the probability of the observed array.

aov_bilitot <- aov(Bilitot~Group, data = metadata_analysis)
aov_bilitot_p <- summary(aov_bilitot)[[1]][1, 4:5]
aov_bilitot_correct <- TukeyHSD(aov_bilitot, ordered = FALSE, conf.level = 0.95)
CHvC_aov_bilitot <- aov_bilitot_correct$Group[1,4]
CHvD_aov_bilitot <- aov_bilitot_correct$Group[2,4]
DvC_aov_bilitot <- aov_bilitot_correct$Group[3,4]

aov_bilidir <- aov(Bilidir~Group, data = metadata_analysis)
aov_bilidir_p <- summary(aov_bilidir)[[1]][1, 4:5]
aov_bilidir_correct <- TukeyHSD(aov_bilidir, ordered = FALSE, conf.level = 0.95)
CHvC_aov_bilidir <- aov_bilidir_correct$Group[1,4]
CHvD_aov_bilidir <- aov_bilidir_correct$Group[2,4]
DvC_aov_bilidir <- aov_bilidir_correct$Group[3,4]

aov_bileacid <- aov(Bileacid~Group, data = metadata_analysis)
aov_bileacid_p <- summary(aov_bileacid)[[1]][1, 4:5]
aov_bileacid_correct<- TukeyHSD(aov_bileacid, ordered = FALSE, conf.level = 0.95)
CHvC_aov_bileacid <- aov_bileacid_correct$Group[1,4]
CHvD_aov_bileacid <- aov_bileacid_correct$Group[2,4]
DvC_aov_bileacid <- aov_bileacid_correct$Group[3,4]


demographic_matrix <- as.data.frame(matrix(ncol = 7, nrow = 5))
colnames(demographic_matrix) <- c("Cholestatic", "Drained", "Control", "Overall p-value", "Cholestatic vs. Control p-value",
                                  "Cholestatic vs. Drained p-value", "Drained vs. Control p-value")
rownames(demographic_matrix) <- c("Gender (M/F)", "Age", 
                                  "Bilirubin total in µmol/L", 
                                  "Bilirubin direct in µmol/L", 
                                  "Bile acids in µmol/L")

demographic_matrix[1,1] <- gender_ratio_chol
demographic_matrix[1,2] <- gender_ratio_drain
demographic_matrix[1,3] <- gender_ratio_control
demographic_matrix[1,4] <- base::format(Fisher_gender$p.value, digits = 2)
demographic_matrix[1,5] <- base::format(Fisher_gender_CHvC$p.value, digits = 2)
demographic_matrix[1,6] <- base::format(Fisher_gender_CHvD$p.value, digits = 2)
demographic_matrix[1,7] <- base::format(Fisher_gender_DvC$p.value, digits = 2)

demographic_matrix[2,1] <- paste(mean_age_chol,"±", std_dev_age_chol)
demographic_matrix[2,2] <- paste(mean_age_drain,"±", std_dev_age_drain)
demographic_matrix[2,3] <- paste(mean_age_control,"±", std_dev_age_control)
demographic_matrix[2,4] <- base::format(aov_age_p$`Pr(>F)`, digits = 2)
demographic_matrix[2,5] <- base::format(CHvC_aov_age, digits = 2)
demographic_matrix[2,6] <- base::format(CHvD_aov_age, digits = 2)
demographic_matrix[2,7] <- base::format(DvC_aov_age, digits = 2)

demographic_matrix[3,1] <- paste(mean_bili_tot_chol,"±", std_dev_bili_tot_chol)
demographic_matrix[3,2] <- paste(mean_bili_tot_drain,"±", std_dev_bili_tot_drain)
demographic_matrix[3,3] <- paste(mean_bili_tot_control,"±", std_dev_bili_tot_control)
demographic_matrix[3,4] <- base::format(aov_bilitot_p$`Pr(>F)`, digits = 2)
demographic_matrix[3,5] <- base::format(CHvC_aov_bilitot, digits = 2, scientific = TRUE)
demographic_matrix[3,6] <- base::format(CHvD_aov_bilitot, digits = 2)
demographic_matrix[3,7] <- base::format(DvC_aov_bilitot, digits = 2)

demographic_matrix[4,1] <- paste(mean_bili_dir_chol,"±", std_dev_bili_dir_chol)
demographic_matrix[4,2] <- paste(mean_bili_dir_drain,"±", std_dev_bili_dir_drain)
demographic_matrix[4,3] <- paste(mean_bili_dir_control,"±", std_dev_bili_dir_control)
demographic_matrix[4,4] <- base::format(aov_bilidir_p$`Pr(>F)`, digits = 2)
demographic_matrix[4,5] <- base::format(CHvC_aov_bilidir, digits = 2, scientific = TRUE)
demographic_matrix[4,6] <- base::format(CHvD_aov_bilidir, digits = 2, scientific = TRUE)
demographic_matrix[4,7] <- base::format(DvC_aov_bilidir, digits = 2)


demographic_matrix[5,1] <- paste(mean_bileacid_chol,"±", std_dev_bileacid_chol)
demographic_matrix[5,2] <- paste(mean_bileacid_drain,"±", std_dev_bileacid_drain)
demographic_matrix[5,3] <- paste(mean_bileacid_control,"±", std_dev_bileacid_control)
demographic_matrix[5,4] <- base::format(aov_bileacid_p$`Pr(>F)`, digits = 2)
demographic_matrix[5,5] <- base::format(CHvC_aov_bileacid, digits = 2)
demographic_matrix[5,6] <- base::format(CHvD_aov_bileacid, digits = 2)
demographic_matrix[5,7] <- base::format(DvC_aov_bileacid, digits = 2)

stars_overall <- stars.pval(as.numeric(demographic_matrix$`Overall p-value`))
stars_CHvC <- stars.pval(as.numeric(demographic_matrix$`Cholestatic vs. Control p-value`))
stars_CHvD <- stars.pval(as.numeric(demographic_matrix$`Cholestatic vs. Drained p-value`))
stars_DvC <- stars.pval(as.numeric(demographic_matrix$`Drained vs. Control p-value`))

demographic_matrix1 <- demographic_matrix

demographic_matrix1$`Overall p-value` <- paste(demographic_matrix$`Overall p-value`,stars_overall)
demographic_matrix1$`Cholestatic vs. Control p-value` <- paste(demographic_matrix$`Cholestatic vs. Control p-value`,stars_CHvC)
demographic_matrix1$`Cholestatic vs. Drained p-value` <- paste(demographic_matrix$`Cholestatic vs. Drained p-value`,stars_CHvD)
demographic_matrix1$`Drained vs. Control p-value` <- paste(demographic_matrix$`Drained vs. Control p-value`,stars_DvC)

formattable(demographic_matrix1, align = c("l", rep("r", NCOL(demographic_matrix1) - 1)))