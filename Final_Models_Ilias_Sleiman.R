# Here lies the best ML models
# This script tests two different datasets considered as blind data, at each step of the predictions of each BEST model, 
# to check if the predictions make sense through Survival Plots. In reality, we don't know the true class of the samples from these TCGA datasets.
# First DS holds 308 samples of TCGA (Splited in half Dead / Alive), 
# and second DS holds 164 samples of TCGA (Splited in half TNBC / NON-TNBC)


library(caret)
library(clusterProfiler)
library(org.Hs.eg.db)
library(edgeR)
library(plyr)
library(skimr)
library(ROCR)
library(caretEnsemble)
library(ggplot2)
library(survival)
library(survminer)
library(dplyr)
library(SummarizedExperiment)

 

TCGA_raw_data <- readRDS("Full_SE_with_308samples_(Vital_Status).rds")      # 308 Blind samples - dead/alive equal
TCGA_raw_data_ilias <- readRDS("exp_whole_cohort_164s.rds")                 # 164 Classified samples - dead/alive NOT equal

a <- read.csv("C:/Users/sleim/OneDrive/Documents/TNBC_PROJECT/French Dataset Results/Gene Prioritization - BIM/GenePrioritization_ALL_ONTOLOGIES.csv")
sign <- as.vector(a$GeneSymbol)
# fortonoume metadata elenas
metadata <- read.csv("C:/Users/sleim/OneDrive/Documents/TNBC_PROJECT/French Dataset Results/French samples_Elenas_Work/samples_heatmap.csv", 
                     header = T,
                     row.names = 1)

y.exp <- edgeR::DGEList(counts = assays(TCGA_raw_data)[[1]])
#TMM normalization
y.exp <- edgeR::normLibSizes(y.exp, method = "TMM")
y.exp2 <- calcNormFactors(y.exp)
logCPM <- edgeR::cpm(y.exp2, prior.count=5, log=TRUE, normalized.lib.sizes = TRUE)

z_score_TCGA=scale(logCPM, center = TRUE, scale = TRUE) 
logCPM_TCGA <- as.data.frame(z_score_TCGA)

# Keep only the common genes in both data frames

logCPM_TCGA_sub <- logCPM_TCGA[sign,]
TCGA_data.t <- as.data.frame(t(logCPM_TCGA_sub))

# Add categorical Features
# Extract the columns from colData of TCGA_raw_data
IRE1_Activity <- TCGA_raw_data@colData@listData$act.level.IRE1
XBP1_Activity <- TCGA_raw_data@colData@listData$act.level.XBP1
RIDD_Activity <- TCGA_raw_data@colData@listData$act.level.RIDD

# Combine the columns with TCGA_data.t
TCGA_data_with_Activity <- cbind(TCGA_data.t, IRE1_Activity, XBP1_Activity, RIDD_Activity)
# Add Class Column as factor in the new df
TCGA_data_with_Activity$Class <- factor(NA)

# set.seed(139)
# 
# preProcValues_TCGA <- preProcess(TCGA_data_with_Activity, method = c("center", "scale"))
# 
# TCGA_data_Scaled <- predict(preProcValues_TCGA,TCGA_data_with_Activity)

# Create dummy variables for the categorical columns
dummies_TCGA <- dummyVars(Class ~ ., data = TCGA_data_with_Activity, levelsOnly = T)

# Transform the data using the dummy variables (DROP CLASS!!!!!!!)
TCGA_data_Activity_EncodedX <- as.data.frame(predict(dummies_TCGA, newdata = TCGA_data_with_Activity))

# Bind "Cancer Abr" column back into combined_df_encoded
TCGA_data_Activity_Encoded <- cbind(TCGA_data_with_Activity$Class, TCGA_data_Activity_EncodedX)
names(TCGA_data_Activity_Encoded)[1] <- "Class" 

# Ftiaxno ligo ta names gt exoun peiraxtei
names(TCGA_data_Activity_Encoded)[2:239] <- names(TCGA_data_with_Activity)[1:238]
names(TCGA_data_Activity_Encoded)[240] <- "IRE1_High_activity_(XBP1+/RIDD+)"
names(TCGA_data_Activity_Encoded)[241] <- "IRE1_Low_activity_(XBP1-/RIDD-)"
names(TCGA_data_Activity_Encoded)[242] = "IRE1_Mid_activity_(XBP1-/RIDD+)"
names(TCGA_data_Activity_Encoded)[243] = "IRE1_Mid_activity_(XBP1+/RIDD-)"
names(TCGA_data_Activity_Encoded)[244] = "XBP1_High_activity"
names(TCGA_data_Activity_Encoded)[245] = "XBP1_Low_activity"
names(TCGA_data_Activity_Encoded)[246] = "RIDD_High_activity"
names(TCGA_data_Activity_Encoded)[247] = "RIDD_Low_activity"

# dokimi montelou
# Model 5 (palio model.comb.50 / 50 features)  RF  -- Vgazei to kalytero survival

set.seed(139)
Model_RF_50 <- readRDS("Model.Top50(rf me TrCtrl + TnGr).rds")
predictions_50 <- predict(Model_RF_50, newdata = TCGA_data_Activity_Encoded)

TCGA_class_counts_50 <- table(predictions_50)

# Print the counts for TNBC and Non-TNBC
print(paste("TNBC:", TCGA_class_counts_50["TNBC"]))
print(paste("Non.TNBC:", TCGA_class_counts_50["Non.TNBC"]))

TCGA_raw_data@colData@listData$Class <- as.character(predictions_50)

function.SurvPlot <- function(exp, clusterCol, filename, title, colors){
  
  TCGAbiolinks::TCGAanalyze_survival(data = SummarizedExperiment::colData(exp),#data = surv.data,
                                     clusterCol = clusterCol,
                                     filename = filename,
                                     height = 8,
                                     width = 12,
                                     dpi = 100,
                                     legend = "Survival Plot",
                                     main = "Kaplan-Meier Overall Survival Curves",
                                     pvalue = TRUE,
                                     risk.table = TRUE,
                                     conf.int= FALSE,
                                     title = title,
                                     color = colors)
}

my_colors <- c("#008000","red")
function.SurvPlot(TCGA_raw_data, "Class", "French Dataset Results/Survival Plots/SurvivalPlot_for_308S_Model_top_50.pdf",
                  "K-M Overall Survival Curves / TCGA-BRCA\n308 TCGA Samples", my_colors)


# Model RF 30 top features  (vgazei akrivos to idio me to TOP 50)
set.seed(139)
Model_RF_30 <- readRDS("Model_RF.FullParam_Top30.rds")
predictions_30 <- predict(Model_RF_30, newdata = TCGA_data_Activity_Encoded)

TCGA_class_counts_30 <- table(predictions_30)

# Print the counts for TNBC and Non-TNBC
print(paste("TNBC:", TCGA_class_counts_30["TNBC"]))
print(paste("Non.TNBC:", TCGA_class_counts_30["Non.TNBC"]))

TCGA_raw_data@colData@listData$Class <- as.character(predictions_30)

function.SurvPlot(TCGA_raw_data, "Class", "French Dataset Results/Survival Plots/SurvivalPlot_for_308S_Model_top_30.pdf",
                  "K-M Overall Survival Curves / TCGA-BRCA\n308 TCGA Samples", my_colors)

#____________________________________________________________________________________________________________________________

# Model 7 (palio model_glmnet.2, tuneLength + trCntr  (BEST SCORE SO FAR))
set.seed(139)
Model_GLM <- readRDS("Model_glm.2(trCtr + TnL.rds")
predictions_GLM <- predict(Model_GLM, newdata = TCGA_data_Activity_Encoded)

TCGA_class_counts_GLM <- table(predictions_GLM)

# Print the counts for TNBC and Non-TNBC
print(paste("TNBC:", TCGA_class_counts_GLM["TNBC"]))
print(paste("Non.TNBC:", TCGA_class_counts_GLM["Non.TNBC"]))

TCGA_raw_data@colData@listData$Class <- as.character(predictions_GLM)

function.SurvPlot(TCGA_raw_data, "Class", "French Dataset Results/Survival Plots/SurvivalPlot_for_308S_Model_GLM.pdf",
                  "K-M Overall Survival Curves / TCGA-BRCA\n308 TCGA Samples", my_colors)

#_________________________________________________________________________________________________________________

# Model 10 (glm top 30, tuneL + trCntr ) 

set.seed(139)
Model_GLM_top30 <- readRDS("Model_glmnet.3(glmnet_me_TrCtrl_&_TnLength).rds")
predictions_GLM_30 <- predict(Model_GLM_top30, newdata = TCGA_data_Activity_Encoded)

TCGA_class_counts_GLM_30 <- table(predictions_GLM_30)

# Print the counts for TNBC and Non-TNBC
print(paste("TNBC:", TCGA_class_counts_GLM_30["TNBC"]))
print(paste("Non.TNBC:", TCGA_class_counts_GLM_30["Non.TNBC"]))

TCGA_raw_data@colData@listData$Class <- as.character(predictions_GLM_30)

function.SurvPlot(TCGA_raw_data, "Class", "French Dataset Results/Survival Plots/SurvivalPlot_for_308S_Model_GLM_top30.pdf")

#_________________________________________________________________________________________________


# Model 4 (palio model.comb.TOP.2 / 150 features)
set.seed(139)
Model_RF_150 <- readRDS("Model.Top150 ( RF me TnGr + trCntr.RDS")
predictions_RF_150 <- predict(Model_RF_150, newdata = TCGA_data_Activity_Encoded)

TCGA_class_counts_RF_150 <- table(predictions_RF_150)

# Print the counts for TNBC and Non-TNBC
print(paste("TNBC:", TCGA_class_counts_RF_150["TNBC"]))
print(paste("Non.TNBC:", TCGA_class_counts_RF_150["Non.TNBC"]))

TCGA_raw_data@colData@listData$Class <- as.character(predictions_RF_150)

function.SurvPlot(TCGA_raw_data, "Class", "French Dataset Results/Survival Plots/SurvivalPlot_for_308S_Model_RF_top150.pdf")

#_________________________________________________________________________________________________________________

# Model 3 (palio model.comb) (kako)

set.seed(139)
Model_RF <- readRDS("Model.Comb(rf me TrCtrl + TnGr).rds")
predictions_TCGA <- predict(Model_RF, newdata = TCGA_data_Activity_Encoded)

TCGA_class_counts <- table(predictions_TCGA)

# Print the counts for TNBC and Non-TNBC
print(paste("TNBC:", TCGA_class_counts["TNBC"]))
print(paste("Non.TNBC:", TCGA_class_counts["Non.TNBC"]))

TCGA_raw_data@colData@listData$Class <- as.character(predictions_TCGA)

function.SurvPlot(TCGA_raw_data, "Class", "French Dataset Results/Survival Plots/SurvivalPlot_for_308S_Model_RF.pdf")
#_________________________________________________________________________________

# Model glm me top 50

set.seed(139)
model_glmnet_50 <- readRDS("Model_glmnet.4(glmnet_me_TrCtrl_&_TnLength_top_50).rds")
predictions_TCGA_glm_50 <- predict(model_glmnet_50, newdata = TCGA_data_Activity_Encoded)

TCGA_class_counts_glm_50 <- table(predictions_TCGA_glm_50)

# Print the counts for TNBC and Non-TNBC
print(paste("TNBC:", TCGA_class_counts_glm_50["TNBC"]))
print(paste("Non.TNBC:", TCGA_class_counts_glm_50["Non.TNBC"]))

TCGA_raw_data@colData@listData$Class <- as.character(predictions_TCGA_glm_50)

function.SurvPlot(TCGA_raw_data, "Class", "French Dataset Results/Survival Plots/SurvivalPlot_for_308S_Model_glm_top_50.pdf",
                  "K-M Overall Survival Curves / TCGA-BRCA\n308 TCGA Samples", my_colors)

###################################################################################################################
################################################################################################################


# Dokimi tou kalyterou modelou sto megalytero SE me olous tous karkinous

exp_BRCA <- readRDS("exp_BRCA_FINAL.rds")


y.BRCA <- edgeR::DGEList(counts = assays(exp_BRCA)[[1]])
#TMM normalization
y.BRCA <- edgeR::normLibSizes(y.BRCA, method = "TMM")
y.eBRCA.2 <- calcNormFactors(y.BRCA)
logCPM_BRCA <- edgeR::cpm(y.eBRCA.2, prior.count=5, log=TRUE, normalized.lib.sizes = TRUE)

z_score_BRCA=scale(logCPM_BRCA, center = TRUE, scale = TRUE) 
logCPM_BRCA_scaled <- as.data.frame(z_score_BRCA)

# Keep only the common genes in both data frames

logCPM_BRCA_scaled_sub <- logCPM_BRCA_scaled[sign,]
BRCA_data.t <- as.data.frame(t(logCPM_BRCA_scaled_sub))

# Add categorical Features
# Extract the columns from colData of TCGA_raw_data
IRE1_Activity <- exp_BRCA@colData@listData$act.level.IRE1
XBP1_Activity <- exp_BRCA@colData@listData$act.level.XBP1
RIDD_Activity <- exp_BRCA@colData@listData$act.level.RIDD

# Combine the columns with TCGA_data.t
BRCA_data_with_Activity <- cbind(BRCA_data.t, IRE1_Activity, XBP1_Activity, RIDD_Activity)
# Add Class Column as factor in the new df
BRCA_data_with_Activity$Class <- factor(NA)


# Create dummy variables for the categorical columns
dummies_BRCA <- dummyVars(Class ~ ., data = BRCA_data_with_Activity, levelsOnly = T)

# Transform the data using the dummy variables (DROP CLASS!!!!!!!)
BRCA_data_Activity_EncodedX <- as.data.frame(predict(dummies_BRCA, newdata = BRCA_data_with_Activity))

# Bind "Cancer Abr" column back into combined_df_encoded
BRCA_data_Activity_Encoded <- cbind(BRCA_data_with_Activity$Class, BRCA_data_Activity_EncodedX)
names(BRCA_data_Activity_Encoded)[1] <- "Class" 

# Ftiaxno ligo ta names gt exoun peiraxtei
names(BRCA_data_Activity_Encoded)[2:239] <- names(BRCA_data_with_Activity)[1:238]
names(BRCA_data_Activity_Encoded)[240] <- "IRE1_High_activity_(XBP1+/RIDD+)"
names(BRCA_data_Activity_Encoded)[241] <- "IRE1_Low_activity_(XBP1-/RIDD-)"
names(BRCA_data_Activity_Encoded)[242] = "IRE1_Mid_activity_(XBP1-/RIDD+)"
names(BRCA_data_Activity_Encoded)[243] = "IRE1_Mid_activity_(XBP1+/RIDD-)"
names(BRCA_data_Activity_Encoded)[244] = "XBP1_High_activity"
names(BRCA_data_Activity_Encoded)[245] = "XBP1_Low_activity"
names(BRCA_data_Activity_Encoded)[246] = "RIDD_High_activity"
names(BRCA_data_Activity_Encoded)[247] = "RIDD_Low_activity"


# dokimi montelou
# Model 5 (palio model.comb.50 / 50 features)  RF  -- Vgazei to kalytero survival K 8A TO XRISIMOPOIISO

set.seed(139)
# Model_RF_50 <- readRDS("Model.Top50(rf me TrCtrl + TnGr).rds")
predictions_BRCA <- predict(Model_RF_50, newdata = BRCA_data_Activity_Encoded)

BRCA_class_counts <- table(predictions_BRCA)

# Print the counts for TNBC and Non-TNBC
print(paste("TNBC:", BRCA_class_counts["TNBC"]))
print(paste("Non.TNBC:", BRCA_class_counts["Non.TNBC"]))

exp_BRCA@colData@listData$Class <- as.character(predictions_BRCA)

function.SurvPlot(exp_BRCA, "Class", "French Dataset Results/Survival Plots/SurvivalPlot_for_BRCA_data_Model_top_50.pdf",
                  "K-M Overall Survival Curves / TCGA-BRCA\n308 TCGA Samples", my_colors)


##########################################################################################################################
###### Dokimi RF top50 me arxiko scale, sta 164 deigmata (ΧΑΛΙΑ)

y.exp <- edgeR::DGEList(counts = assays(TCGA_raw_data_ilias)[[1]])
#TMM normalization
y.exp <- edgeR::normLibSizes(y.exp, method = "TMM")
y.exp2 <- calcNormFactors(y.exp)
logCPM <- edgeR::cpm(y.exp2, prior.count=5, log=TRUE, normalized.lib.sizes = TRUE)

z_score_TCGA=scale(logCPM, center = TRUE, scale = TRUE) 
logCPM_TCGA <- as.data.frame(z_score_TCGA)

# Keep only the common genes in both data frames

logCPM_TCGA_sub <- logCPM_TCGA[sign,]
TCGA_data.t <- as.data.frame(t(logCPM_TCGA_sub))

# Add categorical Features
# Extract the columns from colData of TCGA_raw_data
IRE1_Activity <- TCGA_raw_data_ilias@colData@listData$act.level.IRE1
XBP1_Activity <- TCGA_raw_data_ilias@colData@listData$act.level.XBP1
RIDD_Activity <- TCGA_raw_data_ilias@colData@listData$act.level.RIDD

# Combine the columns with TCGA_data.t
TCGA_data_with_Activity <- cbind(TCGA_data.t, IRE1_Activity, XBP1_Activity, RIDD_Activity)
# Add Class Column as factor in the new df
TCGA_data_with_Activity$Class <- factor(NA)

# set.seed(139)
# 
# preProcValues_TCGA <- preProcess(TCGA_data_with_Activity, method = c("center", "scale"))
# 
# TCGA_data_Scaled <- predict(preProcValues_TCGA,TCGA_data_with_Activity)

# Create dummy variables for the categorical columns
dummies_TCGA <- dummyVars(Class ~ ., data = TCGA_data_with_Activity, levelsOnly = T)

# Transform the data using the dummy variables (DROP CLASS!!!!!!!)
TCGA_data_Activity_EncodedX <- as.data.frame(predict(dummies_TCGA, newdata = TCGA_data_with_Activity))

# Bind "Cancer Abr" column back into combined_df_encoded
TCGA_data_Activity_Encoded <- cbind(TCGA_data_with_Activity$Class, TCGA_data_Activity_EncodedX)
names(TCGA_data_Activity_Encoded)[1] <- "Class" 

# Ftiaxno ligo ta names gt exoun peiraxtei
names(TCGA_data_Activity_Encoded)[2:239] <- names(TCGA_data_with_Activity)[1:238]
names(TCGA_data_Activity_Encoded)[240] <- "IRE1_High_activity_(XBP1+/RIDD+)"
names(TCGA_data_Activity_Encoded)[241] <- "IRE1_Low_activity_(XBP1-/RIDD-)"
names(TCGA_data_Activity_Encoded)[242] = "IRE1_Mid_activity_(XBP1-/RIDD+)"
names(TCGA_data_Activity_Encoded)[243] = "IRE1_Mid_activity_(XBP1+/RIDD-)"
names(TCGA_data_Activity_Encoded)[244] = "XBP1_High_activity"
names(TCGA_data_Activity_Encoded)[245] = "XBP1_Low_activity"
names(TCGA_data_Activity_Encoded)[246] = "RIDD_High_activity"
names(TCGA_data_Activity_Encoded)[247] = "RIDD_Low_activity"

# dokimi montelou
# Model 5 (palio model.comb.50 / 50 features)  RF  -- Vgazei to kalytero survival

set.seed(139)
Model_RF_50 <- readRDS("Model.Top50(rf me TrCtrl + TnGr).rds")
predictions_50 <- predict(Model_RF_50, newdata = TCGA_data_Activity_Encoded)

TCGA_class_counts_50 <- table(predictions_50)

# Print the counts for TNBC and Non-TNBC
print(paste("TNBC:", TCGA_class_counts_50["TNBC"]))
print(paste("Non.TNBC:", TCGA_class_counts_50["Non.TNBC"]))

TCGA_raw_data_ilias@colData@listData$Class <- as.character(predictions_50)

function.SurvPlot <- function(exp, clusterCol, filename, title, colors){
  
  TCGAbiolinks::TCGAanalyze_survival(data = SummarizedExperiment::colData(exp),#data = surv.data,
                                     clusterCol = clusterCol,
                                     filename = filename,
                                     height = 8,
                                     width = 12,
                                     dpi = 100,
                                     legend = "Survival Plot",
                                     main = "Kaplan-Meier Overall Survival Curves",
                                     pvalue = TRUE,
                                     risk.table = TRUE,
                                     conf.int= FALSE,
                                     title = title,
                                     color = colors)
}

my_colors <- c("#008000","red")
function.SurvPlot(TCGA_raw_data_ilias, "Class", "French Dataset Results/Survival Plots/SurvivalPlot_for_164S_Model_top_50.pdf",
                  "K-M Overall Survival Curves / TCGA-BRCA\n308 TCGA Samples", my_colors)

