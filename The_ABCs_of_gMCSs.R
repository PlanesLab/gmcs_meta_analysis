############################################## #####
############# The ABC's of gMCSs ############# #####
############################################## #####
# Main code to compute the perfunctory analysis on gMCS #
####      Basics       ####
rm(list=ls())
cat("\014")
Today <- Sys.Date()

# Author: Danel Olaverri #

#### Load Libraries    ####
library(clusterProfiler)
library(data.table)
library(doSNOW)
library(dplyr)
library(ggbreak) 
library(ggrepel)
library(ggpattern)
library(ggplot2)
library(ggpubr)
library(ggstatsplot)
library(grid)
library(gridExtra)
library(openxlsx)
library(org.Hs.eg.db)
library(patchwork)
library(pheatmap)
library(pROC)
library(reshape2)
library(Rmpi)
library(tibble)

##################################### #####
############# Hart 2015 ############# #####
##################################### #####


# Thresholds to be considered: #
ThresholdVector <- c(0,1,2,2.5,5,10,20)
####    Prepare Data   ####  ####
# This data frame is required for future essentiality analysis, its only purpose is to deem all cells equal. #
Sample_Class_Dummy <- data.frame("Sample_Class" = "Cells")
levels(Sample_Class_Dummy)  <- "Cells"

## Get HumanGEM information ##
HumanGEM_Genes <- as.data.frame(fread("./Data/Genes_HumanGEM_v1_10_0.txt"))
Table_HumanGEM_Genes <- clusterProfiler::bitr(unlist(HumanGEM_Genes), "ENSEMBL", "SYMBOL", org.Hs.eg.db)

## Store function to compute Essential Genes ##
source('./Custom_Functions/Hart/fun_CalculateHartEssentialGenes_gmcsTH_DOM.R')
source('./Custom_Functions/Hart/fun_TweakedHartgMCS_Results.R')
source('./Custom_Functions/ComputeBinaryMatrix.R')

## Load Gene Expression Data (TPM) from Hart ##

Hart_Effect <- read.xlsx("./Data/Hart/Hart2015_BayesFactors.xlsx")[,-c(1,8,9,10)]
Hart_Exp    <- read.table("./Data/Hart/Hart2015_RNAseq.txt", header = TRUE)

## Load gMCSs ##
load('./Data/gMCSs_All_Cases.RData')
gMCS.info <- list(gMCS.info$EssentialTasks_CultureBiomass,
                  gMCS.info$Only_CultureBiomass)
names(gMCS.info) <- c("EssentialTasks_CultureBiomass", "Only_CultureBiomass")

## Prepare Hart Data ##
rownames(Hart_Exp) <- Hart_Exp[,1]
Hart_Exp <- Hart_Exp[,-1]
toENSEMBL <- clusterProfiler::bitr(Hart_Effect[,1], "SYMBOL", "ENSEMBL", org.Hs.eg.db)
colnames(Hart_Effect)[1] <- "SYMBOL"
Hart_Effect_ENSEMBL <- merge(toENSEMBL, Hart_Effect, by = "SYMBOL")
Hart_Effect_ENSEMBL <- Hart_Effect_ENSEMBL[!duplicated(Hart_Effect_ENSEMBL$SYMBOL),]
Hart_Effect_ENSEMBL <- Hart_Effect_ENSEMBL[!duplicated(Hart_Effect_ENSEMBL$ENSEMBL),]
Hart_Effect_SYMBOL  <- Hart_Effect_ENSEMBL[,1]
rownames(Hart_Effect_ENSEMBL) <- Hart_Effect_ENSEMBL[,2]
Hart_Effect_ENSEMBL <- Hart_Effect_ENSEMBL[,-c(1,2)]
colnames(Hart_Effect_ENSEMBL) <- toupper(sapply(strsplit(colnames(Hart_Effect_ENSEMBL), "_"), "[[", 2))
NumberOfgMCS <-  sapply(1:7, 
                        function(x) length(which(rowSums(gMCS.info$EssentialTasks_CultureBiomass$gMCSs.ENSEMBL.mat) == x)))

Effect_Tresholds <-  as.data.frame(rbind(c('DLD1', 3.57), c('GBM', 3.20), c('HCT116', 1.57), c('HELA', 15.47), c('RPE1', 6.84)))
Effect_Tresholds$V2 <- as.numeric(Effect_Tresholds$V2)

####   Compute gmcTH   ####  ####

if("AllThresholdsResults.RDS" %in% dir("./RDSResults/Hart/")){
  ResultList <- readRDS("./RDSResults/Hart/AllThresholdsResults.RDS")
} else {
  if("RoLAllThresholds.RDS" %in% dir("./RDSResults/Hart/")){
    RoLResults <-  readRDS("./RDSResults/Hart/RoLAllThresholds.RDS")
  } else {
    # Initialize list #
    RoLResults <- vector("list", length(ThresholdVector))
    names(RoLResults) <- paste0("Th",ThresholdVector)
    for(Threshold in ThresholdVector){
      if (paste0("RoL_Th",Threshold,".RDS") %in% dir("./RDSResults/Hart/")){
        RoLResults[[paste0("Th",Threshold)]] <- readRDS(paste0("./RDSResults/Hart/RoL_Th",Threshold,".RDS"))
      } else {
        # Compute for each threshold the Essential Genes #
        RoLResults[[paste0("Th",Threshold)]] <- CalculateEssentialGenes_gmcsTH_DOM(gene.exp = Hart_Exp,
                                                                                   gMCS_Analysis = gMCS.info$EssentialTasks_CultureBiomass,
                                                                                   sample.class = Sample_Class_Dummy,
                                                                                   gmcsTH_perc = Threshold/100,
                                                                                   nWorkers = 8,
                                                                                   gMCS_order = 8,
                                                                                   order1 = TRUE)
        saveRDS(RoLResults[[paste0("Th",Threshold)]], file = paste0("./RDSResults/Hart/RoL_Th",Threshold,".RDS"))
      }
    }
    saveRDS(RoLResults, "./RDSResults/Hart/RoLAllThresholds.RDS")
  }
  
  ResultList <- vector("list", length(ThresholdVector))
  names(ResultList) <- paste0("Th",ThresholdVector)
  
  for(Threshold in ThresholdVector){
    Binary_Matrix_Length_List <- ComputeBinaryMatrix(RoLResults[[paste0("Th",Threshold)]],
                                                     gMCS.info$EssentialTasks_CultureBiomass$gMCSs.ENSEMBL.list)
    # Get Results for each threshold #
    ResultsInList <- lapply(1:8, function(x) TweakedHartgMCS_Results(Binary_Matrix = Binary_Matrix_Length_List[[x]],
                                                                     Effect = Hart_Effect_ENSEMBL,
                                                                     Hart_Thresholds = Effect_Tresholds,
                                                                     Table_HumanGEM_Genes = Table_HumanGEM_Genes,
                                                                     gMCS_order =  x,
                                                                     Threshold_Info =  paste0("gMCS_T",Threshold)))
    ResultsInList <- do.call(rbind, ResultsInList)
    ResultsInList$Order <- factor(ResultsInList$Order,
                                  levels=c(paste0("Order_",1:7), "HigherThan_Order_7"))
    ResultsInList <-  ResultsInList %>% rename(TP = "True Positives", FP = "False Positives",
                                               TN = "True Negatives", FN = "False Negatives",
                                               PPV = "Positive Predictive Value",
                                               sensitivity = "Sensitivity", specificity = "Specificity",
                                               accuracy = "Accuracy", MCC = "Matthew's Cor. Coef.")
    # Summarise information and write csv #
    ResultsInList %>% group_by(Order) %>%
      summarise_at(c("True Positives", "False Positives", 
                     "False Negatives", "True Negatives",
                     "Accuracy", "Sensitivity", "Specificity",
                     "Positive Predictive Value", "Matthew's Cor. Coef."), median, na.rm = TRUE) %>%
      write.csv(file = paste0("./Results/Hart/gMCSLengths/",Today,"_LengthMedians_Th",Threshold,".csv"),
                row.names = F)
    ResultList[[paste0("Th",Threshold)]] <- ResultsInList
  }
  saveRDS(ResultList, file = "./RDSResults/Hart/AllThresholdsResults.RDS")
  
  ## Save Not Cumulative  ## #####
  # Not Cumulative data can be used to observe length-specific information #
  if("AllThresholdsResultsNotCumulative.RDS" %in% dir("./RDSResults/Hart/")){
    ResultList_NC <- readRDS("./RDSResults/Hart/AllThresholdsResultsNotCumulative.RDS")
  } else {
    ResultList_NC <- vector("list", length(ThresholdVector))
    names(ResultList_NC) <- paste0("Th",ThresholdVector)
    
    for(Threshold in ThresholdVector){
      Binary_Matrix_Length_List_NotCumulative <- vector("list", 8)
      Order1Genes <- unlist(gMCS.info$EssentialTasks_CultureBiomass$gMCSs.ENSEMBL.list[which(lengths(gMCS.info$EssentialTasks_CultureBiomass$gMCSs.ENSEMBL.list) == 1)])
      
      
      Binary_Matrix_Length_List <- ComputeBinaryMatrix(RoLResults[[paste0("Th",Threshold)]],
                                                       gMCS.info$EssentialTasks_CultureBiomass$gMCSs.ENSEMBL.list,
                                                       FALSE)
      ResultsInList <- lapply(1:8, function(x) TweakedHartgMCS_Results(Binary_Matrix = Binary_Matrix_Length_List[[x]],
                                                                       Effect = Hart_Effect_ENSEMBL,
                                                                       Hart_Thresholds = Effect_Tresholds,
                                                                       Table_HumanGEM_Genes = Table_HumanGEM_Genes,
                                                                       gMCS_order =  x,
                                                                       Threshold_Info =  paste0("gMCS_T",Threshold)))
      ResultsInList <- do.call(rbind, ResultsInList)
      ResultsInList$Order <- factor(ResultsInList$Order,
                                    levels=c(paste0("Order_",1:7), "HigherThan_Order_7"))
      ResultsInList <-  ResultsInList %>% rename(TP = "True Positives", FP = "False Positives",
                                                 TN = "True Negatives", FN = "False Negatives",
                                                 PPV = "Positive Predictive Value",
                                                 sensitivity = "Sensitivity", specificity = "Specificity",
                                                 accuracy = "Accuracy", MCC = "Matthew's Cor. Coef.")
      ResultsInList %>% group_by(Order) %>%
        summarise_at(c("True Positives", "False Positives", 
                       "False Negatives", "True Negatives",
                       "Accuracy", "Sensitivity", "Specificity",
                       "Positive Predictive Value", "Matthew's Cor. Coef."), median, na.rm = TRUE) %>%
        write.csv(file = paste0("./Results/Hart/NotCumulative_gMCSLengths/LengthMedians_Th",Threshold,".csv"),
                  row.names = F)
      ResultList_NC[[paste0("Th",Threshold)]] <- ResultsInList
    }
    saveRDS(ResultList_NC, file = "./RDSResults/Hart/AllThresholdsResultsNotCumulative.RDS")
  }
}


############################### ####
####  LocalT2 - 1244   ####  ####

source("./Custom_Functions/Hart/fun_LocalT2CalculateEssential.R")

if ("Results_localT2.RDS" %in% dir("./RDSResults/Hart/") && "NotCumulativeResults_localT2.RDS" %in% dir("./RDSResults/Hart/")){
  Results_localT2 <- readRDS("./RDSResults/Hart/Results_localT2.RDS")
  NotCumulativeResults_localT2 <- readRDS("./RDSResults/Hart/NotCumulativeResults_localT2.RDS")
} else {
  if ("RoL_gMCS_localT2.RDS" %in% dir("./RDSResults/Hart/")){
    RoL_gMCS_localT2 <- readRDS("./RDSResults/Hart/RoL_gMCS_localT2.RDS")
    
  } else {
    Sample_Class_Dummy <- data.frame("Sample_Class" = c("DLD1","GBM","HCT116","HELA","RPE1"))
    levels(Sample_Class_Dummy)  <- c("DLD1","GBM","HCT116","HELA","RPE1")
    Sample_Cohort_Dummy <- data.frame("Sample_Cohort" = "Hart2015")
    levels(Sample_Cohort_Dummy) <- "Hart2015"
    RoL_gMCS_localT2 <- LocalT2CalculateEssential(gene.exp      = Hart_Exp,
                                                  gMCS.info     = gMCS.info$EssentialTasks_CultureBiomass,
                                                  sample.class  = Sample_Class_Dummy,
                                                  sample.cohort = Sample_Cohort_Dummy,
                                                  localT2_mode  = "all_genes_gMCSs",
                                                  nWorkers      = 8,
                                                  gMCS_order    = 8)
    
    
    saveRDS(object = RoL_gMCS_localT2, file = "./RDSResults/Hart/RoL_gMCS_localT2.RDS")
  }
  #### Cumulative     ####
  BinaryMatrixLengthList <- ComputeBinaryMatrix(RoL_gMCS_localT2,
                                                gMCS.info$EssentialTasks_CultureBiomass$gMCSs.ENSEMBL.list)
  Results_localT2 <- lapply(1:8, function(x) TweakedHartgMCS_Results(Binary_Matrix = BinaryMatrixLengthList[[x]],
                                                                     Effect = Hart_Effect_ENSEMBL,
                                                                     Hart_Thresholds = Effect_Tresholds,
                                                                     Table_HumanGEM_Genes = Table_HumanGEM_Genes,
                                                                     gMCS_order =  x,
                                                                     Threshold_Info =  "localT2"))
  Results_localT2 <- do.call(rbind, Results_localT2)
  Results_localT2$Order <- factor(Results_localT2$Order,levels=c(paste0("Order_",1:7), "HigherThan_Order_7"))
  Results_localT2 <- Results_localT2 %>% dplyr::rename("True Positives"  = TP,
                                                       "False Positives" = FP,
                                                       "True Negatives"  = TN,
                                                       "False Negatives" = FN,
                                                       "Sensitivity" = sensitivity,
                                                       "Specificity" = specificity,
                                                       "Accuracy"    = accuracy,
                                                       "Positive Predictive Value" = PPV,
                                                       "Matthew's Cor. Coef."     = MCC)
  Results_localT2 %>% filter(thMethod == "localT2") %>% group_by(Order) %>%
    summarise_at(c("True Positives", "False Positives", 
                   "False Negatives", "True Negatives",
                   "Accuracy", "Sensitivity", "Specificity",
                   "Positive Predictive Value", "Matthew's Cor. Coef."), median, na.rm = TRUE) %>%
    write.csv(file = "./Results/Hart/localT2/Cumulative_LocalT2.csv", 
              row.names = F)
  
  MeltedLocalT2 <- melt(data = Results_localT2,
                        id.vars = "Order",
                        measure.vars = c("True Positives", "False Positives",
                                         "False Negatives", "True Negatives",
                                         "Accuracy", "Sensitivity", "Specificity",
                                         "Positive Predictive Value", "Matthew's Cor. Coef."))
  saveRDS(Results_localT2, "./RDSResults/Hart/Results_localT2.RDS")
  
  
  #### Not Cumulative ####
  NotCumulativeBinaryMatrixLengthList <- ComputeBinaryMatrix(RoL_gMCS_localT2,
                                                             gMCS.info$EssentialTasks_CultureBiomass$gMCSs.ENSEMBL.list,
                                                             FALSE)
  NotCumulativeResults_localT2 <- lapply(1:8,
                                         function(x) TweakedHartgMCS_Results(Binary_Matrix = NotCumulativeBinaryMatrixLengthList[[x]],
                                                                             Effect = Hart_Effect_ENSEMBL,
                                                                             Hart_Thresholds = Effect_Tresholds,
                                                                             Table_HumanGEM_Genes = Table_HumanGEM_Genes,
                                                                             gMCS_order =  x,
                                                                             Threshold_Info =  "localT2"))
  NotCumulativeResults_localT2 <- do.call(rbind, NotCumulativeResults_localT2)
  NotCumulativeResults_localT2$Order <- factor(NotCumulativeResults_localT2$Order,
                                               levels=c(paste0("Order_",1:7), "HigherThan_Order_7"))
  NotCumulativeResults_localT2 <- NotCumulativeResults_localT2 %>% dplyr::rename("True Positives"  = TP,
                                                                                 "False Positives" = FP,
                                                                                 "True Negatives"  = TN,
                                                                                 "False Negatives" = FN,
                                                                                 "Sensitivity" = sensitivity,
                                                                                 "Specificity" = specificity,
                                                                                 "Accuracy"    = accuracy,
                                                                                 "Positive Predictive Value" = PPV,
                                                                                 "Matthew's Cor. Coef."     = MCC)
  NotCumulativeResults_localT2 %>% filter(thMethod == "localT2") %>% group_by(Order) %>%
    summarise_at(c("True Positives", "False Positives", 
                   "False Negatives", "True Negatives",
                   "Accuracy", "Sensitivity", "Specificity",
                   "Positive Predictive Value", "Matthew's Cor. Coef."), median, na.rm = TRUE) %>%
    write.csv(file = "./Results/Hart/localT2/NotCumulative_LocalT2.csv", 
              row.names = F)
  
  NCMeltedLocalT2 <- melt(data = NotCumulativeResults_localT2,
                          id.vars = "Order",
                          measure.vars = c("True Positives", "False Positives",
                                           "False Negatives", "True Negatives",
                                           "Accuracy", "Sensitivity", "Specificity",
                                           "Positive Predictive Value", "Matthew's Cor. Coef."))
  
  saveRDS(NotCumulativeResults_localT2, "./RDSResults/Hart/NotCumulativeResults_localT2.RDS")
} 
####  LocalT2 - 3650   ####  ####
source("G:/Mi Unidad/Doctorado/R_Functinos/fun_LocalT2CalculateEssential.R")

if ("Results_localT2HumanGEM.RDS" %in% dir("./RDSResults/Hart/") &&
    "NotCumulativeResults_localT2HumanGEM.RDS" %in% dir("./RDSResults/Hart/")){
  Results_localT2HumanGEM <- readRDS("./RDSResults/Hart/Results_localT2HumanGEM.RDS")
  NotCumulativeResults_localT2HumanGEM <- readRDS("./RDSResults/Hart/NotCumulativeResults_localT2HumanGEM.RDS")
} else {
  if ("RoL_gMCS_localT2_HumanGEM.RDS" %in% dir("./RDSResults/Hart/")){
    RoL_gMCS_localT2HumanGEM <- readRDS("./RDSResults/Hart/RoL_gMCS_localT2_HumanGEM.RDS")
  } else {
    Sample_Class_Dummy <- data.frame("Sample_Class" = c("DLD1","GBM","HCT116","HELA","RPE1"))
    levels(Sample_Class_Dummy)  <- c("DLD1","GBM","HCT116","HELA","RPE1")
    Sample_Cohort_Dummy <- data.frame("Sample_Cohort" = "Hart2015")
    levels(Sample_Cohort_Dummy) <- "Hart2015"
    RoL_gMCS_localT2_HumanGEM <- LocalT2CalculateEssential(gene.exp = Hart_Exp,
                                                           gMCS.info = gMCS.info$EssentialTasks_CultureBiomass,
                                                           sample.class = Sample_Class_Dummy,
                                                           sample.cohort = Sample_Cohort_Dummy,
                                                           localT2_mode = "All_HumanGEM",
                                                           nWorkers = 8,
                                                           gMCS_order = 8)
    
  }
  saveRDS(object = RoL_gMCS_localT2_HumanGEM, file = "./RDSResults/Hart/RoL_gMCS_localT2_HumanGEM.RDS")
  
  # Cumulative     ####                     
  BinaryMatrixLengthList <- ComputeBinaryMatrix(RoL_gMCS_localT2_HumanGEM,
                                                gMCS.info$EssentialTasks_CultureBiomass$gMCSs.ENSEMBL.list)
  Results_localT2HumanGEM <- lapply(1:8, function(x) TweakedHartgMCS_Results(Binary_Matrix = BinaryMatrixLengthList[[x]],
                                                                             Effect = Hart_Effect_ENSEMBL,
                                                                             Hart_Thresholds = Effect_Tresholds,
                                                                             Table_HumanGEM_Genes = Table_HumanGEM_Genes,
                                                                             gMCS_order =  x,
                                                                             Threshold_Info =  "localT2HumanGEM"))
  Results_localT2HumanGEM <- do.call(rbind, Results_localT2HumanGEM)
  Results_localT2HumanGEM$Order <- factor(Results_localT2HumanGEM$Order,levels=c(paste0("Order_",1:7), "HigherThan_Order_7"))
  Results_localT2HumanGEM <- Results_localT2HumanGEM %>% dplyr::rename("True Positives"  = TP,
                                                                       "False Positives" = FP,
                                                                       "True Negatives"  = TN,
                                                                       "False Negatives" = FN,
                                                                       "Sensitivity" = sensitivity,
                                                                       "Specificity" = specificity,
                                                                       "Accuracy"    = accuracy,
                                                                       "Positive Predictive Value" = PPV,
                                                                       "Matthew's Cor. Coef."     = MCC)
  saveRDS(Results_localT2HumanGEM, "./RDSResults/Hart/Results_localT2HumanGEM.RDS")
  
  # Not Cumulative ####
  NotCumulativeBinaryMatrixLengthList <- ComputeBinaryMatrix(RoL_gMCS_localT2_HumanGEM,
                                                             gMCS.info$EssentialTasks_CultureBiomass$gMCSs.ENSEMBL.list,
                                                             FALSE)
  NotCumulativeResults_localT2HumanGEM <- lapply(1:8,
                                                 function(x) TweakedHartgMCS_Results(Binary_Matrix = NotCumulativeBinaryMatrixLengthList[[x]],
                                                                                     Effect = Hart_Effect_ENSEMBL,
                                                                                     Hart_Thresholds = Effect_Tresholds,
                                                                                     Table_HumanGEM_Genes = Table_HumanGEM_Genes,
                                                                                     gMCS_order =  x,
                                                                                     Threshold_Info =  "localT2HumanGEM"))
  NotCumulativeResults_localT2HumanGEM <- do.call(rbind, NotCumulativeResults_localT2HumanGEM)
  NotCumulativeResults_localT2HumanGEM$Order <- factor(NotCumulativeResults_localT2HumanGEM$Order,
                                                       levels=c(paste0("Order_",1:7), "HigherThan_Order_7"))
  NotCumulativeResults_localT2HumanGEM <- NotCumulativeResults_localT2HumanGEM %>% dplyr::rename("True Positives"  = TP,
                                                                                                 "False Positives" = FP,
                                                                                                 "True Negatives"  = TN,
                                                                                                 "False Negatives" = FN,
                                                                                                 "Sensitivity" = sensitivity,
                                                                                                 "Specificity" = specificity,
                                                                                                 "Accuracy"    = accuracy,
                                                                                                 "Positive Predictive Value" = PPV,
                                                                                                 "Matthew's Cor. Coef."     = MCC)
  NotCumulativeResults_localT2HumanGEM %>% filter(thMethod == "localT2") %>% group_by(Order) %>%
    summarise_at(c("True Positives", "False Positives", 
                   "False Negatives", "True Negatives",
                   "Accuracy", "Sensitivity", "Specificity",
                   "Positive Predictive Value", "Matthew's Cor. Coef."), median, na.rm = TRUE) %>%
    write.csv(file = "./Results/Hart/localT2/NotCumulative_LocalT2HumanGEM.csv", 
              row.names = F)
  
  saveRDS(NotCumulativeResults_localT2HumanGEM, "./RDSResults/Hart/NotCumulativeResults_localT2HumanGEM.RDS")
}
############################### ####
#         gMCS-to-Many          ####

for(Threshold_Name in c("Th0","Th1","Th2","Th2_5","Th5","Th10","Th20","localT2","localT2_HumanGEM")){
  
  RoL_gMCS <- readRDS(paste0("./RDSResults/Hart/RoL_",Threshold_Name,".RDS"))
  
  # Make gMCS-gene relations          ###
  Results_Order_List_Comparisons <- vector("list",3)
  names(Results_Order_List_Comparisons) <- c("One_Gene_to_One_gMCS_NoEssentials",
                                             "One_Gene_to_Many_gMCS_NoEssentials",
                                             "TrueEssentials")
  
  # Where is each gene Essential? #
  List_gMCS_Essential_Matrix <- RoL_gMCS$list.gMCS.essential.mat
  gMCS_Length <- 8
  Binary_Matrix <- matrix(0, nrow = nrow(RoL_gMCS$mat.essential.gene),
                          ncol = ncol(RoL_gMCS$mat.essential.gene),
                          dimnames = list(rownames(RoL_gMCS$mat.essential.gene),
                                          colnames(RoL_gMCS$mat.essential.gene)))
  NotBinary_Matrix <- Binary_Matrix
  
  for(Essential_Gene in names(List_gMCS_Essential_Matrix)){
    NotBinary_Matrix[Essential_Gene,] <- as.numeric(colSums(List_gMCS_Essential_Matrix[[Essential_Gene]]))
    Binary_Matrix[Essential_Gene,]    <- as.numeric(colSums(List_gMCS_Essential_Matrix[[Essential_Gene]]) > 0)
  }
  
  # Extend data for all HumanGEM genes #
  Essential_Gene_gMCS <-  do.call(rbind, lapply(List_gMCS_Essential_Matrix, colSums))
  Essential_Gene_gMCS_Extended <- matrix(data = 0, 
                                         nrow = length(unique(c(Table_HumanGEM_Genes$ENSEMBL,rownames(Essential_Gene_gMCS)))),
                                         ncol = ncol(Essential_Gene_gMCS),
                                         dimnames = list(unique(c(Table_HumanGEM_Genes$ENSEMBL,rownames(Essential_Gene_gMCS))),
                                                         colnames(Essential_Gene_gMCS)))
  # Essential genes are still Essential #
  Essential_Gene_gMCS_Extended[rownames(Essential_Gene_gMCS_Extended) %in%
                                 rownames(Essential_Gene_gMCS),] <- Essential_Gene_gMCS
  
  # Take Order 1 Genes #
  GetOrder1 <- rownames(Essential_Gene_gMCS_Extended)
  ThisAreOrder1 <- unname(unlist(
    gMCS.info$EssentialTasks_CultureBiomass$gMCSs.ENSEMBL.list[which(lengths(gMCS.info$EssentialTasks_CultureBiomass$gMCSs.ENSEMBL.list)==1)]))
  Order1_Indexes <- which(GetOrder1 %in% ThisAreOrder1)
  Essential_Gene_gMCS_Essentials <- Essential_Gene_gMCS_Extended
  for (Essential_Gene in Order1_Indexes){
    Essential_Gene_gMCS_Essentials[Essential_Gene, which(Essential_Gene_gMCS_Essentials[Essential_Gene,] != 0)] <- NA
  }
  Hart_Effect_ENSEMBL_NoOrder1 <-   Hart_Effect_ENSEMBL[!rownames(Hart_Effect_ENSEMBL) %in% ThisAreOrder1,]
  
  # Make 1 gMCS to 1 gene             ###
  One_to_One_Mat_NoEs <-  NotBinary_Matrix
  One_to_One_Mat_NoEs[One_to_One_Mat_NoEs != 1] <- NA # Caution! Changed! (Previously 0)
  Results_Order_List_Comparisons[[1]] <- One_to_One_Mat_NoEs
  
  
  # Make n gMCS to 1 gene             ###
  Many_to_One_Mat_Ess <-  NotBinary_Matrix
  Many_to_One_Mat_Ess[Many_to_One_Mat_Ess <= 1]  <- NA # Caution! Changed! (Previously 0)
  Many_to_One_Mat_Ess[Many_to_One_Mat_Ess > 1]   <- 1
  Results_Order_List_Comparisons[[2]] <- Many_to_One_Mat_Ess
  
  # All Essential Genes               ###
  Results_Order_List_Comparisons[[3]] <- Binary_Matrix
  
  mean(colSums(One_to_One_Mat_NoEs, na.rm =TRUE));mean(colSums(Many_to_One_Mat_Ess, na.rm =TRUE))
  # Compute Statistics                ###
  
  Comparison_Analysis_Results <- vector("list", length(Results_Order_List_Comparisons))
  k <- 0
  for (Analysis in names(Results_Order_List_Comparisons)){
    try({result_list <- list()
    k <- k+1
    evalMethod = "AllTasks"
    
    if (Analysis == "TrueEssentials"){
      Mat_binary <- matrix(NA, nrow = nrow(Hart_Effect_ENSEMBL), ncol = ncol(Hart_Effect_ENSEMBL))
      for (row in 1:nrow(Effect_Tresholds)){
        Mat_binary[,row] <- Hart_Effect_ENSEMBL[,row] > Effect_Tresholds$V2[row]
      }
      Mat_binary[which(is.na(Mat_binary))] <- 0
      colnames(Mat_binary) <- Effect_Tresholds$V1
      rownames(Mat_binary) <- rownames(Hart_Effect_ENSEMBL)
    } else {
      Mat_binary <- matrix(NA, nrow = nrow(Hart_Effect_ENSEMBL_NoOrder1), ncol = ncol(Hart_Effect_ENSEMBL_NoOrder1))
      for (row in 1:nrow(Effect_Tresholds)){
        Mat_binary[,row] <- Hart_Effect_ENSEMBL_NoOrder1[,row] > Effect_Tresholds$V2[row]
      }
      Mat_binary[which(is.na(Mat_binary))] <- 0
      colnames(Mat_binary) <- Effect_Tresholds$V1
      rownames(Mat_binary) <- rownames(Hart_Effect_ENSEMBL_NoOrder1)
    }
    
    tissues <- intersect(colnames(Results_Order_List_Comparisons[[Analysis]]),
                         colnames(Mat_binary))
    
    # Initialize Variables #
    results = data.frame('cellLine' = tissues,
                         'TP' = NA, 'TN' = NA, 'FP' = NA, 'FN' = NA,
                         "check_sum_genes" = NA, 
                         'accuracy' = NA, 'sensitivity' = NA,
                         'specificity' = NA, 'F1' = NA,
                         'MCC' = NA, 'Penr' = NA, 'logPenr' = NA,
                         'PenrAdj' = NA, 'logPenrAdj' = NA)
    
    for (t in 1:length(tissues)){
      
      # Genes From the Model #
      modelGenes <- unique(Table_HumanGEM_Genes$ENSEMBL) 
      # Essential Genes #
      modelEssential <- rownames(Results_Order_List_Comparisons[[Analysis]])[Results_Order_List_Comparisons[[Analysis]][,tissues[t]]>0]
      modelEssential <- unique(Table_HumanGEM_Genes$ENSEMBL[Table_HumanGEM_Genes$ENSEMBL %in% modelEssential])
      # Non-Essential Genes #
      modelNonEssential <- setdiff(modelGenes, modelEssential)
      
      # expGenes has all the genes that are not NA from Achilles #
      expGenes <-  rownames(Mat_binary)[!is.na(Mat_binary[,tissues[t]])]
      # expEssential has all the genes that are considered essential #
      expEssential <- expGenes[Mat_binary[expGenes,tissues[t]]==1]
      # expNonEssential has all the genes that are NOT essential # 
      expNonEssential <- setdiff(expGenes,expEssential)
      
      results$TP[t] = length(intersect(modelEssential, expEssential))        # True Positives  #
      results$TN[t] = length(intersect(modelNonEssential, expNonEssential))  # True Negatives  #
      results$FP[t] = length(intersect(modelEssential, expNonEssential))     # False Positives #
      results$FN[t] = length(intersect(modelNonEssential, expEssential))     # False Negatives #   
      
      # Intersect Genes from HumanGEM and Genes from ACHILLES #
      pop <- intersect(modelGenes, expGenes)
      # Essential Genes from gMCST5 #
      sample <- intersect(modelEssential, expGenes)
      # Essential Genes from ACHILLES #
      successes <- intersect(expEssential, modelGenes)
      
      results$Penr[t] <- phyper(length(intersect(successes,sample)) - 1, # Number of White Balls Drawn  #
                                length(intersect(successes,pop)),        # Number of White Balls in Urn #
                                length(setdiff(pop, successes)),         # Number of Black Balls in Urn #
                                length(sample),                          # Total Number of Balls Drawn  #  
                                lower.tail = FALSE);
      
    }
    
    
    results$check_sum_genes =  results$TP + results$TN + results$FP + results$FN
    
    # Calculate Metrics
    results$sensitivity = results$TP/(results$TP + results$FN)
    results$specificity = results$TN/(results$TN + results$FP)
    results$accuracy    = (results$TP + results$TN)/(results$TP + results$TN + results$FP + results$FN)
    results$PPV         = results$TP/(results$TP+results$FP)
    results$F1          = 2*results$TP/(2*results$TP + results$FP + results$FN)
    results$MCC         = ((results$TP*results$TN) - (results$FP*results$FN))/(sqrt((results$TP+results$FP)*(results$TP+results$FN))*sqrt((results$TN+results$FP)*(results$TN+results$FN)))  # Matthews correlation coefficient
    results$PenrAdj     = p.adjust(results$Penr,'BH')
    
    results$logPenr     = -log10(results$Penr)
    results$logPenrAdj  = -log10(results$PenrAdj)
    
    results$model <- "Hart2015"
    results$evalMethod <- evalMethod
    
    results$calcMethod <- Analysis
    
    results$PenrAdj = p.adjust(results$Penr,'BH')
    
    results$logPenr = -log10(results$Penr)
    results$logPenrAdj = -log10(results$PenrAdj)
    
    result_list[[evalMethod]] <- results
    
    results <- do.call(rbind, result_list)})
    
    Comparison_Analysis_Results[[k]] <- results
    
  }
  
  gMCS_Statistics <- do.call(rbind, Comparison_Analysis_Results)
  saveRDS(gMCS_Statistics, paste0("./RDSResults/Hart/One_to_Many_Results_",Threshold_Name,".RDS"))
  
}
############################### ####
#         Task-to-Many          ####
for(Threshold_Name in c("Th0","Th1","Th2","Th2_5","Th5","Th10","Th20","localT2","localT2_HumanGEM")){
  RoL_gMCS <- readRDS(paste0("./RDSResults/Hart/RoL_",Threshold_Name,".RDS"))
  
  # Make gMCS-gene relations          ##
  Results_Order_List_Comparisons <- vector("list",3)
  names(Results_Order_List_Comparisons) <- c("One_Gene_to_One_Task_NoEssentials",
                                             "One_Gene_to_Many_Tasks_NoEssentials",
                                             "TrueEssentials")
  
  List_gMCS_Essential_Matrix <- RoL_gMCS$list.gMCS.essential.mat
  gMCS_Length <- 8
  Binary_Matrix <- matrix(0, nrow = nrow(RoL_gMCS$mat.essential.gene),
                          ncol = ncol(RoL_gMCS$mat.essential.gene),
                          dimnames = list(rownames(RoL_gMCS$mat.essential.gene),
                                          colnames(RoL_gMCS$mat.essential.gene)))
  NotBinary_Matrix <- Binary_Matrix
  
  for(Essential_Gene in names(List_gMCS_Essential_Matrix)){
    Binary_Matrix[Essential_Gene,]    <- as.numeric(colSums(List_gMCS_Essential_Matrix[[Essential_Gene]]) > 0)
  }
  
  Essential_Gene_gMCS <-  do.call(rbind, lapply(List_gMCS_Essential_Matrix, colSums))
  # Extend data for all HumanGEM genes #
  Essential_Gene_gMCS_Extended <- matrix(data = 0, 
                                         nrow = length(unique(c(Table_HumanGEM_Genes$ENSEMBL,rownames(Essential_Gene_gMCS)))),
                                         ncol = ncol(Essential_Gene_gMCS),
                                         dimnames = list(unique(c(Table_HumanGEM_Genes$ENSEMBL,rownames(Essential_Gene_gMCS))),
                                                         colnames(Essential_Gene_gMCS)))
  # Essential genes are still Essential #
  for(EssentialGene in rownames(Binary_Matrix)){
    Essential_Gene_gMCS_Extended[rownames(Essential_Gene_gMCS_Extended) %in%
                                   EssentialGene,] <- Binary_Matrix[EssentialGene,]
  }
  
  # Take Order 1 Genes #
  GetOrder1 <- rownames(Essential_Gene_gMCS_Extended)
  ThisAreOrder1 <- unname(unlist(
    gMCS.info$EssentialTasks_CultureBiomass$gMCSs.ENSEMBL.list[which(lengths(gMCS.info$EssentialTasks_CultureBiomass$gMCSs.ENSEMBL.list)==1)]))
  Order1_Indexes <- which(GetOrder1 %in% ThisAreOrder1)
  Essential_Gene_gMCS_Essentials <- Essential_Gene_gMCS_Extended
  for (Essential_Gene in Order1_Indexes){
    Essential_Gene_gMCS_Essentials[Essential_Gene, which(Essential_Gene_gMCS_Essentials[Essential_Gene,] != 0)] <- NA
  }
  Hart_Effect_ENSEMBL_NoOrder1 <-   Hart_Effect_ENSEMBL[!rownames(Hart_Effect_ENSEMBL) %in% ThisAreOrder1,]
  
  # Tasks to Genes                                 #
  Task_DataFrame <- data.frame(matrix(data = NA, 
                                      nrow = length(List_gMCS_Essential_Matrix),
                                      ncol = 2),
                               row.names = names(List_gMCS_Essential_Matrix))
  colnames(Task_DataFrame) <- c("Essential_Gene", "Associated_Tasks")
  for (Essential_Gene in names(List_gMCS_Essential_Matrix)){
    Essentiality_Index <-  which(rowSums(List_gMCS_Essential_Matrix[[Essential_Gene]]) != 0)
    if (length(Essentiality_Index) != 0 ){
      gMCS_Tasks <- names(gMCS.info$EssentialTasks_CultureBiomass$gMCSs.ENSEMBL.txt[Essentiality_Index])
      gMCS_Tasks_All <- unique(sapply(1:length(gMCS_Tasks), function(x) strsplit(gMCS_Tasks[x], "--")[[1]][1]))
      Task_DataFrame[Essential_Gene,] <- c(names(List_gMCS_Essential_Matrix[Essential_Gene]),
                                           length(gMCS_Tasks_All))
    }
  }
  One_Gene_One_Task  <- Task_DataFrame$Essential_Gene[which(Task_DataFrame$Associated_Tasks == 1)]
  One_Gene_One_Task  <- One_Gene_One_Task[!One_Gene_One_Task %in% ThisAreOrder1]
  One_Gene_Many_Task <- Task_DataFrame$Essential_Gene[which(Task_DataFrame$Associated_Tasks != 1)]
  
  # Make 1 task to 1 gene             ##
  One_to_One_Mat_NoEs <-  Essential_Gene_gMCS_Essentials
  One_to_One_Mat_NoEs[rownames(One_to_One_Mat_NoEs) %in% One_Gene_One_Task]  <- 1  
  One_to_One_Mat_NoEs[rownames(One_to_One_Mat_NoEs) %in% One_Gene_Many_Task] <- NA 
  Results_Order_List_Comparisons[[1]] <- One_to_One_Mat_NoEs
  
  
  # Make n tasks to 1 gene            ##
  Task_Analysis_Mat_NoEs <-  Essential_Gene_gMCS_Essentials
  Task_Analysis_Mat_NoEs[rownames(Task_Analysis_Mat_NoEs) %in% One_Gene_One_Task] <- NA 
  Task_Analysis_Mat_NoEs[rownames(Task_Analysis_Mat_NoEs) %in% One_Gene_Many_Task]  <- 1
  Results_Order_List_Comparisons[[2]] <- Task_Analysis_Mat_NoEs
  
  # All Essential Genes               ##
  TrueEssentials <-  Essential_Gene_gMCS_Extended
  TrueEssentials[TrueEssentials != 0]   <- 1
  Results_Order_List_Comparisons[[3]] <- TrueEssentials
  
  # Compute Statistics                ##
  Comparison_Analysis_Results <- vector("list", length(Results_Order_List_Comparisons))
  k <- 0
  for (Analysis in names(Results_Order_List_Comparisons)){
    try({result_list <- list()
    k <- k+1
    evalMethod = "AllTasks"
    
    if (Analysis == "TrueEssentials"){
      Mat_binary <- matrix(NA, nrow = nrow(Hart_Effect_ENSEMBL), ncol = ncol(Hart_Effect_ENSEMBL))
      for (row in 1:nrow(Effect_Tresholds)){
        Mat_binary[,row] <- Hart_Effect_ENSEMBL[,row] > Effect_Tresholds$V2[row]
      }
      Mat_binary[which(is.na(Mat_binary))] <- 0
      colnames(Mat_binary) <- Effect_Tresholds$V1
      rownames(Mat_binary) <- rownames(Hart_Effect_ENSEMBL)
    } else {
      Mat_binary <- matrix(NA, nrow = nrow(Hart_Effect_ENSEMBL_NoOrder1), ncol = ncol(Hart_Effect_ENSEMBL_NoOrder1))
      for (row in 1:nrow(Effect_Tresholds)){
        Mat_binary[,row] <- Hart_Effect_ENSEMBL_NoOrder1[,row] > Effect_Tresholds$V2[row]
      }
      Mat_binary[which(is.na(Mat_binary))] <- 0
      colnames(Mat_binary) <- Effect_Tresholds$V1
      rownames(Mat_binary) <- rownames(Hart_Effect_ENSEMBL_NoOrder1)
    }
    
    tissues <- intersect(colnames(Results_Order_List_Comparisons[[Analysis]]),
                         colnames(Mat_binary))
    
    # Initialize Variables #
    results = data.frame('cellLine' = tissues,
                         'TP' = NA, 'TN' = NA, 'FP' = NA, 'FN' = NA,
                         "check_sum_genes" = NA, 
                         'accuracy' = NA, 'sensitivity' = NA,
                         'specificity' = NA, 'F1' = NA,
                         'MCC' = NA, 'Penr' = NA, 'logPenr' = NA,
                         'PenrAdj' = NA, 'logPenrAdj' = NA)
    
    for (t in 1:length(tissues)){
      
      # Genes From the Model #
      modelGenes <- unique(Table_HumanGEM_Genes$ENSEMBL) 
      # Essential Genes #
      modelEssential <- rownames(Results_Order_List_Comparisons[[Analysis]])[Results_Order_List_Comparisons[[Analysis]][,tissues[t]]>0]
      modelEssential <- unique(Table_HumanGEM_Genes$ENSEMBL[Table_HumanGEM_Genes$ENSEMBL %in% modelEssential])
      # Non-Essential Genes #
      modelNonEssential <- setdiff(modelGenes, modelEssential)
      
      # expGenes has all the genes that are not NA from Achilles #
      expGenes <-  rownames(Mat_binary)[!is.na(Mat_binary[,tissues[t]])]
      # expEssential has all the genes that are considered essential #
      expEssential <- expGenes[Mat_binary[expGenes,tissues[t]]==1]
      # expNonEssential has all the genes that are NOT essential # 
      expNonEssential <- setdiff(expGenes,expEssential)
      
      results$TP[t] = length(intersect(modelEssential, expEssential))        # True Positives  #
      results$TN[t] = length(intersect(modelNonEssential, expNonEssential))  # True Negatives  #
      results$FP[t] = length(intersect(modelEssential, expNonEssential))     # False Positives #
      results$FN[t] = length(intersect(modelNonEssential, expEssential))     # False Negatives #   
      
      # Intersect Genes from HumanGEM and Genes from ACHILLES #
      pop <- intersect(modelGenes, expGenes)
      # Essential Genes from gMCST5 #
      sample <- intersect(modelEssential, expGenes)
      # Essential Genes from ACHILLES #
      successes <- intersect(expEssential, modelGenes)
      
      results$Penr[t] <- phyper(length(intersect(successes,sample)) - 1, # Number of White Balls Drawn  #
                                length(intersect(successes,pop)),        # Number of White Balls in Urn #
                                length(setdiff(pop, successes)),         # Number of Black Balls in Urn #
                                length(sample),                          # Total Number of Balls Drawn  #  
                                lower.tail = FALSE);
      
    }
    
    
    results$check_sum_genes =  results$TP + results$TN + results$FP + results$FN
    
    # Calculate Metrics
    results$sensitivity = results$TP/(results$TP + results$FN)
    results$specificity = results$TN/(results$TN + results$FP)
    results$accuracy    = (results$TP + results$TN)/(results$TP + results$TN + results$FP + results$FN)
    results$PPV         = results$TP/(results$TP+results$FP)
    results$F1          = 2*results$TP/(2*results$TP + results$FP + results$FN)
    results$MCC         = ((results$TP*results$TN) - (results$FP*results$FN))/(sqrt((results$TP+results$FP)*(results$TP+results$FN))*sqrt((results$TN+results$FP)*(results$TN+results$FN)))  # Matthews correlation coefficient
    results$PenrAdj     = p.adjust(results$Penr,'BH')
    
    results$logPenr     = -log10(results$Penr)
    results$logPenrAdj  = -log10(results$PenrAdj)
    
    results$model <- "Hart2015"
    results$evalMethod <- evalMethod
    
    results$calcMethod <- Analysis
    
    results$PenrAdj = p.adjust(results$Penr,'BH')
    
    results$logPenr = -log10(results$Penr)
    results$logPenrAdj = -log10(results$PenrAdj)
    
    result_list[[evalMethod]] <- results
    
    results <- do.call(rbind, result_list)})
    
    Comparison_Analysis_Results[[k]] <- results
    
  }
  
  gMCS_Statistics <- do.call(rbind, Comparison_Analysis_Results)
  saveRDS(gMCS_Statistics, paste0("./RDSResults/Hart/Task_Promiscuity_Results_gMCS",Threshold_Name,".RDS"))
}
############################### ####
#            Biomass            ####
#### Compute gmcTH   #### ####

if("AllThresholdsResults.RDS" %in% dir("./Biomass/RDSResults/Hart")){
  ResultList <- readRDS("./RDSResults/Hart/Biomass/AllThresholdsResults.RDS")
} else {
  if("RoLAllThresholds.RDS" %in% dir("./RDSResults/Hart/Biomass/")){
    RoLResults <-  readRDS("./RDSResults/Hart/Biomass/RoLAllThresholds.RDS")
  } else {
    RoLResults <- vector("list", length(ThresholdVector))
    names(RoLResults) <- paste0("Th",ThresholdVector)
    for(Threshold in ThresholdVector){
      if (paste0("RoL_Th",Threshold,".RDS") %in% dir("./Biomass/RDSResults/Hart")){
        RoLResults[[paste0("Th",Threshold)]] <- readRDS(paste0("./RDSResults/Hart/Biomass/RoL_Th",Threshold,".RDS"))
      } else {
        RoLResults[[paste0("Th",Threshold)]] <- CalculateEssentialGenes_gmcsTH_DOM(gene.exp = Hart_Exp,
                                                                                   gMCS_Analysis = gMCS.info$Only_CultureBiomass,
                                                                                   sample.class = Sample_Class_Dummy,
                                                                                   gmcsTH_perc = Threshold/100,
                                                                                   nWorkers = 8,
                                                                                   gMCS_order = 8,
                                                                                   order1 = TRUE)
        saveRDS(RoLResults[[paste0("Th",Threshold)]], file = paste0("./RDSResults/Hart/Biomass/RoL_Th",Threshold,".RDS"))
      }
    }
    saveRDS(RoLResults, "./RDSResults/Hart/Biomass/RoLAllThresholds.RDS")
  }
  
  ResultList <- vector("list", length(ThresholdVector))
  names(ResultList) <- paste0("Th",ThresholdVector)
  
  for(Threshold in ThresholdVector){
    Binary_Matrix_Length_List <- ComputeBinaryMatrix(RoLResults[[paste0("Th",Threshold)]],
                                                     gMCS.info$Only_CultureBiomass$gMCSs.ENSEMBL.list)
    ResultsInList <- lapply(1:8, function(x) TweakedHartgMCS_Results(Binary_Matrix = Binary_Matrix_Length_List[[x]],
                                                                     Effect = Hart_Effect_ENSEMBL,
                                                                     Hart_Thresholds = Effect_Tresholds,
                                                                     Table_HumanGEM_Genes = Table_HumanGEM_Genes,
                                                                     gMCS_order =  x,
                                                                     Threshold_Info =  paste0("gMCS_T",Threshold)))
    ResultsInList <- do.call(rbind, ResultsInList)
    ResultsInList$Order <- factor(ResultsInList$Order,
                                  levels=c(paste0("Order_",1:7), "HigherThan_Order_7"))
    ResultsInList <-  ResultsInList %>% rename(TP = "True Positives", FP = "False Positives",
                                               TN = "True Negatives", FN = "False Negatives",
                                               PPV = "Positive Predictive Value",
                                               sensitivity = "Sensitivity", specificity = "Specificity",
                                               accuracy = "Accuracy", MCC = "Matthew's Cor. Coef.")
    ResultsInList %>% group_by(Order) %>%
      summarise_at(c("True Positives", "False Positives", 
                     "False Negatives", "True Negatives",
                     "Accuracy", "Sensitivity", "Specificity",
                     "Positive Predictive Value", "Matthew's Cor. Coef."), median, na.rm = TRUE) %>%
      write.csv(file = paste0("./Results/Hart/gMCSLengths/LengthMedians_Th",Threshold,".csv"),
                row.names = F)
    ResultList[[paste0("Th",Threshold)]] <- ResultsInList
  }
  saveRDS(ResultList, file = "./RDSResults/Hart/Biomass/AllThresholdsResults.RDS")
}

####  LocalT2 - 1244 #### ####
if ("Results_localT2.RDS" %in% dir("./RDSResults/Hart/Biomass/") && "NotCumulativeResults_localT2.RDS" %in% dir("./RDSResults/Hart/Biomass/")){
  Results_localT2 <- readRDS("./RDSResults/Hart/Biomass/Results_localT2.RDS")
  NotCumulativeResults_localT2 <- readRDS("./RDSResults/Hart/Biomass/NotCumulativeResults_localT2.RDS")
} else {
  if ("RoL_gMCS_localT2.RDS" %in% dir("./RDSResults/Hart/Biomass/")){
    RoL_gMCS_localT2 <- readRDS("./RDSResults/Hart/Biomass/RoL_gMCS_localT2.RDS")
    
  } else {
    Sample_Class_Dummy <- data.frame("Sample_Class" = c("DLD1","GBM","HCT116","HELA","RPE1"))
    levels(Sample_Class_Dummy)  <- c("DLD1","GBM","HCT116","HELA","RPE1")
    Sample_Cohort_Dummy <- data.frame("Sample_Cohort" = "Hart2015")
    levels(Sample_Cohort_Dummy) <- "Hart2015"
    RoL_gMCS_localT2 <- LocalT2CalculateEssential(gene.exp      = Hart_Exp,
                                                  gMCS.info     = gMCS.info$Only_CultureBiomass,
                                                  sample.class  = Sample_Class_Dummy,
                                                  sample.cohort = Sample_Cohort_Dummy,
                                                  localT2_mode  = "all_genes_gMCSs",
                                                  nWorkers      = 8,
                                                  gMCS_order    = 8)
    
    
    saveRDS(object = RoL_gMCS_localT2, file = "./RDSResults/Hart/Biomass/RoL_gMCS_localT2.RDS")
  }
  # Cumulative     ####
  BinaryMatrixLengthList <- ComputeBinaryMatrix(RoL_gMCS_localT2,
                                                gMCS.info$Only_CultureBiomass$gMCSs.ENSEMBL.list)
  Results_localT2 <- lapply(1:8, function(x) TweakedHartgMCS_Results(Binary_Matrix = BinaryMatrixLengthList[[x]],
                                                                     Effect = Hart_Effect_ENSEMBL,
                                                                     Hart_Thresholds = Effect_Tresholds,
                                                                     Table_HumanGEM_Genes = Table_HumanGEM_Genes,
                                                                     gMCS_order =  x,
                                                                     Threshold_Info =  "localT2"))
  Results_localT2 <- do.call(rbind, Results_localT2)
  Results_localT2$Order <- factor(Results_localT2$Order,levels=c(paste0("Order_",1:7), "HigherThan_Order_7"))
  Results_localT2 <- Results_localT2 %>% dplyr::rename("True Positives"  = TP,
                                                       "False Positives" = FP,
                                                       "True Negatives"  = TN,
                                                       "False Negatives" = FN,
                                                       "Sensitivity" = sensitivity,
                                                       "Specificity" = specificity,
                                                       "Accuracy"    = accuracy,
                                                       "Positive Predictive Value" = PPV,
                                                       "Matthew's Cor. Coef."     = MCC)
  Results_localT2 %>% filter(thMethod == "localT2") %>% group_by(Order) %>%
    summarise_at(c("True Positives", "False Positives", 
                   "False Negatives", "True Negatives",
                   "Accuracy", "Sensitivity", "Specificity",
                   "Positive Predictive Value", "Matthew's Cor. Coef."), median, na.rm = TRUE) %>%
    write.csv(file = "./Results/Biomass/Hart/localT2/Cumulative_LocalT2.csv", 
              row.names = F)
  
  MeltedLocalT2 <- melt(data = Results_localT2,
                        id.vars = "Order",
                        measure.vars = c("True Positives", "False Positives",
                                         "False Negatives", "True Negatives",
                                         "Accuracy", "Sensitivity", "Specificity",
                                         "Positive Predictive Value", "Matthew's Cor. Coef."))
  LocalT2_Length <- ggplot(MeltedLocalT2, aes(x = Order, y = value, fill = Order)) + 
    facet_wrap(~variable,scales = "free") + geom_boxplot() + theme_bw() +  
    # geom_boxplot(width=0.1, fill = NA) +
    labs(title = "LocalT2_1244",
         x = "", y = "",
         fill = "Lengths") + 
    scale_fill_discrete(name = "Maximum gMCS\nLengths (#gMCS)",
                        labels = c(paste0("1 (n = ",NumberOfgMCS[1],")"), paste0("2 (n = ",NumberOfgMCS[2],")"),
                                   paste0("3 (n = ",NumberOfgMCS[3],")"), paste0("4 (n = ",NumberOfgMCS[4],")"),
                                   paste0("5 (n = ",NumberOfgMCS[5],")"), paste0("6 (n = ",NumberOfgMCS[6],")"),
                                   paste0("7 (n = ",NumberOfgMCS[7],")"),
                                   paste0("All (n = ",nrow(gMCS.info$Only_CultureBiomass$gMCSs.ENSEMBL.mat),")")),
                        type = RColorBrewer::brewer.pal(8, "Paired")) +
    theme(plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
          axis.text.x = element_blank(), axis.ticks.x = element_blank(),
          strip.text.x = element_text(size = 18, face = "bold"),
          legend.position = "bottom")
  
  # LocalT2_Length
  saveRDS(Results_localT2, "./RDSResults/Hart/Biomass/Results_localT2.RDS")
  
  
  # Not Cumulative ####
  NotCumulativeBinaryMatrixLengthList <- ComputeBinaryMatrix(RoL_gMCS_localT2,
                                                             gMCS.info$Only_CultureBiomass$gMCSs.ENSEMBL.list,
                                                             FALSE)
  NotCumulativeResults_localT2 <- lapply(1:8,
                                         function(x) TweakedHartgMCS_Results(Binary_Matrix = NotCumulativeBinaryMatrixLengthList[[x]],
                                                                             Effect = Hart_Effect_ENSEMBL,
                                                                             Hart_Thresholds = Effect_Tresholds,
                                                                             Table_HumanGEM_Genes = Table_HumanGEM_Genes,
                                                                             gMCS_order =  x,
                                                                             Threshold_Info =  "localT2"))
  NotCumulativeResults_localT2 <- do.call(rbind, NotCumulativeResults_localT2)
  NotCumulativeResults_localT2$Order <- factor(NotCumulativeResults_localT2$Order,
                                               levels=c(paste0("Order_",1:7), "HigherThan_Order_7"))
  NotCumulativeResults_localT2 <- NotCumulativeResults_localT2 %>% dplyr::rename("True Positives"  = TP,
                                                                                 "False Positives" = FP,
                                                                                 "True Negatives"  = TN,
                                                                                 "False Negatives" = FN,
                                                                                 "Sensitivity" = sensitivity,
                                                                                 "Specificity" = specificity,
                                                                                 "Accuracy"    = accuracy,
                                                                                 "Positive Predictive Value" = PPV,
                                                                                 "Matthew's Cor. Coef."     = MCC)
  NotCumulativeResults_localT2 %>% filter(thMethod == "localT2") %>% group_by(Order) %>%
    summarise_at(c("True Positives", "False Positives", 
                   "False Negatives", "True Negatives",
                   "Accuracy", "Sensitivity", "Specificity",
                   "Positive Predictive Value", "Matthew's Cor. Coef."), median, na.rm = TRUE) %>%
    write.csv(file = "./Results/Biomass/Hart/localT2/NotCumulative_LocalT2.csv", 
              row.names = F)
  
  NCMeltedLocalT2 <- melt(data = NotCumulativeResults_localT2,
                          id.vars = "Order",
                          measure.vars = c("True Positives", "False Positives",
                                           "False Negatives", "True Negatives",
                                           "Accuracy", "Sensitivity", "Specificity",
                                           "Positive Predictive Value", "Matthew's Cor. Coef."))
  NotCumLocalT2_Length <- ggplot(NCMeltedLocalT2, aes(x = Order, y = value, fill = Order)) + 
    facet_wrap(~variable,scales = "free") + geom_boxplot() + theme_bw() +  
    # geom_boxplot(width=0.1, fill = NA) +
    labs(title = "LocalT2_1244",
         x = "", y = "",
         fill = "Lengths") + 
    scale_fill_discrete(name = "Maximum gMCS\nLengths (#gMCS)",
                        labels = c(paste0("1 (n = ",NumberOfgMCS[1],")"), paste0("2 (n = ",NumberOfgMCS[2],")"),
                                   paste0("3 (n = ",NumberOfgMCS[3],")"), paste0("4 (n = ",NumberOfgMCS[4],")"),
                                   paste0("5 (n = ",NumberOfgMCS[5],")"), paste0("6 (n = ",NumberOfgMCS[6],")"),
                                   paste0("7 (n = ",NumberOfgMCS[7],")"),
                                   paste0("All (n = ",nrow(gMCS.info$Only_CultureBiomass$gMCSs.ENSEMBL.mat),")")),
                        type = RColorBrewer::brewer.pal(8, "Paired")) +
    theme(plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
          axis.text.x = element_blank(), axis.ticks.x = element_blank(),
          strip.text.x = element_text(size = 18, face = "bold"),
          legend.position = "bottom")
  
  # NotCumLocalT2_Length
  
  ggsave(paste0("./Results/Hart/Biomass/localT2/NotCumulativeLocalT2.png"),
         plot = NotCumLocalT2_Length, device = "png", width = 16, height = 12, units = "in", dpi = 300)
  
  saveRDS(NotCumulativeResults_localT2, "./RDSResults/Hart/Biomass/NotCumulativeResults_localT2.RDS")
} 
############################### ####
#########   Figures   ######### ####
############################### ####
# Figure 1A - HART ####
toAdd <- do.call(rbind,readRDS("./RDSResults/Hart/AllThresholdsResults.RDS"))
toKeep <- c("gMCS_T0","gMCS_T1","gMCS_T2","gMCS_T2.5","gMCS_T5","gMCS_T10","gMCS_T20")
toRemove <- setdiff(unique(toAdd$thMethod), toKeep)

for(element in toRemove){toAdd <- toAdd[which(toAdd$thMethod != element),]}
toAdd$thMethod <- factor(toAdd$thMethod, levels=toKeep)
ThresholdColors <- viridis::magma(8)

All_Results_toDF     <- rbind(toAdd,
                 readRDS("./RDSResults/Hart/Results_localT2.RDS"),
                 readRDS("./RDSResults/Hart/Results_localT2HumanGEM.RDS"))

All_Results_toDF$thMethod <- factor(All_Results_toDF$thMethod, levels = c("gMCS_T0","gMCS_T1","gMCS_T2","gMCS_T2.5",
                                                "gMCS_T5","gMCS_T10","gMCS_T20",
                                                "NA",
                                                "localT2", "localT2HumanGEM"))
All_Results_toDF <- All_Results_toDF %>% filter(Order == "HigherThan_Order_7")
All_Results_toDF <- reshape2::melt(All_Results_toDF, id.vars = c("thMethod"),
                      measure.vars =  c("True Positives", "False Positives", "Positive Predicted Value"))

# Prepare TP Plot #
TruePositive_tmp          <- All_Results_toDF %>% filter(variable == "True Positives")
if(FALSE){
  TruePositive_tmp %>% kruskal.test(value~thMethod, data=.)
  # t.test   ####
  
  pValThDF <- as.data.frame(matrix(1,nrow = length(unique(TruePositive_tmp$thMethod)),
                                   ncol = length(unique(TruePositive_tmp$thMethod))))
  colnames(pValThDF) = rownames(pValThDF) = unique(TruePositive_tmp$thMethod)
  
  for(i in unique(TruePositive_tmp$thMethod)){
    for(j in setdiff(unique(TruePositive_tmp$thMethod),i)){
      pValThDF[i,j] <- round(t.test(TruePositive_tmp$value[TruePositive_tmp$thMethod == i],
                                    TruePositive_tmp$value[TruePositive_tmp$thMethod == j],
                                    paired = TRUE)$p.value, 5)
      
    }
  }
  write.csv(x = pValThDF, file = paste0("./Results/Hart/Threshold/",Today,"_Threshold_TP_pvals_ttest.csv"),
            quote = FALSE, row.names = TRUE)
  # Wilcoxon ####
  pValThDF <- as.data.frame(matrix(1,nrow = length(unique(TruePositive_tmp$thMethod)),
                                   ncol = length(unique(TruePositive_tmp$thMethod))))
  colnames(pValThDF) = rownames(pValThDF) = unique(TruePositive_tmp$thMethod)
  
  for(i in unique(TruePositive_tmp$thMethod)){
    for(j in setdiff(unique(TruePositive_tmp$thMethod),i)){
      pValThDF[i,j] <- round(wilcox.test(TruePositive_tmp$value[TruePositive_tmp$thMethod == i],
                                         TruePositive_tmp$value[TruePositive_tmp$thMethod == j],
                                         paired = TRUE)$p.value, 5)
      
    }
  }
  write.csv(x = pValThDF, file = paste0("./Results/Hart/Threshold/",Today,"_Threshold_TP_pvals_wilcox.csv"),
            quote = FALSE, row.names = TRUE)
  # Continue ####
  for(i in unique(TruePositive_tmp$thMethod)){
    shtest <- shapiro.test(TruePositive_tmp$value[TruePositive_tmp$thMethod == i])
    varthr <- var(TruePositive_tmp$value[TruePositive_tmp$thMethod == i])
    if(i != "localT2" & i != "localT2HumanGEM"){
      cat(paste0("Threshold: ", sapply(strsplit(i,"_T"), "[[", 2), "\t Shapiro-Test: ", shtest$p.value, "\t Variance: ", varthr, "\n"))
    } else {
      cat(paste0("Threshold: ", i, "\t Shapiro-Test: ", shtest$p.value, "\t Variance: ", varthr, "\n"))
    }
  }
  
}

TruePositive_tmp_mean <- TruePositive_tmp %>% 
  group_by(thMethod) %>% 
  summarize(Mean = mean(value)) %>%
  ungroup()


TruePositive_tmp <- rbind(TruePositive_tmp, list("NA","True Positives",NA))
change_lab <- "True Positives"; names(change_lab) <- "True Positives"
TPlot <- ggplot(TruePositive_tmp, aes(x = thMethod, y = value, fill = thMethod)) + 
  geom_boxplot(aes(fill = thMethod, alpha = 0.4)) + 
  scale_fill_manual(values = c("gMCS_T0"   = ThresholdColors[1],"gMCS_T1"   = ThresholdColors[2],
                               "gMCS_T2"   = ThresholdColors[3],"gMCS_T2.5" = ThresholdColors[4],
                               "gMCS_T5"   = ThresholdColors[5],"gMCS_T10"  = ThresholdColors[6],
                               "gMCS_T20"  = ThresholdColors[7],
                               "NA" = "#000000",
                               "localT2" = "#B2DF8A", "l ocalT2HumanGEM" = "#6ba35b")) +
  geom_point(aes(color = thMethod, size = 2)) + 
  facet_wrap("variable",labeller = labeller(variable = change_lab))  +
  scale_color_manual(values = c("gMCS_T0"   = ThresholdColors[1],"gMCS_T1"   = ThresholdColors[2],
                                "gMCS_T2"   = ThresholdColors[3],"gMCS_T2.5" = ThresholdColors[4],
                                "gMCS_T5"   = ThresholdColors[5],"gMCS_T10"  = ThresholdColors[6],
                                "gMCS_T20"  = ThresholdColors[7],
                                "NA" = "#000000",
                                "localT2" = "#B2DF8A", "localT2HumanGEM" = "#6ba35b")) +
  theme_bw() + xlab("") + ylab("") + theme(legend.position = "none") + 
  geom_point(data    = TruePositive_tmp_mean,
             mapping = aes(x = thMethod, y = Mean, size = 1.5),
             color   = "dark red") + 
  geom_line(data = TruePositive_tmp_mean[grepl("gMCS", TruePositive_tmp_mean$thMethod),], 
            mapping = aes(x = thMethod, y = Mean, group = 1, size = 1.25)) +
  geom_line(data = TruePositive_tmp_mean[!grepl("gMCS", TruePositive_tmp_mean$thMethod),], 
            mapping = aes(x = thMethod, y = Mean, group = 1, size = 1.25)) +
  geom_vline(xintercept = 8, size = 2, color = "dark red") +
  scale_x_discrete(labels=c("gMCS_T0"   = "0%",
                            "gMCS_T1"   = "1%",
                            "gMCS_T2"   = "2%",
                            "gMCS_T2.5" = "2.5%",
                            "gMCS_T5"   = "5%",
                            "gMCS_T10"  = "10%",
                            "gMCS_T20"  = "20%",
                            "NA" = "",
                            "localT2"   = "LocalT2",
                            "localT2HumanGEM"  = "LocalT2-H1"),
                   limits = c ("gMCS_T0","gMCS_T1","gMCS_T2","gMCS_T2.5",
                               "gMCS_T5","gMCS_T10","gMCS_T20", "NA",
                               "localT2","localT2HumanGEM")) +
  scale_y_continuous(limits = c(50,NA), breaks = c(seq(from = 50, to = 150, by = 10))) +
  theme(strip.text = element_text(face  = "bold", size = 28),
        strip.background = element_rect(fill = "#afc3dc"),
        axis.text.x = element_text(angle = 30, hjust = 1, vjust = 1, size = 20),
        axis.text.y = element_text(size = 23)) 
# Prepare FP Plot #
FalsePositive_tmp          <- All_Results_toDF %>% filter(variable == "False Positives")
# Get P values #
if(FALSE){
  FalsePositive_tmp %>% kruskal.test(value~thMethod, data=.)
  # t.test   ####
  pValThDF <- as.data.frame(matrix(1,nrow = length(unique(FalsePositive_tmp$thMethod)),
                                   ncol = length(unique(FalsePositive_tmp$thMethod))))
  colnames(pValThDF) = rownames(pValThDF) = unique(FalsePositive_tmp$thMethod)
  
  for(i in unique(FalsePositive_tmp$thMethod)){
    for(j in setdiff(unique(FalsePositive_tmp$thMethod),i)){
      pValThDF[i,j] <- round(wilcox.test(FalsePositive_tmp$value[FalsePositive_tmp$thMethod == i],
                                         FalsePositive_tmp$value[FalsePositive_tmp$thMethod == j],
                                         paired = TRUE)$p.value,5)
      
    }
  }
  write.csv(x = pValThDF, file = paste0("./Results/Hart/Threshold/",Today,"_Threshold_FP_ttest.csv"),
            quote = FALSE, row.names = TRUE)
  # Wilcoxon ####
  pValThDF <- as.data.frame(matrix(1,nrow = length(unique(FalsePositive_tmp$thMethod)),
                                   ncol = length(unique(FalsePositive_tmp$thMethod))))
  colnames(pValThDF) = rownames(pValThDF) = unique(FalsePositive_tmp$thMethod)
  
  for(i in unique(FalsePositive_tmp$thMethod)){
    for(j in setdiff(unique(FalsePositive_tmp$thMethod),i)){
      pValThDF[i,j] <- round(wilcox.test(FalsePositive_tmp$value[FalsePositive_tmp$thMethod == i],
                                         FalsePositive_tmp$value[FalsePositive_tmp$thMethod == j],
                                         paired = TRUE)$p.value,5)
      
    }
  }
  write.csv(x = pValThDF, file = paste0("./Results/Hart/Threshold/",Today,"_Threshold_FP_wilcoxon.csv"),
            quote = FALSE, row.names = TRUE)
  # Continue ####
  for(i in unique(FalsePositive_tmp$thMethod)){
    shtest <- shapiro.test(FalsePositive_tmp$value[FalsePositive_tmp$thMethod == i])
    varthr <- var(FalsePositive_tmp$value[FalsePositive_tmp$thMethod == i])
    if(i != "localT2" & i != "localT2HumanGEM"){
      cat(paste0("Threshold: ", sapply(strsplit(i,"_T"), "[[", 2), "\t Shapiro-Test: ", shtest$p.value, "\t Variance: ", varthr, "\n"))
    } else {
      cat(paste0("Threshold: ", i, "\t Shapiro-Test: ", shtest$p.value, "\t Variance: ", varthr, "\n"))
    }
  }
  
}
FalsePositive_tmp_mean <- FalsePositive_tmp %>% 
  group_by(thMethod) %>% 
  summarize(Mean = mean(value)) %>%
  ungroup()

FalsePositive_tmp <- rbind(FalsePositive_tmp, list("NA","False Positives",NA))
change_lab <- "False Positives"; names(change_lab) <- "False Positives"
FPlot <- ggplot(FalsePositive_tmp, aes(x = thMethod, y = value, fill = thMethod)) + 
  geom_boxplot(aes(fill = thMethod, alpha = 0.4)) + 
  scale_fill_manual(values = c("gMCS_T0"   = ThresholdColors[1],"gMCS_T1"   = ThresholdColors[2],
                               "gMCS_T2"   = ThresholdColors[3],"gMCS_T2.5" = ThresholdColors[4],
                               "gMCS_T5"   = ThresholdColors[5],"gMCS_T10"  = ThresholdColors[6],
                               "gMCS_T20"  = ThresholdColors[7],
                               "NA" = "#000000",
                               "localT2" = "#B2DF8A", "localT2HumanGEM" = "#6ba35b")) +
  geom_point(aes(color = thMethod, size = 2)) + 
  facet_wrap("variable",labeller = labeller(variable = change_lab))  +
  scale_color_manual(values = c("gMCS_T0"   = ThresholdColors[1],"gMCS_T1"   = ThresholdColors[2],
                                "gMCS_T2"   = ThresholdColors[3],"gMCS_T2.5" = ThresholdColors[4],
                                "gMCS_T5"   = ThresholdColors[5],"gMCS_T10"  = ThresholdColors[6],
                                "gMCS_T20"  = ThresholdColors[7],
                                "NA" = "#000000",
                                "localT2" = "#B2DF8A", "localT2HumanGEM" = "#6ba35b")) +
  theme_bw() + xlab("") + ylab("") +
  theme(legend.position = "none") + 
  geom_point(data    = FalsePositive_tmp_mean,
             mapping = aes(x = thMethod, y = Mean, size = 1.5),
             color   = "dark red") + 
  geom_line(data = FalsePositive_tmp_mean[grepl("gMCS", FalsePositive_tmp_mean$thMethod),], 
            mapping = aes(x = thMethod, y = Mean, group = 1, size = 1.25)) +
  geom_line(data = FalsePositive_tmp_mean[!grepl("gMCS", FalsePositive_tmp_mean$thMethod),], 
            mapping = aes(x = thMethod, y = Mean, group = 1, size = 1.25)) +
  geom_vline(xintercept = 8, size = 2, color = "dark red") +
  scale_x_discrete(labels=c("gMCS_T0"   = "0%",
                            "gMCS_T1"   = "1%",
                            "gMCS_T2"   = "2%",
                            "gMCS_T2.5" = "2.5%",
                            "gMCS_T5"   = "5%",
                            "gMCS_T10"  = "10%",
                            "gMCS_T20"  = "20%",
                            "NA" = "",
                            "localT2"   = "LocalT2",
                            "localT2HumanGEM"  = "LocalT2-H1"),
                   limits = c ("gMCS_T0","gMCS_T1","gMCS_T2","gMCS_T2.5",
                               "gMCS_T5","gMCS_T10","gMCS_T20", "NA",
                               "localT2","localT2HumanGEM")) +
  scale_y_continuous(limits = c(100,NA), breaks = c(seq(from = 100, to = 250, by = 15))) +
  theme(strip.text   = element_text(face  = "bold", size = 28),
        strip.background = element_rect(fill = "#afc3dc"),
        axis.text.x  = element_text(angle = 30, hjust = 1, vjust = 1, size = 20),
        axis.text.y  = element_text(size = 23)) 
FPlot
# Prepare PPV Plot #
PositivePredictiveValue_tmp          <- All_Results_toDF %>% filter(variable == "Positive Predicted Value")
# Get P values #
PositivePredictiveValue_tmp %>% kruskal.test(value~thMethod, data=.)
if(FALSE){
  # t.test   ####
  pValThDF <- as.data.frame(matrix(1,nrow = length(unique(PositivePredictiveValue_tmp$thMethod)),
                                   ncol = length(unique(PositivePredictiveValue_tmp$thMethod))))
  colnames(pValThDF) = rownames(pValThDF) = unique(PositivePredictiveValue_tmp$thMethod)
  
  for(i in unique(PositivePredictiveValue_tmp$thMethod)){
    for(j in setdiff(unique(PositivePredictiveValue_tmp$thMethod),i)){
      pValThDF[i,j] <- round(t.test(PositivePredictiveValue_tmp$value[PositivePredictiveValue_tmp$thMethod == i],
                                    PositivePredictiveValue_tmp$value[PositivePredictiveValue_tmp$thMethod == j],
                                    paired = TRUE)$p.value,5)
      
    }
  }
  write.csv(x = pValThDF, file = paste0("./Results/Hart/Threshold/",Today,"_Threshold_PPV_ttest.csv"),
            quote = FALSE, row.names = TRUE)
  # Wilcoxon ####
  pValThDF <- as.data.frame(matrix(1,nrow = length(unique(PositivePredictiveValue_tmp$thMethod)),
                                   ncol = length(unique(PositivePredictiveValue_tmp$thMethod))))
  colnames(pValThDF) = rownames(pValThDF) = unique(PositivePredictiveValue_tmp$thMethod)
  
  for(i in unique(PositivePredictiveValue_tmp$thMethod)){
    for(j in setdiff(unique(PositivePredictiveValue_tmp$thMethod),i)){
      pValThDF[i,j] <- round(wilcox.test(PositivePredictiveValue_tmp$value[PositivePredictiveValue_tmp$thMethod == i],
                                         PositivePredictiveValue_tmp$value[PositivePredictiveValue_tmp$thMethod == j],
                                         paired = TRUE)$p.value,5)
      
    }
  }
  write.csv(x = pValThDF, file = paste0("./Results/Hart/Threshold/",Today,"_Threshold_PPV_wilcoxon.csv"),
            quote = FALSE, row.names = TRUE)
  # Continue ####
  for(i in unique(PositivePredictiveValue_tmp$thMethod)){
    shtest <- shapiro.test(PositivePredictiveValue_tmp$value[PositivePredictiveValue_tmp$thMethod == i])
    varthr <- var(PositivePredictiveValue_tmp$value[PositivePredictiveValue_tmp$thMethod == i])
    if(i != "localT2" & i != "localT2HumanGEM"){
      cat(paste0("Threshold: ", sapply(strsplit(i,"_T"), "[[", 2), "\t Shapiro-Test: ", shtest$p.value, "\t Variance: ", varthr, "\n"))
    } else {
      cat(paste0("Threshold: ", i, "\t Shapiro-Test: ", shtest$p.value, "\t Variance: ", varthr, "\n"))
    }
  }
  
}
PositivePredictiveValue_tmp_mean <- PositivePredictiveValue_tmp %>% 
  group_by(thMethod) %>% 
  summarize(Mean = mean(value)) %>%
  ungroup()

PositivePredictiveValue_tmp <- rbind(PositivePredictiveValue_tmp, list("NA","Positive Predicted Value",NA))
change_lab <- "Positive Predictive Value"; names(change_lab) <- "Positive Predicted Value"
PPVlot <- ggplot(PositivePredictiveValue_tmp, aes(x = thMethod, y = value, fill = thMethod)) + 
  geom_boxplot(aes(fill = thMethod, alpha = 0.4)) + 
  scale_fill_manual(values = c("gMCS_T0"   = ThresholdColors[1],"gMCS_T1"   = ThresholdColors[2],
                               "gMCS_T2"   = ThresholdColors[3],"gMCS_T2.5" = ThresholdColors[4],
                               "gMCS_T5"   = ThresholdColors[5],"gMCS_T10"  = ThresholdColors[6],
                               "gMCS_T20"  = ThresholdColors[7],
                               "NA" = "#000000",
                               "localT2" = "#B2DF8A", "localT2HumanGEM" = "#6ba35b")) +
  geom_point(aes(color = thMethod, size = 2)) + 
  facet_wrap("variable",labeller = labeller(variable = change_lab))  +
  scale_color_manual(values = c("gMCS_T0"   = ThresholdColors[1],"gMCS_T1"   = ThresholdColors[2],
                                "gMCS_T2"   = ThresholdColors[3],"gMCS_T2.5" = ThresholdColors[4],
                                "gMCS_T5"   = ThresholdColors[5],"gMCS_T10"  = ThresholdColors[6],
                                "gMCS_T20"  = ThresholdColors[7],
                                "NA" = "#000000",
                                "localT2" = "#B2DF8A", "localT2HumanGEM" = "#6ba35b")) +
  theme_bw() + xlab("") + ylab("") +
  theme(legend.position = "none") + 
  geom_point(data    = PositivePredictiveValue_tmp_mean,
             mapping = aes(x = thMethod, y = Mean, size = 1.5),
             color   = "dark red") + 
  geom_line(data = PositivePredictiveValue_tmp_mean[grepl("gMCS", PositivePredictiveValue_tmp_mean$thMethod),], 
            mapping = aes(x = thMethod, y = Mean, group = 1, size = 1.25)) +
  geom_line(data = PositivePredictiveValue_tmp_mean[!grepl("gMCS", PositivePredictiveValue_tmp_mean$thMethod),], 
            mapping = aes(x = thMethod, y = Mean, group = 1, size = 1.25)) +
  geom_vline(xintercept = 8, size = 2, color = "dark red") +
  scale_x_discrete(labels=c("gMCS_T0"   = "0%",
                            "gMCS_T1"   = "1%",
                            "gMCS_T2"   = "2%",
                            "gMCS_T2.5" = "2.5%",
                            "gMCS_T5"   = "5%",
                            "gMCS_T10"  = "10%",
                            "gMCS_T20"  = "20%",
                            "NA" = "",
                            "localT2"   = "LocalT2",
                            "localT2HumanGEM"  = "LocalT2-H1"),
                   limits = c ("gMCS_T0","gMCS_T1","gMCS_T2","gMCS_T2.5",
                               "gMCS_T5","gMCS_T10","gMCS_T20", "NA",
                               "localT2","localT2HumanGEM")) +
  scale_y_continuous(limits = c(0,0.6), breaks = c(seq(from = 0, to = 1, by = 0.1))) +
  theme(strip.text  = element_text(face  = "bold", size = 28),
        strip.background = element_rect(fill = "#afc3dc"),
        axis.text.x = element_text(angle = 30, hjust = 1, vjust = 1, size = 20),
        axis.text.y = element_text(size = 23)) 
PPVlot
# Add All Three Together #
plotlistThreshold <-  list(TPlot, FPlot, PPVlot)
Boxplotting_Thresholds <- ggarrange(plotlist = plotlistThreshold,
                                    ncol = 3, nrow = 1)

# Boxplotting_Thresholds <- annotate_figure(p = Boxplotting_Thresholds,top = text_grob("Threshold Analysis", face = "bold", size = 18))
Boxplotting_Thresholds
ggsave(filename = paste0("./Results/Hart/",Today,"_Hart_Thresholds.png"),
       plot = Boxplotting_Thresholds, device = "png", width = 21, height = 8,
       units = "in",dpi = 300)

# Figure 1C - HART ####
ThresholdColors <- viridis::magma(8)
Results <- do.call(rbind, readRDS("./RDSResults/Hart/AllThresholdsResults.RDS"))
index <-  Results$thMethod == "gMCS_T0" |Results$thMethod == "gMCS_T1" | Results$thMethod == "gMCS_T2" | Results$thMethod == "gMCS_T2.5" | Results$thMethod == "gMCS_T5" | Results$thMethod == "gMCS_T10" | Results$thMethod == "gMCS_T20"
Results_gMCST  <- Results[index,]
tmp     <- rbind(Results_gMCST,
                 readRDS("./RDSResults/Hart/Results_localT2.RDS"),
                 readRDS("./RDSResults/Hart/Results_localT2HumanGEM.RDS"))
tmp <- reshape2::melt(tmp, id.vars = c("Order", "thMethod"),
                      measure.vars =  c("True Positives", "False Positives", "Positive Predicted Value"))

# Prepare TP Plot #
TruePositive_tmp          <- tmp %>% filter(variable == "True Positives")
pValThDF <- as.data.frame(matrix(1,nrow = length(unique(TruePositive_tmp$thMethod)),
                                 ncol = length(unique(TruePositive_tmp$thMethod))))
colnames(pValThDF) = rownames(pValThDF) = unique(TruePositive_tmp$thMethod)

for(ordering in unique(TruePositive_tmp$Order)){
  for(i in c("gMCS_T2", "gMCS_T5", "localT2")){
    for(j in setdiff(unique(TruePositive_tmp$thMethod),i)){
      tmp_TP <- TruePositive_tmp[TruePositive_tmp$Order == ordering,]
      pValThDF[i,j] <- wilcox.test(tmp_TP$value[tmp_TP$thMethod == i],
                                   tmp_TP$value[tmp_TP$thMethod == j],
                                   paired = TRUE)$p.value
      
    }
  }
}

TruePositive_tmp_mean <- TruePositive_tmp %>% 
  group_by(Order, thMethod) %>% 
  summarize(Mean = mean(value)) %>%
  ungroup()
if(is.na(sum(TruePositive_tmp_mean$Mean))){
  TruePositive_tmp_mean$Mean[which(is.na(TruePositive_tmp_mean$Mean))] <-  0
}
TruePositive_tmp_mean$variable <-  "True_Positives"
change_lab <- "True Positives"; names(change_lab) <- "True_Positives"
TruePositive_tmp_mean$thMethod <- factor(TruePositive_tmp_mean$thMethod, 
                                         levels = c("gMCS_T0", "gMCS_T1", "gMCS_T2", "gMCS_T2.5", "gMCS_T5", 
                                                    "gMCS_T10", "gMCS_T20", "localT2", "localT2HumanGEM"))
TPlot <- ggplot(TruePositive_tmp_mean, aes(x = Order, y = Mean, color = thMethod, group = thMethod)) + 
  geom_point(size = 4.5, alpha = 0.75) + geom_line(size = 1.5, alpha = 0.75) + 
  facet_wrap("variable",labeller = labeller(variable = change_lab))  +
  scale_color_manual(values = c("gMCS_T0"   = ThresholdColors[1],"gMCS_T1"   = ThresholdColors[2],
                                "gMCS_T2"   = ThresholdColors[3],"gMCS_T2.5" = ThresholdColors[4],
                                "gMCS_T5"   = ThresholdColors[5],"gMCS_T10"  = ThresholdColors[6],
                                "gMCS_T20"  = ThresholdColors[7],
                                "NA" = "#000000",
                                "localT2" = "#B2DF8A", "localT2HumanGEM" = "#6ba35b"),
                     labels = c("0%","1%","2%","2.5%","5%",
                                "10%","20%","localT2","localT2-H1")) +
  theme_bw() + xlab("") + ylab("") +
  guides(color=guide_legend(title="Threshold Method")) +
  theme(legend.position = "bottom") + 
  scale_x_discrete(labels=c("Order_1" = "1",
                            "Order_2" = "2",
                            "Order_3" = "3",
                            "Order_4" = "4",
                            "Order_5" = "5",
                            "Order_6" = "6",
                            "Order_7" = "7",
                            "HigherThan_Order_7" = ">7")) +
  scale_y_continuous(limits = c(78,105), breaks = c(seq(from = 78, to = 105, by = 3))) +
  theme(strip.text = element_text(face  = "bold", size = 28),
        strip.background = element_rect(fill = "#afc3dc"),
        axis.text.x = element_text(size = 23),
        axis.text.y = element_text(size = 23),
        legend.key.size = unit(3, 'cm'),
        legend.text = element_text(size = 18),
        legend.title    = element_text(size = 14))
TPlot
# Prepare FP Plot #
FalsePositive_tmp          <- tmp %>% filter(variable == "False Positives")

FalsePositive_tmp_mean <- FalsePositive_tmp %>% 
  group_by(Order, thMethod) %>% 
  summarize(Mean = mean(value)) %>%
  ungroup()

if(is.na(sum(FalsePositive_tmp_mean$Mean))){
  FalsePositive_tmp_mean$Mean[which(is.na(FalsePositive_tmp_mean$Mean))] <-  0
}
FalsePositive_tmp_mean$variable <-  "False_Positives"
FalsePositive_tmp_mean$thMethod <- factor(FalsePositive_tmp_mean$thMethod, 
                                          levels = c("gMCS_T0", "gMCS_T1", "gMCS_T2", "gMCS_T2.5", "gMCS_T5", 
                                                     "gMCS_T10", "gMCS_T20", "localT2", "localT2HumanGEM"))
change_lab <- "False Positives"; names(change_lab) <- "False_Positives"
FPlot <- ggplot(FalsePositive_tmp_mean, aes(x = Order, y = Mean, color = thMethod, group = thMethod)) + 
  geom_point(size = 4.5, alpha = 0.75) + geom_line(size = 1.5, alpha = 0.75) + 
  facet_wrap("variable",labeller = labeller(variable = change_lab))  +
  scale_color_manual(values = c("gMCS_T0"   = ThresholdColors[1],"gMCS_T1"   = ThresholdColors[2],
                                "gMCS_T2"   = ThresholdColors[3],"gMCS_T2.5" = ThresholdColors[4],
                                "gMCS_T5"   = ThresholdColors[5],"gMCS_T10"  = ThresholdColors[6],
                                "gMCS_T20"  = ThresholdColors[7],
                                "NA" = "#000000",
                                "localT2" = "#B2DF8A", "localT2HumanGEM" = "#6ba35b"),
                     labels = c("0%","1%","2%","2.5%","5%",
                                "10%","20%","localT2","localT2-H1")) +
  theme_bw() + xlab("") + ylab("") + 
  guides(color=guide_legend(title="Threshold Method")) +
  theme(legend.position = "bottom") + 
  scale_x_discrete(labels=c("Order_1" = "1",
                            "Order_2" = "2",
                            "Order_3" = "3",
                            "Order_4" = "4",
                            "Order_5" = "5",
                            "Order_6" = "6",
                            "Order_7" = "7",
                            "HigherThan_Order_7" = ">7")) +
  scale_y_continuous(limits = c(100,184), breaks = c(seq(from = 100, to = 185, by = 15))) +
  theme(strip.text = element_text(face  = "bold", size = 28), 
        strip.background = element_rect(fill = "#afc3dc"),
        axis.text.x = element_text(size = 23),
        axis.text.y = element_text(size = 23),
        legend.key.size = unit(3, 'cm'),
        legend.text     = element_text(size = 18),
        legend.title    = element_text(size = 14))
FPlot
# Prepare PPV Plot #
PositivePredictiveValue_tmp          <- tmp %>% filter(variable == "Positive Predicted Value")
PositivePredictiveValue_tmp_mean <- PositivePredictiveValue_tmp %>% 
  group_by(Order, thMethod) %>% 
  summarize(Mean = mean(value)) %>%
  ungroup()
if(is.na(sum(PositivePredictiveValue_tmp_mean$Mean))){
  PositivePredictiveValue_tmp_mean$Mean[which(is.na(PositivePredictiveValue_tmp_mean$Mean))] <-  0
}
PositivePredictiveValue_tmp_mean$variable <-  "Positive Predicted Value"
PositivePredictiveValue_tmp_mean$thMethod <- factor(PositivePredictiveValue_tmp_mean$thMethod, 
                                                    levels = c("gMCS_T0", "gMCS_T1", "gMCS_T2", "gMCS_T2.5", "gMCS_T5", 
                                                               "gMCS_T10", "gMCS_T20", "localT2", "localT2HumanGEM"))
change_lab <- "Positive Predictive Value"; names(change_lab) <- "Positive Predicted Value"
PPVlot <- ggplot(PositivePredictiveValue_tmp_mean, aes(x = Order, y = Mean, color = thMethod, group = thMethod)) + 
  geom_point(size = 4.5, alpha = 0.75) + geom_line(size = 1.5, alpha = 0.75) + 
  facet_wrap("variable",labeller = labeller(variable = change_lab))  +
  scale_color_manual(values = c("gMCS_T0"   = ThresholdColors[1],"gMCS_T1"   = ThresholdColors[2],
                                "gMCS_T2"   = ThresholdColors[3],"gMCS_T2.5" = ThresholdColors[4],
                                "gMCS_T5"   = ThresholdColors[5],"gMCS_T10"  = ThresholdColors[6],
                                "gMCS_T20"  = ThresholdColors[7],
                                "localT2" = "#B2DF8A", "localT2HumanGEM" = "#6ba35b"),
                     labels=c("0%","1%","2%","2.5%","5%",
                              "10%","20%","localT2","localT2-H1")) +
  theme_bw() + xlab("") + ylab("") + 
  guides(color=guide_legend(title="Threshold Method")) +
  theme(legend.position = "bottom") + 
  scale_x_discrete(labels=c("Order_1" = "1",
                            "Order_2" = "2",
                            "Order_3" = "3",
                            "Order_4" = "4",
                            "Order_5" = "5",
                            "Order_6" = "6",
                            "Order_7" = "7",
                            "HigherThan_Order_7" = ">7")) +
  scale_y_continuous(limits = c(0.36,0.44), breaks = c(seq(from = 0.35, to = 0.45, by = 0.01))) +
  theme(strip.text = element_text(face  = "bold", size = 28), 
        strip.background = element_rect(fill = "#afc3dc"),
        axis.text.x = element_text(size = 23),
        axis.text.y = element_text(size = 23),
        legend.key.size = unit(3, 'cm'),
        legend.text = element_text(size = 18),
        legend.title    = element_text(size = 14)) 

PPVlot
# Add All Three Together #
LinelistLength <-  list(TPlot, FPlot, PPVlot)
Lineplotting_Length <- ggarrange(plotlist = LinelistLength,
                                 ncol = 3, nrow = 1,
                                 common.legend =  TRUE)
Lineplotting_Length_Unlegend <- ggarrange(plotlist = LinelistLength,
                                          ncol = 3, nrow = 1,
                                          common.legend =  TRUE, legend = "none")
ggsave(filename = paste0("./Results/Hart/",Today,"_Hart_LineLengths.png"),
       plot = Lineplotting_Length, device = "png", width = 21, height = 8,units = "in",dpi = 300)
ggsave(filename = paste0("./Results/Hart/",Today,"_Hart_LineLengths_NoLegend.png"),
       plot = Lineplotting_Length_Unlegend, device = "png", width = 21, height = 8,units = "in",dpi = 300)



# Figure 2A - HART ####

toAdd <- do.call(rbind,readRDS("./RDSResults/Hart/AllThresholdsResultsNotCumulative.RDS"))
toKeep <- c("gMCS_T1","gMCS_T2","gMCS_T2.5","gMCS_T3.5","gMCS_T5","gMCS_T10","gMCS_T20")
toRemove <- setdiff(unique(toAdd$thMethod), toKeep)
myColors <- viridis::viridis(10)

for(element in toRemove){toAdd <- toAdd[which(toAdd$thMethod != element),]}
toAdd$thMethod <- factor(toAdd$thMethod, levels=toKeep)

tmp     <- rbind(toAdd,
                 readRDS("./RDSResults/Hart/Results_localT2.RDS"),
                 readRDS("./RDSResults/Hart/Results_localT2HumanGEM.RDS"))

tmp <-  tmp %>% filter(thMethod == "gMCS_T2")

tmp <- reshape2::melt(tmp, id.vars = c("Order"),
                      measure.vars =  c("True Positives", "False Positives", "Positive Predicted Value"))

# Prepare TP Plot #
TruePositive_tmp          <- tmp %>% filter(variable == "True Positives")
TruePositive_tmp_mean <- TruePositive_tmp %>% 
  group_by(Order) %>% 
  summarize(Mean = mean(value)) %>%
  ungroup()
if(is.na(sum(TruePositive_tmp_mean$Mean))){
  TruePositive_tmp_mean$Mean[which(is.na(TruePositive_tmp_mean$Mean))] <-  0
}
TruePositive_tmp_mean2 <- TruePositive_tmp_mean
colnames(TruePositive_tmp_mean2)[2] <- "value"
change_lab <- "True Positives"; names(change_lab) <- "True Positives"

TPlot <- ggplot(TruePositive_tmp, aes(x = Order, y = value, fill = Order)) + 
  geom_boxplot(aes(fill = Order, alpha = 0.4)) + 
  scale_fill_manual(values = c("Order_1" = myColors[1],
                               "Order_2" = myColors[2],
                               "Order_3" = myColors[3],
                               "Order_4" = myColors[4],
                               "Order_5" = myColors[5],
                               "Order_6" = myColors[6],
                               "Order_7" = myColors[7],
                               "HigherThan_Order_7" = myColors[8])) + 
  geom_point(aes(color = Order, size = 2)) + 
  facet_wrap("variable",labeller = labeller(variable = change_lab))  +
  scale_color_manual(values = c("Order_1" = myColors[1],
                                "Order_2" = myColors[2],
                                "Order_3" = myColors[3],
                                "Order_4" = myColors[4],
                                "Order_5" = myColors[5],
                                "Order_6" = myColors[6],
                                "Order_7" = myColors[7],
                                "HigherThan_Order_7" = myColors[8])) + 
  scale_y_continuous(limits = c(0,NA), breaks = c(seq(from = 0, to = 100, by = 10))) +
  theme_bw() + xlab("") + ylab("") +
  theme(legend.position = "none") + 
  geom_point(data    = TruePositive_tmp_mean,
             mapping = aes(x = Order, y = Mean, size = 1.5),
             color   = "dark red") + 
  geom_line(data = TruePositive_tmp_mean, 
            mapping = aes(x = Order, y = Mean, group = 1, size = 1.25)) +
  scale_x_discrete(labels=c("Order_1" = "1",
                            "Order_2" = "2",
                            "Order_3" = "3",
                            "Order_4" = "4",
                            "Order_5" = "5",
                            "Order_6" = "6",
                            "Order_7" = "7",
                            "HigherThan_Order_7" = ">7")) +
  theme(strip.text  = element_text(face  = "bold", size = 28), 
        strip.background = element_rect(fill = "#afc3dc"),
        axis.text.x = element_text(size = 23),
        axis.text.y = element_text(size = 23))

TPlot
# Prepare FP Plot #
FalsePositive_tmp          <- tmp %>% filter(variable == "False Positives")
FalsePositive_tmp_mean <- FalsePositive_tmp %>% 
  group_by(Order) %>% 
  summarize(Mean = mean(value)) %>%
  ungroup()
if(is.na(sum(FalsePositive_tmp_mean$Mean))){
  FalsePositive_tmp_mean$Mean[which(is.na(FalsePositive_tmp_mean$Mean))] <-  0
}

change_lab <- "False Positives"; names(change_lab) <- "False Positives"
FPlot <- ggplot(FalsePositive_tmp, aes(x = Order, y = value, fill = Order)) + 
  geom_boxplot(aes(fill = Order, alpha = 0.4)) + 
  scale_fill_manual(values = c("Order_1" = myColors[1],
                               "Order_2" = myColors[2],
                               "Order_3" = myColors[3],
                               "Order_4" = myColors[4],
                               "Order_5" = myColors[5],
                               "Order_6" = myColors[6],
                               "Order_7" = myColors[7],
                               "HigherThan_Order_7" = myColors[8])) + 
  geom_point(aes(color = Order, size = 2)) + 
  facet_wrap("variable",labeller = labeller(variable = change_lab))  +
  scale_color_manual(values = c("Order_1" = myColors[1],
                                "Order_2" = myColors[2],
                                "Order_3" = myColors[3],
                                "Order_4" = myColors[4],
                                "Order_5" = myColors[5],
                                "Order_6" = myColors[6],
                                "Order_7" = myColors[7],
                                "HigherThan_Order_7" = myColors[8])) + 
  theme_bw() + xlab("") + ylab("") +
  theme(legend.position = "none") + 
  scale_y_continuous(limits = c(0,NA), breaks = seq(from = 0, to = 150, by = 15)) +
  geom_point(data    = FalsePositive_tmp_mean,
             mapping = aes(x = Order, y = Mean, size = 1.5),
             color   = "dark red") + 
  geom_line(data = FalsePositive_tmp_mean, 
            mapping = aes(x = Order, y = Mean, group = 1, size = 1.25)) +
  scale_x_discrete(labels=c("Order_1" = "1",
                            "Order_2" = "2",
                            "Order_3" = "3",
                            "Order_4" = "4",
                            "Order_5" = "5",
                            "Order_6" = "6",
                            "Order_7" = "7",
                            "HigherThan_Order_7" = ">7")) +
  theme(strip.text   = element_text(face  = "bold", size = 28),
        strip.background = element_rect(fill = "#afc3dc"),
        axis.text.x  = element_text(size = 23),
        axis.text.y  = element_text(size = 23)) 
FPlot
# Prepare PPV Plot #
PositivePredictiveValue_tmp          <- tmp %>% filter(variable == "Positive Predicted Value")
PositivePredictiveValue_tmp_mean <- PositivePredictiveValue_tmp %>% 
  group_by(Order) %>% 
  summarize(Mean = mean(value, na.rm = TRUE)) %>%
  ungroup()
if(is.nan(sum(PositivePredictiveValue_tmp_mean$Mean))){
  PositivePredictiveValue_tmp_mean$Mean[which(is.nan(PositivePredictiveValue_tmp_mean$Mean))] <-  0
}
change_lab <- "Positive Predictive Value"; names(change_lab) <- "Positive Predicted Value"
PPVlot <- ggplot(PositivePredictiveValue_tmp, aes(x = Order, y = value, fill = Order)) + 
  geom_boxplot(aes(fill = Order, alpha = 0.4)) + 
  scale_fill_manual(values = c("Order_1" = myColors[1],
                               "Order_2" = myColors[2],
                               "Order_3" = myColors[3],
                               "Order_4" = myColors[4],
                               "Order_5" = myColors[5],
                               "Order_6" = myColors[6],
                               "Order_7" = myColors[7],
                               "HigherThan_Order_7" = myColors[8])) + 
  geom_point(aes(color = Order, size = 2)) + 
  scale_y_continuous(limits = c(0,NA), breaks = seq(from = 0, to = 1, by = 0.1)) +
  facet_wrap("variable",labeller = labeller(variable = change_lab))  +
  scale_color_manual(values = c("Order_1" = myColors[1],
                                "Order_2" = myColors[2],
                                "Order_3" = myColors[3],
                                "Order_4" = myColors[4],
                                "Order_5" = myColors[5],
                                "Order_6" = myColors[6],
                                "Order_7" = myColors[7],
                                "HigherThan_Order_7" = myColors[8])) + 
  theme_bw() + xlab("") + ylab("") +
  theme(legend.position = "none") + 
  geom_point(data    = PositivePredictiveValue_tmp_mean,
             mapping = aes(x = Order, y = Mean, size = 1.5),
             color   = "dark red") + 
  geom_line(data = PositivePredictiveValue_tmp_mean, 
            mapping = aes(x = Order, y = Mean, group = 1, size = 1.25)) +
  scale_x_discrete(labels=c("Order_1" = "1",
                            "Order_2" = "2",
                            "Order_3" = "3",
                            "Order_4" = "4",
                            "Order_5" = "5",
                            "Order_6" = "6",
                            "Order_7" = "7",
                            "HigherThan_Order_7" = ">7")) +
  theme(strip.text = element_text(face  = "bold", size = 28), 
        strip.background = element_rect(fill = "#afc3dc"),
        axis.text.x = element_text(size = 23),
        axis.text.y = element_text(size = 23)) 
PPVlot
# Add All Three Together #
plotlistLength <-  list(TPlot, FPlot, PPVlot)
Boxplotting_Length <- ggarrange(plotlist = plotlistLength,
                                ncol = 3, nrow = 1)
# Boxplotting_Length <- annotate_figure(p = Boxplotting_Length,top = text_grob("gMCST5 Length Analysis", face = "bold", size = 18))
Boxplotting_Length
ggsave(filename = paste0("./Results/Hart/",Today,"_Hart_Lengths_gMCSth2.png"),
       plot = Boxplotting_Length, device = "png", width = 21, height = 8,units = "in",dpi = 300)


# Figure 2B - HART ####

toAdd <- do.call(rbind,readRDS("./RDSResults/Hart/AllThresholdsResultsNotCumulative.RDS"))
toKeep <- c("gMCS_T1","gMCS_T2","gMCS_T2.5","gMCS_T3.5","gMCS_T5","gMCS_T10","gMCS_T20")
toRemove <- setdiff(unique(toAdd$thMethod), toKeep)
myColors <- viridis::viridis(10)

for(element in toRemove){toAdd <- toAdd[which(toAdd$thMethod != element),]}
toAdd$thMethod <- factor(toAdd$thMethod, levels=toKeep)

tmp     <- rbind(toAdd,
                 readRDS("./RDSResults/Hart/Results_localT2.RDS"),
                 readRDS("./RDSResults/Hart/Results_localT2HumanGEM.RDS"))

tmp <-  tmp %>% filter(thMethod == "gMCS_T5")

tmp <- reshape2::melt(tmp, id.vars = c("Order"),
                      measure.vars =  c("True Positives", "False Positives", "Positive Predicted Value"))

# Prepare TP Plot #
TruePositive_tmp          <- tmp %>% filter(variable == "True Positives")
TruePositive_tmp_mean <- TruePositive_tmp %>% 
  group_by(Order) %>% 
  summarize(Mean = mean(value)) %>%
  ungroup()
if(is.na(sum(TruePositive_tmp_mean$Mean))){
  TruePositive_tmp_mean$Mean[which(is.na(TruePositive_tmp_mean$Mean))] <-  0
}
TruePositive_tmp_mean2 <- TruePositive_tmp_mean
colnames(TruePositive_tmp_mean2)[2] <- "value"
change_lab <- "True Positives"; names(change_lab) <- "True Positives"

TPlot <- ggplot(TruePositive_tmp, aes(x = Order, y = value, fill = Order)) + 
  geom_boxplot(aes(fill = Order, alpha = 0.4)) + 
  scale_fill_manual(values = c("Order_1" = myColors[1],
                               "Order_2" = myColors[2],
                               "Order_3" = myColors[3],
                               "Order_4" = myColors[4],
                               "Order_5" = myColors[5],
                               "Order_6" = myColors[6],
                               "Order_7" = myColors[7],
                               "HigherThan_Order_7" = myColors[8])) + 
  geom_point(aes(color = Order, size = 2)) + 
  facet_wrap("variable",labeller = labeller(variable = change_lab))  +
  scale_color_manual(values = c("Order_1" = myColors[1],
                                "Order_2" = myColors[2],
                                "Order_3" = myColors[3],
                                "Order_4" = myColors[4],
                                "Order_5" = myColors[5],
                                "Order_6" = myColors[6],
                                "Order_7" = myColors[7],
                                "HigherThan_Order_7" = myColors[8])) + 
  scale_y_continuous(limits = c(0,NA), breaks = c(seq(from = 0, to = 100, by = 10))) +
  theme_bw() + xlab("") + ylab("") +
  theme(legend.position = "none") + 
  geom_point(data    = TruePositive_tmp_mean,
             mapping = aes(x = Order, y = Mean, size = 1.5),
             color   = "dark red") + 
  geom_line(data = TruePositive_tmp_mean, 
            mapping = aes(x = Order, y = Mean, group = 1, size = 1.25)) +
  scale_x_discrete(labels=c("Order_1" = "1",
                            "Order_2" = "2",
                            "Order_3" = "3",
                            "Order_4" = "4",
                            "Order_5" = "5",
                            "Order_6" = "6",
                            "Order_7" = "7",
                            "HigherThan_Order_7" = ">7")) +
  theme(strip.text  = element_text(face  = "bold", size = 28), 
        strip.background = element_rect(fill = "#afc3dc"),
        axis.text.x = element_text(size = 23),
        axis.text.y = element_text(size = 23))

TPlot
# Prepare FP Plot #
FalsePositive_tmp          <- tmp %>% filter(variable == "False Positives")
FalsePositive_tmp_mean <- FalsePositive_tmp %>% 
  group_by(Order) %>% 
  summarize(Mean = mean(value)) %>%
  ungroup()
if(is.na(sum(FalsePositive_tmp_mean$Mean))){
  FalsePositive_tmp_mean$Mean[which(is.na(FalsePositive_tmp_mean$Mean))] <-  0
}

change_lab <- "False Positives"; names(change_lab) <- "False Positives"
FPlot <- ggplot(FalsePositive_tmp, aes(x = Order, y = value, fill = Order)) + 
  geom_boxplot(aes(fill = Order, alpha = 0.4)) + 
  scale_fill_manual(values = c("Order_1" = myColors[1],
                               "Order_2" = myColors[2],
                               "Order_3" = myColors[3],
                               "Order_4" = myColors[4],
                               "Order_5" = myColors[5],
                               "Order_6" = myColors[6],
                               "Order_7" = myColors[7],
                               "HigherThan_Order_7" = myColors[8])) + 
  geom_point(aes(color = Order, size = 2)) + 
  facet_wrap("variable",labeller = labeller(variable = change_lab))  +
  scale_color_manual(values = c("Order_1" = myColors[1],
                                "Order_2" = myColors[2],
                                "Order_3" = myColors[3],
                                "Order_4" = myColors[4],
                                "Order_5" = myColors[5],
                                "Order_6" = myColors[6],
                                "Order_7" = myColors[7],
                                "HigherThan_Order_7" = myColors[8])) + 
  theme_bw() + xlab("") + ylab("") +
  theme(legend.position = "none") + 
  scale_y_continuous(limits = c(0,NA), breaks = seq(from = 0, to = 150, by = 15)) +
  geom_point(data    = FalsePositive_tmp_mean,
             mapping = aes(x = Order, y = Mean, size = 1.5),
             color   = "dark red") + 
  geom_line(data = FalsePositive_tmp_mean, 
            mapping = aes(x = Order, y = Mean, group = 1, size = 1.25)) +
  scale_x_discrete(labels=c("Order_1" = "1",
                            "Order_2" = "2",
                            "Order_3" = "3",
                            "Order_4" = "4",
                            "Order_5" = "5",
                            "Order_6" = "6",
                            "Order_7" = "7",
                            "HigherThan_Order_7" = ">7")) +
  theme(strip.text   = element_text(face  = "bold", size = 28),
        strip.background = element_rect(fill = "#afc3dc"),
        axis.text.x  = element_text(size = 23),
        axis.text.y  = element_text(size = 23)) 
FPlot
# Prepare PPV Plot #
PositivePredictiveValue_tmp          <- tmp %>% filter(variable == "Positive Predicted Value")
PositivePredictiveValue_tmp_mean <- PositivePredictiveValue_tmp %>% 
  group_by(Order) %>% 
  summarize(Mean = mean(value)) %>%
  ungroup()
if(is.na(sum(PositivePredictiveValue_tmp_mean$Mean))){
  PositivePredictiveValue_tmp_mean$Mean[which(is.na(PositivePredictiveValue_tmp_mean$Mean))] <-  0
}
change_lab <- "Positive Predictive Value"; names(change_lab) <- "Positive Predicted Value"
PPVlot <- ggplot(PositivePredictiveValue_tmp, aes(x = Order, y = value, fill = Order)) + 
  geom_boxplot(aes(fill = Order, alpha = 0.4)) + 
  scale_fill_manual(values = c("Order_1" = myColors[1],
                               "Order_2" = myColors[2],
                               "Order_3" = myColors[3],
                               "Order_4" = myColors[4],
                               "Order_5" = myColors[5],
                               "Order_6" = myColors[6],
                               "Order_7" = myColors[7],
                               "HigherThan_Order_7" = myColors[8])) + 
  geom_point(aes(color = Order, size = 2)) + 
  scale_y_continuous(limits = c(0,NA), breaks = seq(from = 0, to = 0.7, by = 0.1)) +
  facet_wrap("variable",labeller = labeller(variable = change_lab))  +
  scale_color_manual(values = c("Order_1" = myColors[1],
                                "Order_2" = myColors[2],
                                "Order_3" = myColors[3],
                                "Order_4" = myColors[4],
                                "Order_5" = myColors[5],
                                "Order_6" = myColors[6],
                                "Order_7" = myColors[7],
                                "HigherThan_Order_7" = myColors[8])) + 
  theme_bw() + xlab("") + ylab("") +
  theme(legend.position = "none") + 
  geom_point(data    = PositivePredictiveValue_tmp_mean,
             mapping = aes(x = Order, y = Mean, size = 1.5),
             color   = "dark red") + 
  geom_line(data = PositivePredictiveValue_tmp_mean, 
            mapping = aes(x = Order, y = Mean, group = 1, size = 1.25)) +
  scale_x_discrete(labels=c("Order_1" = "1",
                            "Order_2" = "2",
                            "Order_3" = "3",
                            "Order_4" = "4",
                            "Order_5" = "5",
                            "Order_6" = "6",
                            "Order_7" = "7",
                            "HigherThan_Order_7" = ">7")) +
  theme(strip.text = element_text(face  = "bold", size = 28), 
        strip.background = element_rect(fill = "#afc3dc"),
        axis.text.x = element_text(size = 23),
        axis.text.y = element_text(size = 23)) 
PPVlot
# Add All Three Together #
plotlistLength <-  list(TPlot, FPlot, PPVlot)
Boxplotting_Length <- ggarrange(plotlist = plotlistLength,
                                ncol = 3, nrow = 1)
# Boxplotting_Length <- annotate_figure(p = Boxplotting_Length,top = text_grob("gMCST5 Length Analysis", face = "bold", size = 18))
Boxplotting_Length
ggsave(filename = paste0("./Results/Hart/",Today,"_Hart_Lengths_gMCSth5.png"),
       plot = Boxplotting_Length, device = "png", width = 21, height = 8,units = "in",dpi = 300)


# Figure 2C - HART ####

tmp     <- readRDS("./RDSResults/Hart/NotCumulativeResults_localT2.RDS")

tmp$Order <- factor(tmp$Order,levels=c(paste0("Order_",1:7), "HigherThan_Order_7"))
tmp <- reshape2::melt(tmp, id.vars = c("Order"),
                      measure.vars =  c("True Positives", "False Positives", "Positive Predicted Value"))

# Prepare TP Plot #
TruePositive_tmp          <- tmp %>% filter(variable == "True Positives")
TruePositive_tmp_mean <- TruePositive_tmp %>% 
  group_by(Order) %>% 
  summarize(Mean = mean(value)) %>%
  ungroup()
if(is.na(sum(FalsePositive_tmp_mean$Mean))){
  FalsePositive_tmp_mean$Mean[which(is.na(FalsePositive_tmp_mean$Mean))] <-  0
}

change_lab <- "True Positives"; names(change_lab) <- "True Positives"
TPlot <- ggplot(TruePositive_tmp, aes(x = Order, y = value, fill = Order)) + 
  geom_boxplot(aes(fill = Order, alpha = 0.4)) + 
  scale_fill_manual(values = c("Order_1" = myColors[1],
                               "Order_2" = myColors[2],
                               "Order_3" = myColors[3],
                               "Order_4" = myColors[4],
                               "Order_5" = myColors[5],
                               "Order_6" = myColors[6],
                               "Order_7" = myColors[7],
                               "HigherThan_Order_7" = myColors[8])) + 
  geom_point(aes(color = Order, size = 2)) + 
  facet_wrap("variable",labeller = labeller(variable = change_lab))  +
  scale_color_manual(values = c("Order_1" = myColors[1],
                                "Order_2" = myColors[2],
                                "Order_3" = myColors[3],
                                "Order_4" = myColors[4],
                                "Order_5" = myColors[5],
                                "Order_6" = myColors[6],
                                "Order_7" = myColors[7],
                                "HigherThan_Order_7" = myColors[8])) + 
  theme_bw() + xlab("") + ylab("") +
  theme(legend.position = "none") + 
  geom_point(data    = TruePositive_tmp_mean,
             mapping = aes(x = Order, y = Mean, size = 1.5),
             color   = "dark red") + 
  geom_line(data = TruePositive_tmp_mean, 
            mapping = aes(x = Order, y = Mean, group = 1, size = 1.25)) +
  scale_x_discrete(labels=c("Order_1" = "1",
                            "Order_2" = "2",
                            "Order_3" = "3",
                            "Order_4" = "4",
                            "Order_5" = "5",
                            "Order_6" = "6",
                            "Order_7" = "7",
                            "HigherThan_Order_7" = ">7")) +
  scale_y_continuous(limits = c(0,NA), breaks = c(seq(from = 0, to = 100, by = 10))) +
  theme(strip.text = element_text(face  = "bold", size = 28), 
        strip.background = element_rect(fill = "#afc3dc"),
        axis.text.x = element_text(size = 23),
        axis.text.y = element_text(size = 23)) 
TPlot
# Prepare FP Plot #
FalsePositive_tmp          <- tmp %>% filter(variable == "False Positives")
FalsePositive_tmp_mean <- FalsePositive_tmp %>% 
  group_by(Order) %>% 
  summarize(Mean = mean(value)) %>%
  ungroup()
if(is.na(sum(FalsePositive_tmp_mean$Mean))){
  FalsePositive_tmp_mean$Mean[which(is.na(FalsePositive_tmp_mean$Mean))] <-  0
}

change_lab <- "False Positives"; names(change_lab) <- "False Positives"
FPlot <- ggplot(FalsePositive_tmp, aes(x = Order, y = value, fill = Order)) + 
  geom_boxplot(aes(fill = Order, alpha = 0.4)) + 
  scale_fill_manual(values = c("Order_1" = myColors[1],
                               "Order_2" = myColors[2],
                               "Order_3" = myColors[3],
                               "Order_4" = myColors[4],
                               "Order_5" = myColors[5],
                               "Order_6" = myColors[6],
                               "Order_7" = myColors[7],
                               "HigherThan_Order_7" = myColors[8])) + 
  geom_point(aes(color = Order, size = 2)) + 
  facet_wrap("variable",labeller = labeller(variable = change_lab))  +
  scale_color_manual(values = c("Order_1" = myColors[1],
                                "Order_2" = myColors[2],
                                "Order_3" = myColors[3],
                                "Order_4" = myColors[4],
                                "Order_5" = myColors[5],
                                "Order_6" = myColors[6],
                                "Order_7" = myColors[7],
                                "HigherThan_Order_7" = myColors[8])) + 
  theme_bw() + xlab("") + ylab("") +
  theme(legend.position = "none") + 
  geom_point(data    = FalsePositive_tmp_mean,
             mapping = aes(x = Order, y = Mean, size = 1.5),
             color   = "dark red") + 
  geom_line(data = FalsePositive_tmp_mean, 
            mapping = aes(x = Order, y = Mean, group = 1, size = 1.25)) +
  scale_x_discrete(labels=c("Order_1" = "1",
                            "Order_2" = "2",
                            "Order_3" = "3",
                            "Order_4" = "4",
                            "Order_5" = "5",
                            "Order_6" = "6",
                            "Order_7" = "7",
                            "HigherThan_Order_7" = ">7")) +
  scale_y_continuous(limits = c(0,NA), breaks = c(seq(from = 0, to = 140, by = 15))) +
  theme(strip.text   = element_text(face  = "bold", size = 28), 
        strip.background = element_rect(fill = "#afc3dc"),
        axis.text.x  = element_text(size = 23),
        axis.text.y  = element_text(size = 23)) 
FPlot
# Prepare PPV Plot #
PositivePredictiveValue_tmp          <- tmp %>% filter(variable == "Positive Predicted Value")
PositivePredictiveValue_tmp_mean <- PositivePredictiveValue_tmp %>% 
  group_by(Order) %>% 
  summarize(Mean = mean(value)) %>%
  ungroup()
if(is.na(sum(PositivePredictiveValue_tmp_mean$Mean))){
  PositivePredictiveValue_tmp_mean$Mean[which(is.na(PositivePredictiveValue_tmp_mean$Mean))] <-  0
}
change_lab <- "Positive Predictive Value"; names(change_lab) <- "Positive Predicted Value"
PPVlot <- ggplot(PositivePredictiveValue_tmp, aes(x = Order, y = value, fill = Order)) + 
  geom_boxplot(aes(fill = Order, alpha = 0.4)) + 
  scale_fill_manual(values = c("Order_1" = myColors[1],
                               "Order_2" = myColors[2],
                               "Order_3" = myColors[3],
                               "Order_4" = myColors[4],
                               "Order_5" = myColors[5],
                               "Order_6" = myColors[6],
                               "Order_7" = myColors[7],
                               "HigherThan_Order_7" = myColors[8])) + 
  geom_point(aes(color = Order, size = 2)) + 
  facet_wrap("variable",labeller = labeller(variable = change_lab))  +
  scale_color_manual(values = c("Order_1" = myColors[1],
                                "Order_2" = myColors[2],
                                "Order_3" = myColors[3],
                                "Order_4" = myColors[4],
                                "Order_5" = myColors[5],
                                "Order_6" = myColors[6],
                                "Order_7" = myColors[7],
                                "HigherThan_Order_7" = myColors[8])) + 
  theme_bw() + xlab("") + ylab("") +
  theme(legend.position = "none") + 
  geom_point(data    = PositivePredictiveValue_tmp_mean,
             mapping = aes(x = Order, y = Mean, size = 1.5),
             color   = "dark red") + 
  geom_line(data = PositivePredictiveValue_tmp_mean, 
            mapping = aes(x = Order, y = Mean, group = 1, size = 1.25)) +
  scale_x_discrete(labels=c("Order_1" = "1",
                            "Order_2" = "2",
                            "Order_3" = "3",
                            "Order_4" = "4",
                            "Order_5" = "5",
                            "Order_6" = "6",
                            "Order_7" = "7",
                            "HigherThan_Order_7" = ">7")) +
  scale_y_continuous(limits = c(0,0.667), breaks = c(seq(from = 0, to = 0.7, by = 0.10))) +
  theme(strip.text = element_text(face  = "bold", size = 28), 
        strip.background = element_rect(fill = "#afc3dc"),
        axis.text.x = element_text(size = 23),
        axis.text.y = element_text(size = 23)) 
PPVlot
# Add All Three Together #
plotlistLocalT2 <-  list(TPlot, FPlot, PPVlot)
Boxplotting_LocalT2 <- ggarrange(plotlist = plotlistLocalT2,
                                 ncol = 3, nrow = 1)
Boxplotting_LocalT2
# Boxplotting_LocalT2 <- annotate_figure(p = Boxplotting_LocalT2,top = text_grob("LocalT2 Length Analysis", face = "bold", size = 18))
ggsave(filename = paste0("./Results/Hart/",Today,"_Hart_Lengths_LocalT2.png"),
       plot = Boxplotting_LocalT2, device = "png", width = 21, height = 8,units = "in",dpi = 300)


# Figure 3A - HART ####
gMCS_to_Many_ResultsT2 <- readRDS("./RDSResults/Hart/One_to_Many_Results_gMCST2.RDS")
gMCS_to_Many_ResultsT2$thMethod <- "gMCSth2"
gMCS_to_Many_ResultsT5 <- readRDS("./RDSResults/Hart/One_to_Many_Results_gMCST5.RDS")
gMCS_to_Many_ResultsT5$thMethod <- "gMCSth5"
gMCS_to_Many_ResultslocalT2 <- readRDS("./RDSResults/Hart/One_to_Many_Results_localT2.RDS")
gMCS_to_Many_ResultslocalT2$thMethod <- "localT2"

AllgMCS_toMany <- rbind(gMCS_to_Many_ResultsT2,gMCS_to_Many_ResultsT5,gMCS_to_Many_ResultslocalT2)

colnames(AllgMCS_toMany)[grep("TP", colnames(AllgMCS_toMany))]  <- "True Positives"
colnames(AllgMCS_toMany)[grep("FP", colnames(AllgMCS_toMany))]  <- "False Positives"
colnames(AllgMCS_toMany)[grep("PPV",colnames(AllgMCS_toMany))] <- "Positive Predictive Value"
myColors <- viridis::viridis(4)

AllgMCS_toMany <- reshape2::melt(AllgMCS_toMany, id.vars = c("calcMethod", "thMethod"),
                                 measure.vars =  c("True Positives", "False Positives", "Positive Predictive Value"))
AllgMCS_toMany <- as_tibble(AllgMCS_toMany)
AllgMCS_toMany <- AllgMCS_toMany[AllgMCS_toMany$calcMethod != "TrueEssentials",]
AllgMCS_toMany$calcMethod <- factor(AllgMCS_toMany$calcMethod, 
                                    levels = c("One_Gene_to_One_gMCS_NoEssentials", "One_Gene_to_Many_gMCS_NoEssentials"))
## Prepare TP Plot ##
TruePositive_AllgMCS_toMany          <- AllgMCS_toMany %>% filter(variable == "True Positives")

# Compute p values #
if(FALSE){
  # t.test   ####
  for(i in c("gMCSth2", "gMCSth5", "localT2")){
    curTPDF <- TruePositive_AllgMCS_toMany[TruePositive_AllgMCS_toMany$thMethod == i,]
    print(round(t.test(curTPDF$value[curTPDF$calcMethod == "One_Gene_to_One_gMCS_NoEssentials"],
                       curTPDF$value[curTPDF$calcMethod == "One_Gene_to_Many_gMCS_NoEssentials"],
                       paired = TRUE)$p.value, 5))
  }
  
  # Wilcoxon ####
  for(i in c("gMCSth2", "gMCSth5", "localT2")){
    curTPDF <- TruePositive_AllgMCS_toMany[TruePositive_AllgMCS_toMany$thMethod == i,]
    print(round(wilcox.test(curTPDF$value[curTPDF$calcMethod == "One_Gene_to_One_gMCS_NoEssentials"],
                            curTPDF$value[curTPDF$calcMethod == "One_Gene_to_Many_gMCS_NoEssentials"],
                            paired = TRUE)$p.value, 5))
  }
  
  # Other2   ####
  for(i in c("gMCSth2", "gMCSth5", "localT2")){
    for(j in unique(TruePositive_AllgMCS_toMany$calcMethod)){
      shtest <- 0
      shtest <- try(shapiro.test(TruePositive_AllgMCS_toMany$value[TruePositive_AllgMCS_toMany$thMethod == i & TruePositive_AllgMCS_toMany$calcMethod == j])$p.value)
      varthr <- var(TruePositive_AllgMCS_toMany$value[TruePositive_AllgMCS_toMany$thMethod == i & TruePositive_AllgMCS_toMany$calcMethod == j])
      if(i != "localT2" & i != "localT2HumanGEM"){
        cat(paste0("Threshold: ", sapply(strsplit(i,"th"), "[[", 2), "\tMethod: ",j, "\t Shapiro-Test: ", shtest, "\t Variance: ", varthr, "\n"))
      } else {
        cat(paste0("Threshold: ", i, "\tMethod: ",j,"\tShapiro-Test: ", shtest, "\t Variance: ", varthr, "\n"))
      }
    }
  }
}
# Continue #
for(i in unique(TruePositive_AllgMCS_toMany$thMethod)){
  tmp <- TruePositive_AllgMCS_toMany[TruePositive_AllgMCS_toMany$thMethod == i,]
  tmp$value <- as.numeric(tmp$value)
  print(wilcox.test(value~calcMethod, data = tmp) )
}
TruePositive_AllgMCS_toMany %>% wilcox.test(value~thMethod, data=.)
TruePositive_AllgMCS_toMany_mean <- TruePositive_AllgMCS_toMany %>% 
  group_by(calcMethod) %>% 
  summarize(Mean = mean(value)) %>%
  ungroup()
if(is.na(sum(TruePositive_AllgMCS_toMany_mean$Mean))){
  TruePositive_AllgMCS_toMany_mean$Mean[which(is.na(TruePositive_AllgMCS_toMany_mean$Mean))] <-  0
}


change_lab <- "True Positives"; names(change_lab) <- "True Positives"
TPlot <- ggplot(TruePositive_AllgMCS_toMany, aes(x = calcMethod, y = value, fill = thMethod)) + 
  geom_boxplot(alpha = 0.6) + 
  geom_point(position=position_dodge(width=0.75),aes(group=thMethod, color = thMethod, size = 5, alpha = 0.5)) + 
  facet_wrap("variable",labeller = labeller(variable = change_lab))  +
  scale_fill_manual(values = c("gMCSth2" = viridis::magma(8)[3],
                               "gMCSth5" = viridis::magma(8)[5],
                               "localT2" = "#B2DF8A"), 
                    labels = c("gMCSth2"   = "2%","gMCSth5"   = "5%","localT2"   = "LocalT2"),
                    name = "Threshold Method") + 
  scale_color_manual(values = c("gMCSth2" = viridis::magma(8)[3],
                                "gMCSth5" = viridis::magma(8)[5],
                                "localT2" = "#B2DF8A"), 
                     labels = c("gMCSth2"   = "2%","gMCSth5"   = "5%","localT2"   = "LocalT2"),
                     guide = "none") + 
  scale_alpha(guide = "none") + scale_size(guide = "none") +
  theme_bw() + xlab("") + ylab("") +
  geom_vline(xintercept = 1.5, size = 2, color = "dark red") +
  theme(legend.position = "bottom") + 
  scale_x_discrete(labels=c("One_Gene_to_Many_gMCS_NoEssentials" = "One-to-Many",
                            "One_Gene_to_One_gMCS_NoEssentials" = "One-to-One")) +
  scale_y_continuous(limits = c(0,NA), breaks = c(seq(from = 0, to = 30, by = 5))) +
  theme(strip.text  = element_text(face  = "bold", size = 28), 
        strip.background = element_rect(fill = "#afc3dc"),
        axis.text.x = element_text(size = 23),
        axis.text.y = element_text(size = 23)) 

## Prepare FP Plot ##
FalsePositive_AllgMCS_toMany          <- AllgMCS_toMany %>% filter(variable == "False Positives")
# Compute pvalues #
if(FALSE){
  # t.test   ####
  for(i in c("gMCSth2", "gMCSth5", "localT2")){
    curFPDF <- FalsePositive_AllgMCS_toMany[FalsePositive_AllgMCS_toMany$thMethod == i,]
    print(round(t.test(curFPDF$value[curFPDF$calcMethod == "One_Gene_to_One_gMCS_NoEssentials"],
                       curFPDF$value[curFPDF$calcMethod == "One_Gene_to_Many_gMCS_NoEssentials"],
                       paired = TRUE)$p.value, 5))
  }
  
  # Wilcoxon ####
  for(i in c("gMCSth2", "gMCSth5", "localT2")){
    curFPDF <- FalsePositive_AllgMCS_toMany[FalsePositive_AllgMCS_toMany$thMethod == i,]
    print(round(wilcox.test(curFPDF$value[curFPDF$calcMethod == "One_Gene_to_One_gMCS_NoEssentials"],
                            curFPDF$value[curFPDF$calcMethod == "One_Gene_to_Many_gMCS_NoEssentials"],
                            paired = TRUE)$p.value, 5))
  }
  
  # Other2 ####
  for(i in c("gMCSth2", "gMCSth5", "localT2")){
    for(j in unique(FalsePositive_AllgMCS_toMany$calcMethod)){
      shtest <- 0
      shtest <- try(shapiro.test(FalsePositive_AllgMCS_toMany$value[FalsePositive_AllgMCS_toMany$thMethod == i & FalsePositive_AllgMCS_toMany$calcMethod == j])$p.value)
      varthr <- var(FalsePositive_AllgMCS_toMany$value[FalsePositive_AllgMCS_toMany$thMethod == i & FalsePositive_AllgMCS_toMany$calcMethod == j])
      if(i != "localT2" & i != "localT2HumanGEM"){
        cat(paste0("Threshold: ", sapply(strsplit(i,"th"), "[[", 2), "\tMethod: ",j, "\t Shapiro-Test: ", shtest, "\t Variance: ", varthr, "\n"))
      } else {
        cat(paste0("Threshold: ", i, "\tMethod: ",j,"\tShapiro-Test: ", shtest, "\t Variance: ", varthr, "\n"))
      }
    }
  }
}
# Continue #
FalsePositive_AllgMCS_toMany_mean <- FalsePositive_AllgMCS_toMany %>% 
  group_by(calcMethod) %>% 
  summarize(Mean = mean(value)) %>%
  ungroup()
if(is.na(sum(FalsePositive_AllgMCS_toMany_mean$Mean))){
  FalsePositive_AllgMCS_toMany_mean$Mean[which(is.na(FalsePositive_AllgMCS_toMany_mean$Mean))] <-  0
}

change_lab <- "False Positives"; names(change_lab) <- "False Positives"
FPlot <- ggplot(FalsePositive_AllgMCS_toMany, aes(x = calcMethod, y = value, fill = thMethod)) + 
  geom_boxplot(alpha = 0.6) + 
  geom_point(position=position_dodge(width=0.75),aes(group=thMethod, color = thMethod, size = 5, alpha = 0.5)) + 
  facet_wrap("variable",labeller = labeller(variable = change_lab))  +
  scale_fill_manual(values = c("gMCSth2" = viridis::magma(8)[3],
                               "gMCSth5" = viridis::magma(8)[5],
                               "localT2" = "#B2DF8A"), 
                    labels = c("gMCSth2"   = "2%","gMCSth5"   = "5%","localT2"   = "LocalT2"),
                    name = "Threshold Method") + 
  scale_color_manual(values = c("gMCSth2" = viridis::magma(8)[3],
                                "gMCSth5" = viridis::magma(8)[5],
                                "localT2" = "#B2DF8A"), 
                     labels = c("gMCSth2"   = "2%","gMCSth5"   = "5%","localT2"   = "LocalT2"),
                     guide = "none") + 
  scale_alpha(guide = "none") + scale_size(guide = "none") +
  theme_bw() + xlab("") + ylab("") +
  geom_vline(xintercept = 1.5, size = 2, color = "dark red") +
  theme(legend.position = "bottom") + 
  scale_x_discrete(labels=c("One_Gene_to_Many_gMCS_NoEssentials" = "One-to-Many",
                            "One_Gene_to_One_gMCS_NoEssentials" = "One-to-One")) +
  scale_y_continuous(limits = c(0,NA), breaks = c(seq(from = 0, to = 75, by = 15))) +
  theme(strip.text  = element_text(face  = "bold", size = 28), 
        strip.background = element_rect(fill = "#afc3dc"),
        axis.text.x = element_text(size = 23),
        axis.text.y = element_text(size = 23)) 
FPlot
## Prepare PPV Plot ##
PositivePredictiveValue_AllgMCS_toMany          <- AllgMCS_toMany %>% filter(variable == "Positive Predictive Value")
# Compute pvalues #
if(FALSE){
  # t.test   ####
  for(i in c("gMCSth2", "gMCSth5", "localT2")){
    curPPVDF <- PositivePredictiveValue_AllgMCS_toMany[PositivePredictiveValue_AllgMCS_toMany$thMethod == i,]
    print(round(t.test(curPPVDF$value[curPPVDF$calcMethod == "One_Gene_to_One_gMCS_NoEssentials"],
                       curPPVDF$value[curPPVDF$calcMethod == "One_Gene_to_Many_gMCS_NoEssentials"],
                       paired = TRUE)$p.value, 5))
  }
  # Wilcoxon ####
  for(i in c("gMCSth2", "gMCSth5", "localT2")){
    curPPVDF <- PositivePredictiveValue_AllgMCS_toMany[PositivePredictiveValue_AllgMCS_toMany$thMethod == i,]
    print(round(wilcox.test(curPPVDF$value[curPPVDF$calcMethod == "One_Gene_to_One_gMCS_NoEssentials"],
                            curPPVDF$value[curPPVDF$calcMethod == "One_Gene_to_Many_gMCS_NoEssentials"],
                            paired = TRUE)$p.value, 5))
  }
  # Other2 ####
  for(i in c("gMCSth2", "gMCSth5", "localT2")){
    for(j in unique(PositivePredictiveValue_AllgMCS_toMany$calcMethod)){
      shtest <- 0
      shtest <- try(shapiro.test(PositivePredictiveValue_AllgMCS_toMany$value[PositivePredictiveValue_AllgMCS_toMany$thMethod == i & PositivePredictiveValue_AllgMCS_toMany$calcMethod == j])$p.value)
      varthr <- var(PositivePredictiveValue_AllgMCS_toMany$value[PositivePredictiveValue_AllgMCS_toMany$thMethod == i & PositivePredictiveValue_AllgMCS_toMany$calcMethod == j])
      if(i != "localT2" & i != "localT2HumanGEM"){
        cat(paste0("Threshold: ", sapply(strsplit(i,"th"), "[[", 2), "\tMethod: ",j, "\t Shapiro-Test: ", shtest, "\t Variance: ", varthr, "\n"))
      } else {
        cat(paste0("Threshold: ", i, "\tMethod: ",j,"\tShapiro-Test: ", shtest, "\t Variance: ", varthr, "\n"))
      }
    }
  }
}
# Continue ##
for(i in unique(PositivePredictiveValue_AllgMCS_toMany$thMethod)){
  tmp <- PositivePredictiveValue_AllgMCS_toMany[PositivePredictiveValue_AllgMCS_toMany$thMethod == i,]
  tmp$value <- as.numeric(tmp$value)
  print(wilcox.test(value~calcMethod, data = tmp) )
}
PositivePredictiveValue_AllgMCS_toMany_mean <- PositivePredictiveValue_AllgMCS_toMany %>% 
  group_by(calcMethod) %>% 
  summarize(Mean = mean(value)) %>%
  ungroup()
if(is.na(sum(PositivePredictiveValue_AllgMCS_toMany_mean$Mean))){
  PositivePredictiveValue_AllgMCS_toMany_mean$Mean[which(is.na(PositivePredictiveValue_AllgMCS_toMany_mean$Mean))] <-  0
}
change_lab <- "Positive Predictive Value"; names(change_lab) <- "Positive Predictive Value"
PPVlot <-ggplot(PositivePredictiveValue_AllgMCS_toMany, aes(x = calcMethod, y = value, fill = thMethod)) + 
  geom_boxplot(alpha = 0.6) + 
  geom_point(position=position_dodge(width=0.75),aes(group=thMethod, color = thMethod, size = 5, alpha = 0.5)) + 
  facet_wrap("variable",labeller = labeller(variable = change_lab))  +
  scale_fill_manual(values = c("gMCSth2" = viridis::magma(8)[3],
                               "gMCSth5" = viridis::magma(8)[5],
                               "localT2" = "#B2DF8A"), 
                    labels = c("gMCSth2"   = "2%","gMCSth5"   = "5%","localT2"   = "LocalT2"),
                    name = "Threshold Method") + 
  scale_color_manual(values = c("gMCSth2" = viridis::magma(8)[3],
                                "gMCSth5" = viridis::magma(8)[5],
                                "localT2" = "#B2DF8A"), 
                     labels = c("gMCSth2"   = "2%","gMCSth5"   = "5%","localT2"   = "LocalT2"),
                     guide = "none") + 
  scale_alpha(guide = "none") + scale_size(guide = "none") +
  geom_vline(xintercept = 1.5, size = 2, color = "dark red") +
  theme_bw() + xlab("") + ylab("") +
  theme(legend.position = "bottom") + 
  scale_x_discrete(labels=c("One_Gene_to_Many_gMCS_NoEssentials" = "One-to-Many",
                            "One_Gene_to_One_gMCS_NoEssentials" = "One-to-One")) +
  scale_y_continuous(limits = c(0,NA), breaks = c(seq(from = 0, to = 1, by = 0.1))) +
  theme(strip.text  = element_text(face  = "bold", size = 28), 
        strip.background = element_rect(fill = "#afc3dc"),
        axis.text.x = element_text(size = 23),
        axis.text.y = element_text(size = 23)) 
PPVlot
# Add All Three Together #
plotlistgMCS <-  list(TPlot, FPlot, PPVlot)
Boxplotting_gMCS <- ggarrange(plotlist = plotlistgMCS,
                              ncol = 3, nrow = 1, common.legend = TRUE)
Boxplotting_gMCS_NoLegend <- ggarrange(plotlist = plotlistgMCS,
                                       ncol = 3, nrow = 1, legend = "none")
Boxplotting_gMCS
ggsave(filename = paste0("./Results/Hart/",Today,"_Hart_gMCS2Many_Figure.png"),
       plot = Boxplotting_gMCS, device = "png", width = 21, height = 8,units = "in",dpi = 300)
ggsave(filename = paste0("./Results/Hart/",Today,"_Hart_gMCS2Many_Figure_NoLegend.png"),
       plot = Boxplotting_gMCS_NoLegend, device = "png", width = 21, height = 8,units = "in",dpi = 300)

# Figure 3C - HART ####
Task_to_Many_ResultsT2 <- readRDS("./RDSResults/Hart/Task_Promiscuity_Results_gMCST2.RDS")
Task_to_Many_ResultsT2$thMethod <- "gMCSth2"
Task_to_Many_ResultsT5 <- readRDS("./RDSResults/Hart/Task_Promiscuity_Results_gMCST5.RDS")
Task_to_Many_ResultsT5$thMethod <- "gMCSth5"
Task_to_Many_ResultslocalT2 <- readRDS("./RDSResults/Hart/Task_Promiscuity_Results_localT2.RDS")
Task_to_Many_ResultslocalT2$thMethod <- "localT2"

AllTasks_toMany <- rbind(Task_to_Many_ResultsT2,Task_to_Many_ResultsT5,Task_to_Many_ResultslocalT2)

colnames(AllTasks_toMany)[grep("TP", colnames(AllTasks_toMany))]  <- "True Positives"
colnames(AllTasks_toMany)[grep("FP", colnames(AllTasks_toMany))]  <- "False Positives"
colnames(AllTasks_toMany)[grep("PPV",colnames(AllTasks_toMany))] <- "Positive Predictive Value"
myColors <- viridis::viridis(4)

AllTasks_toMany <- reshape2::melt(AllTasks_toMany, id.vars = c("calcMethod", "thMethod"),
                                  measure.vars =  c("True Positives", "False Positives", "Positive Predictive Value"))
AllTasks_toMany <- as_tibble(AllTasks_toMany)
AllTasks_toMany <- AllTasks_toMany[AllTasks_toMany$calcMethod != "TrueEssentials",]
AllTasks_toMany$calcMethod <- factor(AllTasks_toMany$calcMethod, 
                                     levels = c("One_Gene_to_One_Task_NoEssentials", "One_Gene_to_Many_Tasks_NoEssentials"))
## Prepare TP Plot ##
TruePositive_AllTasks_toMany          <- AllTasks_toMany %>% filter(variable == "True Positives")
# Compute pvalues #
if(FALSE){
  # t.test ####
  for(i in c("gMCSth2", "gMCSth5", "localT2")){
    curTPDF <- TruePositive_AllTasks_toMany[TruePositive_AllTasks_toMany$thMethod == i,]
    print(round(t.test(curTPDF$value[curTPDF$calcMethod == "One_Gene_to_One_Task_NoEssentials"],
                       curTPDF$value[curTPDF$calcMethod == "One_Gene_to_Many_Tasks_NoEssentials"],
                       paired = TRUE)$p.value, 5))
  }
  
  # wilcox ####
  for(i in c("gMCSth2", "gMCSth5", "localT2")){
    curTPDF <- TruePositive_AllTasks_toMany[TruePositive_AllTasks_toMany$thMethod == i,]
    print(round(wilcox.test(curTPDF$value[curTPDF$calcMethod == "One_Gene_to_One_Task_NoEssentials"],
                            curTPDF$value[curTPDF$calcMethod == "One_Gene_to_Many_Tasks_NoEssentials"],
                            paired = TRUE)$p.value, 5))
  }
  
  # Other2 ####
  for(i in c("gMCSth2", "gMCSth5", "localT2")){
    for(j in unique(TruePositive_AllTasks_toMany$calcMethod)){
      shtest <- 0
      shtest <- try(shapiro.test(TruePositive_AllTasks_toMany$value[TruePositive_AllTasks_toMany$thMethod == i & TruePositive_AllTasks_toMany$calcMethod == j])$p.value)
      varthr <- var(TruePositive_AllTasks_toMany$value[TruePositive_AllTasks_toMany$thMethod == i & TruePositive_AllTasks_toMany$calcMethod == j])
      if(i != "localT2" & i != "localT2HumanGEM"){
        cat(paste0("Threshold: ", sapply(strsplit(i,"th"), "[[", 2), "\tMethod: ",j, "\t Shapiro-Test: ", shtest, "\t Variance: ", varthr, "\n"))
      } else {
        cat(paste0("Threshold: ", i, "\tMethod: ",j,"\tShapiro-Test: ", shtest, "\t Variance: ", varthr, "\n"))
      }
    }
  }
}
# Continue #
TruePositive_AllTasks_toMany_mean <- TruePositive_AllTasks_toMany %>% 
  group_by(calcMethod) %>% 
  summarize(Mean = mean(value)) %>%
  ungroup()
if(is.na(sum(TruePositive_AllTasks_toMany_mean$Mean))){
  TruePositive_AllTasks_toMany_mean$Mean[which(is.na(TruePositive_AllTasks_toMany_mean$Mean))] <-  0
}


change_lab <- "True Positives"; names(change_lab) <- "True Positives"
TPlot <- ggplot(TruePositive_AllTasks_toMany, aes(x = calcMethod, y = value, fill = thMethod)) + 
  geom_boxplot(alpha = 0.6) + 
  geom_point(position=position_dodge(width=0.75),aes(group=thMethod, color = thMethod, size = 5, alpha = 0.5)) + 
  facet_wrap("variable",labeller = labeller(variable = change_lab))  +
  scale_fill_manual(values = c("gMCSth2" = viridis::magma(8)[3],
                               "gMCSth5" = viridis::magma(8)[5],
                               "localT2" = "#B2DF8A"), 
                    labels = c("gMCSth2"   = "2%","gMCSth5"   = "5%","localT2"   = "LocalT2"),
                    name = "Threshold Method") + 
  scale_color_manual(values = c("gMCSth2" = viridis::magma(8)[3],
                                "gMCSth5" = viridis::magma(8)[5],
                                "localT2" = "#B2DF8A"), 
                     labels = c("gMCSth2"   = "2%","gMCSth5"   = "5%","localT2"   = "LocalT2"),
                     guide = "none") + 
  scale_alpha(guide = "none") + scale_size(guide = "none") +
  geom_vline(xintercept = 1.5, size = 2, color = "dark red") +
  theme_bw() + xlab("") + ylab("") +
  theme(legend.position = "bottom") + 
  scale_x_discrete(labels=c("One_Gene_to_Many_Tasks_NoEssentials" = "Miscellaneous",
                            "One_Gene_to_One_Task_NoEssentials" = "Focused")) +
  scale_y_continuous(limits = c(0,NA), breaks = c(seq(from = 0, to = 45, by = 5))) +
  theme(strip.text  = element_text(face  = "bold", size = 28), 
        strip.background = element_rect(fill = "#afc3dc"),
        axis.text.x = element_text(size = 23),
        axis.text.y = element_text(size = 23)) 
TPlot
## Prepare FP Plot ##
FalsePositive_AllTasks_toMany          <- AllTasks_toMany %>% filter(variable == "False Positives")
if(FALSE){
  # t.test ####
  for(i in c("gMCSth2", "gMCSth5", "localT2")){
    curFPDF <- FalsePositive_AllTasks_toMany[FalsePositive_AllTasks_toMany$thMethod == i,]
    print(round(t.test(curFPDF$value[curFPDF$calcMethod == "One_Gene_to_One_Task_NoEssentials"],
                       curFPDF$value[curFPDF$calcMethod == "One_Gene_to_Many_Tasks_NoEssentials"],
                       paired = TRUE)$p.value, 5))
  }
  
  # Wilcox ####
  for(i in c("gMCSth2", "gMCSth5", "localT2")){
    curFPDF <- FalsePositive_AllTasks_toMany[FalsePositive_AllTasks_toMany$thMethod == i,]
    print(round(wilcox.test(curFPDF$value[curFPDF$calcMethod == "One_Gene_to_One_Task_NoEssentials"],
                            curFPDF$value[curFPDF$calcMethod == "One_Gene_to_Many_Tasks_NoEssentials"],
                            paired = TRUE)$p.value, 5))
  }
  # Other2 ####
  for(i in c("gMCSth2", "gMCSth5", "localT2")){
    for(j in unique(FalsePositive_AllTasks_toMany$calcMethod)){
      shtest <- 0
      shtest <- try(shapiro.test(FalsePositive_AllTasks_toMany$value[FalsePositive_AllTasks_toMany$thMethod == i & FalsePositive_AllTasks_toMany$calcMethod == j])$p.value)
      varthr <- var(FalsePositive_AllTasks_toMany$value[FalsePositive_AllTasks_toMany$thMethod == i & FalsePositive_AllTasks_toMany$calcMethod == j])
      if(i != "localT2" & i != "localT2HumanGEM"){
        cat(paste0("Threshold: ", sapply(strsplit(i,"th"), "[[", 2), "\tMethod: ",j, "\t Shapiro-Test: ", shtest, "\t Variance: ", varthr, "\n"))
      } else {
        cat(paste0("Threshold: ", i, "\tMethod: ",j,"\tShapiro-Test: ", shtest, "\t Variance: ", varthr, "\n"))
      }
    }
  }
}
# Continue #
FalsePositive_AllTasks_toMany_mean <- FalsePositive_AllTasks_toMany %>% 
  group_by(calcMethod) %>% 
  summarize(Mean = mean(value)) %>%
  ungroup()
if(is.na(sum(FalsePositive_AllTasks_toMany_mean$Mean))){
  FalsePositive_AllTasks_toMany_mean$Mean[which(is.na(FalsePositive_AllTasks_toMany_mean$Mean))] <-  0
}

change_lab <- "False Positives"; names(change_lab) <- "False Positives"
FPlot <- ggplot(FalsePositive_AllTasks_toMany, aes(x = calcMethod, y = value, fill = thMethod)) + 
  geom_boxplot(alpha = 0.6) + 
  geom_point(position=position_dodge(width=0.75),aes(group=thMethod, color = thMethod, size = 5, alpha = 0.5)) + 
  facet_wrap("variable",labeller = labeller(variable = change_lab))  +
  scale_fill_manual(values = c("gMCSth2" = viridis::magma(8)[3],
                               "gMCSth5" = viridis::magma(8)[5],
                               "localT2" = "#B2DF8A"), 
                    labels = c("gMCSth2"   = "2%","gMCSth5"   = "5%","localT2"   = "LocalT2"),
                    name = "Threshold Method") + 
  scale_color_manual(values = c("gMCSth2" = viridis::magma(8)[3],
                                "gMCSth5" = viridis::magma(8)[5],
                                "localT2" = "#B2DF8A"), 
                     labels = c("gMCSth2"   = "2%","gMCSth5"   = "5%","localT2"   = "LocalT2"),
                     guide = "none") + 
  scale_alpha(guide = "none") + scale_size(guide = "none") +
  geom_vline(xintercept = 1.5, size = 2, color = "dark red") +
  theme_bw() + xlab("") + ylab("") +
  theme(legend.position = "bottom") + 
  scale_x_discrete(labels=c("One_Gene_to_Many_Tasks_NoEssentials" = "Miscellaneous",
                            "One_Gene_to_One_Task_NoEssentials" = "Focused")) +
  scale_y_continuous(limits = c(0,NA), breaks = c(seq(from = 0, to = 115, by = 15))) +
  theme(strip.text  = element_text(face  = "bold", size = 28), 
        strip.background = element_rect(fill = "#afc3dc"),
        axis.text.x = element_text(size = 23),
        axis.text.y = element_text(size = 23)) 
FPlot
## Prepare PPV Plot ##
PositivePredictiveValue_AllTasks_toMany          <- AllTasks_toMany %>% filter(variable == "Positive Predictive Value")
# Compute pvalues #
if(FALSE){
  # t.test ####
  for(i in c("gMCSth2", "gMCSth5", "localT2")){
    curPPVDF <- PositivePredictiveValue_AllTasks_toMany[PositivePredictiveValue_AllTasks_toMany$thMethod == i,]
    print(round(t.test(curPPVDF$value[curPPVDF$calcMethod == "One_Gene_to_One_Task_NoEssentials"],
                       curPPVDF$value[curPPVDF$calcMethod == "One_Gene_to_Many_Tasks_NoEssentials"],
                       paired = TRUE)$p.value, 5))
  }
  
  # Wilcox ####
  for(i in c("gMCSth2", "gMCSth5", "localT2")){
    curPPVDF <- PositivePredictiveValue_AllTasks_toMany[PositivePredictiveValue_AllTasks_toMany$thMethod == i,]
    print(round(wilcox.test(curPPVDF$value[curPPVDF$calcMethod == "One_Gene_to_One_Task_NoEssentials"],
                            curPPVDF$value[curPPVDF$calcMethod == "One_Gene_to_Many_Tasks_NoEssentials"],
                            paired = TRUE)$p.value, 5))
  }
  
  # Other2 ####
  for(i in c("gMCSth2", "gMCSth5", "localT2")){
    for(j in unique(PositivePredictiveValue_AllTasks_toMany$calcMethod)){
      shtest <- 0
      shtest <- try(shapiro.test(PositivePredictiveValue_AllTasks_toMany$value[PositivePredictiveValue_AllTasks_toMany$thMethod == i & PositivePredictiveValue_AllTasks_toMany$calcMethod == j])$p.value)
      varthr <- var(PositivePredictiveValue_AllTasks_toMany$value[PositivePredictiveValue_AllTasks_toMany$thMethod == i & PositivePredictiveValue_AllTasks_toMany$calcMethod == j])
      if(i != "localT2" & i != "localT2HumanGEM"){
        cat(paste0("Threshold: ", sapply(strsplit(i,"th"), "[[", 2), "\tMethod: ",j, "\t Shapiro-Test: ", shtest, "\t Variance: ", varthr, "\n"))
      } else {
        cat(paste0("Threshold: ", i, "\tMethod: ",j,"\tShapiro-Test: ", shtest, "\t Variance: ", varthr, "\n"))
      }
    }
  }
}
# Continue #
PositivePredictiveValue_AllTasks_toMany_mean <- PositivePredictiveValue_AllTasks_toMany %>% 
  group_by(calcMethod) %>% 
  summarize(Mean = mean(value)) %>%
  ungroup()
if(is.na(sum(PositivePredictiveValue_AllTasks_toMany_mean$Mean))){
  PositivePredictiveValue_AllTasks_toMany_mean$Mean[which(is.na(PositivePredictiveValue_AllTasks_toMany_mean$Mean))] <-  0
}
change_lab <- "Positive Predictive Value"; names(change_lab) <- "Positive Predictive Value"
PPVlot <-ggplot(PositivePredictiveValue_AllTasks_toMany, aes(x = calcMethod, y = value, fill = thMethod)) + 
  geom_boxplot(alpha = 0.6) + 
  geom_point(position=position_dodge(width=0.75),aes(group=thMethod, color = thMethod, size = 5, alpha = 0.5)) + 
  facet_wrap("variable",labeller = labeller(variable = change_lab))  +
  scale_fill_manual(values = c("gMCSth2" = viridis::magma(8)[3],
                               "gMCSth5" = viridis::magma(8)[5],
                               "localT2" = "#B2DF8A"), 
                    labels = c("gMCSth2"   = "2%","gMCSth5"   = "5%","localT2"   = "LocalT2"),
                    name = "Threshold Method") + 
  scale_color_manual(values = c("gMCSth2" = viridis::magma(8)[3],
                                "gMCSth5" = viridis::magma(8)[5],
                                "localT2" = "#B2DF8A"), 
                     labels = c("gMCSth2"   = "2%","gMCSth5"   = "5%","localT2"   = "LocalT2"),
                     guide = "none") + 
  scale_alpha(guide = "none") + scale_size(guide = "none") +
  geom_vline(xintercept = 1.5, size = 2, color = "dark red") +
  theme_bw() + xlab("") + ylab("") +
  theme(legend.position = "bottom") + 
  scale_x_discrete(labels=c("One_Gene_to_Many_Tasks_NoEssentials" = "Miscellaneous",
                            "One_Gene_to_One_Task_NoEssentials" = "Focused")) +
  scale_y_continuous(limits = c(0,NA), breaks = c(seq(from = 0, to = 1, by = 0.1))) +
  theme(strip.text  = element_text(face  = "bold", size = 28), 
        strip.background = element_rect(fill = "#afc3dc"),
        axis.text.x = element_text(size = 23),
        axis.text.y = element_text(size = 23)) 
PPVlot
# Add All Three Together #
plotlistgMCS <-  list(TPlot, FPlot, PPVlot)
Boxplotting_gMCS <- ggarrange(plotlist = plotlistgMCS,
                              ncol = 3, nrow = 1, common.legend = TRUE)
Boxplotting_gMCS_NoLegend <- ggarrange(plotlist = plotlistgMCS,
                                       ncol = 3, nrow = 1, legend = "none")
Boxplotting_gMCS
ggsave(filename = paste0("./Results/Hart/",Today,"_Hart_Task_Figure.png"),
       plot = Boxplotting_gMCS, device = "png", width = 21, height = 8,units = "in",dpi = 300)
ggsave(filename = paste0("./Results/Hart/",Today,"_Hart_Task_Figure_NoLegend.png"),
       plot = Boxplotting_gMCS_NoLegend, device = "png", width = 21, height = 8,units = "in",dpi = 300)
# Figure 4A - HART ####
ThresholdColors <- viridis::magma(8)
BiomassLines  <- do.call(rbind, readRDS("./RDSResults/Hart/Biomass/AllThresholdsResults.RDS"))
BiomassLinesth     <- rbind(BiomassLines,
                            readRDS("./RDSResults/Hart/Biomass/Results_localT2.RDS"))
BiomassLinesth <- reshape2::melt(BiomassLinesth, id.vars = c("Order", "thMethod"),
                                 measure.vars =  c("True Positives", "False Positives", "Positive Predictive Value"))

BiomassLinesth$evalMethod <- "Biomass"
AllTasksLines <- do.call(rbind, readRDS("./RDSResults/Hart/AllThresholdsResults.RDS"))
AllTasksLines     <- rbind(AllTasksLines,
                           readRDS("./RDSResults/Hart/Results_localT2.RDS"),
                           readRDS("./RDSResults/Hart/Results_localT2HumanGEM.RDS"))
AllTasksLinesth <- reshape2::melt(AllTasksLines, id.vars = c("Order", "thMethod"),
                                  measure.vars =  c("True Positives", "False Positives", "Positive Predicted Value"))
AllTasksLinesth$evalMethod <- "AllTasks"

BothLinesth <- rbind(BiomassLinesth, AllTasksLinesth)
BothLinesth <- BothLinesth[BothLinesth$thMethod == "gMCS_T2" | BothLinesth$thMethod == "gMCS_T5" | BothLinesth$thMethod == "localT2",]
BothLinesth$eval_th <- paste0(BothLinesth$thMethod,"--",BothLinesth$evalMethod)
BothLinesth$Order <- as.character(BothLinesth$Order)
BothLinesth$Order[BothLinesth$thMethod == "gMCS_T5"] <- paste0(BothLinesth$Order[BothLinesth$thMethod == "gMCS_T5"], "_15")
BothLinesth$Order[BothLinesth$thMethod == "localT2"] <- paste0(BothLinesth$Order[BothLinesth$thMethod == "localT2"], "_30")
BothLinesth$Order <- factor(BothLinesth$Order,
                            levels = c("Order_1","Order_1_15", "Order_1_30", 
                                       "Order_2","Order_2_15", "Order_2_30", 
                                       "Order_3","Order_3_15", "Order_3_30", 
                                       "Order_4","Order_4_15", "Order_4_30", 
                                       "Order_5","Order_5_15", "Order_5_30", 
                                       "Order_6","Order_6_15", "Order_6_30", 
                                       "Order_7","Order_7_15", "Order_7_30", 
                                       "HigherThan_Order_7",
                                       "HigherThan_Order_7_15", 
                                       "HigherThan_Order_7_30"))
# Prepare TP Plot #
TruePositive_BothLinesth          <- BothLinesth %>% filter(variable == "True Positives")

TruePositive_BothLinesth_mean <- TruePositive_BothLinesth %>% 
  group_by(Order, eval_th, evalMethod, thMethod) %>% 
  summarize(Mean = mean(value)) %>%
  ungroup()
if(is.na(sum(TruePositive_BothLinesth_mean$Mean))){
  TruePositive_BothLinesth_mean$Mean[which(is.na(TruePositive_BothLinesth_mean$Mean))] <-  0
}
TruePositive_BothLinesth_mean$variable <-  "True_Positives"
change_lab <- "True Positives"; names(change_lab) <- "True_Positives"
TruePositive_BothLinesth_mean$evalMethod <- factor(TruePositive_BothLinesth_mean$evalMethod,
                                                   levels = c("Biomass", "AllTasks"))
for(thresholding in unique(TruePositive_BothLinesth_mean$thMethod)){
  for(ordering in unique(TruePositive_BothLinesth_mean$Order)){
    
    Total <- TruePositive_BothLinesth_mean$Mean[TruePositive_BothLinesth_mean$Order == ordering &
                                                  TruePositive_BothLinesth_mean$thMethod == thresholding &
                                                  TruePositive_BothLinesth_mean$evalMethod == "AllTasks"]
    Substract <- TruePositive_BothLinesth_mean$Mean[TruePositive_BothLinesth_mean$Order == ordering &
                                                      TruePositive_BothLinesth_mean$thMethod == thresholding &
                                                      TruePositive_BothLinesth_mean$evalMethod != "AllTasks"]
    TruePositive_BothLinesth_mean$Mean[TruePositive_BothLinesth_mean$Order == ordering &
                                         TruePositive_BothLinesth_mean$thMethod == thresholding &
                                         TruePositive_BothLinesth_mean$evalMethod == "AllTasks"] <- Total - Substract
    
  }  
}


TPlot <- ggplot(TruePositive_BothLinesth_mean, aes(x = Order, y = Mean,
                                                   alpha = evalMethod, color = evalMethod,
                                                   group = eval_th, fill = thMethod)) + 
  geom_bar(stat = "identity") + 
  facet_wrap("variable",labeller = labeller(variable = change_lab))  +
  scale_color_manual(values = c("AllTasks"   = "black", "Biomass"   = "grey20")) +
  scale_alpha_manual(values = c("AllTasks" = 0.2, "Biomass" = 0.8)) +
  theme_bw() + xlab("") + ylab("") +
  guides(color=guide_legend(title="Analysis")) +
  theme(legend.position = "bottom") + 
  scale_x_discrete(labels=c("Order_1" = "",
                            "Order_1_15" = "1",
                            "Order_1_30" = "",
                            "Order_1_45" = "",
                            "Order_2" = "",
                            "Order_2_15" = "2",
                            "Order_2_30" = "",
                            "Order_2_45" = "",
                            "Order_3" = "",
                            "Order_3_15" = "3",
                            "Order_3_30" = "",
                            "Order_3_45" = "",
                            "Order_4" = "",
                            "Order_4_15" = "4",
                            "Order_4_30" = "",
                            "Order_4_45" = "",
                            "Order_5" = "",
                            "Order_5_15" = "5",
                            "Order_5_30" = "",
                            "Order_5_45" = "",
                            "Order_6" = "",
                            "Order_6_15" = "6",
                            "Order_6_30" = "",
                            "Order_6_45" = "",
                            "Order_7" = "",
                            "Order_7_15" = "7",
                            "Order_7_30" = "",
                            "Order_7_45" = "",
                            "HigherThan_Order_7" = "",
                            "HigherThan_Order_7_15" = ">7",
                            "HigherThan_Order_7_30" = ""
  ),
  limits = c("Order_1","Order_1_15","Order_1_30",
             "Order_1_45",
             "Order_2","Order_2_15","Order_2_30",
             "Order_2_45",
             "Order_3","Order_3_15","Order_3_30",
             "Order_3_45",
             "Order_4","Order_4_15","Order_4_30",
             "Order_4_45",
             "Order_5","Order_5_15","Order_5_30",
             "Order_5_45",
             "Order_6","Order_6_15","Order_6_30",
             "Order_6_45",
             "Order_7","Order_7_15","Order_7_30",
             "Order_7_45",
             "HigherThan_Order_7","HigherThan_Order_7_15","HigherThan_Order_7_30")) +
  scale_fill_manual(values = c("gMCS_T2" = viridis::magma(8)[3],
                               "gMCS_T5" = viridis::magma(8)[5],
                               "localT2" = "#B2DF8A"), 
                    labels = c("gMCS_T2"   = "2%","gMCS_T5"   = "5%","localT2"   = "LocalT2"),
                    name = "Threshold Method") + 
  scale_y_continuous(limits = c(0,NA), breaks = c(seq(from = 0, to = 250, by = 15))) +
  theme(strip.text = element_text(face  = "bold", size = 50),
        strip.background = element_rect(fill = "#afc3dc"),
        axis.text.x = element_text(size = 45),
        axis.ticks.x = element_blank(),
        axis.text.y = element_text(size = 45),
        legend.position = "none",
        legend.text = element_text(size = 28),
        legend.title    = element_text(size = 28))
TPlot
ggsave(filename = paste0("./Results/Hart/",Today,"_Hart_Bio_vs_Tasks_BarPlot_TP_Biomass.png"),
       plot = TPlot, device = "png", width = 21, height = 12,units = "in",dpi = 300)
# Figure 4B - HART ####
# Prepare FP Plot #
FalsePositive_BothLinesth          <- BothLinesth %>% filter(variable == "False Positives")

FalsePositive_BothLinesth_mean <- FalsePositive_BothLinesth %>% 
  group_by(Order, eval_th, evalMethod, thMethod) %>% 
  summarize(Mean = mean(value)) %>%
  ungroup()

if(is.na(sum(FalsePositive_BothLinesth_mean$Mean))){
  FalsePositive_BothLinesth_mean$Mean[which(is.na(FalsePositive_BothLinesth_mean$Mean))] <-  0
}
FalsePositive_BothLinesth_mean$variable <-  "False_Positives"
change_lab <- "False Positives"; names(change_lab) <- "False_Positives"
FalsePositive_BothLinesth_mean$evalMethod <- factor(FalsePositive_BothLinesth_mean$evalMethod,
                                                    levels = c("Biomass", "AllTasks"))
for(thresholding in unique(FalsePositive_BothLinesth_mean$thMethod)){
  for(ordering in unique(FalsePositive_BothLinesth_mean$Order)){
    
    Total <- FalsePositive_BothLinesth_mean$Mean[FalsePositive_BothLinesth_mean$Order == ordering &
                                                   FalsePositive_BothLinesth_mean$thMethod == thresholding &
                                                   FalsePositive_BothLinesth_mean$evalMethod == "AllTasks"]
    Substract <- FalsePositive_BothLinesth_mean$Mean[FalsePositive_BothLinesth_mean$Order == ordering &
                                                       FalsePositive_BothLinesth_mean$thMethod == thresholding &
                                                       FalsePositive_BothLinesth_mean$evalMethod != "AllTasks"]
    FalsePositive_BothLinesth_mean$Mean[FalsePositive_BothLinesth_mean$Order == ordering &
                                          FalsePositive_BothLinesth_mean$thMethod == thresholding &
                                          FalsePositive_BothLinesth_mean$evalMethod == "AllTasks"] <- Total - Substract
    
  }  
}

FPlot <- ggplot(FalsePositive_BothLinesth_mean, aes(x = Order, y = Mean,
                                                    alpha = evalMethod, color = evalMethod,
                                                    group = eval_th, fill = thMethod)) + 
  geom_bar(stat = "identity") + 
  facet_wrap("variable",labeller = labeller(variable = change_lab))  +
  scale_color_manual(values = c("AllTasks"   = "black", "Biomass"   = "grey20")) +
  scale_alpha_manual(values = c("AllTasks" = 0.2, "Biomass" = 0.8)) +
  theme_bw() + xlab("") + ylab("") +
  guides(color=guide_legend(title="Analysis")) +
  theme(legend.position = "bottom") + 
  scale_x_discrete(labels=c("Order_1" = "",
                            "Order_1_15" = "1",
                            "Order_1_30" = "",
                            "Order_1_45" = "",
                            "Order_2" = "",
                            "Order_2_15" = "2",
                            "Order_2_30" = "",
                            "Order_2_45" = "",
                            "Order_3" = "",
                            "Order_3_15" = "3",
                            "Order_3_30" = "",
                            "Order_3_45" = "",
                            "Order_4" = "",
                            "Order_4_15" = "4",
                            "Order_4_30" = "",
                            "Order_4_45" = "",
                            "Order_5" = "",
                            "Order_5_15" = "5",
                            "Order_5_30" = "",
                            "Order_5_45" = "",
                            "Order_6" = "",
                            "Order_6_15" = "6",
                            "Order_6_30" = "",
                            "Order_6_45" = "",
                            "Order_7" = "",
                            "Order_7_15" = "7",
                            "Order_7_30" = "",
                            "Order_7_45" = "",
                            "HigherThan_Order_7" = "",
                            "HigherThan_Order_7_15" = ">7",
                            "HigherThan_Order_7_30" = ""
  ),
  limits = c("Order_1","Order_1_15","Order_1_30",
             "Order_1_45",
             "Order_2","Order_2_15","Order_2_30",
             "Order_2_45",
             "Order_3","Order_3_15","Order_3_30",
             "Order_3_45",
             "Order_4","Order_4_15","Order_4_30",
             "Order_4_45",
             "Order_5","Order_5_15","Order_5_30",
             "Order_5_45",
             "Order_6","Order_6_15","Order_6_30",
             "Order_6_45",
             "Order_7","Order_7_15","Order_7_30",
             "Order_7_45",
             "HigherThan_Order_7","HigherThan_Order_7_15","HigherThan_Order_7_30")) +
  scale_fill_manual(values = c("gMCS_T2" = viridis::magma(8)[3],
                               "gMCS_T5" = viridis::magma(8)[5],
                               "localT2" = "#B2DF8A"), 
                    labels = c("gMCS_T2"   = "2%","gMCS_T5"   = "5%","localT2"   = "LocalT2"),
                    name = "Threshold Method") + 
  scale_y_continuous(limits = c(0,NA), breaks = c(seq(from = 0, to = 500, by = 25))) +
  theme(strip.text = element_text(face  = "bold", size = 50),
        strip.background = element_rect(fill = "#afc3dc"),
        axis.text.x = element_text(size = 45),
        axis.ticks.x = element_blank(),
        axis.text.y = element_text(size = 45),
        legend.position = "none",
        legend.text = element_text(size = 18),
        legend.title    = element_text(size = 14))
FPlot
ggsave(filename = paste0("./Results/Hart/",Today,"_Hart_Bio_vs_Tasks_BarPlot_FP_NoLegend_Biomass.png"),
       plot = FPlot, device = "png", width = 21, height = 12,units = "in",dpi = 300)
# Supp - PPV       ####
PositivePredictedValue_BothLinesth          <- BothLinesth %>% filter(variable == "Positive Predicted Value" | variable == "Positive Predictive Value")
PositivePredictedValue_BothLinesth_mean <- PositivePredictedValue_BothLinesth %>% 
  group_by(Order, eval_th, evalMethod, thMethod) %>% 
  summarize(Mean = mean(value)) %>%
  ungroup()
if(is.na(sum(PositivePredictedValue_BothLinesth_mean$Mean))){
  PositivePredictedValue_BothLinesth_mean$Mean[which(is.na(PositivePredictedValue_BothLinesth_mean$Mean))] <-  0
}
PositivePredictedValue_BothLinesth_mean$variable <-  "Positive Predicted Value"
change_lab <- "Positive Predicted Value"; names(change_lab) <- "Positive Predicted Value"
PositivePredictedValue_BothLinesth_mean$evalMethod <- factor(PositivePredictedValue_BothLinesth_mean$evalMethod,
                                                             levels = c("AllTasks", "Biomass"))
PPVlot <- ggplot(PositivePredictedValue_BothLinesth_mean, aes(x = Order, y = Mean,
                                                              alpha = evalMethod, color = evalMethod,
                                                              group = eval_th, fill = thMethod)) + 
  geom_bar(stat = "identity", position = position_dodge()) + 
  facet_wrap("variable",labeller = labeller(variable = change_lab))  +
  scale_color_manual(values = c("gMCS_T2"   = ThresholdColors[3], "gMCS_T5"   = ThresholdColors[5],
                                "localT2" = "#B2DF8A")) +
  scale_alpha_manual(values = c("AllTasks" = 0.6, "Biomass" = 0.8)) +
  theme_bw() + xlab("") + ylab("") +
  theme(legend.position = "bottom") + 
  scale_x_discrete(labels=c("Order_1" = "",
                            "Order_1_15" = "1",
                            "Order_1_30" = "",
                            "Order_1_45" = "",
                            "Order_2" = "",
                            "Order_2_15" = "2",
                            "Order_2_30" = "",
                            "Order_2_45" = "",
                            "Order_3" = "",
                            "Order_3_15" = "3",
                            "Order_3_30" = "",
                            "Order_3_45" = "",
                            "Order_4" = "",
                            "Order_4_15" = "4",
                            "Order_4_30" = "",
                            "Order_4_45" = "",
                            "Order_5" = "",
                            "Order_5_15" = "5",
                            "Order_5_30" = "",
                            "Order_5_45" = "",
                            "Order_6" = "",
                            "Order_6_15" = "6",
                            "Order_6_30" = "",
                            "Order_6_45" = "",
                            "Order_7" = "",
                            "Order_7_15" = "7",
                            "Order_7_30" = "",
                            "Order_7_45" = "",
                            "HigherThan_Order_7" = "",
                            "HigherThan_Order_7_15" = ">7",
                            "HigherThan_Order_7_30" = ""
  ),
  limits = c("Order_1","Order_1_15","Order_1_30",
             "Order_1_45",
             "Order_2","Order_2_15","Order_2_30",
             "Order_2_45",
             "Order_3","Order_3_15","Order_3_30",
             "Order_3_45",
             "Order_4","Order_4_15","Order_4_30",
             "Order_4_45",
             "Order_5","Order_5_15","Order_5_30",
             "Order_5_45",
             "Order_6","Order_6_15","Order_6_30",
             "Order_6_45",
             "Order_7","Order_7_15","Order_7_30",
             "Order_7_45",
             "HigherThan_Order_7","HigherThan_Order_7_15","HigherThan_Order_7_30")) +
  scale_fill_manual(values = c("gMCS_T2" = viridis::magma(8)[3],
                               "gMCS_T5" = viridis::magma(8)[5],
                               "localT2" = "#B2DF8A"), 
                    labels = c("gMCS_T2"   = "2%","gMCS_T5"   = "5%","localT2"   = "LocalT2"),
                    name = "Threshold Method") + 
  scale_y_continuous(limits = c(0,NA), breaks = c(seq(from = 0, to = 1, by = 0.1))) +
  theme(strip.text = element_text(face  = "bold", size = 23),
        strip.background = element_rect(fill = "#afc3dc"),
        axis.text.x = element_text(size = 23),
        axis.ticks.x = element_blank(),
        axis.text.y = element_text(size = 23),
        legend.key.size = unit(3, 'cm'),
        legend.text = element_text(size = 18),
        legend.title    = element_text(size = 14))
PPVlot
ggsave(filename = paste0("./Results/Hart/",Today,"_Hart_Bio_vs_Tasks_BarPlot_PPV_NoLegend_Biomass.png"),
       plot = PPVlot, device = "png", width = 21, height = 12,units = "in",dpi = 300)
# Add All Three Together #
LinelistBioVSAll <-  list(TPlot, FPlot, PPVlot)
LineplottingBioVSAll <- ggarrange(plotlist = LinelistBioVSAll,
                                  ncol = 3, nrow = 1,
                                  common.legend =  TRUE, legend = "none")
LineplottingBioVSAll_NoLegend <- ggarrange(plotlist = LinelistBioVSAll,
                                           ncol = 3, nrow = 1,
                                           common.legend =  TRUE)

ggsave(filename = paste0("./Results/Hart/",Today,"_Hart_Bio_vs_Tasks_BarPlot_Figure_Biomass.png"),
       plot = LineplottingBioVSAll, device = "png", width = 21, height = 8,units = "in",dpi = 300)
ggsave(filename = paste0("./Results/Hart//",Today,"_Hart_Bio_vs_Tasks_BarPlot_Figure_NoLegend_Biomass.png"),
       plot = LineplottingBioVSAll_NoLegend, device = "png", width = 21, height = 8,units = "in",dpi = 300)
# Write Table #
toSave <- rbind(TruePositive_BothLinesth_mean, FalsePositive_BothLinesth_mean, PositivePredictedValue_BothLinesth_mean)
toSave$Order[toSave$Order == "Order_1_15" | toSave$Order == "Order_1_30"] <-  "Order_1"
toSave$Order[toSave$Order == "Order_2_15" | toSave$Order == "Order_2_30"] <-  "Order_2"
toSave$Order[toSave$Order == "Order_3_15" | toSave$Order == "Order_3_30"] <-  "Order_3"
toSave$Order[toSave$Order == "Order_4_15" | toSave$Order == "Order_4_30"] <-  "Order_4"
toSave$Order[toSave$Order == "Order_5_15" | toSave$Order == "Order_5_30"] <-  "Order_5"
toSave$Order[toSave$Order == "Order_6_15" | toSave$Order == "Order_6_30"] <-  "Order_6"
toSave$Order[toSave$Order == "Order_7_15" | toSave$Order == "Order_7_30"] <-  "Order_7"
toSave$Order[toSave$Order == "HigherThan_Order_7_15" | toSave$Order == "HigherThan_Order_7_30"] <-  "HigherThan_Order_7"
toSave$variable[toSave$variable == "Positive Predicted Value"] <- "Positive Predictive Value"
write.csv(toSave, file = "./Results/Hart/Bio_vs_Tasks/Hart_Biomass_AllTasks_Mean_Results_Biomass.csv",quote = FALSE, row.names = FALSE)
##################################### #####
############ DepMap 21Q4 ############ #####
##################################### #####
ThresholdVector <- c(0,1,2,2.5,5,10,20)

####    Prepare Data   ####  ####
# This data frame is required for future essentiality analysis, its only purpose is to deem all cells equal. #
Sample_Class_Dummy <- data.frame("Sample_Class" = "Cells")
levels(Sample_Class_Dummy)  <- "Cells"

## Get HumanGEM information ##
HumanGEM_Genes <- as.data.frame(fread("./Data/Genes_HumanGEM_v1_10_0.txt"))
Table_HumanGEM_Genes <- clusterProfiler::bitr(unlist(HumanGEM_Genes), "ENSEMBL", "SYMBOL", org.Hs.eg.db)

## Store function to compute Essential Genes ##
source('./Custom_Functions/DepMap/fun_All_Tasks_CalculateEssentialGenes_gmcsTH_DOM.R')
source('./Custom_Functions/DepMap/fun_TweakedDepMapgMCS_Results.R')
source('./Custom_Functions/ComputeBinaryMatrix.R')

# Load Gene Expression Data (TPM) from CCLE #
CCLE_Exp <- readRDS("./Data/DepMap/CCLE_Expression_TPM.RDS")
dictionary_CCLE <- as.data.frame(fread("./Data/DepMap/sample_info_21Q4.csv"))
rownames(dictionary_CCLE) <- dictionary_CCLE[,1]
CCLE_Effect <- as.data.frame(fread("./Data/DepMap/Achilles_gene_effect_21Q4.csv"))

## Load gMCSs ##
load('./Data/gMCSs_All_Cases.RData')
gMCS.info <- list(gMCS.info.raw$EssentialTasks_CultureBiomass,
                  gMCS.info.raw$Only_CultureMedium)
names(gMCS.info) <- c("EssentialTasks_CultureBiomass", "Only_CultureMedium")
rm(gMCS.info)
NumberOfgMCS <-  sapply(1:7, function(x) length(which(rowSums(gMCS.info$EssentialTasks_CultureBiomass$gMCSs.ENSEMBL.mat) == x)))

## Load Achilles ##
rownames(CCLE_Exp) <- CCLE_Exp[,1]
CCLE_Exp <- CCLE_Exp[,-1]
rownames(CCLE_Effect) <- CCLE_Effect[,1]
CCLE_Effect <- CCLE_Effect[,-1]

Achilles_Gene_Symbol <- sapply(strsplit(colnames(CCLE_Effect), " "), "[[", 1)
colnames(CCLE_Effect) <- Achilles_Gene_Symbol

CCLE_Effect <- t(CCLE_Effect)
CCLE_Exp    <- t(CCLE_Exp)
CCLE_Exp    <- CCLE_Exp[,intersect(colnames(CCLE_Exp),colnames(CCLE_Effect))]
CCLE_Effect <- CCLE_Effect[,intersect(colnames(CCLE_Exp),colnames(CCLE_Effect))]

gene_Names_Full <-  strsplit(rownames(CCLE_Exp), "[.]") # Separate SYMBOL & ENSEMBL
rownames(CCLE_Exp) <- sapply(gene_Names_Full, function(x) x[length(x)]) # Take Only ENSEMBL

CCLE_Effect_ENSEMBL <- as.data.frame(CCLE_Effect)
toENSEMBL           <- clusterProfiler::bitr(rownames(CCLE_Effect_ENSEMBL), "SYMBOL", "ENSEMBL", org.Hs.eg.db)
toENSEMBL <- toENSEMBL[!duplicated(toENSEMBL$ENSEMBL),]
CCLE_Effect_ENSEMBL$SYMBOL <- rownames(CCLE_Effect_ENSEMBL)
CCLE_Effect_ENSEMBL <- merge(CCLE_Effect_ENSEMBL, toENSEMBL, by = "SYMBOL")
rownames(CCLE_Effect_ENSEMBL) <- CCLE_Effect_ENSEMBL$ENSEMBL
CCLE_Effect_ENSEMBL <- CCLE_Effect_ENSEMBL[,-c(which(colnames(CCLE_Effect_ENSEMBL) == "SYMBOL"),
                                               which(colnames(CCLE_Effect_ENSEMBL) == "ENSEMBL"))]

####   Compute gmcTH   ####  ####
if("AllThresholdsResults.RDS" %in% dir("./RDSResults/DepMap")){
  ResultList <- readRDS("./RDSResults/DepMap/AllThresholdsResults.RDS")
} else {
  if("RoLAllThresholds.RDS" %in% dir("./RDSResults/DepMap/")){
    RoLResults <-  readRDS("./RDSResults/DepMap/RoLAllThresholds.RDS")
  } else {
    RoLResults <- vector("list", length(ThresholdVector))
    names(RoLResults) <- paste0("Th",ThresholdVector)
    for(Threshold in ThresholdVector){
      if (paste0("RoL_Th",Threshold,".RDS") %in% dir("./RDSResults/DepMap")){
        RoLResults[[paste0("Th",Threshold)]] <- readRDS(paste0("./RDSResults/DepMap/RoL_Th",Threshold,".RDS"))
      } else {
        RoLResults[[paste0("Th",Threshold)]] <- CalculateEssentialGenes_gmcsTH_DOM(gene.exp = CCLE_Exp,
                                                                                   gMCS_Analysis = gMCS.info$EssentialTasks_CultureBiomass,
                                                                                   sample.class = Sample_Class_Dummy,
                                                                                   gmcsTH_perc = Threshold/100,
                                                                                   nWorkers = 8,
                                                                                   gMCS_order = 8,
                                                                                   order1 = TRUE)
        saveRDS(RoLResults[[paste0("Th",Threshold)]], file = paste0("./RDSResults/DepMap/RoL_Th",Threshold,".RDS"))
      }
    }
    saveRDS(RoLResults, "./RDSResults/DepMap/RoLAllThresholds.RDS")
  }
  
  ResultList <- vector("list", length(ThresholdVector))
  names(ResultList) <- paste0("Th",ThresholdVector)
  
  for(Threshold in ThresholdVector){
    Binary_Matrix_Length_List <- ComputeBinaryMatrix(RoLResults[[paste0("Th",Threshold)]],
                                                     gMCS.info$EssentialTasks_CultureBiomass$gMCSs.ENSEMBL.list)
    ResultsInList <- lapply(1:8, function(x) TweakedDepMapgMCS_Results(Binary_Matrix = Binary_Matrix_Length_List[[x]],
                                                                       Effect = CCLE_Effect,
                                                                       Table_HumanGEM_Genes = Table_HumanGEM_Genes,
                                                                       gMCS_order =  x,
                                                                       Threshold_Info =  paste0("gMCS_T",Threshold)))
    ResultsInList <- do.call(rbind, ResultsInList)
    ResultsInList$Order <- factor(ResultsInList$Order,
                                  levels=c(paste0("Order_",1:7), "HigherThan_Order_7"))
    ResultsInList <-  ResultsInList %>% rename(TP = "True Positives", FP = "False Positives",
                                               TN = "True Negatives", FN = "False Negatives",
                                               PPV = "Positive Predicted Value",
                                               sensitivity = "Sensitivity", specificity = "Specificity",
                                               accuracy = "Accuracy", MCC = "Matthew's Cor. Coef.")
    ResultsInList %>% group_by(Order) %>%
      summarise_at(c("True Positives", "False Positives", 
                     "False Negatives", "True Negatives",
                     "Accuracy", "Sensitivity", "Specificity",
                     "Positive Predicted Value", "Matthew's Cor. Coef."), median, na.rm = TRUE) %>%
      write.csv(file = paste0(".Results/DepMap/gMCSLengths/LengthMedians_Th",Threshold,".csv"),
                row.names = F)
    ResultList[[paste0("Th",Threshold)]] <- ResultsInList
  }
  saveRDS(ResultList, file = "./RDSResults/DepMap/AllThresholdsResults.RDS")
}
## Not Cumulative ##
if("AllThresholdsResultsNotCumulative.RDS" %in% dir("./RDSResults/DepMap")){
  ResultList_NC <- readRDS("./RDSResults/DepMap/AllThresholdsResultsNotCumulative.RDS")
} else {
  ResultList_NC <- vector("list", length(ThresholdVector))
  names(ResultList_NC) <- paste0("Th",ThresholdVector)
  
  for(Threshold in ThresholdVector){
    Binary_Matrix_Length_List_NotCumulative <- vector("list", 8)
    Order1Genes <- unlist(gMCS.info$EssentialTasks_CultureBiomass$gMCSs.ENSEMBL.list[which(lengths(gMCS.info$EssentialTasks_CultureBiomass$gMCSs.ENSEMBL.list) == 1)])
    
    
    Binary_Matrix_Length_List <- ComputeBinaryMatrix(RoLResults[[paste0("Th",Threshold)]],
                                                     gMCS.info$EssentialTasks_CultureBiomass$gMCSs.ENSEMBL.list,
                                                     FALSE)
    ResultsInList <- lapply(1:8, function(x) TweakedDepMapgMCS_Results(Binary_Matrix = Binary_Matrix_Length_List[[x]],
                                                                       Effect = CCLE_Effect,
                                                                       Table_HumanGEM_Genes = Table_HumanGEM_Genes,
                                                                       gMCS_order =  x,
                                                                       Threshold_Info =  paste0("gMCS_T",Threshold)))
    ResultsInList <- do.call(rbind, ResultsInList)
    ResultsInList$Order <- factor(ResultsInList$Order,
                                  levels=c(paste0("Order_",1:7), "HigherThan_Order_7"))
    ResultsInList <-  ResultsInList %>% rename(TP = "True Positives", FP = "False Positives",
                                               TN = "True Negatives", FN = "False Negatives",
                                               PPV = "Positive Predicted Value",
                                               sensitivity = "Sensitivity", specificity = "Specificity",
                                               accuracy = "Accuracy", MCC = "Matthew's Cor. Coef.")
    ResultsInList %>% group_by(Order) %>%
      summarise_at(c("True Positives", "False Positives", 
                     "False Negatives", "True Negatives",
                     "Accuracy", "Sensitivity", "Specificity",
                     "Positive Predicted Value", "Matthew's Cor. Coef."), median, na.rm = TRUE) %>%
      write.csv(file = paste0(".Results/DepMap/NotCumulative/LengthMedians_Th",Threshold,".csv"),
                row.names = F)
    ResultList_NC[[paste0("Th",Threshold)]] <- ResultsInList
  }
  saveRDS(ResultList_NC, file = "./RDSResults/DepMap/AllThresholdsResultsNotCumulative.RDS")
}
############################### ####
####  LocalT2 - 1244   #### ####
source("./CustomFunctions/DepMap/fun_DepMapLocalT2CalculateEssential.R")
if ("Results_localT2.RDS" %in% dir("./RDSResults/DepMap/") && "NotCumulativeResults_localT2.RDS" %in% dir("./RDSResults/DepMap/")){
  Results_localT2 <- readRDS("./RDSResults/DepMap/Results_localT2.RDS")
  NotCumulativeResults_localT2 <- readRDS("./RDSResults/DepMap/NotCumulativeResults_localT2.RDS")
} else {
  if ("RoL_gMCS_localT2.RDS" %in% dir("./RDSResults/DepMap/")){
    RoL_gMCS_localT2 <- readRDS("./RDSResults/DepMap/RoL_gMCS_localT2.RDS")
    
  } else {
    Sample_Class_Dummy <- data.frame("Sample_Class" = colnames(CCLE_Exp))
    levels(Sample_Class_Dummy)  <- colnames(CCLE_Exp)
    Sample_Cohort_Dummy <- data.frame("Sample_Cohort" = "DepMap")
    levels(Sample_Cohort_Dummy) <- "DepMap"
    CCLE_Exp <- matrix(as.numeric(CCLE_Exp), nrow = nrow(CCLE_Exp), ncol = ncol(CCLE_Exp),
                       dimnames = list(rownames(CCLE_Exp), colnames(CCLE_Exp)))
    
    RoL_gMCS_localT2 <- LocalT2CalculateEssential(gene.exp      = CCLE_Exp,
                                                  gMCS.info     = gMCS.info$EssentialTasks_CultureBiomass,
                                                  sample.class  = Sample_Class_Dummy,
                                                  sample.cohort = Sample_Cohort_Dummy,
                                                  localT2_mode  = "all_genes_gMCSs",
                                                  nWorkers      = 8,
                                                  gMCS_order    = 8)
    
    
    saveRDS(object = RoL_gMCS_localT2, file = "./RDSResults/DepMap/RoL_gMCS_localT2.RDS")
  }
  # Cumulative     ####
  BinaryMatrixLengthList <- ComputeBinaryMatrix(RoL_gMCS_localT2,
                                                gMCS.info$EssentialTasks_CultureBiomass$gMCSs.ENSEMBL.list)
  Results_localT2 <- lapply(1:8, function(x) TweakedDepMapgMCS_Results(Binary_Matrix = BinaryMatrixLengthList[[x]],
                                                                       Effect = CCLE_Effect,
                                                                       Table_HumanGEM_Genes = Table_HumanGEM_Genes,
                                                                       gMCS_order =  x,
                                                                       Threshold_Info =  "localT2"))
  Results_localT2 <- do.call(rbind, Results_localT2)
  Results_localT2$Order <- factor(Results_localT2$Order,levels=c(paste0("Order_",1:7), "HigherThan_Order_7"))
  Results_localT2 <- Results_localT2 %>% dplyr::rename("True Positives"  = TP,
                                                       "False Positives" = FP,
                                                       "True Negatives"  = TN,
                                                       "False Negatives" = FN,
                                                       "Sensitivity" = sensitivity,
                                                       "Specificity" = specificity,
                                                       "Accuracy"    = accuracy,
                                                       "Positive Predicted Value" = PPV,
                                                       "Matthew's Cor. Coef."     = MCC)
  Results_localT2 %>% filter(thMethod == "localT2") %>% group_by(Order) %>%
    summarise_at(c("True Positives", "False Positives", 
                   "False Negatives", "True Negatives",
                   "Accuracy", "Sensitivity", "Specificity",
                   "Positive Predicted Value", "Matthew's Cor. Coef."), median, na.rm = TRUE) %>%
    write.csv(file = ".Results/DepMap/localT2/Cumulative_LocalT2.csv", 
              row.names = F)
  
  MeltedLocalT2 <- melt(data = Results_localT2,
                        id.vars = "Order",
                        measure.vars = c("True Positives", "False Positives",
                                         "False Negatives", "True Negatives",
                                         "Accuracy", "Sensitivity", "Specificity",
                                         "Positive Predicted Value", "Matthew's Cor. Coef."))
  LocalT2_Length <- ggplot(MeltedLocalT2, aes(x = Order, y = value, fill = Order)) + 
    facet_wrap(~variable,scales = "free") + geom_boxplot() + theme_bw() +  
    # geom_boxplot(width=0.1, fill = NA) +
    labs(title = "LocalT2_1244",
         x = "", y = "",
         fill = "Lengths") + 
    scale_fill_discrete(name = "Maximum gMCS\nLengths (#gMCS)",
                        labels = c(paste0("1 (n = ",NumberOfgMCS[1],")"), paste0("2 (n = ",NumberOfgMCS[2],")"),
                                   paste0("3 (n = ",NumberOfgMCS[3],")"), paste0("4 (n = ",NumberOfgMCS[4],")"),
                                   paste0("5 (n = ",NumberOfgMCS[5],")"), paste0("6 (n = ",NumberOfgMCS[6],")"),
                                   paste0("7 (n = ",NumberOfgMCS[7],")"),
                                   paste0("All (n = ",nrow(gMCS.info$EssentialTasks_CultureBiomass$gMCSs.ENSEMBL.mat),")")),
                        type = RColorBrewer::brewer.pal(8, "Paired")) +
    theme(plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
          axis.text.x = element_blank(), axis.ticks.x = element_blank(),
          strip.text.x = element_text(size = 18, face = "bold"),
          legend.position = "bottom")
  
  # LocalT2_Length
  
  ggsave(paste0(".Results/DepMap/localT2/CumulativeLocalT2.png"),
         plot = LocalT2_Length, device = "png", width = 16, height = 12, units = "in", dpi = 300)
  
  saveRDS(Results_localT2, "./RDSResults/DepMap/Results_localT2.RDS")
  
  
  # Not Cumulative ####
  NotCumulativeBinaryMatrixLengthList <- ComputeBinaryMatrix(RoL_gMCS_localT2,
                                                             gMCS.info$EssentialTasks_CultureBiomass$gMCSs.ENSEMBL.list,
                                                             FALSE)
  NotCumulativeResults_localT2 <- lapply(1:8, function(x) TweakedDepMapgMCS_Results(Binary_Matrix = NotCumulativeBinaryMatrixLengthList[[x]],
                                                                                    Effect = CCLE_Effect,
                                                                                    Table_HumanGEM_Genes = Table_HumanGEM_Genes,
                                                                                    gMCS_order =  x,
                                                                                    Threshold_Info =  "localT2"))
  NotCumulativeResults_localT2 <- do.call(rbind, NotCumulativeResults_localT2)
  NotCumulativeResults_localT2$Order <- factor(NotCumulativeResults_localT2$Order,
                                               levels=c(paste0("Order_",1:7), "HigherThan_Order_7"))
  NotCumulativeResults_localT2 <- NotCumulativeResults_localT2 %>% dplyr::rename("True Positives"  = TP,
                                                                                 "False Positives" = FP,
                                                                                 "True Negatives"  = TN,
                                                                                 "False Negatives" = FN,
                                                                                 "Sensitivity" = sensitivity,
                                                                                 "Specificity" = specificity,
                                                                                 "Accuracy"    = accuracy,
                                                                                 "Positive Predicted Value" = PPV,
                                                                                 "Matthew's Cor. Coef."     = MCC)
  NotCumulativeResults_localT2 %>% filter(thMethod == "localT2") %>% group_by(Order) %>%
    summarise_at(c("True Positives", "False Positives", 
                   "False Negatives", "True Negatives",
                   "Accuracy", "Sensitivity", "Specificity",
                   "Positive Predicted Value", "Matthew's Cor. Coef."), median, na.rm = TRUE) %>%
    write.csv(file = ".Results/DepMap/localT2/NotCumulative_LocalT2.csv", 
              row.names = F)
  
  NCMeltedLocalT2 <- melt(data = NotCumulativeResults_localT2,
                          id.vars = "Order",
                          measure.vars = c("True Positives", "False Positives",
                                           "False Negatives", "True Negatives",
                                           "Accuracy", "Sensitivity", "Specificity",
                                           "Positive Predicted Value", "Matthew's Cor. Coef."))
  NotCumLocalT2_Length <- ggplot(NCMeltedLocalT2, aes(x = Order, y = value, fill = Order)) + 
    facet_wrap(~variable,scales = "free") + geom_boxplot() + theme_bw() +  
    # geom_boxplot(width=0.1, fill = NA) +
    labs(title = "LocalT2_1244",
         x = "", y = "",
         fill = "Lengths") + 
    scale_fill_discrete(name = "Maximum gMCS\nLengths (#gMCS)",
                        labels = c(paste0("1 (n = ",NumberOfgMCS[1],")"), paste0("2 (n = ",NumberOfgMCS[2],")"),
                                   paste0("3 (n = ",NumberOfgMCS[3],")"), paste0("4 (n = ",NumberOfgMCS[4],")"),
                                   paste0("5 (n = ",NumberOfgMCS[5],")"), paste0("6 (n = ",NumberOfgMCS[6],")"),
                                   paste0("7 (n = ",NumberOfgMCS[7],")"),
                                   paste0("All (n = ",nrow(gMCS.info$EssentialTasks_CultureBiomass$gMCSs.ENSEMBL.mat),")")),
                        type = RColorBrewer::brewer.pal(8, "Paired")) +
    theme(plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
          axis.text.x = element_blank(), axis.ticks.x = element_blank(),
          strip.text.x = element_text(size = 18, face = "bold"),
          legend.position = "bottom")
  
  # NotCumLocalT2_Length
  
  ggsave(paste0(".Results/DepMap/localT2/NotCumulativeLocalT2.png"),
         plot = NotCumLocalT2_Length, device = "png", width = 16, height = 12, units = "in", dpi = 300)
  
  saveRDS(NotCumulativeResults_localT2, "./RDSResults/DepMap/NotCumulativeResults_localT2.RDS")
} 
####  LocalT2 - 3650   #### ####
if ("Results_localT2HumanGEM.RDS" %in% dir("./RDSResults/DepMap/") &&
    "NotCumulativeResults_localT2HumanGEM.RDS" %in% dir("./RDSResults/DepMap/")){
  Results_localT2HumanGEM <- readRDS("./RDSResults/DepMap/Results_localT2HumanGEM.RDS")
  NotCumulativeResults_localT2HumanGEM <- readRDS("./RDSResults/DepMap/NotCumulativeResults_localT2HumanGEM.RDS")
} else {
  if ("RoL_gMCS_localT2_HumanGEM.RDS" %in% dir("./RDSResults/DepMap/")){
    RoL_gMCS_localT2HumanGEM <- readRDS("./RDSResults/DepMap/RoL_gMCS_localT2_HumanGEM.RDS")
  } else {
    Sample_Class_Dummy <- data.frame("Sample_Class" = colnames(CCLE_Exp))
    levels(Sample_Class_Dummy)  <- colnames(CCLE_Exp)
    Sample_Cohort_Dummy <- data.frame("Sample_Cohort" = "DepMap")
    levels(Sample_Cohort_Dummy) <- "DepMap"
    CCLE_Exp <- matrix(as.numeric(CCLE_Exp), nrow = nrow(CCLE_Exp), ncol = ncol(CCLE_Exp),
                       dimnames = list(rownames(CCLE_Exp), colnames(CCLE_Exp)))
    RoL_gMCS_localT2_HumanGEM <- LocalT2CalculateEssential(gene.exp = CCLE_Exp,
                                                           gMCS.info = gMCS.info$EssentialTasks_CultureBiomass,
                                                           sample.class = Sample_Class_Dummy,
                                                           sample.cohort = Sample_Cohort_Dummy,
                                                           localT2_mode = "All_HumanGEM",
                                                           nWorkers = 8,
                                                           gMCS_order = 8)
    
  }
  saveRDS(object = RoL_gMCS_localT2_HumanGEM, file = "./RDSResults/DepMap/RoL_gMCS_localT2_HumanGEM.RDS")
  
  # Cumulative     ####                     
  BinaryMatrixLengthList <- ComputeBinaryMatrix(RoL_gMCS_localT2_HumanGEM,
                                                gMCS.info$EssentialTasks_CultureBiomass$gMCSs.ENSEMBL.list)
  Results_localT2HumanGEM <- lapply(1:8, function(x) TweakedDepMapgMCS_Results(Binary_Matrix = BinaryMatrixLengthList[[x]],
                                                                               Effect = CCLE_Effect,
                                                                               Table_HumanGEM_Genes = Table_HumanGEM_Genes,
                                                                               gMCS_order =  x,
                                                                               Threshold_Info =  "localT2HumanGEM"))
  Results_localT2HumanGEM <- do.call(rbind, Results_localT2HumanGEM)
  Results_localT2HumanGEM$Order <- factor(Results_localT2HumanGEM$Order,levels=c(paste0("Order_",1:7), "HigherThan_Order_7"))
  Results_localT2HumanGEM <- Results_localT2HumanGEM %>% dplyr::rename("True Positives"  = TP,
                                                                       "False Positives" = FP,
                                                                       "True Negatives"  = TN,
                                                                       "False Negatives" = FN,
                                                                       "Sensitivity" = sensitivity,
                                                                       "Specificity" = specificity,
                                                                       "Accuracy"    = accuracy,
                                                                       "Positive Predicted Value" = PPV,
                                                                       "Matthew's Cor. Coef."     = MCC)
  Results_localT2HumanGEM %>% filter(thMethod == "localT2") %>% group_by(Order) %>%
    summarise_at(c("True Positives", "False Positives", 
                   "False Negatives", "True Negatives",
                   "Accuracy", "Sensitivity", "Specificity",
                   "Positive Predicted Value", "Matthew's Cor. Coef."), median, na.rm = TRUE) %>%
    write.csv(file = ".Results/DepMap/localT2/Cumulative_LocalT2HumanGEM.csv", 
              row.names = F)
  
  MeltedLocalT2 <- melt(data = Results_localT2HumanGEM,
                        id.vars = "Order",
                        measure.vars = c("True Positives", "False Positives",
                                         "False Negatives", "True Negatives",
                                         "Accuracy", "Sensitivity", "Specificity",
                                         "Positive Predicted Value", "Matthew's Cor. Coef."))
  
  LocalT2_Length <- ggplot(MeltedLocalT2, aes(x = Order, y = value, fill = Order)) + 
    facet_wrap(~variable,scales = "free") + geom_boxplot() + theme_bw() +  
    # geom_boxplot(width=0.1, fill = NA) +
    labs(title = "LocalT2_3650",
         x = "", y = "",
         fill = "Lengths") + 
    scale_fill_discrete(name = "Maximum gMCS\nLengths (#gMCS)",
                        labels = c(paste0("1 (n = ",NumberOfgMCS[1],")"), paste0("2 (n = ",NumberOfgMCS[2],")"),
                                   paste0("3 (n = ",NumberOfgMCS[3],")"), paste0("4 (n = ",NumberOfgMCS[4],")"),
                                   paste0("5 (n = ",NumberOfgMCS[5],")"), paste0("6 (n = ",NumberOfgMCS[6],")"),
                                   paste0("7 (n = ",NumberOfgMCS[7],")"),
                                   paste0("All (n = ",nrow(gMCS.info$EssentialTasks_CultureBiomass$gMCSs.ENSEMBL.mat),")")),
                        type = RColorBrewer::brewer.pal(8, "Paired")) +
    theme(plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
          axis.text.x = element_blank(), axis.ticks.x = element_blank(),
          strip.text.x = element_text(size = 18, face = "bold"),
          legend.position = "bottom")
  
  # LocalT2_Length
  
  ggsave(paste0(".Results/DepMap/localT2/CumulativeLocalT2HumanGEM.png"),
         plot = LocalT2_Length, device = "png", width = 16, height = 12, units = "in", dpi = 300)
  
  saveRDS(Results_localT2HumanGEM, "./RDSResults/DepMap/Results_localT2HumanGEM.RDS")
  
  # Not Cumulative ####
  NotCumulativeBinaryMatrixLengthList <- ComputeBinaryMatrix(RoL_gMCS_localT2_HumanGEM,
                                                             gMCS.info$EssentialTasks_CultureBiomass$gMCSs.ENSEMBL.list,
                                                             FALSE)
  NotCumulativeResults_localT2HumanGEM <- lapply(1:8,
                                                 function(x) TweakedDepMapgMCS_Results(Binary_Matrix = NotCumulativeBinaryMatrixLengthList[[x]],
                                                                                       Effect = CCLE_Effect,
                                                                                       Table_HumanGEM_Genes = Table_HumanGEM_Genes,
                                                                                       gMCS_order =  x,
                                                                                       Threshold_Info =  "localT2HumanGEM"))
  NotCumulativeResults_localT2HumanGEM <- do.call(rbind, NotCumulativeResults_localT2HumanGEM)
  NotCumulativeResults_localT2HumanGEM$Order <- factor(NotCumulativeResults_localT2HumanGEM$Order,
                                                       levels=c(paste0("Order_",1:7), "HigherThan_Order_7"))
  NotCumulativeResults_localT2HumanGEM <- NotCumulativeResults_localT2HumanGEM %>% dplyr::rename("True Positives"  = TP,
                                                                                                 "False Positives" = FP,
                                                                                                 "True Negatives"  = TN,
                                                                                                 "False Negatives" = FN,
                                                                                                 "Sensitivity" = sensitivity,
                                                                                                 "Specificity" = specificity,
                                                                                                 "Accuracy"    = accuracy,
                                                                                                 "Positive Predicted Value" = PPV,
                                                                                                 "Matthew's Cor. Coef."     = MCC)
  NotCumulativeResults_localT2HumanGEM %>% filter(thMethod == "localT2") %>% group_by(Order) %>%
    summarise_at(c("True Positives", "False Positives", 
                   "False Negatives", "True Negatives",
                   "Accuracy", "Sensitivity", "Specificity",
                   "Positive Predicted Value", "Matthew's Cor. Coef."), median, na.rm = TRUE) %>%
    write.csv(file = ".Results/DepMap/localT2/NotCumulative_LocalT2HumanGEM.csv", 
              row.names = F)
  NotCumulativeMeltedLocalT2 <- melt(data = NotCumulativeResults_localT2HumanGEM,
                                     id.vars = "Order",
                                     measure.vars = c("True Positives", "False Positives",
                                                      "False Negatives", "True Negatives",
                                                      "Accuracy", "Sensitivity", "Specificity",
                                                      "Positive Predicted Value", "Matthew's Cor. Coef."))
  NotCumulativeLocalT2_Length <- ggplot(NotCumulativeMeltedLocalT2, aes(x = Order, y = value, fill = Order)) + 
    facet_wrap(~variable,scales = "free") + geom_boxplot() + theme_bw() +  
    # geom_boxplot(width=0.1, fill = NA) +
    labs(title = "LocalT2_3650",
         x = "", y = "",
         fill = "Lengths") + 
    scale_fill_discrete(name = "Maximum gMCS\nLengths (#gMCS)",
                        labels = c(paste0("1 (n = ",NumberOfgMCS[1],")"), paste0("2 (n = ",NumberOfgMCS[2],")"),
                                   paste0("3 (n = ",NumberOfgMCS[3],")"), paste0("4 (n = ",NumberOfgMCS[4],")"),
                                   paste0("5 (n = ",NumberOfgMCS[5],")"), paste0("6 (n = ",NumberOfgMCS[6],")"),
                                   paste0("7 (n = ",NumberOfgMCS[7],")"),
                                   paste0("All (n = ",nrow(gMCS.info$EssentialTasks_CultureBiomass$gMCSs.ENSEMBL.mat),")")),
                        type = RColorBrewer::brewer.pal(8, "Paired")) +
    theme(plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
          axis.text.x = element_blank(), axis.ticks.x = element_blank(),
          strip.text.x = element_text(size = 18, face = "bold"),
          legend.position = "bottom")
  
  # NotCumulativeLocalT2_Length
  
  ggsave(paste0(".Results/DepMap/localT2/NotCumulativeLocalT2HumanGEM.png"),
         plot = NotCumulativeLocalT2_Length, device = "png", width = 16, height = 12, units = "in", dpi = 300)
  
  saveRDS(NotCumulativeResults_localT2HumanGEM, "./RDSResults/DepMap/NotCumulativeResults_localT2HumanGEM.RDS")
}
############################### ####
#         gMCS-to-Many          ####
for(Threshold_Name in c("Th0","Th1","Th2","Th2_5","Th5","Th10","Th20","localT2","localT2_HumanGEM")){
  
  RoL_gMCS <- readRDS(paste0("./RDSResults/DepMap/RoL_",Threshold_Name,".RDS"))
  
  # Make gMCS-gene relations          ##
  Results_Order_List_Comparisons <- vector("list",3)
  names(Results_Order_List_Comparisons) <- c("One_Gene_to_One_g_NoEssentials",
                                             "One_Gene_to_Many_gMCSs_NoEssentials",
                                             "TrueEssentials")
  
  # Where is each gene Essential? #
  List_gMCS_Essential_Matrix <- RoL_gMCS$list.gMCS.essential.mat
  gMCS_Length <- 8
  Binary_Matrix <- matrix(0, nrow = nrow(RoL_gMCS$mat.essential.gene),
                          ncol = ncol(RoL_gMCS$mat.essential.gene),
                          dimnames = list(rownames(RoL_gMCS$mat.essential.gene),
                                          colnames(RoL_gMCS$mat.essential.gene)))
  NotBinary_Matrix <- Binary_Matrix
  
  for(Essential_Gene in names(List_gMCS_Essential_Matrix)){
    NotBinary_Matrix[Essential_Gene,] <- as.numeric(colSums(List_gMCS_Essential_Matrix[[Essential_Gene]]))
    Binary_Matrix[Essential_Gene,]    <- as.numeric(colSums(List_gMCS_Essential_Matrix[[Essential_Gene]]) > 0)
  }
  Essential_Gene_gMCS <-  do.call(rbind, lapply(List_gMCS_Essential_Matrix, colSums))
  # Extend data for all HumanGEM genes #
  Essential_Gene_gMCS_Extended <- matrix(data = 0, 
                                         nrow = length(unique(c(Table_HumanGEM_Genes$ENSEMBL,rownames(Essential_Gene_gMCS)))),
                                         ncol = ncol(Essential_Gene_gMCS),
                                         dimnames = list(unique(c(Table_HumanGEM_Genes$ENSEMBL,rownames(Essential_Gene_gMCS))),
                                                         colnames(Essential_Gene_gMCS)))
  # Essential genes are still Essential #
  Essential_Gene_gMCS_Extended[rownames(Essential_Gene_gMCS_Extended) %in%
                                 rownames(Essential_Gene_gMCS),] <- Essential_Gene_gMCS
  
  # Take Order 1 Genes #
  GetOrder1 <- rownames(Essential_Gene_gMCS_Extended)
  ThisAreOrder1 <- unname(unlist(
    gMCS.info$EssentialgMCSs_CultureBiomass$gMCSs.ENSEMBL.list[which(lengths(gMCS.info$EssentialgMCSs_CultureBiomass$gMCSs.ENSEMBL.list)==1)]))
  Order1_Indexes <- which(GetOrder1 %in% ThisAreOrder1)
  Essential_Gene_gMCS_Essentials <- Essential_Gene_gMCS_Extended
  for (Essential_Gene in Order1_Indexes){
    Essential_Gene_gMCS_Essentials[Essential_Gene, which(Essential_Gene_gMCS_Essentials[Essential_Gene,] != 0)] <- NA
  }
  
  CCLE_Effect_NoOrder1 <- CCLE_Effect_ENSEMBL[!rownames(CCLE_Effect_ENSEMBL) %in% ThisAreOrder1,]
  
  # Make 1 gMCS to 1 gene             ##
  One_to_One_Mat_NoEs <-  NotBinary_Matrix
  One_to_One_Mat_NoEs[One_to_One_Mat_NoEs != 1] <- NA # Caution! Changed! (Previously 0)
  Results_Order_List_Comparisons[[1]] <- One_to_One_Mat_NoEs
  
  
  # Make n gMCS to 1 gene             ##
  Many_to_One_Mat_Ess <-  NotBinary_Matrix
  Many_to_One_Mat_Ess[Many_to_One_Mat_Ess <= 1]  <- NA # Caution! Changed! (Previously 0)
  Many_to_One_Mat_Ess[Many_to_One_Mat_Ess > 1]   <- 1
  Results_Order_List_Comparisons[[2]] <- Many_to_One_Mat_Ess
  
  # All Essential Genes               ##
  Results_Order_List_Comparisons[[3]] <- Binary_Matrix
  mean(colSums(One_to_One_Mat_NoEs, na.rm =TRUE));mean(colSums(Many_to_One_Mat_Ess, na.rm =TRUE))
  # Compute Statistics                ##
  
  Comparison_Analysis_Results <- vector("list", length(Results_Order_List_Comparisons))
  k <- 0
  for (Analysis in names(Results_Order_List_Comparisons)){
    try({result_list <- list()
    k <- k+1
    evalMethod = "All_gMCSs"
    
    if (Analysis == "TrueEssentials"){
      Mat_binary <- CCLE_Effect_ENSEMBL < (-0.6)
      Mat_binary[which(is.na(Mat_binary))] <- 0
    } else {
      Mat_binary <- CCLE_Effect_NoOrder1 < (-0.6)
      Mat_binary[which(is.na(Mat_binary))] <- 0
      
    }
    
    tissues <- intersect(colnames(Results_Order_List_Comparisons[[Analysis]]),
                         colnames(Mat_binary))
    
    # Initialize Variables #
    results = data.frame('cellLine' = tissues,
                         'TP' = NA, 'TN' = NA, 'FP' = NA, 'FN' = NA,
                         "check_sum_genes" = NA, 
                         'accuracy' = NA, 'sensitivity' = NA,
                         'specificity' = NA, 'F1' = NA,
                         'MCC' = NA, 'Penr' = NA, 'logPenr' = NA,
                         'PenrAdj' = NA, 'logPenrAdj' = NA)
    
    for (t in 1:length(tissues)){
      
      # Genes From the Model #
      modelGenes <- unique(Table_HumanGEM_Genes$ENSEMBL) 
      # Essential Genes #
      modelEssential <- rownames(Results_Order_List_Comparisons[[Analysis]])[Results_Order_List_Comparisons[[Analysis]][,tissues[t]]>0]
      modelEssential <- unique(Table_HumanGEM_Genes$ENSEMBL[Table_HumanGEM_Genes$ENSEMBL %in% modelEssential])
      # Non-Essential Genes #
      modelNonEssential <- setdiff(modelGenes, modelEssential)
      
      # expGenes has all the genes that are not NA from Achilles #
      expGenes <-  rownames(Mat_binary)[!is.na(Mat_binary[,tissues[t]])]
      # expEssential has all the genes that are considered essential #
      expEssential <- expGenes[Mat_binary[expGenes,tissues[t]]==1]
      # expNonEssential has all the genes that are NOT essential # 
      expNonEssential <- setdiff(expGenes,expEssential)
      
      results$TP[t] = length(intersect(modelEssential, expEssential))        # True Positives  #
      results$TN[t] = length(intersect(modelNonEssential, expNonEssential))  # True Negatives  #
      results$FP[t] = length(intersect(modelEssential, expNonEssential))     # False Positives #
      results$FN[t] = length(intersect(modelNonEssential, expEssential))     # False Negatives #   
      
      # Intersect Genes from HumanGEM and Genes from ACHILLES #
      pop <- intersect(modelGenes, expGenes)
      # Essential Genes from gMCST0 #
      sample <- intersect(modelEssential, expGenes)
      # Essential Genes from ACHILLES #
      successes <- intersect(expEssential, modelGenes)
      
      results$Penr[t] <- phyper(length(intersect(successes,sample)) - 1, # Number of White Balls Drawn  #
                                length(intersect(successes,pop)),        # Number of White Balls in Urn #
                                length(setdiff(pop, successes)),         # Number of Black Balls in Urn #
                                length(sample),                          # Total Number of Balls Drawn  #  
                                lower.tail = FALSE);
      
    }
    
    
    results$check_sum_genes =  results$TP + results$TN + results$FP + results$FN
    
    # Calculate Metrics
    results$sensitivity = results$TP/(results$TP + results$FN)
    results$specificity = results$TN/(results$TN + results$FP)
    results$accuracy    = (results$TP + results$TN)/(results$TP + results$TN + results$FP + results$FN)
    results$PPV         = results$TP/(results$TP+results$FP)
    results$F1          = 2*results$TP/(2*results$TP + results$FP + results$FN)
    results$MCC         = ((results$TP*results$TN) - (results$FP*results$FN))/(sqrt((results$TP+results$FP)*(results$TP+results$FN))*sqrt((results$TN+results$FP)*(results$TN+results$FN)))  # Matthews correlation coefficient
    results$PenrAdj     = p.adjust(results$Penr,'BH')
    
    results$logPenr     = -log10(results$Penr)
    results$logPenrAdj  = -log10(results$PenrAdj)
    
    results$model <- "DepMap"
    results$evalMethod <- evalMethod
    
    results$calcMethod <- Analysis
    
    results$PenrAdj = p.adjust(results$Penr,'BH')
    
    results$logPenr = -log10(results$Penr)
    results$logPenrAdj = -log10(results$PenrAdj)
    
    result_list[[evalMethod]] <- results
    
    results <- do.call(rbind, result_list)})
    
    Comparison_Analysis_Results[[k]] <- results
    
  }
  
  gMCS_Statistics <- do.call(rbind, Comparison_Analysis_Results)
  saveRDS(gMCS_Statistics, paste0("./RDSResults/DepMap/",Today,"_One_to_Many_Results_gMCS",ThresholdName,".RDS"))
}
############################### ####
#         Task-to-Many          ####
for(Threshold_Name in c("Th0","Th1","Th2","Th2_5","Th5","Th10","Th20","localT2","localT2_HumanGEM")){
  RoL_gMCS <- readRDS(paste0("./RDSResults/DepMap/RoL_",Threshold_Name,".RDS"))
  
  
  # Make gMCS-gene relations          ##
  Results_Order_List_Comparisons <- vector("list",3)
  names(Results_Order_List_Comparisons) <- c("One_Gene_to_One_Task_NoEssentials",
                                             "One_Gene_to_Many_Tasks_NoEssentials",
                                             "TrueEssentials")
  List_gMCS_Essential_Matrix <- RoL_gMCS$list.gMCS.essential.mat
  
  gMCS_Length <- 8
  Binary_Matrix <- matrix(0, nrow = nrow(RoL_gMCS$mat.essential.gene),
                          ncol = ncol(RoL_gMCS$mat.essential.gene),
                          dimnames = list(rownames(RoL_gMCS$mat.essential.gene),
                                          colnames(RoL_gMCS$mat.essential.gene)))
  NotBinary_Matrix <- Binary_Matrix
  
  for(Essential_Gene in names(List_gMCS_Essential_Matrix)){
    Binary_Matrix[Essential_Gene,]    <- as.numeric(colSums(List_gMCS_Essential_Matrix[[Essential_Gene]]) > 0)
  }
  Essential_Gene_gMCS <-  do.call(rbind, lapply(List_gMCS_Essential_Matrix, colSums))
  # Extend data for all HumanGEM genes #
  Essential_Gene_gMCS_Extended <- matrix(data = 0, 
                                         nrow = length(unique(c(Table_HumanGEM_Genes$ENSEMBL,rownames(Essential_Gene_gMCS)))),
                                         ncol = ncol(Essential_Gene_gMCS),
                                         dimnames = list(unique(c(Table_HumanGEM_Genes$ENSEMBL,rownames(Essential_Gene_gMCS))),
                                                         colnames(Essential_Gene_gMCS)))
  # Essential genes are still Essential #
  for(EssentialGene in rownames(Binary_Matrix)){
    Essential_Gene_gMCS_Extended[rownames(Essential_Gene_gMCS_Extended) %in%
                                   EssentialGene,] <- Binary_Matrix[EssentialGene,]
  }
  
  # Take Order 1 Genes #
  GetOrder1 <- rownames(Essential_Gene_gMCS_Extended)
  ThisAreOrder1 <- unname(unlist(
    gMCS.info$EssentialTasks_CultureBiomass$gMCSs.ENSEMBL.list[which(lengths(gMCS.info$EssentialTasks_CultureBiomass$gMCSs.ENSEMBL.list)==1)]))
  Order1_Indexes <- which(GetOrder1 %in% ThisAreOrder1)
  Essential_Gene_gMCS_Essentials <- Essential_Gene_gMCS_Extended
  for (Essential_Gene in Order1_Indexes){
    Essential_Gene_gMCS_Essentials[Essential_Gene, which(Essential_Gene_gMCS_Essentials[Essential_Gene,] != 0)] <- NA
  }
  CCLE_Effect_ENSEMBL_NoOrder1 <-   CCLE_Effect_ENSEMBL[!rownames(CCLE_Effect_ENSEMBL) %in% ThisAreOrder1,]
  
  # Tasks to Genes                                 #
  Task_DataFrame <- data.frame(matrix(data = NA, 
                                      nrow = length(List_gMCS_Essential_Matrix),
                                      ncol = 2),
                               row.names = names(List_gMCS_Essential_Matrix))
  colnames(Task_DataFrame) <- c("Essential_Gene", "Associated_Tasks")
  for (Essential_Gene in names(List_gMCS_Essential_Matrix)){
    Essentiality_Index <-  which(rowSums(List_gMCS_Essential_Matrix[[Essential_Gene]]) != 0)
    if (length(Essentiality_Index) != 0 ){
      gMCS_Tasks <- names(gMCS.info$EssentialTasks_CultureBiomass$gMCSs.ENSEMBL.txt[Essentiality_Index])
      gMCS_Tasks_All <- unique(sapply(1:length(gMCS_Tasks), function(x) strsplit(gMCS_Tasks[x], "--")[[1]][1]))
      Task_DataFrame[Essential_Gene,] <- c(names(List_gMCS_Essential_Matrix[Essential_Gene]),
                                           length(gMCS_Tasks_All))
    }
  }
  One_Gene_One_Task  <- Task_DataFrame$Essential_Gene[which(Task_DataFrame$Associated_Tasks == 1)]
  One_Gene_One_Task  <- One_Gene_One_Task[!One_Gene_One_Task %in% ThisAreOrder1]
  One_Gene_Many_Task <- Task_DataFrame$Essential_Gene[which(Task_DataFrame$Associated_Tasks != 1)]
  length(One_Gene_One_Task);length(One_Gene_Many_Task);length(ThisAreOrder1)
  # Make 1 task to 1 gene             ##
  One_to_One_Mat_NoEs <-  Essential_Gene_gMCS_Essentials
  One_to_One_Mat_NoEs[rownames(One_to_One_Mat_NoEs) %in% One_Gene_One_Task]  <- 1  
  One_to_One_Mat_NoEs[rownames(One_to_One_Mat_NoEs) %in% One_Gene_Many_Task] <- NA 
  Results_Order_List_Comparisons[[1]] <- One_to_One_Mat_NoEs
  
  
  # Make n tasks to 1 gene            ##
  Task_Analysis_Mat_NoEs <-  Essential_Gene_gMCS_Essentials
  Task_Analysis_Mat_NoEs[rownames(Task_Analysis_Mat_NoEs) %in% One_Gene_One_Task] <- NA 
  Task_Analysis_Mat_NoEs[rownames(Task_Analysis_Mat_NoEs) %in% One_Gene_Many_Task]  <- 1
  Results_Order_List_Comparisons[[2]] <- Task_Analysis_Mat_NoEs
  
  # All Essential Genes               ##
  TrueEssentials <-  Essential_Gene_gMCS_Extended
  TrueEssentials[TrueEssentials != 0]   <- 1
  Results_Order_List_Comparisons[[3]] <- TrueEssentials
  
  # Compute Statistics                ##
  Comparison_Analysis_Results <- vector("list", length(Results_Order_List_Comparisons))
  k <- 0
  for (Analysis in names(Results_Order_List_Comparisons)){
    try({result_list <- list()
    k <- k+1
    evalMethod = "AllTasks"
    
    if (Analysis == "TrueEssentials"){
      Mat_binary <- CCLE_Effect_ENSEMBL < (-0.6)
      Mat_binary[which(is.na(Mat_binary))] <- 0
    } else {
      Mat_binary <- CCLE_Effect_ENSEMBL_NoOrder1 < (-0.6)
      Mat_binary[which(is.na(Mat_binary))] <- 0
    }
    
    tissues <- intersect(colnames(Results_Order_List_Comparisons[[Analysis]]),
                         colnames(Mat_binary))
    
    # Initialize Variables #
    results = data.frame('cellLine' = tissues,
                         'TP' = NA, 'TN' = NA, 'FP' = NA, 'FN' = NA,
                         "check_sum_genes" = NA, 
                         'accuracy' = NA, 'sensitivity' = NA,
                         'specificity' = NA, 'F1' = NA,
                         'MCC' = NA, 'Penr' = NA, 'logPenr' = NA,
                         'PenrAdj' = NA, 'logPenrAdj' = NA)
    
    for (t in 1:length(tissues)){
      
      # Genes From the Model #
      modelGenes <- unique(Table_HumanGEM_Genes$ENSEMBL) 
      # Essential Genes #
      modelEssential <- rownames(Results_Order_List_Comparisons[[Analysis]])[Results_Order_List_Comparisons[[Analysis]][,tissues[t]]>0]
      modelEssential <- unique(Table_HumanGEM_Genes$ENSEMBL[Table_HumanGEM_Genes$ENSEMBL %in% modelEssential])
      # Non-Essential Genes #
      modelNonEssential <- setdiff(modelGenes, modelEssential)
      
      # expGenes has all the genes that are not NA from Achilles #
      expGenes <-  rownames(Mat_binary)[!is.na(Mat_binary[,tissues[t]])]
      # expEssential has all the genes that are considered essential #
      expEssential <- expGenes[Mat_binary[expGenes,tissues[t]]==1]
      # expNonEssential has all the genes that are NOT essential # 
      expNonEssential <- setdiff(expGenes,expEssential)
      
      results$TP[t] = length(intersect(modelEssential, expEssential))        # True Positives  #
      results$TN[t] = length(intersect(modelNonEssential, expNonEssential))  # True Negatives  #
      results$FP[t] = length(intersect(modelEssential, expNonEssential))     # False Positives #
      results$FN[t] = length(intersect(modelNonEssential, expEssential))     # False Negatives #   
      
      # Intersect Genes from HumanGEM and Genes from ACHILLES #
      pop <- intersect(modelGenes, expGenes)
      # Essential Genes from gMCST5 #
      sample <- intersect(modelEssential, expGenes)
      # Essential Genes from ACHILLES #
      successes <- intersect(expEssential, modelGenes)
      
      results$Penr[t] <- phyper(length(intersect(successes,sample)) - 1, # Number of White Balls Drawn  #
                                length(intersect(successes,pop)),        # Number of White Balls in Urn #
                                length(setdiff(pop, successes)),         # Number of Black Balls in Urn #
                                length(sample),                          # Total Number of Balls Drawn  #  
                                lower.tail = FALSE);
      
    }
    
    
    results$check_sum_genes =  results$TP + results$TN + results$FP + results$FN
    
    # Calculate Metrics
    results$sensitivity = results$TP/(results$TP + results$FN)
    results$specificity = results$TN/(results$TN + results$FP)
    results$accuracy    = (results$TP + results$TN)/(results$TP + results$TN + results$FP + results$FN)
    results$PPV         = results$TP/(results$TP+results$FP)
    results$F1          = 2*results$TP/(2*results$TP + results$FP + results$FN)
    results$MCC         = ((results$TP*results$TN) - (results$FP*results$FN))/(sqrt((results$TP+results$FP)*(results$TP+results$FN))*sqrt((results$TN+results$FP)*(results$TN+results$FN)))  # Matthews correlation coefficient
    results$PenrAdj     = p.adjust(results$Penr,'BH')
    
    results$logPenr     = -log10(results$Penr)
    results$logPenrAdj  = -log10(results$PenrAdj)
    
    results$model <- "CCLE"
    results$evalMethod <- evalMethod
    
    results$calcMethod <- Analysis
    
    results$PenrAdj = p.adjust(results$Penr,'BH')
    
    results$logPenr = -log10(results$Penr)
    results$logPenrAdj = -log10(results$PenrAdj)
    
    result_list[[evalMethod]] <- results
    
    results <- do.call(rbind, result_list)})
    
    Comparison_Analysis_Results[[k]] <- results
    
  }
  
  gMCS_Statistics <- do.call(rbind, Comparison_Analysis_Results)
  saveRDS(gMCS_Statistics, paste0("./RDSResults/DepMap/",Today,"_Task_Promiscuity_Results_gMCS",ThresholdName,".RDS"))
}
############################### ####
#            Biomass            ####
#### Compute gmcTH   #### ####
ThresholdVector <- c(0,1,2,2.5,5,10,20)
Sample_Class_Dummy <- data.frame("Sample_Class" = "Cells")
levels(Sample_Class_Dummy)  <- "Cells"

if("AllThresholdsResults.RDS" %in% dir("./Biomass/RDSResults/DepMap")){
  ResultList <- readRDS("./Biomass/RDSResults/DepMap/AllThresholdsResults.RDS")
} else {
  if("RoLAllThresholds.RDS" %in% dir("./Biomass/RDSResults/DepMap/")){
    RoLResults <-  readRDS("./Biomass/RDSResults/DepMap/RoLAllThresholds.RDS")
  } else {
    RoLResults <- vector("list", length(ThresholdVector))
    names(RoLResults) <- paste0("Th",ThresholdVector)
    for(Threshold in ThresholdVector){
      if (paste0("RoL_Th",Threshold,".RDS") %in% dir("./Biomass/RDSResults/DepMap")){
        RoLResults[[paste0("Th",Threshold)]] <- readRDS(paste0("./Biomass/RDSResults/DepMap/RoL_Th",Threshold,".RDS"))
      } else {
        RoLResults[[paste0("Th",Threshold)]] <- CalculateEssentialGenes_gmcsTH_DOM(gene.exp = CCLE_Exp,
                                                                                   gMCS_Analysis = gMCS.info$Only_CultureBiomass,
                                                                                   sample.class = Sample_Class_Dummy,
                                                                                   gmcsTH_perc = Threshold/100,
                                                                                   nWorkers = 8,
                                                                                   gMCS_order = 8,
                                                                                   order1 = TRUE)
        saveRDS(RoLResults[[paste0("Th",Threshold)]], file = paste0("./Biomass/RDSResults/DepMap/RoL_Th",Threshold,".RDS"))
      }
    }
    saveRDS(RoLResults, "./Biomass/RDSResults/DepMap/RoLAllThresholds.RDS")
  }
  
  ResultList <- vector("list", length(ThresholdVector))
  names(ResultList) <- paste0("Th",ThresholdVector)
  
  for(Threshold in ThresholdVector){
    Binary_Matrix_Length_List <- ComputeBinaryMatrix(RoLResults[[paste0("Th",Threshold)]],
                                                     gMCS.info$Only_CultureBiomass$gMCSs.ENSEMBL.list)
    ResultsInList <- lapply(1:8, function(x) TweakedDepMapgMCS_Results(Binary_Matrix = Binary_Matrix_Length_List[[x]],
                                                                       Effect = CCLE_Effect,
                                                                       Table_HumanGEM_Genes = Table_HumanGEM_Genes,
                                                                       gMCS_order =  x,
                                                                       Threshold_Info =  paste0("gMCS_T",Threshold)))
    ResultsInList <- do.call(rbind, ResultsInList)
    ResultsInList$Order <- factor(ResultsInList$Order,
                                  levels=c(paste0("Order_",1:7), "HigherThan_Order_7"))
    ResultsInList <-  ResultsInList %>% rename(TP = "True Positives", FP = "False Positives",
                                               TN = "True Negatives", FN = "False Negatives",
                                               PPV = "Positive Predicted Value",
                                               sensitivity = "Sensitivity", specificity = "Specificity",
                                               accuracy = "Accuracy", MCC = "Matthew's Cor. Coef.")
    ResultsInList %>% group_by(Order) %>%
      summarise_at(c("True Positives", "False Positives", 
                     "False Negatives", "True Negatives",
                     "Accuracy", "Sensitivity", "Specificity",
                     "Positive Predicted Value", "Matthew's Cor. Coef."), median, na.rm = TRUE) %>%
      write.csv(file = paste0(".Results/DepMap/gMCSLengths/BIOLengthMedians_Th",Threshold,".csv"),
                row.names = F)
    ResultList[[paste0("Th",Threshold)]] <- ResultsInList
  }
  saveRDS(ResultList, file = "./Biomass/RDSResults/DepMap/AllThresholdsResults.RDS")
}
####  LocalT2 - 1244 #### ####
source("./CustomFunctions/DepMap/fun_DepMapLocalT2CalculateEssential.R")
if ("Results_localT2.RDS" %in% dir("./Biomass/RDSResults/DepMap/") && "NotCumulativeResults_localT2.RDS" %in% dir("./Biomass/RDSResults/DepMap/")){
  Results_localT2 <- readRDS("./RDSResults/DepMap/Results_localT2.RDS")
  NotCumulativeResults_localT2 <- readRDS("./Biomass/RDSResults/DepMap/NotCumulativeResults_localT2.RDS")
} else {
  if ("RoL_gMCS_localT2.RDS" %in% dir("./Biomass/RDSResults/DepMap/")){
    RoL_gMCS_localT2 <- readRDS("./Biomass/RDSResults/DepMap/RoL_gMCS_localT2.RDS")
    
  } else {
    Sample_Class_Dummy <- data.frame("Sample_Class" = colnames(CCLE_Exp))
    levels(Sample_Class_Dummy)  <- colnames(CCLE_Exp)
    Sample_Cohort_Dummy <- data.frame("Sample_Cohort" = "DepMap")
    levels(Sample_Cohort_Dummy) <- "DepMap"
    CCLE_Exp <- matrix(as.numeric(CCLE_Exp), nrow = nrow(CCLE_Exp), ncol = ncol(CCLE_Exp),
                       dimnames = list(rownames(CCLE_Exp), colnames(CCLE_Exp)))
    
    RoL_gMCS_localT2 <- LocalT2CalculateEssential(gene.exp      = CCLE_Exp,
                                                  gMCS.info     = gMCS.info$Only_CultureBiomass,
                                                  sample.class  = Sample_Class_Dummy,
                                                  sample.cohort = Sample_Cohort_Dummy,
                                                  localT2_mode  = "all_genes_gMCSs",
                                                  nWorkers      = 8,
                                                  gMCS_order    = 8)
    
    
    saveRDS(object = RoL_gMCS_localT2, file = "./Biomass/RDSResults/DepMap/RoL_gMCS_localT2.RDS")
  }
  # Cumulative     ####
  BinaryMatrixLengthList <- ComputeBinaryMatrix(RoL_gMCS_localT2,
                                                gMCS.info$Only_CultureBiomass$gMCSs.ENSEMBL.list)
  Results_localT2 <- lapply(1:8, function(x) TweakedDepMapgMCS_Results(Binary_Matrix = BinaryMatrixLengthList[[x]],
                                                                       Effect = CCLE_Effect,
                                                                       Table_HumanGEM_Genes = Table_HumanGEM_Genes,
                                                                       gMCS_order =  x,
                                                                       Threshold_Info =  "localT2"))
  Results_localT2 <- do.call(rbind, Results_localT2)
  Results_localT2$Order <- factor(Results_localT2$Order,levels=c(paste0("Order_",1:7), "HigherThan_Order_7"))
  Results_localT2 <- Results_localT2 %>% dplyr::rename("True Positives"  = TP,
                                                       "False Positives" = FP,
                                                       "True Negatives"  = TN,
                                                       "False Negatives" = FN,
                                                       "Sensitivity" = sensitivity,
                                                       "Specificity" = specificity,
                                                       "Accuracy"    = accuracy,
                                                       "Positive Predicted Value" = PPV,
                                                       "Matthew's Cor. Coef."     = MCC)
  Results_localT2 %>% filter(thMethod == "localT2") %>% group_by(Order) %>%
    summarise_at(c("True Positives", "False Positives", 
                   "False Negatives", "True Negatives",
                   "Accuracy", "Sensitivity", "Specificity",
                   "Positive Predicted Value", "Matthew's Cor. Coef."), median, na.rm = TRUE) %>%
    write.csv(file = "./BiomassResults/DepMap/localT2/Cumulative_LocalT2.csv", 
              row.names = F)
  
  MeltedLocalT2 <- melt(data = Results_localT2,
                        id.vars = "Order",
                        measure.vars = c("True Positives", "False Positives",
                                         "False Negatives", "True Negatives",
                                         "Accuracy", "Sensitivity", "Specificity",
                                         "Positive Predicted Value", "Matthew's Cor. Coef."))
  LocalT2_Length <- ggplot(MeltedLocalT2, aes(x = Order, y = value, fill = Order)) + 
    facet_wrap(~variable,scales = "free") + geom_boxplot() + theme_bw() +  
    # geom_boxplot(width=0.1, fill = NA) +
    labs(title = "LocalT2_1244",
         x = "", y = "",
         fill = "Lengths") + 
    scale_fill_discrete(name = "Maximum gMCS\nLengths (#gMCS)",
                        labels = c(paste0("1 (n = ",NumberOfgMCS[1],")"), paste0("2 (n = ",NumberOfgMCS[2],")"),
                                   paste0("3 (n = ",NumberOfgMCS[3],")"), paste0("4 (n = ",NumberOfgMCS[4],")"),
                                   paste0("5 (n = ",NumberOfgMCS[5],")"), paste0("6 (n = ",NumberOfgMCS[6],")"),
                                   paste0("7 (n = ",NumberOfgMCS[7],")"),
                                   paste0("All (n = ",nrow(gMCS.info$Only_CultureBiomass$gMCSs.ENSEMBL.mat),")")),
                        type = RColorBrewer::brewer.pal(8, "Paired")) +
    theme(plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
          axis.text.x = element_blank(), axis.ticks.x = element_blank(),
          strip.text.x = element_text(size = 18, face = "bold"),
          legend.position = "bottom")
  
  # LocalT2_Length
  
  ggsave(paste0("./BiomassResults/DepMap/localT2/CumulativeLocalT2.png"),
         plot = LocalT2_Length, device = "png", width = 16, height = 12, units = "in", dpi = 300)
  
  saveRDS(Results_localT2, "./Biomass/RDSResults/DepMap/Results_localT2.RDS")
  
  
  # Not Cumulative ####
  NotCumulativeBinaryMatrixLengthList <- ComputeBinaryMatrix(RoL_gMCS_localT2,
                                                             gMCS.info$Only_CultureBiomass$gMCSs.ENSEMBL.list,
                                                             FALSE)
  NotCumulativeResults_localT2 <- lapply(1:8, function(x) TweakedDepMapgMCS_Results(Binary_Matrix = NotCumulativeBinaryMatrixLengthList[[x]],
                                                                                    Effect = CCLE_Effect,
                                                                                    Table_HumanGEM_Genes = Table_HumanGEM_Genes,
                                                                                    gMCS_order =  x,
                                                                                    Threshold_Info =  "localT2"))
  NotCumulativeResults_localT2 <- do.call(rbind, NotCumulativeResults_localT2)
  NotCumulativeResults_localT2$Order <- factor(NotCumulativeResults_localT2$Order,
                                               levels=c(paste0("Order_",1:7), "HigherThan_Order_7"))
  NotCumulativeResults_localT2 <- NotCumulativeResults_localT2 %>% dplyr::rename("True Positives"  = TP,
                                                                                 "False Positives" = FP,
                                                                                 "True Negatives"  = TN,
                                                                                 "False Negatives" = FN,
                                                                                 "Sensitivity" = sensitivity,
                                                                                 "Specificity" = specificity,
                                                                                 "Accuracy"    = accuracy,
                                                                                 "Positive Predicted Value" = PPV,
                                                                                 "Matthew's Cor. Coef."     = MCC)
  NotCumulativeResults_localT2 %>% filter(thMethod == "localT2") %>% group_by(Order) %>%
    summarise_at(c("True Positives", "False Positives", 
                   "False Negatives", "True Negatives",
                   "Accuracy", "Sensitivity", "Specificity",
                   "Positive Predicted Value", "Matthew's Cor. Coef."), median, na.rm = TRUE) %>%
    write.csv(file = "./BiomassResults/DepMap/localT2/NotCumulative_LocalT2.csv", 
              row.names = F)
  
  NCMeltedLocalT2 <- melt(data = NotCumulativeResults_localT2,
                          id.vars = "Order",
                          measure.vars = c("True Positives", "False Positives",
                                           "False Negatives", "True Negatives",
                                           "Accuracy", "Sensitivity", "Specificity",
                                           "Positive Predicted Value", "Matthew's Cor. Coef."))
  NotCumLocalT2_Length <- ggplot(NCMeltedLocalT2, aes(x = Order, y = value, fill = Order)) + 
    facet_wrap(~variable,scales = "free") + geom_boxplot() + theme_bw() +  
    # geom_boxplot(width=0.1, fill = NA) +
    labs(title = "LocalT2_1244",
         x = "", y = "",
         fill = "Lengths") + 
    scale_fill_discrete(name = "Maximum gMCS\nLengths (#gMCS)",
                        labels = c(paste0("1 (n = ",NumberOfgMCS[1],")"), paste0("2 (n = ",NumberOfgMCS[2],")"),
                                   paste0("3 (n = ",NumberOfgMCS[3],")"), paste0("4 (n = ",NumberOfgMCS[4],")"),
                                   paste0("5 (n = ",NumberOfgMCS[5],")"), paste0("6 (n = ",NumberOfgMCS[6],")"),
                                   paste0("7 (n = ",NumberOfgMCS[7],")"),
                                   paste0("All (n = ",nrow(gMCS.info$Only_CultureBiomass$gMCSs.ENSEMBL.mat),")")),
                        type = RColorBrewer::brewer.pal(8, "Paired")) +
    theme(plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
          axis.text.x = element_blank(), axis.ticks.x = element_blank(),
          strip.text.x = element_text(size = 18, face = "bold"),
          legend.position = "bottom")
  
  # NotCumLocalT2_Length
  
  ggsave(paste0("./BiomassResults/DepMap/localT2/NotCumulativeLocalT2.png"),
         plot = NotCumLocalT2_Length, device = "png", width = 16, height = 12, units = "in", dpi = 300)
  
  saveRDS(NotCumulativeResults_localT2, "./Biomass/RDSResults/DepMap/NotCumulativeResults_localT2.RDS")
} 
############################### ####
#########   Figures   ######### ####
############################### ####
# Figure 1B - DepMap ####
toAdd <- do.call(rbind,readRDS("./RDSResults/DepMap/AllThresholdsResults.RDS"))
toKeep <- c("gMCS_T0","gMCS_T1","gMCS_T2","gMCS_T2.5","gMCS_T3.5","gMCS_T5","gMCS_T10","gMCS_T20")
toRemove <- setdiff(unique(toAdd$thMethod), toKeep)

for(element in toRemove){toAdd <- toAdd[which(toAdd$thMethod != element),]}
toAdd$thMethod <- factor(toAdd$thMethod, levels=toKeep)
ThresholdColors <- viridis::magma(8)

tmp     <- rbind(toAdd,
                 readRDS("./RDSResults/DepMap/Results_localT2.RDS"),
                 readRDS("./RDSResults/DepMap/Results_localT2HumanGEM.RDS"))

tmp$thMethod <- factor(tmp$thMethod, levels = c("gMCS_T0", "gMCS_T1", "gMCS_T2",
                                                "gMCS_T2.5", "gMCS_T5", "gMCS_T10",
                                                "gMCS_T20", "NA",
                                                "localT2", "localT2HumanGEM"))
tmp <- tmp %>% filter(Order == "HigherThan_Order_7")
tmp <- reshape2::melt(tmp, id.vars = c("thMethod"),
                      measure.vars =  c("True Positives", "False Positives", "Positive Predicted Value"))

## Prepare TP Plot ##
TruePositive_tmp          <- tmp %>% filter(variable == "True Positives")

# Compute pvalues #
if(FALSE){
  TruePositive_tmp %>% kruskal.test(value~thMethod, data=.)
  # wilcox   ####
  pValThDF <- as.data.frame(matrix(1,nrow = length(unique(TruePositive_tmp$thMethod)),
                                   ncol = length(unique(TruePositive_tmp$thMethod))))
  colnames(pValThDF) = rownames(pValThDF) = unique(TruePositive_tmp$thMethod)
  
  for(i in unique(TruePositive_tmp$thMethod)){
    for(j in setdiff(unique(TruePositive_tmp$thMethod),i)){
      pValThDF[i,j] <- round(wilcox.test(TruePositive_tmp$value[TruePositive_tmp$thMethod == i],
                                         TruePositive_tmp$value[TruePositive_tmp$thMethod == j],
                                         paired = TRUE)$p.value, 5)
      
    }
  }
  write.csv(x = pValThDF, file = paste0(".Results/DepMap/Threshold/",Today,"_Threshold_TP_wilcox_DepMap.csv"),
            quote = FALSE, row.names = TRUE)
  # Continue ####
  for(i in unique(TruePositive_tmp$thMethod)){
    shtest <- shapiro.test(TruePositive_tmp$value[TruePositive_tmp$thMethod == i])
    varthr <- var(TruePositive_tmp$value[TruePositive_tmp$thMethod == i])
    if(i != "localT2" & i != "localT2HumanGEM"){
      cat(paste0("Threshold: ", sapply(strsplit(i,"_T"), "[[", 2), "\t Shapiro-Test: ", shtest$p.value, "\t Variance: ", varthr, "\n"))
    } else {
      cat(paste0("Threshold: ", i, "\t Shapiro-Test: ", shtest$p.value, "\t Variance: ", varthr, "\n"))
    }
  }
  
}
# Continue #
TruePositive_tmp_mean <- TruePositive_tmp %>% 
  group_by(thMethod) %>% 
  summarize(Mean = mean(value)) %>%
  ungroup()

TruePositive_tmp <- rbind(TruePositive_tmp, list("NA","True Positives",NA))
change_lab <- "True Positives"; names(change_lab) <- "True Positives"
TPlot <- ggplot(TruePositive_tmp, aes(x = thMethod, y = value, fill = thMethod)) + 
  geom_boxplot(alpha = 0.6) + 
  facet_wrap("variable",labeller = labeller(variable = change_lab))  +
  scale_fill_manual(values = c("gMCS_T0"   = ThresholdColors[1], "gMCS_T1"   = ThresholdColors[2],
                               "gMCS_T2"   = ThresholdColors[3], "gMCS_T2.5" = ThresholdColors[4],
                               "gMCS_T5"   = ThresholdColors[5],"gMCS_T10"  = ThresholdColors[6],
                               "gMCS_T20"  = ThresholdColors[7],
                               "NA" = "#000000",
                               "localT2" = "#B2DF8A", "localT2HumanGEM" = "#6ba35b")) +
  theme_bw() + xlab("") + ylab("") + 
  theme(legend.position = "none") + 
  geom_point(data    = TruePositive_tmp_mean,
             mapping = aes(x = thMethod, y = Mean), size  = 2.5,
             color   = "dark red") + 
  geom_line(data = TruePositive_tmp_mean[grepl("gMCS", TruePositive_tmp_mean$thMethod),], 
            mapping = aes(x = thMethod, y = Mean, group = 1), size = 0.9) +
  geom_line(data = TruePositive_tmp_mean[!grepl("gMCS", TruePositive_tmp_mean$thMethod),], 
            mapping = aes(x = thMethod, y = Mean, group = 1), size = 0.9) +
  geom_vline(xintercept = 8, size = 2, color = "dark red") +
  scale_x_discrete(labels=c("gMCS_T0"   = "0%",
                            "gMCS_T1"   = "1%",
                            "gMCS_T2"   = "2%",
                            "gMCS_T2.5" = "2.5%",
                            "gMCS_T5"   = "5%",
                            "gMCS_T10"  = "10%",
                            "gMCS_T20"  = "20%",
                            "NA" = "",
                            "localT2"   = "LocalT2",
                            "localT2HumanGEM"  = "LocalT2-H1"),
                   limits = c("gMCS_T0","gMCS_T1","gMCS_T2","gMCS_T2.5",
                              "gMCS_T5","gMCS_T10","gMCS_T20","gMCS_T20",
                              "localT2","localT2HumanGEM")) +
  scale_y_continuous(limits = c(60,180), breaks = c(seq(from = 60, to = 180, by = 20))) +
  theme(strip.text = element_text(face  = "bold", size = 28),
        strip.background = element_rect(fill="#b09c7b"),
        axis.text.x = element_text(angle = 30, hjust = 1, vjust = 1, size = 20),
        axis.text.y = element_text(size = 23)) 
TPlot
## Prepare FP Plot ##
FalsePositive_tmp          <- tmp %>% filter(variable == "False Positives")
# Compute pvalues #
if(FALSE){
  FalsePositive_tmp %>% kruskal.test(value~thMethod, data=.)
  # wilcox ####
  pValThDF <- as.data.frame(matrix(1,nrow = length(unique(FalsePositive_tmp$thMethod)),
                                   ncol = length(unique(FalsePositive_tmp$thMethod))))
  colnames(pValThDF) = rownames(pValThDF) = unique(FalsePositive_tmp$thMethod)
  
  for(i in unique(FalsePositive_tmp$thMethod)){
    for(j in setdiff(unique(FalsePositive_tmp$thMethod),i)){
      pValThDF[i,j] <- round(wilcox.test(FalsePositive_tmp$value[FalsePositive_tmp$thMethod == i],
                                         FalsePositive_tmp$value[FalsePositive_tmp$thMethod == j],
                                         paired = TRUE)$p.value, 5)
      
    }
  }
  write.csv(x = pValThDF, file = paste0(".Results/DepMap/Threshold/",Today,"_Threshold_FP_wilcox_DepMap.csv"),
            quote = FALSE, row.names = TRUE)
  # Continue ####
  for(i in unique(FalsePositive_tmp$thMethod)){
    shtest <- shapiro.test(FalsePositive_tmp$value[FalsePositive_tmp$thMethod == i])
    varthr <- var(FalsePositive_tmp$value[FalsePositive_tmp$thMethod == i])
    if(i != "localT2" & i != "localT2HumanGEM"){
      cat(paste0("Threshold: ", sapply(strsplit(i,"_T"), "[[", 2), "\t Shapiro-Test: ", shtest$p.value, "\t Variance: ", varthr, "\n"))
    } else {
      cat(paste0("Threshold: ", i, "\t Shapiro-Test: ", shtest$p.value, "\t Variance: ", varthr, "\n"))
    }
  }
  
}
# Continue #
FalsePositive_tmp_mean <- FalsePositive_tmp %>% 
  group_by(thMethod) %>% 
  summarize(Mean = mean(value)) %>%
  ungroup()

FalsePositive_tmp <- rbind(FalsePositive_tmp, list("NA","False Positives",NA))
change_lab <- "False Positives"; names(change_lab) <- "False Positives"
FPlot <- ggplot(FalsePositive_tmp, aes(x = thMethod, y = value, fill = thMethod)) + 
  geom_boxplot(alpha = 0.6) +  
  facet_wrap("variable",labeller = labeller(variable = change_lab))  +
  scale_fill_manual(values = c("gMCS_T0"   = ThresholdColors[1], "gMCS_T1"   = ThresholdColors[2],
                               "gMCS_T2"   = ThresholdColors[3], "gMCS_T2.5" = ThresholdColors[4],
                               "gMCS_T5"   = ThresholdColors[5],"gMCS_T10"  = ThresholdColors[6],
                               "gMCS_T20"  = ThresholdColors[7],
                               "NA" = "#000000",
                               "localT2" = "#B2DF8A", "localT2HumanGEM" = "#6ba35b")) +
  theme_bw() + xlab("") + ylab("") +
  theme(legend.position = "none") + 
  geom_point(data    = FalsePositive_tmp_mean,
             mapping = aes(x = thMethod, y = Mean),size  = 2.5,
             color   = "dark red") + 
  geom_line(data = FalsePositive_tmp_mean[grepl("gMCS", FalsePositive_tmp_mean$thMethod),], 
            mapping = aes(x = thMethod, y = Mean, group = 1), size = 0.9) +
  geom_line(data = FalsePositive_tmp_mean[!grepl("gMCS", FalsePositive_tmp_mean$thMethod),], 
            mapping = aes(x = thMethod, y = Mean, group = 1), size = 0.9) +
  geom_vline(xintercept = 8, size = 2, color = "dark red") +
  scale_x_discrete(labels=c("gMCS_T0"   = "0%",
                            "gMCS_T1"   = "1%",
                            "gMCS_T2"   = "2%",
                            "gMCS_T2.5" = "2.5%",
                            "gMCS_T5"   = "5%",
                            "gMCS_T10"  = "10%",
                            "gMCS_T20"  = "20%",
                            "NA" = "",
                            "localT2"   = "LocalT2",
                            "localT2HumanGEM"  = "LocalT2-H1"),
                   limits = c("gMCS_T0","gMCS_T1","gMCS_T2","gMCS_T2.5",
                              "gMCS_T5","gMCS_T10","gMCS_T20","gMCS_T20",
                              "localT2","localT2HumanGEM")) +
  scale_y_continuous(limits = c(100,275), breaks = c(seq(from = 100, to = 275, by = 25))) +
  theme(strip.text = element_text(face  = "bold", size = 28),
        strip.background = element_rect(fill="#b09c7b"),
        axis.text.x = element_text(angle = 30, hjust = 1, vjust = 1, size = 20),
        axis.text.y = element_text(size = 23)) 
FPlot
## Prepare PPV Plot ##
PositivePredictedValue_tmp          <- tmp %>% filter(variable == "Positive Predicted Value")
# Compute pvalues #
if(FALSE){
  PositivePredictiveValue_tmp %>% kruskal.test(value~thMethod, data=.)
  # Wilcox ####
  pValThDF <- as.data.frame(matrix(1,nrow = length(unique(PositivePredictiveValue_tmp$thMethod)),
                                   ncol = length(unique(PositivePredictiveValue_tmp$thMethod))))
  colnames(pValThDF) = rownames(pValThDF) = unique(PositivePredictiveValue_tmp$thMethod)
  
  for(i in unique(PositivePredictiveValue_tmp$thMethod)){
    for(j in setdiff(unique(PositivePredictiveValue_tmp$thMethod),i)){
      pValThDF[i,j] <- round(wilcox.test(PositivePredictiveValue_tmp$value[PositivePredictiveValue_tmp$thMethod == i],
                                         PositivePredictiveValue_tmp$value[PositivePredictiveValue_tmp$thMethod == j],
                                         paired = TRUE)$p.value, 5)
      
    }
  }
  write.csv(x = pValThDF, file = paste0(".Results/DepMap/Threshold/",Today,"_Threshold_PPV_wilcoxon_DepMap.csv"),
            quote = FALSE, row.names = TRUE)
  # Continue #
  for(i in unique(PositivePredictiveValue_tmp$thMethod)){
    shtest <- shapiro.test(PositivePredictiveValue_tmp$value[PositivePredictiveValue_tmp$thMethod == i])
    varthr <- var(PositivePredictiveValue_tmp$value[PositivePredictiveValue_tmp$thMethod == i])
    if(i != "localT2" & i != "localT2HumanGEM"){
      cat(paste0("Threshold: ", sapply(strsplit(i,"_T"), "[[", 2), "\t Shapiro-Test: ", shtest$p.value, "\t Variance: ", varthr, "\n"))
    } else {
      cat(paste0("Threshold: ", i, "\t Shapiro-Test: ", shtest$p.value, "\t Variance: ", varthr, "\n"))
    }
  }
  
}
# Continue #
PositivePredictedValue_tmp_mean <- PositivePredictedValue_tmp %>% 
  group_by(thMethod) %>% 
  summarize(Mean = mean(value)) %>%
  ungroup()

PositivePredictedValue_tmp <- rbind(PositivePredictedValue_tmp, list("NA","Positive Predicted Value",NA))
change_lab <- "Positive Predicted Value"; names(change_lab) <- "Positive Predicted Value"
PPVlot <- ggplot(PositivePredictedValue_tmp, aes(x = thMethod, y = value, fill = thMethod)) + 
  geom_boxplot(alpha = 0.6) + 
  facet_wrap("variable",labeller = labeller(variable = change_lab))  +
  scale_fill_manual(values = c("gMCS_T0"   = ThresholdColors[1], "gMCS_T1"   = ThresholdColors[2],
                               "gMCS_T2"   = ThresholdColors[3], "gMCS_T2.5" = ThresholdColors[4],
                               "gMCS_T5"   = ThresholdColors[5],"gMCS_T10"  = ThresholdColors[6],
                               "gMCS_T20"  = ThresholdColors[7],
                               "NA" = "#000000",
                               "localT2" = "#B2DF8A", "localT2HumanGEM" = "#6ba35b")) +
  theme_bw() + xlab("") + ylab("") +
  theme(legend.position = "none") + 
  geom_point(data    = PositivePredictedValue_tmp_mean,
             mapping = aes(x = thMethod, y = Mean), size = 2.5,
             color   = "dark red") + 
  geom_line(data = PositivePredictedValue_tmp_mean[grepl("gMCS", PositivePredictedValue_tmp_mean$thMethod),], 
            mapping = aes(x = thMethod, y = Mean, group = 1), size = 0.9) +
  geom_line(data = PositivePredictedValue_tmp_mean[!grepl("gMCS", PositivePredictedValue_tmp_mean$thMethod),], 
            mapping = aes(x = thMethod, y = Mean, group = 1), size = 0.9) +
  geom_vline(xintercept = 8, size = 2, color = "dark red") +
  scale_x_discrete(labels=c("gMCS_T0"   = "0%",
                            "gMCS_T1"   = "1%",
                            "gMCS_T2"   = "2%",
                            "gMCS_T2.5" = "2.5%",
                            "gMCS_T5"   = "5%",
                            "gMCS_T10"  = "10%",
                            "gMCS_T20"  = "20%",
                            "NA" = "",
                            "localT2"   = "LocalT2",
                            "localT2HumanGEM"  = "LocalT2-H1"),
                   limits = c("gMCS_T0","gMCS_T1","gMCS_T2","gMCS_T2.5",
                              "gMCS_T5","gMCS_T10","gMCS_T20","NA",
                              "localT2","localT2HumanGEM")) +
  scale_y_continuous(limits = c(0.2,0.65), breaks = c(seq(from = 0.2, to = 0.65, by = 0.1))) +
  theme(strip.text  = element_text(face  = "bold", size = 28),
        strip.background = element_rect(fill="#b09c7b"),
        axis.text.x = element_text(angle = 30, hjust = 1, vjust = 1, size = 20),
        axis.text.y = element_text(size = 23)) 
PPVlot
# Add All Three Together #
plotlistThreshold <-  list(TPlot, FPlot, PPVlot)
Boxplotting_Thresholds <- ggarrange(plotlist = plotlistThreshold,
                                    ncol = 3, nrow = 1)

# Boxplotting_Thresholds <- annotate_figure(p = Boxplotting_Thresholds,
#                                           top = text_grob("Threshold Analysis", face = "bold", size = 18))
Boxplotting_Thresholds
ggsave(filename = paste0(".Results/DepMap/Threshold/",Today,"_DepMap_Threshold.png"),
       plot = Boxplotting_Thresholds, device = "png", width = 21, height = 8,units = "in",dpi = 300)


# Figure 1D - DepMap ####
ThresholdColors <- viridis::magma(8)
Results <- do.call(rbind, readRDS("./RDSResults/DepMap/AllThresholdsResults.RDS"))
index <-  Results$thMethod == "gMCS_T0" | Results$thMethod == "gMCS_T1" |Results$thMethod == "gMCS_T2" | Results$thMethod == "gMCS_T2.5" | Results$thMethod == "gMCS_T5" | Results$thMethod == "gMCS_T10" | Results$thMethod == "gMCS_T20"
Results_gMCST  <- Results[index,]
tmp     <- rbind(Results_gMCST,
                 readRDS("./RDSResults/DepMap/Results_localT2.RDS"),
                 readRDS("./RDSResults/DepMap/Results_localT2HumanGEM.RDS"))
tmp <- reshape2::melt(tmp, id.vars = c("Order", "thMethod"),
                      measure.vars =  c("True Positives", "False Positives", "Positive Predicted Value"))

# Prepare TP Plot #
TruePositive_tmp          <- tmp %>% filter(variable == "True Positives")
TruePositive_tmp_mean <- TruePositive_tmp %>% 
  group_by(Order, thMethod) %>% 
  summarize(Mean = mean(value)) %>%
  ungroup()
if(is.na(sum(FalsePositive_tmp_mean$Mean))){
  FalsePositive_tmp_mean$Mean[which(is.na(FalsePositive_tmp_mean$Mean))] <-  0
}
TruePositive_tmp_mean$variable <-  "True_Positives"
change_lab <- "True Positives"; names(change_lab) <- "True_Positives"
tmptmp <- TruePositive_tmp_mean %>% filter(thMethod == "gMCS_T5")
TPlot <- ggplot(TruePositive_tmp_mean, aes(x = Order, y = Mean, color = thMethod, group = thMethod)) + 
  geom_point(size = 4.5, alpha = 0.75) + geom_line(size = 1.5, alpha = 0.75) + 
  facet_wrap("variable",labeller = labeller(variable = change_lab))  +
  scale_color_manual(values = c("gMCS_T0"   = ThresholdColors[1],"gMCS_T1"   = ThresholdColors[2],
                                "gMCS_T2"   = ThresholdColors[3],"gMCS_T2.5" = ThresholdColors[4],
                                "gMCS_T5"   = ThresholdColors[5],"gMCS_T10"  = ThresholdColors[6],
                                "gMCS_T20"  = ThresholdColors[7],
                                "localT2" = "#B2DF8A", "localT2HumanGEM" = "#6ba35b")) +
  theme_bw() + xlab("") + ylab("") + 
  theme(legend.position = "bottom") + 
  scale_x_discrete(labels=c("Order_1" = "1",
                            "Order_2" = "2",
                            "Order_3" = "3",
                            "Order_4" = "4",
                            "Order_5" = "5",
                            "Order_6" = "6",
                            "Order_7" = "7",
                            "HigherThan_Order_7" = ">7")) +
  scale_y_continuous(limits = c(99,119), breaks = c(seq(from = 95, to = 120, by = 2))) +
  theme(strip.text = element_text(face  = "bold", size = 28), 
        strip.background = element_rect(fill="#b09c7b"), axis.text.x = element_text(size = 23),
        axis.text.y = element_text(size = 23))
TPlot
# Prepare FP Plot #
FalsePositive_tmp          <- tmp %>% filter(variable == "False Positives")

FalsePositive_tmp_mean <- FalsePositive_tmp %>% 
  group_by(Order, thMethod) %>% 
  summarize(Mean = mean(value)) %>%
  ungroup()

if(is.na(sum(FalsePositive_tmp_mean$Mean))){
  FalsePositive_tmp_mean$Mean[which(is.na(FalsePositive_tmp_mean$Mean))] <-  0
}
FalsePositive_tmp_mean$variable <-  "False_Positives"
change_lab <- "False Positives"; names(change_lab) <- "False_Positives"
FPlot <- ggplot(FalsePositive_tmp_mean, aes(x = Order, y = Mean, color = thMethod, group = thMethod)) + 
  geom_point(size = 4.5, alpha = 0.75) + geom_line(size = 1.5, alpha = 0.75) + 
  facet_wrap("variable",labeller = labeller(variable = change_lab))  +
  scale_color_manual(values = c("gMCS_T0"   = ThresholdColors[1],"gMCS_T1"   = ThresholdColors[2],
                                "gMCS_T2"   = ThresholdColors[3],"gMCS_T2.5" = ThresholdColors[4],
                                "gMCS_T5"   = ThresholdColors[5],"gMCS_T10"  = ThresholdColors[6],
                                "gMCS_T20"  = ThresholdColors[7],
                                "localT2" = "#B2DF8A", "localT2HumanGEM" = "#6ba35b")) +
  theme_bw() + xlab("") + ylab("") + 
  theme(legend.position = "bottom") + 
  scale_x_discrete(labels=c("Order_1" = "1",
                            "Order_2" = "2",
                            "Order_3" = "3",
                            "Order_4" = "4",
                            "Order_5" = "5",
                            "Order_6" = "6",
                            "Order_7" = "7",
                            "HigherThan_Order_7" = ">7")) +
  scale_y_continuous(limits = c(125,225), breaks = c(seq(from = 125, to = 225, by = 10))) +
  theme(strip.text = element_text(face  = "bold", size = 28), 
        strip.background = element_rect(fill="#b09c7b"), axis.text.x = element_text(size = 23),
        axis.text.y = element_text(size = 23))
FPlot
# Prepare PPV Plot #
PositivePredictedValue_tmp          <- tmp %>% filter(variable == "Positive Predicted Value")
PositivePredictedValue_tmp_mean <- PositivePredictedValue_tmp %>% 
  group_by(Order, thMethod) %>% 
  summarize(Mean = mean(value)) %>%
  ungroup()
if(is.na(sum(PositivePredictedValue_tmp_mean$Mean))){
  PositivePredictedValue_tmp_mean$Mean[which(is.na(PositivePredictedValue_tmp_mean$Mean))] <-  0
}
PositivePredictedValue_tmp_mean$variable <-  "Positive Predicted Value"
change_lab <- "Positive Predicted Value"; names(change_lab) <- "Positive Predicted Value"
PPVlot <- ggplot(PositivePredictedValue_tmp_mean, aes(x = Order, y = Mean, color = thMethod, group = thMethod)) + 
  geom_point(size = 4.5, alpha = 0.75) + geom_line(size = 1.5, alpha = 0.75) + 
  facet_wrap("variable",labeller = labeller(variable = change_lab))  +
  scale_color_manual(values = c("gMCS_T0"   = ThresholdColors[1],"gMCS_T1"   = ThresholdColors[2],
                                "gMCS_T2"   = ThresholdColors[3],"gMCS_T2.5" = ThresholdColors[4],
                                "gMCS_T5"   = ThresholdColors[5],"gMCS_T10"  = ThresholdColors[6],
                                "gMCS_T20"  = ThresholdColors[7],
                                "localT2" = "#B2DF8A", "localT2HumanGEM" = "#6ba35b")) +
  theme_bw() + xlab("") + ylab("") + 
  theme(legend.position = "bottom") + 
  scale_x_discrete(labels=c("Order_1" = "1",
                            "Order_2" = "2",
                            "Order_3" = "3",
                            "Order_4" = "4",
                            "Order_5" = "5",
                            "Order_6" = "6",
                            "Order_7" = "7",
                            "HigherThan_Order_7" = ">7")) +
  scale_y_continuous(limits = c(0.35,0.44), breaks = c(seq(from = 0.35, to = 0.44, by = 0.01))) +
  theme(strip.text = element_text(face  = "bold", size = 28), 
        strip.background = element_rect(fill="#b09c7b"), axis.text.x = element_text(size = 23),
        axis.text.y = element_text(size = 23)) 
PPVlot
# Add All Three Together #
LinelistLength <-  list(TPlot, FPlot, PPVlot)
Lineplotting_Length <- ggarrange(plotlist = LinelistLength,
                                 ncol = 3, nrow = 1, common.legend =  TRUE, legend = "bottom")
Lineplotting_Length_Unlegend <- ggarrange(plotlist = LinelistLength,
                                          ncol = 3, nrow = 1, common.legend =  TRUE, legend = "none")
# Lineplotting_Length <- annotate_figure(p = Lineplotting_Length,
#                                        top = text_grob("Lengths in different threshold methods", face = "bold", size = 18))
Lineplotting_Length
ggsave(filename = paste0(".Results/DepMap/Threshold/",Today,"_DepMap_LineLengths.png"),
       plot = Lineplotting_Length, device = "png", width = 21, height = 8,units = "in",dpi = 300)
ggsave(filename = paste0(".Results/DepMap/Threshold/",Today,"_DepMap_LineLengths_NoLegend.png"),
       plot = Lineplotting_Length_Unlegend, device = "png", width = 21, height = 8,units = "in",dpi = 300)

# Figure 2D - DepMap ####

toAdd <- do.call(rbind,readRDS("./RDSResults/DepMap/AllThresholdsResultsNotCumulative.RDS"))
toKeep <- c("gMCS_T0","gMCS_T1","gMCS_T2","gMCS_T2.5","gMCS_T5","gMCS_T10","gMCS_T20")
toRemove <- setdiff(unique(toAdd$thMethod), toKeep)
myColors <- viridis::viridis(10)

for(element in toRemove){toAdd <- toAdd[which(toAdd$thMethod != element),]}
toAdd$thMethod <- factor(toAdd$thMethod, levels=toKeep)

tmp     <- rbind(toAdd,
                 readRDS("./RDSResults/DepMap/Results_localT2.RDS"),
                 readRDS("./RDSResults/DepMap/Results_localT2HumanGEM.RDS"))
tmp <-  tmp %>% filter(thMethod == "gMCS_T5")


tmp <- reshape2::melt(tmp, id.vars = c("Order"),
                      measure.vars =  c("True Positives", "False Positives", "Positive Predicted Value"))

# Prepare TP Plot #
TruePositive_tmp          <- tmp %>% filter(variable == "True Positives")
TruePositive_tmp_mean <- TruePositive_tmp %>% 
  group_by(Order) %>% 
  summarize(Mean = mean(value)) %>%
  ungroup()
if(is.na(sum(TruePositive_tmp_mean$Mean))){
  TruePositive_tmp_mean$Mean[which(is.na(TruePositive_tmp_mean$Mean))] <-  0
}


change_lab <- "True Positives"; names(change_lab) <- "True Positives"
TPlot <- ggplot(TruePositive_tmp, aes(x = Order, y = value, fill = Order)) + 
  geom_boxplot(alpha = 0.4) + 
  geom_point(data    = TruePositive_tmp_mean,
             mapping = aes(x = Order, y = Mean), size = 2.5,
             color   = "dark red") + 
  geom_line(data = TruePositive_tmp_mean, 
            mapping = aes(x = Order, y = Mean, group = 1), size = 0.9) +
  
  facet_wrap("variable",labeller = labeller(variable = change_lab))  +
  scale_fill_manual(values = c("Order_1" = myColors[1],
                               "Order_2" = myColors[2],
                               "Order_3" = myColors[3],
                               "Order_4" = myColors[4],
                               "Order_5" = myColors[5],
                               "Order_6" = myColors[6],
                               "Order_7" = myColors[7],
                               "HigherThan_Order_7" = myColors[8])) + 
  theme_bw() + xlab("") + ylab("") +
  theme(legend.position = "none") + 
  scale_x_discrete(labels=c("Order_1" = "1",
                            "Order_2" = "2",
                            "Order_3" = "3",
                            "Order_4" = "4",
                            "Order_5" = "5",
                            "Order_6" = "6",
                            "Order_7" = "7",
                            "HigherThan_Order_7" = ">7")) +
  scale_y_continuous(limits = c(0,NA), breaks = c(seq(from = 0, to = 150, by = 15))) +
  theme(strip.text  = element_text(face  = "bold", size = 28),
        strip.background = element_rect(fill="#b09c7b"), axis.text.x = element_text(size = 23),
        axis.text.y = element_text(size = 23)) 
TPlot
# Prepare FP Plot #
FalsePositive_tmp          <- tmp %>% filter(variable == "False Positives")
FalsePositive_tmp_mean <- FalsePositive_tmp %>% 
  group_by(Order) %>% 
  summarize(Mean = mean(value)) %>%
  ungroup()
if(is.na(sum(FalsePositive_tmp_mean$Mean))){
  FalsePositive_tmp_mean$Mean[which(is.na(FalsePositive_tmp_mean$Mean))] <-  0
}

change_lab <- "False Positives"; names(change_lab) <- "False Positives"
FPlot <- ggplot(FalsePositive_tmp, aes(x = Order, y = value, fill = Order)) + 
  geom_boxplot(alpha = 0.6) + 
  facet_wrap("variable",labeller = labeller(variable = change_lab))  +
  scale_fill_manual(values = c("Order_1" = myColors[1],
                               "Order_2" = myColors[2],
                               "Order_3" = myColors[3],
                               "Order_4" = myColors[4],
                               "Order_5" = myColors[5],
                               "Order_6" = myColors[6],
                               "Order_7" = myColors[7],
                               "HigherThan_Order_7" = myColors[8])) + 
  theme_bw() + xlab("Length") + ylab("") +
  theme(legend.position = "none") + 
  geom_point(data    = FalsePositive_tmp_mean,
             mapping = aes(x = Order, y = Mean), size = 2.5,
             color   = "dark red") + 
  geom_line(data = FalsePositive_tmp_mean, 
            mapping = aes(x = Order, y = Mean, group = 1), size  = 0.9) +
  scale_x_discrete(labels=c("Order_1" = "1",
                            "Order_2" = "2",
                            "Order_3" = "3",
                            "Order_4" = "4",
                            "Order_5" = "5",
                            "Order_6" = "6",
                            "Order_7" = "7",
                            "HigherThan_Order_7" = ">7")) +
  scale_y_continuous(limits = c(0,NA), breaks = c(seq(from = 0, to = 150, by = 15))) +
  theme(strip.text = element_text(face  = "bold", size = 28),
        strip.background = element_rect(fill="#b09c7b"), axis.text.x = element_text(size = 23),
        axis.text.y = element_text(size = 23)) 
FPlot
# Prepare PPV Plot #
PositivePredictedValue_tmp          <- tmp %>% filter(variable == "Positive Predicted Value")
PositivePredictedValue_tmp_mean <- PositivePredictedValue_tmp %>% 
  group_by(Order) %>% 
  summarize(Mean = mean(value)) %>%
  ungroup()
if(is.na(sum(PositivePredictedValue_tmp_mean$Mean))){
  PositivePredictedValue_tmp_mean$Mean[which(is.na(PositivePredictedValue_tmp_mean$Mean))] <-  0
}
change_lab <- "Positive Predicted Value"; names(change_lab) <- "Positive Predicted Value"
PPVlot <- ggplot(PositivePredictedValue_tmp, aes(x = Order, y = value, fill = Order)) + 
  geom_boxplot(alpha = 0.6) +  
  facet_wrap("variable",labeller = labeller(variable = change_lab))  +
  scale_fill_manual(values = c("Order_1" = myColors[1],
                               "Order_2" = myColors[2],
                               "Order_3" = myColors[3],
                               "Order_4" = myColors[4],
                               "Order_5" = myColors[5],
                               "Order_6" = myColors[6],
                               "Order_7" = myColors[7],
                               "HigherThan_Order_7" = myColors[8])) + 
  geom_point(data    = PositivePredictedValue_tmp_mean,
             mapping = aes(x = Order, y = Mean), size  = 2.5,
             color   = "dark red") + 
  geom_line(data = PositivePredictedValue_tmp_mean, 
            mapping = aes(x = Order, y = Mean, group = 1), size  = 0.9) +
  theme_bw() + xlab("") + ylab("") +
  theme(legend.position = "none") +
  scale_x_discrete(labels=c("Order_1" = "1",
                            "Order_2" = "2",
                            "Order_3" = "3",
                            "Order_4" = "4",
                            "Order_5" = "5",
                            "Order_6" = "6",
                            "Order_7" = "7",
                            "HigherThan_Order_7" = ">7")) +
  scale_y_continuous(limits = c(0,NA), breaks = c(seq(from = 0, to = 1, by = 0.1))) +
  theme(strip.text = element_text(face  = "bold", size = 28), 
        strip.background = element_rect(fill="#b09c7b"), axis.text.x = element_text(size = 23),
        axis.text.y = element_text(size = 23)) 
PPVlot
# Add All Three Together #
plotlistLength <-  list(TPlot, FPlot, PPVlot)
Boxplotting_Length <- ggarrange(plotlist = plotlistLength,
                                ncol = 3, nrow = 1)
# Boxplotting_Length <- annotate_figure(p = Boxplotting_Length,
#                                       top = text_grob("gMCST5 Length Analysis", face = "bold", size = 18))
Boxplotting_Length
ggsave(filename = paste0(".Results/DepMap/Length/",Today,"_DepMap_Lengths_gMCSTH5.png"),
       plot = Boxplotting_Length, device = "png", width = 21, height = 8,units = "in",dpi = 300)
# Figure 3B - DepMap ####
gMCS_to_Many_ResultsT2 <- readRDS("./RDSResults/DepMap/2023-04-25_One_to_Many_Results_gMCST2.RDS")
gMCS_to_Many_ResultsT2$thMethod <- "gMCSth2"
gMCS_to_Many_ResultsT5 <- readRDS("./RDSResults/DepMap/2023-04-25_One_to_Many_Results_gMCST5.RDS")
gMCS_to_Many_ResultsT5$thMethod <- "gMCSth5"
gMCS_to_Many_ResultslocalT2 <- readRDS("./RDSResults/DepMap/2023-04-25_One_to_Many_Results_gMCSlocalT2.RDS")
gMCS_to_Many_ResultslocalT2$thMethod <- "localT2"

AllgMCS_toMany <- rbind(gMCS_to_Many_ResultsT2,gMCS_to_Many_ResultsT5,gMCS_to_Many_ResultslocalT2)

colnames(AllgMCS_toMany)[grep("TP", colnames(AllgMCS_toMany))]  <- "True Positives"
colnames(AllgMCS_toMany)[grep("FP", colnames(AllgMCS_toMany))]  <- "False Positives"
colnames(AllgMCS_toMany)[grep("PPV",colnames(AllgMCS_toMany))] <- "Positive Predictive Value"
myColors <- viridis::viridis(4)

AllgMCS_toMany <- reshape2::melt(AllgMCS_toMany, id.vars = c("calcMethod", "thMethod"),
                                 measure.vars =  c("True Positives", "False Positives", "Positive Predictive Value"))
AllgMCS_toMany <- as_tibble(AllgMCS_toMany)
AllgMCS_toMany <- AllgMCS_toMany[AllgMCS_toMany$calcMethod != "TrueEssentials",]
# CAREFUL, THE PRELOADED DATA HAS THESE FACTORS WRONG, IT IS LATER CHANGED, BUT WE ARE STILL WORKING WITH gMCSs
AllgMCS_toMany$calcMethod <- factor(AllgMCS_toMany$calcMethod, 
                                    levels = c("One_Gene_to_One_Task_NoEssentials", "One_Gene_to_Many_Tasks_NoEssentials"))
# Prepare TP Plot #
TruePositive_AllgMCS_toMany          <- AllgMCS_toMany %>% filter(variable == "True Positives")
for(i in unique(TruePositive_AllgMCS_toMany$thMethod)){
  tmp <- TruePositive_AllgMCS_toMany[TruePositive_AllgMCS_toMany$thMethod == i,]
  tmp$value <- as.numeric(tmp$value)
  print(wilcox.test(value~calcMethod, data = tmp) )
}
TruePositive_AllgMCS_toMany_mean <- TruePositive_AllgMCS_toMany %>% 
  group_by(calcMethod) %>% 
  summarize(Mean = mean(value)) %>%
  ungroup()
if(is.na(sum(TruePositive_AllgMCS_toMany_mean$Mean))){
  TruePositive_AllgMCS_toMany_mean$Mean[which(is.na(TruePositive_AllgMCS_toMany_mean$Mean))] <-  0
}


change_lab <- "True Positives"; names(change_lab) <- "True Positives"
TPlot <- ggplot(TruePositive_AllgMCS_toMany, aes(x = calcMethod, y = value, fill = thMethod)) + 
  geom_boxplot(alpha = 0.6) + 
  facet_wrap("variable",labeller = labeller(variable = change_lab))  +
  scale_fill_manual(values = c("gMCSth2" = viridis::magma(8)[3],
                               "gMCSth5" = viridis::magma(8)[5],
                               "localT2" = "#B2DF8A"), 
                    labels = c("gMCSth2"   = "2%","gMCSth5"   = "5%","localT2"   = "LocalT2"),
                    name = "Threshold Method") + 
  scale_color_manual(values = c("gMCSth2" = viridis::magma(8)[3],
                                "gMCSth5" = viridis::magma(8)[5],
                                "localT2" = "#B2DF8A"), 
                     labels = c("gMCSth2"   = "2%","gMCSth5"   = "5%","localT2"   = "LocalT2"),
                     guide = "none") + 
  scale_alpha(guide = "none") + scale_size(guide = "none") +
  theme_bw() + xlab("") + ylab("") +
  geom_vline(xintercept = 1.5, size = 2, color = "dark red") +
  theme(legend.position = "bottom") + 
  scale_x_discrete(labels=c("One_Gene_to_Many_Tasks_NoEssentials" = "One-to-Many",
                            "One_Gene_to_One_Task_NoEssentials" = "One-to-One")) +
  scale_y_continuous(limits = c(0,NA), breaks = c(seq(from = 0, to = 45, by = 5))) +
  theme(strip.text  = element_text(face  = "bold", size = 28), 
        strip.background = element_rect(fill = "#b09c7b"),
        axis.text.x = element_text(size = 23),
        axis.text.y = element_text(size = 23)) 
TPlot
# Prepare FP Plot #
FalsePositive_AllgMCS_toMany          <- AllgMCS_toMany %>% filter(variable == "False Positives")
for(i in unique(FalsePositive_AllgMCS_toMany$thMethod)){
  tmp <- FalsePositive_AllgMCS_toMany[FalsePositive_AllgMCS_toMany$thMethod == i,]
  tmp$value <- as.numeric(tmp$value)
  print(wilcox.test(value~calcMethod, data = tmp) )
}
FalsePositive_AllgMCS_toMany_mean <- FalsePositive_AllgMCS_toMany %>% 
  group_by(calcMethod) %>% 
  summarize(Mean = mean(value)) %>%
  ungroup()
if(is.na(sum(FalsePositive_AllgMCS_toMany_mean$Mean))){
  FalsePositive_AllgMCS_toMany_mean$Mean[which(is.na(FalsePositive_AllgMCS_toMany_mean$Mean))] <-  0
}

change_lab <- "False Positives"; names(change_lab) <- "False Positives"
FPlot <- ggplot(FalsePositive_AllgMCS_toMany, aes(x = calcMethod, y = value, fill = thMethod)) + 
  geom_boxplot(alpha = 0.6) + 
  facet_wrap("variable",labeller = labeller(variable = change_lab))  +
  scale_fill_manual(values = c("gMCSth2" = viridis::magma(8)[3],
                               "gMCSth5" = viridis::magma(8)[5],
                               "localT2" = "#B2DF8A"), 
                    labels = c("gMCSth2"   = "2%","gMCSth5"   = "5%","localT2"   = "LocalT2"),
                    name = "Threshold Method") + 
  scale_color_manual(values = c("gMCSth2" = viridis::magma(8)[3],
                                "gMCSth5" = viridis::magma(8)[5],
                                "localT2" = "#B2DF8A"), 
                     labels = c("gMCSth2"   = "2%","gMCSth5"   = "5%","localT2"   = "LocalT2"),
                     guide = "none") + 
  scale_alpha(guide = "none") + scale_size(guide = "none") +
  theme_bw() + xlab("") + ylab("") +
  geom_vline(xintercept = 1.5, size = 2, color = "dark red") +
  theme(legend.position = "bottom") + 
  scale_x_discrete(labels=c("One_Gene_to_Many_Tasks_NoEssentials" = "One-to-Many",
                            "One_Gene_to_One_Task_NoEssentials" = "One-to-One")) +
  scale_y_continuous(limits = c(0,NA), breaks = c(seq(from = 0, to = 100, by = 10))) +
  theme(strip.text  = element_text(face  = "bold", size = 28), 
        strip.background = element_rect(fill = "#b09c7b"),
        axis.text.x = element_text(size = 23),
        axis.text.y = element_text(size = 23)) 
FPlot
# Prepare PPV Plot #
PositivePredictiveValue_AllgMCS_toMany          <- AllgMCS_toMany %>% filter(variable == "Positive Predictive Value")
for(i in unique(PositivePredictiveValue_AllgMCS_toMany$thMethod)){
  tmp <- PositivePredictiveValue_AllgMCS_toMany[PositivePredictiveValue_AllgMCS_toMany$thMethod == i,]
  tmp$value <- as.numeric(tmp$value)
  print(wilcox.test(value~calcMethod, data = tmp) )
}
PositivePredictiveValue_AllgMCS_toMany_mean <- PositivePredictiveValue_AllgMCS_toMany %>% 
  group_by(calcMethod) %>% 
  summarize(Mean = mean(value)) %>%
  ungroup()
if(is.na(sum(PositivePredictiveValue_AllgMCS_toMany_mean$Mean))){
  PositivePredictiveValue_AllgMCS_toMany_mean$Mean[which(is.na(PositivePredictiveValue_AllgMCS_toMany_mean$Mean))] <-  0
}
change_lab <- "Positive Predictive Value"; names(change_lab) <- "Positive Predictive Value"
PPVlot <-ggplot(PositivePredictiveValue_AllgMCS_toMany, aes(x = calcMethod, y = value, fill = thMethod)) + 
  geom_boxplot(alpha = 0.6) + 
  facet_wrap("variable",labeller = labeller(variable = change_lab))  +
  scale_fill_manual(values = c("gMCSth2" = viridis::magma(8)[3],
                               "gMCSth5" = viridis::magma(8)[5],
                               "localT2" = "#B2DF8A"), 
                    labels = c("gMCSth2"   = "2%","gMCSth5"   = "5%","localT2"   = "LocalT2"),
                    name = "Threshold Method") + 
  scale_color_manual(values = c("gMCSth2" = viridis::magma(8)[3],
                                "gMCSth5" = viridis::magma(8)[5],
                                "localT2" = "#B2DF8A"), 
                     labels = c("gMCSth2"   = "2%","gMCSth5"   = "5%","localT2"   = "LocalT2"),
                     guide = "none") + 
  scale_alpha(guide = "none") + scale_size(guide = "none") +
  theme_bw() + xlab("") + ylab("") +
  geom_vline(xintercept = 1.5, size = 2, color = "dark red") +
  theme(legend.position = "bottom") + 
  scale_x_discrete(labels=c("One_Gene_to_Many_Tasks_NoEssentials" = "One-to-Many",
                            "One_Gene_to_One_Task_NoEssentials" = "One-to-One")) +
  scale_y_continuous(limits = c(0,NA), breaks = c(seq(from = 0, to = 1, by = 0.1))) +
  theme(strip.text  = element_text(face  = "bold", size = 28), 
        strip.background = element_rect(fill = "#b09c7b"),
        axis.text.x = element_text(size = 23),
        axis.text.y = element_text(size = 23)) 
PPVlot
# Add All Three Together #
plotlistgMCS <-  list(TPlot, FPlot, PPVlot)
Boxplotting_gMCS <- ggarrange(plotlist = plotlistgMCS,
                              ncol = 3, nrow = 1, common.legend = TRUE)
Boxplotting_gMCS_Nolegend <- ggarrange(plotlist = plotlistgMCS,
                                       ncol = 3, nrow = 1, common.legend = TRUE, legend = "none")
Boxplotting_gMCS
ggsave(filename = paste0(".Results/DepMap/OneToMany/",Today,"_DepMap_gMCS2Many_Figure.png"),
       plot = Boxplotting_gMCS, device = "png", width = 21, height = 8,units = "in",dpi = 300)
ggsave(filename = paste0(".Results/DepMap/OneToMany/",Today,"_DepMap_gMCS2Many_Figure_NOLegend.png"),
       plot = Boxplotting_gMCS_Nolegend, device = "png", width = 21, height = 8,units = "in",dpi = 300)
# Figure 3D - DepMap ####
Task_to_Many_ResultsT2 <- readRDS("./RDSResults/DepMap/2023-04-25_Task_Promiscuity_Results_gMCST2.RDS")
Task_to_Many_ResultsT2$thMethod <- "gMCSth2"
Task_to_Many_ResultsT5 <- readRDS("./RDSResults/DepMap/2023-04-25_Task_Promiscuity_Results_gMCST5.RDS")
Task_to_Many_ResultsT5$thMethod <- "gMCSth5"
Task_to_Many_ResultslocalT2 <- readRDS("./RDSResults/DepMap/2023-04-25_Task_Promiscuity_Results_localT2.RDS")
Task_to_Many_ResultslocalT2$thMethod <- "localT2"

AllTasks_toMany <- rbind(Task_to_Many_ResultsT2,Task_to_Many_ResultsT5,Task_to_Many_ResultslocalT2)

colnames(AllTasks_toMany)[grep("TP", colnames(AllTasks_toMany))]  <- "True Positives"
colnames(AllTasks_toMany)[grep("FP", colnames(AllTasks_toMany))]  <- "False Positives"
colnames(AllTasks_toMany)[grep("PPV",colnames(AllTasks_toMany))] <- "Positive Predictive Value"
myColors <- viridis::viridis(4)

AllTasks_toMany <- reshape2::melt(AllTasks_toMany, id.vars = c("calcMethod", "thMethod"),
                                  measure.vars =  c("True Positives", "False Positives", "Positive Predictive Value"))
AllTasks_toMany <- as_tibble(AllTasks_toMany)
AllTasks_toMany <- AllTasks_toMany[AllTasks_toMany$calcMethod != "TrueEssentials",]
AllTasks_toMany$calcMethod <- factor(AllTasks_toMany$calcMethod, 
                                     levels = c("One_Gene_to_One_Task_NoEssentials", "One_Gene_to_Many_Tasks_NoEssentials"))
# Prepare TP Plot #
TruePositive_AllTasks_toMany          <- AllTasks_toMany %>% filter(variable == "True Positives")
for(i in unique(TruePositive_AllTasks_toMany$thMethod)){
  tmp <- TruePositive_AllTasks_toMany[TruePositive_AllTasks_toMany$thMethod == i,]
  tmp$value <- as.numeric(tmp$value)
  print(wilcox.test(value~calcMethod, data = tmp) )
}
TruePositive_AllTasks_toMany_mean <- TruePositive_AllTasks_toMany %>% 
  group_by(calcMethod) %>% 
  summarize(Mean = mean(value)) %>%
  ungroup()
if(is.na(sum(TruePositive_AllTasks_toMany_mean$Mean))){
  TruePositive_AllTasks_toMany_mean$Mean[which(is.na(TruePositive_AllTasks_toMany_mean$Mean))] <-  0
}


change_lab <- "True Positives"; names(change_lab) <- "True Positives"
TPlot <- ggplot(TruePositive_AllTasks_toMany, aes(x = calcMethod, y = value,fill = thMethod)) + 
  geom_boxplot(alpha = 0.6) + 
  facet_wrap("variable",labeller = labeller(variable = change_lab))  +
  scale_fill_manual(values = c("gMCSth2" = viridis::magma(8)[3],
                               "gMCSth5" = viridis::magma(8)[5],
                               "localT2" = "#B2DF8A"), 
                    labels = c("gMCSth2"   = "2%","gMCSth5"   = "5%","localT2"   = "LocalT2"),
                    name = "Threshold Method") + 
  scale_color_manual(values = c("gMCSth2" = viridis::magma(8)[3],
                                "gMCSth5" = viridis::magma(8)[5],
                                "localT2" = "#B2DF8A"), 
                     labels = c("gMCSth2"   = "2%","gMCSth5"   = "5%","localT2"   = "LocalT2"),
                     guide = "none") + 
  scale_alpha(guide = "none") + scale_size(guide = "none") +
  geom_vline(xintercept = 1.5, size = 2, color = "dark red") +
  theme_bw() + xlab("") + ylab("") +
  theme(legend.position = "bottom") + 
  scale_x_discrete(labels=c("One_Gene_to_Many_Tasks_NoEssentials" = "Miscellaneous",
                            "One_Gene_to_One_Task_NoEssentials" = "Focused")) +
  scale_y_continuous(limits = c(0,NA), breaks = c(seq(from = 0, to = 75, by = 5))) +
  theme(strip.text  = element_text(face  = "bold", size = 28), 
        strip.background = element_rect(fill = "#b09c7b"),
        axis.text.x = element_text(size = 23),
        axis.text.y = element_text(size = 23)) 
TPlot
# Prepare FP Plot #
FalsePositive_AllTasks_toMany          <- AllTasks_toMany %>% filter(variable == "False Positives")
for(i in unique(FalsePositive_AllTasks_toMany$thMethod)){
  tmp <- FalsePositive_AllTasks_toMany[FalsePositive_AllTasks_toMany$thMethod == i,]
  tmp$value <- as.numeric(tmp$value)
  print(wilcox.test(value~calcMethod, data = tmp) )
}
FalsePositive_AllTasks_toMany_mean <- FalsePositive_AllTasks_toMany %>% 
  group_by(calcMethod) %>% 
  summarize(Mean = mean(value)) %>%
  ungroup()
if(is.na(sum(FalsePositive_AllTasks_toMany_mean$Mean))){
  FalsePositive_AllTasks_toMany_mean$Mean[which(is.na(FalsePositive_AllTasks_toMany_mean$Mean))] <-  0
}

change_lab <- "False Positives"; names(change_lab) <- "False Positives"
FPlot <- ggplot(FalsePositive_AllTasks_toMany, aes(x = calcMethod, y = value, fill = thMethod)) + 
  geom_boxplot(alpha = 0.6) + 
  facet_wrap("variable",labeller = labeller(variable = change_lab))  +
  scale_fill_manual(values = c("gMCSth2" = viridis::magma(8)[3],
                               "gMCSth5" = viridis::magma(8)[5],
                               "localT2" = "#B2DF8A"), 
                    labels = c("gMCSth2"   = "2%","gMCSth5"   = "5%","localT2"   = "LocalT2"),
                    name = "Threshold Method") + 
  scale_color_manual(values = c("gMCSth2" = viridis::magma(8)[3],
                                "gMCSth5" = viridis::magma(8)[5],
                                "localT2" = "#B2DF8A"), 
                     labels = c("gMCSth2"   = "2%","gMCSth5"   = "5%","localT2"   = "LocalT2"),
                     guide = "none") + 
  scale_alpha(guide = "none") + scale_size(guide = "none") +
  geom_vline(xintercept = 1.5, size = 2, color = "dark red") +
  theme_bw() + xlab("") + ylab("") +
  theme(legend.position = "bottom") + 
  scale_x_discrete(labels=c("One_Gene_to_Many_Tasks_NoEssentials" = "Miscellaneous",
                            "One_Gene_to_One_Task_NoEssentials" = "Focused")) +
  scale_y_continuous(limits = c(0,NA), breaks = c(seq(from = 0, to = 375, by = 25))) +
  theme(strip.text  = element_text(face  = "bold", size = 28), 
        strip.background = element_rect(fill = "#b09c7b"),
        axis.text.x = element_text(size = 23),
        axis.text.y = element_text(size = 23)) 
FPlot
# Prepare PPV Plot #
PositivePredictiveValue_AllTasks_toMany          <- AllTasks_toMany %>% filter(variable == "Positive Predictive Value")
for(i in unique(PositivePredictiveValue_AllTasks_toMany$thMethod)){
  tmp <- PositivePredictiveValue_AllTasks_toMany[PositivePredictiveValue_AllTasks_toMany$thMethod == i,]
  tmp$value <- as.numeric(tmp$value)
  print(wilcox.test(value~calcMethod, data = tmp) )
}
PositivePredictiveValue_AllTasks_toMany_mean <- PositivePredictiveValue_AllTasks_toMany %>% 
  group_by(calcMethod) %>% 
  summarize(Mean = mean(value)) %>%
  ungroup()
if(is.na(sum(PositivePredictiveValue_AllTasks_toMany_mean$Mean))){
  PositivePredictiveValue_AllTasks_toMany_mean$Mean[which(is.na(PositivePredictiveValue_AllTasks_toMany_mean$Mean))] <-  0
}
change_lab <- "Positive Predictive Value"; names(change_lab) <- "Positive Predictive Value"
PPVlot <-ggplot(PositivePredictiveValue_AllTasks_toMany, aes(x = calcMethod, y = value,fill = thMethod)) + 
  geom_boxplot(alpha = 0.6) + 
  facet_wrap("variable",labeller = labeller(variable = change_lab))  +
  scale_fill_manual(values = c("gMCSth2" = viridis::magma(8)[3],
                               "gMCSth5" = viridis::magma(8)[5],
                               "localT2" = "#B2DF8A"), 
                    labels = c("gMCSth2"   = "2%","gMCSth5"   = "5%","localT2"   = "LocalT2"),
                    name = "Threshold Method") + 
  scale_color_manual(values = c("gMCSth2" = viridis::magma(8)[3],
                                "gMCSth5" = viridis::magma(8)[5],
                                "localT2" = "#B2DF8A"), 
                     labels = c("gMCSth2"   = "2%","gMCSth5"   = "5%","localT2"   = "LocalT2"),
                     guide = "none") + 
  scale_alpha(guide = "none") + scale_size(guide = "none") +
  geom_vline(xintercept = 1.5, size = 2, color = "dark red") +
  theme_bw() + xlab("") + ylab("") +
  theme(legend.position = "bottom") + 
  scale_x_discrete(labels=c("One_Gene_to_Many_Tasks_NoEssentials" = "Miscellaneous",
                            "One_Gene_to_One_Task_NoEssentials" = "Focused")) +
  scale_y_continuous(limits = c(0,NA), breaks = c(seq(from = 0, to = 1, by = 0.1))) +
  theme(strip.text  = element_text(face  = "bold", size = 28), 
        strip.background = element_rect(fill = "#b09c7b"),
        axis.text.x = element_text(size = 23),
        axis.text.y = element_text(size = 23)) 
PPVlot
# Add All Three Together #
plotlistgMCS <-  list(TPlot, FPlot, PPVlot)
Boxplotting_gMCS <- ggarrange(plotlist = plotlistgMCS,
                              ncol = 3, nrow = 1, common.legend = TRUE)
Boxplotting_gMCS_NoLegend <- ggarrange(plotlist = plotlistgMCS,
                                       ncol = 3, nrow = 1, common.legend = TRUE, legend = "none")
Boxplotting_gMCS
ggsave(filename = paste0(".Results/DepMap/OneToMany/",Today,"_DepMap_Task_Figure.png"),
       plot = Boxplotting_gMCS, device = "png", width = 21, height = 8,units = "in",dpi = 300)
ggsave(filename = paste0(".Results/DepMap/OneToMany/",Today,"_DepMap_Task_Figure_NoLegend.png"),
       plot = Boxplotting_gMCS_NoLegend, device = "png", width = 21, height = 8,units = "in",dpi = 300)
# Figure 4C - DepMap ####
ThresholdColors <- viridis::magma(8)
BiomassLines  <- do.call(rbind, readRDS("./Biomass/RDSResults/DepMap/AllThresholdsResults.RDS"))
BiomassLinesth     <- rbind(BiomassLines,
                            readRDS("./Biomass/RDSResults/DepMap/Results_localT2.RDS"),
                            readRDS("./Biomass/RDSResults/DepMap/Results_localT2HumanGEM.RDS"))
BiomassLinesth <- reshape2::melt(BiomassLinesth, id.vars = c("Order", "thMethod"),
                                 measure.vars =  c("True Positives", "False Positives", "Positive Predicted Value"))

BiomassLinesth$evalMethod <- "Biomass"
AllTasksLines <- do.call(rbind, readRDS("./RDSResults/DepMap/AllThresholdsResults.RDS"))
AllTasksLines     <- rbind(AllTasksLines,
                           readRDS("./RDSResults/DepMap/Results_localT2.RDS"),
                           readRDS("./RDSResults/DepMap/Results_localT2HumanGEM.RDS"))
AllTasksLinesth <- reshape2::melt(AllTasksLines, id.vars = c("Order", "thMethod"),
                                  measure.vars =  c("True Positives", "False Positives", "Positive Predicted Value"))
AllTasksLinesth$evalMethod <- "AllTasks"

BothLinesth <- rbind(BiomassLinesth, AllTasksLinesth)
BothLinesth <- BothLinesth[BothLinesth$thMethod == "gMCS_T2" | BothLinesth$thMethod == "gMCS_T5" | BothLinesth$thMethod == "localT2",]
BothLinesth$eval_th <- paste0(BothLinesth$thMethod,"--",BothLinesth$evalMethod)
BothLinesth$Order <- as.character(BothLinesth$Order)
BothLinesth$Order[BothLinesth$thMethod == "gMCS_T5"] <- paste0(BothLinesth$Order[BothLinesth$thMethod == "gMCS_T5"], "_15")
BothLinesth$Order[BothLinesth$thMethod == "localT2"] <- paste0(BothLinesth$Order[BothLinesth$thMethod == "localT2"], "_30")
BothLinesth$Order <- factor(BothLinesth$Order,
                            levels = c("Order_1","Order_1_15", "Order_1_30", 
                                       "Order_2","Order_2_15", "Order_2_30", 
                                       "Order_3","Order_3_15", "Order_3_30", 
                                       "Order_4","Order_4_15", "Order_4_30", 
                                       "Order_5","Order_5_15", "Order_5_30", 
                                       "Order_6","Order_6_15", "Order_6_30", 
                                       "Order_7","Order_7_15", "Order_7_30", 
                                       "HigherThan_Order_7",
                                       "HigherThan_Order_7_15", 
                                       "HigherThan_Order_7_30"))
# Prepare TP Plot #
TruePositive_BothLinesth          <- BothLinesth %>% filter(variable == "True Positives")
TruePositive_BothLinesth_mean <- TruePositive_BothLinesth %>% 
  group_by(Order, eval_th, evalMethod, thMethod) %>% 
  summarize(Mean = mean(value)) %>%
  ungroup()
if(is.na(sum(TruePositive_BothLinesth_mean$Mean))){
  TruePositive_BothLinesth_mean$Mean[which(is.na(TruePositive_BothLinesth_mean$Mean))] <-  0
}
TruePositive_BothLinesth_mean$variable <-  "True_Positives"
change_lab <- "True Positives"; names(change_lab) <- "True_Positives"
TruePositive_BothLinesth_mean$evalMethod <- factor(TruePositive_BothLinesth_mean$evalMethod,
                                                   levels = c("Biomass", "AllTasks"))
for(thresholding in unique(TruePositive_BothLinesth_mean$thMethod)){
  for(ordering in unique(TruePositive_BothLinesth_mean$Order)){
    
    Total <- TruePositive_BothLinesth_mean$Mean[TruePositive_BothLinesth_mean$Order == ordering &
                                                  TruePositive_BothLinesth_mean$thMethod == thresholding &
                                                  TruePositive_BothLinesth_mean$evalMethod == "AllTasks"]
    Substract <- TruePositive_BothLinesth_mean$Mean[TruePositive_BothLinesth_mean$Order == ordering &
                                                      TruePositive_BothLinesth_mean$thMethod == thresholding &
                                                      TruePositive_BothLinesth_mean$evalMethod != "AllTasks"]
    TruePositive_BothLinesth_mean$Mean[TruePositive_BothLinesth_mean$Order == ordering &
                                         TruePositive_BothLinesth_mean$thMethod == thresholding &
                                         TruePositive_BothLinesth_mean$evalMethod == "AllTasks"] <- Total - Substract
    
  }  
}


TPlot <- ggplot(TruePositive_BothLinesth_mean, aes(x = Order, y = Mean,
                                                   alpha = evalMethod, color = evalMethod,
                                                   group = eval_th, fill = thMethod)) + 
  geom_bar(stat = "identity") + 
  facet_wrap("variable",labeller = labeller(variable = change_lab))  +
  scale_color_manual(values = c("AllTasks"   = "black", "Biomass"   = "grey20")) +
  scale_alpha_manual(values = c("AllTasks" = 0.2, "Biomass" = 0.8)) +
  theme_bw() + xlab("") + ylab("") +
  guides(color=guide_legend(title="Analysis")) +
  theme(legend.position = "bottom") + 
  scale_x_discrete(labels=c("Order_1" = "",
                            "Order_1_15" = "1",
                            "Order_1_30" = "",
                            "Order_1_45" = "",
                            "Order_2" = "",
                            "Order_2_15" = "2",
                            "Order_2_30" = "",
                            "Order_2_45" = "",
                            "Order_3" = "",
                            "Order_3_15" = "3",
                            "Order_3_30" = "",
                            "Order_3_45" = "",
                            "Order_4" = "",
                            "Order_4_15" = "4",
                            "Order_4_30" = "",
                            "Order_4_45" = "",
                            "Order_5" = "",
                            "Order_5_15" = "5",
                            "Order_5_30" = "",
                            "Order_5_45" = "",
                            "Order_6" = "",
                            "Order_6_15" = "6",
                            "Order_6_30" = "",
                            "Order_6_45" = "",
                            "Order_7" = "",
                            "Order_7_15" = "7",
                            "Order_7_30" = "",
                            "Order_7_45" = "",
                            "HigherThan_Order_7" = "",
                            "HigherThan_Order_7_15" = ">7",
                            "HigherThan_Order_7_30" = ""
  ),
  limits = c("Order_1","Order_1_15","Order_1_30",
             "Order_1_45",
             "Order_2","Order_2_15","Order_2_30",
             "Order_2_45",
             "Order_3","Order_3_15","Order_3_30",
             "Order_3_45",
             "Order_4","Order_4_15","Order_4_30",
             "Order_4_45",
             "Order_5","Order_5_15","Order_5_30",
             "Order_5_45",
             "Order_6","Order_6_15","Order_6_30",
             "Order_6_45",
             "Order_7","Order_7_15","Order_7_30",
             "Order_7_45",
             "HigherThan_Order_7","HigherThan_Order_7_15","HigherThan_Order_7_30")) +
  scale_fill_manual(values = c("gMCS_T2" = viridis::magma(8)[3],
                               "gMCS_T5" = viridis::magma(8)[5],
                               "localT2" = "#B2DF8A"), 
                    labels = c("gMCS_T2"   = "2%","gMCS_T5"   = "5%","localT2"   = "LocalT2"),
                    name = "Threshold Method") + 
  scale_y_continuous(limits = c(0,NA), breaks = c(seq(from = 0, to = 250, by = 15))) +
  theme(strip.text = element_text(face  = "bold", size = 50),
        strip.background = element_rect(fill = "#b09c7b"),
        axis.text.x = element_text(size = 45),
        axis.ticks.x = element_blank(),
        axis.text.y = element_text(size = 45),
        legend.position = "none",
        legend.text = element_text(size = 18),
        legend.title    = element_text(size = 14))
TPlot
ggsave(filename = paste0("./BiomassResults/DepMap/Bio_vs_Tasks/",Today,"_DepMap_Bio_vs_Tasks_BarPlot_TP.png"),
       plot = TPlot, device = "png", width = 21, height = 12,units = "in",dpi = 300)

# Figure 4D - DepMap ####
FalsePositive_BothLinesth          <- BothLinesth %>% filter(variable == "False Positives")

FalsePositive_BothLinesth_mean <- FalsePositive_BothLinesth %>% 
  group_by(Order, eval_th, evalMethod, thMethod) %>% 
  summarize(Mean = mean(value)) %>%
  ungroup()

if(is.na(sum(FalsePositive_BothLinesth_mean$Mean))){
  FalsePositive_BothLinesth_mean$Mean[which(is.na(FalsePositive_BothLinesth_mean$Mean))] <-  0
}
FalsePositive_BothLinesth_mean$variable <-  "False_Positives"
change_lab <- "False Positives"; names(change_lab) <- "False_Positives"
FalsePositive_BothLinesth_mean$evalMethod <- factor(FalsePositive_BothLinesth_mean$evalMethod,
                                                    levels = c("Biomass", "AllTasks"))
for(thresholding in unique(FalsePositive_BothLinesth_mean$thMethod)){
  for(ordering in unique(FalsePositive_BothLinesth_mean$Order)){
    
    Total <- FalsePositive_BothLinesth_mean$Mean[FalsePositive_BothLinesth_mean$Order == ordering &
                                                   FalsePositive_BothLinesth_mean$thMethod == thresholding &
                                                   FalsePositive_BothLinesth_mean$evalMethod == "AllTasks"]
    Substract <- FalsePositive_BothLinesth_mean$Mean[FalsePositive_BothLinesth_mean$Order == ordering &
                                                       FalsePositive_BothLinesth_mean$thMethod == thresholding &
                                                       FalsePositive_BothLinesth_mean$evalMethod != "AllTasks"]
    FalsePositive_BothLinesth_mean$Mean[FalsePositive_BothLinesth_mean$Order == ordering &
                                          FalsePositive_BothLinesth_mean$thMethod == thresholding &
                                          FalsePositive_BothLinesth_mean$evalMethod == "AllTasks"] <- Total - Substract
    
  }  
}

FPlot <- ggplot(FalsePositive_BothLinesth_mean, aes(x = Order, y = Mean,
                                                    alpha = evalMethod, color = evalMethod,
                                                    group = eval_th, fill = thMethod)) + 
  geom_bar(stat = "identity") + 
  facet_wrap("variable",labeller = labeller(variable = change_lab))  +
  scale_color_manual(values = c("AllTasks"   = "black", "Biomass"   = "grey20")) +
  scale_alpha_manual(values = c("AllTasks" = 0.2, "Biomass" = 0.8)) +
  theme_bw() + xlab("") + ylab("") +
  guides(color=guide_legend(title="Analysis")) +
  theme(legend.position = "bottom") + 
  scale_x_discrete(labels=c("Order_1" = "",
                            "Order_1_15" = "1",
                            "Order_1_30" = "",
                            "Order_1_45" = "",
                            "Order_2" = "",
                            "Order_2_15" = "2",
                            "Order_2_30" = "",
                            "Order_2_45" = "",
                            "Order_3" = "",
                            "Order_3_15" = "3",
                            "Order_3_30" = "",
                            "Order_3_45" = "",
                            "Order_4" = "",
                            "Order_4_15" = "4",
                            "Order_4_30" = "",
                            "Order_4_45" = "",
                            "Order_5" = "",
                            "Order_5_15" = "5",
                            "Order_5_30" = "",
                            "Order_5_45" = "",
                            "Order_6" = "",
                            "Order_6_15" = "6",
                            "Order_6_30" = "",
                            "Order_6_45" = "",
                            "Order_7" = "",
                            "Order_7_15" = "7",
                            "Order_7_30" = "",
                            "Order_7_45" = "",
                            "HigherThan_Order_7" = "",
                            "HigherThan_Order_7_15" = ">7",
                            "HigherThan_Order_7_30" = ""
  ),
  limits = c("Order_1","Order_1_15","Order_1_30",
             "Order_1_45",
             "Order_2","Order_2_15","Order_2_30",
             "Order_2_45",
             "Order_3","Order_3_15","Order_3_30",
             "Order_3_45",
             "Order_4","Order_4_15","Order_4_30",
             "Order_4_45",
             "Order_5","Order_5_15","Order_5_30",
             "Order_5_45",
             "Order_6","Order_6_15","Order_6_30",
             "Order_6_45",
             "Order_7","Order_7_15","Order_7_30",
             "Order_7_45",
             "HigherThan_Order_7","HigherThan_Order_7_15","HigherThan_Order_7_30")) +
  scale_fill_manual(values = c("gMCS_T2" = viridis::magma(8)[3],
                               "gMCS_T5" = viridis::magma(8)[5],
                               "localT2" = "#B2DF8A"), 
                    labels = c("gMCS_T2"   = "2%","gMCS_T5"   = "5%","localT2"   = "LocalT2"),
                    name = "Threshold Method") + 
  scale_y_continuous(limits = c(0,NA), breaks = c(seq(from = 0, to = 500, by = 25))) +
  theme(strip.text = element_text(face  = "bold", size = 50),
        strip.background = element_rect(fill = "#b09c7b"),
        axis.text.x = element_text(size = 45),
        axis.ticks.x = element_blank(),
        axis.text.y = element_text(size = 45),
        legend.position = "none",
        legend.text = element_text(size = 18),
        legend.title    = element_text(size = 14))
FPlot
ggsave(filename = paste0("./BiomassResults/DepMap/Bio_vs_Tasks/",Today,"_DepMap_Bio_vs_Tasks_BarPlot_FP_NoLegend.png"),
       plot = FPlot, device = "png", width = 21, height = 12,units = "in",dpi = 300)

# Supp - PPV ####
PositivePredictedValue_BothLinesth          <- BothLinesth %>% filter(variable == "Positive Predicted Value")
PositivePredictedValue_BothLinesth_mean <- PositivePredictedValue_BothLinesth %>% 
  group_by(Order, eval_th, evalMethod, thMethod) %>% 
  summarize(Mean = mean(value)) %>%
  ungroup()
if(is.na(sum(PositivePredictedValue_BothLinesth_mean$Mean))){
  PositivePredictedValue_BothLinesth_mean$Mean[which(is.na(PositivePredictedValue_BothLinesth_mean$Mean))] <-  0
}
PositivePredictedValue_BothLinesth_mean$variable <-  "Positive Predicted Value"
change_lab <- "Positive Predicted Value"; names(change_lab) <- "Positive Predicted Value"
PositivePredictedValue_BothLinesth_mean$evalMethod <- factor(PositivePredictedValue_BothLinesth_mean$evalMethod,
                                                             levels = c("AllTasks", "Biomass"))
PPVlot <- ggplot(PositivePredictedValue_BothLinesth_mean, aes(x = Order, y = Mean,
                                                              alpha = evalMethod, color = evalMethod,
                                                              group = eval_th, fill = thMethod)) + 
  geom_bar(stat = "identity", position = position_dodge()) + 
  facet_wrap("variable",labeller = labeller(variable = change_lab))  +
  scale_color_manual(values = c("gMCS_T2"   = ThresholdColors[3], "gMCS_T5"   = ThresholdColors[5],
                                "localT2" = "#B2DF8A")) +
  scale_alpha_manual(values = c("AllTasks" = 0.6, "Biomass" = 0.8)) +
  theme_bw() + xlab("") + ylab("") +
  theme(legend.position = "bottom") + 
  scale_x_discrete(labels=c("Order_1" = "",
                            "Order_1_15" = "1",
                            "Order_1_30" = "",
                            "Order_1_45" = "",
                            "Order_2" = "",
                            "Order_2_15" = "2",
                            "Order_2_30" = "",
                            "Order_2_45" = "",
                            "Order_3" = "",
                            "Order_3_15" = "3",
                            "Order_3_30" = "",
                            "Order_3_45" = "",
                            "Order_4" = "",
                            "Order_4_15" = "4",
                            "Order_4_30" = "",
                            "Order_4_45" = "",
                            "Order_5" = "",
                            "Order_5_15" = "5",
                            "Order_5_30" = "",
                            "Order_5_45" = "",
                            "Order_6" = "",
                            "Order_6_15" = "6",
                            "Order_6_30" = "",
                            "Order_6_45" = "",
                            "Order_7" = "",
                            "Order_7_15" = "7",
                            "Order_7_30" = "",
                            "Order_7_45" = "",
                            "HigherThan_Order_7" = "",
                            "HigherThan_Order_7_15" = ">7",
                            "HigherThan_Order_7_30" = ""
  ),
  limits = c("Order_1","Order_1_15","Order_1_30",
             "Order_1_45",
             "Order_2","Order_2_15","Order_2_30",
             "Order_2_45",
             "Order_3","Order_3_15","Order_3_30",
             "Order_3_45",
             "Order_4","Order_4_15","Order_4_30",
             "Order_4_45",
             "Order_5","Order_5_15","Order_5_30",
             "Order_5_45",
             "Order_6","Order_6_15","Order_6_30",
             "Order_6_45",
             "Order_7","Order_7_15","Order_7_30",
             "Order_7_45",
             "HigherThan_Order_7","HigherThan_Order_7_15","HigherThan_Order_7_30")) +
  scale_fill_manual(values = c("gMCS_T2" = viridis::magma(8)[3],
                               "gMCS_T5" = viridis::magma(8)[5],
                               "localT2" = "#B2DF8A"), 
                    labels = c("gMCS_T2"   = "2%","gMCS_T5"   = "5%","localT2"   = "LocalT2"),
                    name = "Threshold Method") + 
  scale_y_continuous(limits = c(0,NA), breaks = c(seq(from = 0, to = 1, by = 0.1))) +
  theme(strip.text = element_text(face  = "bold", size = 23),
        strip.background = element_rect(fill = "#b09c7b"),
        axis.text.x = element_text(size = 23),
        axis.ticks.x = element_blank(),
        axis.text.y = element_text(size = 23),
        legend.key.size = unit(3, 'cm'),
        legend.text = element_text(size = 18),
        legend.title    = element_text(size = 14))
PPVlot
# Add All Three Together #
LinelistBioVSAll <-  list(TPlot, FPlot, PPVlot)
LineplottingBioVSAll <- ggarrange(plotlist = LinelistBioVSAll,
                                  ncol = 3, nrow = 1,
                                  common.legend =  TRUE, legend = "none")
LineplottingBioVSAll_NoLegend <- ggarrange(plotlist = LinelistBioVSAll,
                                           ncol = 3, nrow = 1,
                                           common.legend =  TRUE)

# ggsave(filename = paste0("./BiomassResults/DepMap/Bio_vs_Tasks/",Today,"_DepMap_Bio_vs_Tasks_BarPlot_Figure.png"),
#        plot = LineplottingBioVSAll, device = "png", width = 21, height = 8,units = "in",dpi = 300)
# ggsave(filename = paste0("./BiomassResults/DepMap/Bio_vs_Tasks/",Today,"_DepMap_Bio_vs_Tasks_BarPlot_Figure_NoLegend.png"),
#        plot = LineplottingBioVSAll_NoLegend, device = "png", width = 21, height = 8,units = "in",dpi = 300)
# Write Table #
toSave <- rbind(TruePositive_BothLinesth_mean, FalsePositive_BothLinesth_mean, PositivePredictedValue_BothLinesth_mean)
toSave$Order[toSave$Order == "Order_1_15" | toSave$Order == "Order_1_30"] <-  "Order_1"
toSave$Order[toSave$Order == "Order_2_15" | toSave$Order == "Order_2_30"] <-  "Order_2"
toSave$Order[toSave$Order == "Order_3_15" | toSave$Order == "Order_3_30"] <-  "Order_3"
toSave$Order[toSave$Order == "Order_4_15" | toSave$Order == "Order_4_30"] <-  "Order_4"
toSave$Order[toSave$Order == "Order_5_15" | toSave$Order == "Order_5_30"] <-  "Order_5"
toSave$Order[toSave$Order == "Order_6_15" | toSave$Order == "Order_6_30"] <-  "Order_6"
toSave$Order[toSave$Order == "Order_7_15" | toSave$Order == "Order_7_30"] <-  "Order_7"
toSave$Order[toSave$Order == "HigherThan_Order_7_15" | toSave$Order == "HigherThan_Order_7_30"] <-  "HigherThan_Order_7"
toSave$variable[toSave$variable == "Positive Predicted Value"] <- "Positive Predictive Value"
write.csv(toSave, file = "./DepMap_Biomass_AllTasks_Mean_Results.csv",quote = FALSE, row.names = FALSE)
