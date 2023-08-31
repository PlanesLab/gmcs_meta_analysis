TweakedDepMapgMCS_Results <- function(Binary_Matrix, 
                                      Effect,
                                      Table_HumanGEM_Genes,
                                      gMCS_order,
                                      Threshold_Info,
                                      evalMethod = "All_Tasks")
  # Input - DepMap #
  #' @param Binary_Matrix, binary matrix, 1 if the gene is essential 0 if it is not, per cell line.
  #' @param Effect, information regarding essentiality values based on Achilles scores.
  #' @param Table_HumanGEM_Genes, used to translate genes from SYMBOL to ENSEMBL or vice versa.
  #' @param gMCS_order, which is the order of gMCS we are working on, default = 8
  #' @param Threshold_Info, threshold information for statistics purposes.
  #' @param evalMethod, mode of localT2: only used genes (~1244) or all the Human1 universe?
  
  # Output: #
  #' @param results, a dataframe containing: TP/TN/FP/FN/Accuracy/Sensitivity/Specificity/F1 score/ Matthew's Correlation Coefficient
  #'  per cell line.

  
  try({
    result_list <- list()
    evalMethod = evalMethod
    # browser()
    ExpMatBinary <- Effect < (-0.6)
    ExpMatBinary[is.na(Effect)] <- 0 
    cellLines <- intersect(colnames(Binary_Matrix),
                         colnames(ExpMatBinary))
    
    # Initialize Variables #
    results = data.frame('cellLine' = cellLines,
                         'TP' = NA, 'TN' = NA, 'FP' = NA, 'FN' = NA,
                         "check_sum_genes" = NA, 
                         'accuracy' = NA, 'sensitivity' = NA,
                         'specificity' = NA, 'F1' = NA,
                         'MCC' = NA, 'Penr' = NA, 'logPenr' = NA,
                         'PenrAdj' = NA, 'logPenrAdj' = NA)
    
    for (t in 1:length(cellLines)){
      
      # Genes From the Model #
      modelGenes       <- unique(Table_HumanGEM_Genes$ENSEMBL) 
      modelGenesSYMBOL <- unique(Table_HumanGEM_Genes$SYMBOL) 
      # Essential Genes #
      modelEssential        <- rownames(Binary_Matrix)[Binary_Matrix[,cellLines[t]]>0]
      modelEssentialENSEMBL <- unique(Table_HumanGEM_Genes$ENSEMBL[Table_HumanGEM_Genes$ENSEMBL %in% modelEssential])
      modelEssentialSYMBOL  <- unique(Table_HumanGEM_Genes$SYMBOL[Table_HumanGEM_Genes$ENSEMBL %in% modelEssential])
      # Non-Essential Genes #
      modelNonEssential       <- setdiff(modelGenes, modelEssential)
      modelNonEssentialSYMBOL <- unique(Table_HumanGEM_Genes$SYMBOL[Table_HumanGEM_Genes$ENSEMBL %in% modelNonEssential])
      
      # expGenes has all the genes that are not NA from Achilles #
      expGenes <-  rownames(ExpMatBinary)[!is.na(ExpMatBinary[,cellLines[t]])]
      # expEssential has all the genes that are considered essential #
      expEssential <- expGenes[ExpMatBinary[expGenes,cellLines[t]]==1]
      # expNonEssential has all the genes that are NOT essential # 
      expNonEssential <- setdiff(expGenes,expEssential)
      
      
      results$TP[t] = length(intersect(modelEssentialSYMBOL, expEssential))        # True Positives  #
      results$TN[t] = length(intersect(modelNonEssentialSYMBOL, expNonEssential))  # True Negatives  #
      results$FP[t] = length(intersect(modelEssentialSYMBOL, expNonEssential))     # False Positives #
      results$FN[t] = length(intersect(modelNonEssentialSYMBOL, expEssential))     # False Negatives #   
      
      # Intersect Genes from HumanGEM and Genes from ACHILLES #
      pop <- intersect(modelGenesSYMBOL, expGenes)
      # Essential Genes from gMCST5 #
      sample <- intersect(modelEssentialSYMBOL, expGenes)
      # Essential Genes from ACHILLES #
      successes <- intersect(expEssential, modelGenesSYMBOL)
      
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
    
    results$model <- "HumanGEM"
    results$evalMethod <- evalMethod
    
    results$thMethod <- Threshold_Info
    if(gMCS_order > 7){
      results$Order    <- paste0("HigherThan_Order_7")
    } else {
      results$Order    <- paste0("Order_", gMCS_order)
    }
    
    results$PenrAdj = p.adjust(results$Penr,'BH')
    
    results$logPenr = -log10(results$Penr)
    results$logPenrAdj = -log10(results$PenrAdj)
    
    result_list[[evalMethod]] <- results
    
    results <- do.call(rbind, result_list)
    return(results)})