ComputeBinaryMatrix <- function(gMCSAnalysis, gMCSList, Cumulative = TRUE,
                                max_order = 8, HigherOrders = TRUE, isNew = FALSE){
  
  # Input #
  #' @param gMCSAnalysis, list regarding the gMCS Analysis, usually named RoL_[Threshold_Name]
  #' @param gMCSList, list regarding core gMCS structural information.
  #' @param Cumulative, to compute cumulative or not cumulatively, useful for length independent analyses.
  #' @param max_order, max order of gMCSs to take into account.
  #' @param HigherOrders, should higher thresholds should be considered?
  #' @param isNew, to check whether the binary matrix has already been computed
  
  # Output: #
  #' @param Binary_Matrix_Length_List, A list of binary matrices for each length marked in max order.
  
  if(isNew){
    List_gMCS_Essential_Matrix <- gMCSAnalysis$list.gene.of.gMCS.by.sample
  } else {
    List_gMCS_Essential_Matrix <- gMCSAnalysis$list.gMCS.essential.mat
  }
  Binary_Matrix_Length_List <- vector("list", max_order)
  
  if(Cumulative){
    for (gMCS_Length in 1:max_order){
      Binary_Matrix <- matrix(0, nrow = nrow(gMCSAnalysis$mat.essential.gene),
                              ncol = ncol(gMCSAnalysis$mat.essential.gene),
                              dimnames = list(rownames(gMCSAnalysis$mat.essential.gene),
                                              colnames(gMCSAnalysis$mat.essential.gene)))
      
      for(Essential_Gene in names(List_gMCS_Essential_Matrix)){
        if(gMCS_Length == max_order){
          Binary_Matrix[Essential_Gene,] <- as.numeric(colSums(List_gMCS_Essential_Matrix[[Essential_Gene]]) > 0)
        } else {
          gMCS_to_Take <- which(lengths(gMCSList) <= gMCS_Length)
          if(length(gMCS_to_Take)!=1){
            Binary_Matrix[Essential_Gene,] <- as.numeric(colSums(List_gMCS_Essential_Matrix[[Essential_Gene]][gMCS_to_Take,]) > 0)
            
          }else{
            Binary_Matrix[Essential_Gene,] <- as.numeric((List_gMCS_Essential_Matrix[[Essential_Gene]][gMCS_to_Take,]) > 0)
          }
          
        }
        
      }
      Binary_Matrix_Length_List[[gMCS_Length]] <- Binary_Matrix
    }
  } else { # Not Cumulative Case #
    for (gMCS_Length in 1:max_order){
      Binary_Matrix <- matrix(0, nrow = nrow(gMCSAnalysis$mat.essential.gene),
                              ncol = ncol(gMCSAnalysis$mat.essential.gene),
                              dimnames = list(rownames(gMCSAnalysis$mat.essential.gene),
                                              colnames(gMCSAnalysis$mat.essential.gene)))
      
      if(gMCS_Length == max_order){
        if(!HigherOrders){ # 
          gMCS_to_Take <- which(lengths(gMCSList) == gMCS_Length)
        } else { gMCS_to_Take <- which(lengths(gMCSList) >= gMCS_Length)}
      }
      else { 
        gMCS_to_Take <- which(lengths(gMCSList) == gMCS_Length)
      }
      
      for(Essential_Gene in names(List_gMCS_Essential_Matrix)){
        Binary_Matrix[Essential_Gene,] <- as.numeric(colSums(List_gMCS_Essential_Matrix[[Essential_Gene]][gMCS_to_Take,]) > 0)
      }
      
      Binary_Matrix_Length_List[[gMCS_Length]] <- Binary_Matrix
    }
  }
  return(Binary_Matrix_Length_List)
}