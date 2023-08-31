CalculateEssentialGenes_gmcsTH_DOM <-  function(gene.exp, 
                                                gMCS_Analysis, 
                                                sample.class, 
                                                gmcsTH_perc = 0.05, 
                                                nWorkers = 4,
                                                gMCS_order = 1,
                                                order1 = TRUE){
  ## FUNCTION TO COMPUTE gMCSTh - HART ###
  # Input #
  #' @param gene.exp, gene expression matrix, in this case, CCLE TPM
  #' @param gMCS.info, information about gMCS
  #' @param sample.class, this is not important right now, will not be touched
  #' @param gmcsTH_perc, percentage to consider highly expressed, default = 0.05
  #' @param nWorkers, to parallelize, default = 4
  #' @param gMCS_order, which is the order of gMCS we are working on, default = 1
    
  # Output - list with the following elements: #
  #' @param gene.ratio # [All Genes, Samples], Ratio of Expression/Expression Threshold
  #' @param num.essential.gene, Subset of Genes from HumanGEM, nomenclature info, percentage of samples affected
  #' @param ratio.essential.gene, Same as before but ratio per Sample Class levels
  #' @param mat.essential.gene, Matrix of [HumanGEM Subset, Samples], how many gMCS predict the gene to be essential in a given Cell Line
  #' @param list.gMCS.essential, Gene table, nomenclature, associated gMCS and percentage of samples in which the gene is essential
  #' @param list.gMCS.essential.mat, List made of essential genes, in which a sparse matrix is in with [gMCS x CellLine]
  
  if(gMCS_order > 7){cat(paste0("Computing for gMCS orders higher than ", gMCS_order-1,"."))
  } else{cat(paste0("Computing for gMCS of order ", gMCS_order,"."))}
  
  # browser()
  threshold_logFC =  1e-3
  
  # Prepare cluster #
  
  cl <- parallel::makeCluster(nWorkers)
  registerDoSNOW(cl)
  
  # Extract info #
  gMCSs.ENSEMBL.mat    <- gMCS_Analysis$gMCSs.ENSEMBL.mat    # Matrix that associates gMCS and genes
  if(order1){init <-  1}else{init=2}
  if(gMCS_order <= 7){
    gMCSs.ENSEMBL.mat    <- sapply(init:gMCS_order, 
                                   function(x) gMCSs.ENSEMBL.mat[which(rowSums(gMCSs.ENSEMBL.mat) == x),])
    gMCSs.ENSEMBL.mat    <- do.call(rbind, gMCSs.ENSEMBL.mat)
  }else{
    if(!order1){
      gMCSs.ENSEMBL.mat <- gMCSs.ENSEMBL.mat[which(rowSums(gMCSs.ENSEMBL.mat) > 2),]
    }
  }
  
  genes.gMCSs.ENSEMBL  <- gMCS_Analysis$genes.gMCSs.ENSEMBL  # List of genes (ENSEMBL)
  gMCSs.ENSEMBL.list   <- gMCS_Analysis$gMCSs.ENSEMBL.list   # gMCS List
  if(gMCS_order <= 7){
    if(order1){
      gMCSs.ENSEMBL.list   <- gMCSs.ENSEMBL.list[which(lengths(gMCSs.ENSEMBL.list) <= gMCS_order)]
    }else{
      gMCSs.ENSEMBL.list   <- gMCSs.ENSEMBL.list[lengths(gMCSs.ENSEMBL.list) <= gMCS_order & 
                                                   lengths(gMCSs.ENSEMBL.list) > 1]
      
    }
  }else{
    if(!order1){
      gMCSs.ENSEMBL.list <- gMCSs.ENSEMBL.list[!lengths(gMCSs.ENSEMBL.list)==1]
      }
  }
  gMCSs.ENSEMBL.txt    <- gMCS_Analysis$gMCSs.ENSEMBL.txt    # Genes per gMCS (text)
  table.genes.HumanGEM <- gMCS_Analysis$table.genes.HumanGEM # Genes - ENSEMBL & SYMBOL
  
  # Calculate first and second most expressed gene in each gMCS in each sample ####
  
  # Iniatialize Variables ####
  gene.exp[is.na(gene.exp)] <- 0 # Eliminate nans
  dim(gene.exp)
  
  nn <- dim(gMCSs.ENSEMBL.mat)[1] # Number of gMCS
  
  gene.first.exp = gene.first.ENSEMBL = # Initialize gMCS per cell line matrices
    gene.second.exp = gene.second.ENSEMBL =
    gene.first.exp.2 = gene.first.ENSEMBL =
    gene.second.exp.2 = gene.second.ENSEMBL =
    matrix(NA,nrow = nn, ncol = dim(gene.exp)[2])
  
  colnames(gene.first.exp) = colnames(gene.first.ENSEMBL) = # Put colnames to all variables
    colnames(gene.second.exp) = colnames(gene.second.ENSEMBL) = 
    colnames(gene.exp)
  
  index_gmcs_length_1 <- as.numeric(which(apply(gMCSs.ENSEMBL.mat,1,sum)==1)) # Store gMCS of length 1
  
  k <- list()
  
  gene.exp.gMCSs <- gene.exp[intersect(rownames(gene.exp), genes.gMCSs.ENSEMBL),] # Save all genes
  rownames(gene.exp.gMCSs) <- intersect(rownames(gene.exp), genes.gMCSs.ENSEMBL)  # Name the rows
  gene.exp.gMCSs[is.na(gene.exp.gMCSs)] <- 0       # NAs are 0 now  
  aux_zeros_mat <- matrix(0, nrow = 3, 
                          ncol = dim(gene.exp)[2]) # Make auxiliary matrix to get gMCS of different lengths
  
  cat(paste("Processing the", nn, "gMCSs for", ncol(gene.exp), "samples!"))
  
  # system.time({
  
  # Obtain Most Expressed Genes ####
  obtain_1_2_3_genes_gMCS <- function(i = 1){
    gmcs <- gMCSs.ENSEMBL.list[[i]]
    if (length(gmcs)==1){
      gMCSs.ENSEMBL.exp <- rbind(matrix(gene.exp[gmcs,], nrow = 1), aux_zeros_mat)
      gMCSs.ENSEMBL.exp <- apply(gMCSs.ENSEMBL.exp,2,unlist)
      colnames(gMCSs.ENSEMBL.exp) <- colnames(gene.exp)
    } else {
      gMCSs.ENSEMBL.exp <- rbind(as.matrix(gene.exp[gmcs,]), aux_zeros_mat)
    }
    gmcs <- c(gmcs, rep("", nrow(aux_zeros_mat)))
    rownames(gMCSs.ENSEMBL.exp) <- gmcs
    
    # Save expression of most expressed genes
    # Which one is the most expressed gene in the set per Cell Line?
    aux <- apply(gMCSs.ENSEMBL.exp,2,function(x){sort(x, decreasing = T)[1:2]})
    # Save which gene is the most expressed
    # Which one is the most expressed gene in the set per Cell Line?
    aux2 <- apply(gMCSs.ENSEMBL.exp,2,function(x){order(x, decreasing = T)[1:2]})
    
    return(list(first.exp = aux[1,],
                second.exp = aux[2,],
                first.ENSEMBL = gmcs[aux2[1,]],
                second.ENSEMBL = gmcs[aux2[2,]]))
    
  }
  if(order1){
    results <- foreach(i=1:nn) %dopar% obtain_1_2_3_genes_gMCS(i)
  }else{
    results <- foreach(i=1:nn) %dopar% obtain_1_2_3_genes_gMCS(i)
  }
  
  # Expression for the MOST expressed genes in all gMCS across all samples #
  gene.first.exp <- do.call(rbind,lapply(results, function(x){matrix(x$first.exp, nrow = 1)}))
  # Expression for the SECOND MOST expressed genes in all gMCS across all samples #
  gene.second.exp <- do.call(rbind,lapply(results, function(x){matrix(x$second.exp, nrow = 1)}))
  
  # Gene Names for the MOST expressed genes in all gMCS across all samples #
  gene.first.ENSEMBL <- do.call(rbind,lapply(results, function(x){matrix(x$first.ENSEMBL, nrow = 1)}))
  # Gene Names for the SECOND MOST expressed genes in all gMCS across all samples #
  gene.second.ENSEMBL <- do.call(rbind,lapply(results, function(x){matrix(x$second.ENSEMBL, nrow = 1)}))
  
  # Rename columns #
  colnames(gene.first.exp) = colnames(gene.first.ENSEMBL) = 
    colnames(gene.second.exp) = colnames(gene.second.ENSEMBL) = 
    colnames(gene.exp)
  
  # Remove Not Expressed Genes #
  
  for (sample in 1:dim(gene.exp)[2]) {
    idx_not_expressed <- which(gene.first.exp[,sample]<1e-3)
    gene.first.exp[idx_not_expressed,sample] <- 0
    gene.first.ENSEMBL[idx_not_expressed,sample] <- ""
    
    idx_not_expressed <- which(gene.second.exp[,sample]<1e-3)
    gene.second.exp[idx_not_expressed,sample] <- 0
    gene.second.ENSEMBL[idx_not_expressed,sample] <- ""
    
    # Filter for gMCS with length 1 #
    gene.second.exp[index_gmcs_length_1,sample] <- 0
    gene.second.ENSEMBL[index_gmcs_length_1,sample] <- ""
  }
  
  gene.second.exp[index_gmcs_length_1,] <- 0
  gene.second.ENSEMBL[index_gmcs_length_1,] <- ""
  
  # Calculate first and second most expressed gene RATIO in each gMCSs in each sample ####
  #browser()
  
  # system.time({
  ratio_threshold_list <- list()
  
  ratio_threshold <- rep(NA,dim(gene.first.exp)[2])
  names(ratio_threshold) <- colnames(gene.first.exp)
  
  # Compute Ratio Thresholds #
  
  for (i in 1:dim(gene.first.exp)[2]){
    
    aux <- as.list(gene.first.exp[,i])
    
    for (f in 1:length(aux)){
      names(aux[[f]]) <- gene.first.ENSEMBL[f,i]
    }
    
    names(aux) <- gene.first.ENSEMBL[,i]
    aux <- unlist(unique(aux))
    ratio_threshold[i] <- quantile(aux[aux>0], gmcsTH_perc) # Ratio per Cell Line
  }
  
  gene.first.ratio <- gene.first.exp
  gene.second.ratio <- gene.second.exp
  gene.ratio <- gene.exp
  
  # If the value i in sample j in the matrix gene.ratio > 1 it means that
  # the gene is expressed more than the threshold, in the other hand,
  # if the value i in sample j has a value of < 1, it implies that is not
  # properly expressed and should be discarded.
  
  if (dim(gene.first.exp)[1]>0){
    for (i in 1:dim(gene.first.ratio)[2]){
      
      gene.first.ratio[,i] <- gene.first.exp[,i]/ratio_threshold[i]
      gene.second.ratio[,i] <- gene.second.exp[,i]/ratio_threshold[i]
      gene.ratio[,i] <- gene.exp[,i]/ratio_threshold[i]
    }
  }
  
  # Calculate essential genes by sample ####
  # browser()
  gene.first.log2ratio <- log2(gene.first.ratio)
  gene.second.log2ratio <- log2(gene.second.ratio)
  
  genes.ENSEMBL.essential <- gMCSs.ENSEMBL.txt[!grepl('--', gMCSs.ENSEMBL.txt)]
  
  ################################################################
  ####    if the gMCS is only one gene, make it essential     ####
  ################################################################
  
  gene.first.log2ratio[rowSums(gMCSs.ENSEMBL.mat)==1,] <- (+1)
  gene.second.log2ratio[rowSums(gMCSs.ENSEMBL.mat)==1,] <- (-1)
  
  ################################################################
  ####                Single essential genes                  ####
  ################################################################
  
  # Enumerate all the genes in the gMCS #
  genes.gMCSs.ENSEMBL.12 <- unique(c(unique(gene.first.ENSEMBL), unique(gene.second.ENSEMBL)))
  genes.gMCSs.ENSEMBL.12 <- genes.gMCSs.ENSEMBL.12[genes.gMCSs.ENSEMBL.12!=""]
  length(genes.gMCSs.ENSEMBL.12)
  
  # browser() 
  # Essential Genes per Sample #
  essential.single.gMCS <- gene.first.log2ratio > (+threshold_logFC) & gene.second.log2ratio < (-threshold_logFC)
  
  # Initialize Number of Essential Gene per Cell Line Matrix #
  num.essential.gene <- as.data.frame(matrix(0,nrow = length(genes.gMCSs.ENSEMBL),
                                             ncol = length(levels(sample.class))+2))
  colnames(num.essential.gene) <- c("gen","num.gMCS",levels(sample.class))
  rownames(num.essential.gene) <- genes.gMCSs.ENSEMBL
  num.essential.gene$gen <- genes.gMCSs.ENSEMBL
  
  # Initialize Genes per Cell Line Matrix #
  mat.essential.gene <- matrix(0,nrow = length(genes.gMCSs.ENSEMBL),
                               ncol = dim(gene.exp)[2])
  colnames(mat.essential.gene) <- colnames(gene.exp)
  rownames(mat.essential.gene) <- genes.gMCSs.ENSEMBL
  
  # Define function to obtain essential genes ####
  
  obtain_essential_single_genes_gMCS <- function(i){
    gene <- genes.gMCSs.ENSEMBL[i]
    return(colSums(gene.first.ENSEMBL == gene & essential.single.gMCS))
  }
  
  obtain_essential_single_gMCS_simple_class <- function(i){
    gene <- genes.gMCSs.ENSEMBL.12[i] # Take gene Information
    aux1 <- gene.first.ENSEMBL == gene & essential.single.gMCS # Where is the gene Essential?
    aux2 <- apply(aux1*1, 1, sum) # Make a numeric vector of cell lines affected per gMCS
    aux3 <- aux2/ncol(aux1) # Make it range from 0 to 1
    aux4 <- order(aux3, decreasing = T) # Order it from most affected to least
    aux3 <- matrix(aux3, ncol = 1) 
    aux5 <- data.frame(gene = gene,  # Make dataframe with important first
                       task = "combined", gMCS = aux4,aux3[aux4])
    colnames(aux5)[-c(1:3)] <- levels(sample.class)
    aux5 <- aux5[aux5[,4]>0,]
    
    return(aux5)
  }
  
  obtain_essential_single_gMCS_multiple_class <- function(i){
    # This function is similar to last one, but taking into account sample types
    gene = genes.gMCSs.ENSEMBL.12[i]
    aux1 <- gene.first.ENSEMBL == gene & essential.single.gMCS
    
    if (sum(aux1)==0){
      return(NA)
    } else {
      
      aux2 <- t(rowsum(t((aux1)*1), sample.class))
      rownames(aux2) <- 1:dim(aux2)[1]
      aux3 <- aux2[order(rowSums(aux2),decreasing = T),]
      if (sum(rowSums(aux3)>0)>1){
        aux3 <- aux3[rowSums(aux3)>0,]
      } else {
        aux3 <- aux3[1:4,]
      }
      
      aux5 <- data.frame(gen = gen, task = "combined", gMCS = rownames(aux3),aux3)
      aux5 <- aux5[rowSums(aux5[,-c(1:3)])>0,]
      
      for (col in colnames(aux3)){
        aux5[,col] <- aux5[,col]/sum(sample.class==col)
      }   
      aux5 <- aux5[order(apply(aux5[,-c(1:3)],1,mean),decreasing = T),]
      
      return(aux5)
    }
  }
  
  nn <- length(genes.gMCSs.ENSEMBL)
  
  # Run Essential Genes #  
  results <- foreach(i=1:nn) %dopar% obtain_essential_single_genes_gMCS(i) 
  
  # Run essential gMCS for each gen
  if (length(levels(sample.class))>1){
    list.gMCS.essential <- foreach(i=1:nn) %dopar% obtain_essential_single_gMCS_multiple_class(i)
  } else {
    list.gMCS.essential <- foreach(i=1:nn) %dopar% obtain_essential_single_gMCS_simple_class(i) 
  }
  
  mat.essential.gene <- do.call(rbind,results)
  rownames(mat.essential.gene) <- genes.gMCSs.ENSEMBL
  
  list.gMCS.essential <- do.call(rbind, list.gMCS.essential)
  list.gMCS.essential <- na.omit(list.gMCS.essential)
  
  # Compute Ratio, divide number of samples in which the gene n1 
  # is expressed by all the samples in the analysis #
  num.essential.gene$gen <- genes.gMCSs.ENSEMBL
  num.essential.gene$num.gMCS <- apply(gMCSs.ENSEMBL.mat[,genes.gMCSs.ENSEMBL],2,sum)
  if (length(levels(sample.class))>1) {
    num.essential.gene[,levels(sample.class)] <- t(apply(mat.essential.gene, 1,
                                                         function(x){rowsum((x>0)*1, sample.class)}))
  } else {
    num.essential.gene[,levels(sample.class)] <- apply(mat.essential.gene>0, 1, sum)
  }
  
  ratio.essential.gene <- num.essential.gene
  for (col in levels(sample.class)){
    ratio.essential.gene[,col] <- ratio.essential.gene[,col]/sum(sample.class==col)
  }
  
  
  # Set essential genes (gMCS of only 1 gene) #
  for (col in levels(sample.class)){
    num.essential.gene[genes.ENSEMBL.essential, col] <- sum(sample.class==col)
    ratio.essential.gene[genes.ENSEMBL.essential, col] <- 1
  }
  mat.essential.gene[genes.ENSEMBL.essential, ] <- 1
  
  
  # Save Essential gMCSs by gene by sample (aux info)
  # library(doSNOW)
  # cl <- makeCluster(nWorkers); registerDoSNOW(cl)
  
  cat(paste("Calculating essential gMCS for each gene for ",ncol(gene.exp), "samples (part2):"))
  pb <- txtProgressBar(max = sum(rowSums(mat.essential.gene)>0), style = 3)
  progress <- function(i) {setTxtProgressBar(pb, i)}
  opts <- list(progress = progress)
  list.gMCS.essential.mat <- function(g){Matrix::Matrix(gene.first.ENSEMBL==g & essential.single.gMCS, sparse = T)}
  # list.gMCS.essential.mat <- foreach(g = genes.gMCSs.ENSEMBL[rowSums(mat.essential.gene)>0],
  list.gMCS.essential.mat <- foreach(g = genes.gMCSs.ENSEMBL,
                                     .options.snow = opts) %dopar% list.gMCS.essential.mat(g) 
  names(list.gMCS.essential.mat) <- genes.gMCSs.ENSEMBL
  # names(list.gMCS.essential.mat) <- genes.gMCSs.ENSEMBL[rowSums(mat.essential.gene)>0]
  for(EssentialGene in genes.ENSEMBL.essential){
    list.gMCS.essential.mat[[EssentialGene]][which(unlist(gMCS_Analysis$gMCSs.ENSEMBL.list) == EssentialGene)[1],] <- rep(TRUE,5)
  }
  
  # stopCluster(cl = cl)
  
  
  # close cluster
  stopCluster(cl = cl)
  
  # browser()
  # Add the gene information (symbol, entrez, isEssential)
  colnames(num.essential.gene)[colnames(num.essential.gene)=="gen"] <- "ENSEMBL"
  colnames(ratio.essential.gene)[colnames(ratio.essential.gene)=="gen"] <- "ENSEMBL"
  colnames(list.gMCS.essential)[colnames(list.gMCS.essential)=="gen"] <- "ENSEMBL"
  
  num.essential.gene   <- merge(table.genes.HumanGEM, num.essential.gene)
  ratio.essential.gene <- merge(table.genes.HumanGEM, ratio.essential.gene)
  list.gMCS.essential  <- merge(table.genes.HumanGEM, list.gMCS.essential)
  
  # Export the results 
  ResultsEssentiality <- list(gene.ratio = gene.ratio,                     # [All Genes, Samples] Ratio of Expression/Expression Threshold
                              num.essential.gene   = num.essential.gene,   # Subset of Genes from HumanGEM, nomenclature info, percentage of samples affected
                              ratio.essential.gene = ratio.essential.gene, # Same as before but ratio per Sample Class levels
                              mat.essential.gene   = mat.essential.gene,   # Matrix of [HumanGEM Subset, Samples], how many gMCS predict the gene to be essential in a given Cell Line
                              list.gMCS.essential  = list.gMCS.essential,  # Gene table, nomenclature, associated gMCS and percentage of samples in which the gene is essential
                              list.gMCS.essential.mat = list.gMCS.essential.mat) # List made of essential genes, in which a sparse matrix is in with [gMCS x CellLine]
  
  
  
}
