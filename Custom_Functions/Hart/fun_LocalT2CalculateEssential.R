LocalT2CalculateEssential <- function(gene.exp, gMCS.info,
                                      sample.class, sample.cohort,
                                      gMCS_order = 8,
                                      localT2_mode, nWorkers){
  
  # browser()
  ## FUNCTION TO COMPUTE LocalT2 - Hart ###
  # Input #
  #' @param gene.exp, gene expression matrix, in this case, CCLE TPM
  #' @param gMCS.info, information about gMCS
  #' @param sample.class, this is not important right now, will not be touched
  #' @param sample.cohort, needed for localT2 assignments.
  #' @param localT2_mode, mode of localT2: only used genes (~1244) or all the Human1 universe?
  #' @param nWorkers, to parallelize, default = 4
  #' @param gMCS_order, which is the order of gMCS we are working on, default = 8
  
  # Output - list with the following elements: #
  #' @param localT2, value for each gene.
  #' @param gene.ON.OFF, data frame of genes that are ON or OFF.
  #' @param num.essential.gene, Subset of Genes from HumanGEM, nomenclature info, percentage of samples affected
  #' @param ratio.essential.gene, Same as before but ratio per Sample Class levels
  #' @param mat.essential.gene, Matrix of [HumanGEM Subset, Samples], how many gMCS predict the gene to be essential in a given Cell Line
  #' @param list.gMCS.essential, Gene table, nomenclature, associated gMCS and percentage of samples in which the gene is essential
  #' @param list.gMCS.essential.mat, List made of essential genes, in which a sparse matrix is in with [gMCS x CellLine]
  
  
  if(gMCS_order > 7){print(paste0("Computing for gMCS orders higher than ", gMCS_order-1,"."))
  } else{print(paste0("Computing for gMCS of order ", gMCS_order,"."))}
  
  # prepare cluster
  cl <- parallel::makeCluster(nWorkers)
  registerDoSNOW(cl)
  
  
  # extract info
  gMCSs.ENSEMBL.mat <- gMCS.info$gMCSs.ENSEMBL.mat # Matrix that associates gMCS and genes
  if(gMCS_order <= 7){
    gMCSs.ENSEMBL.mat    <- sapply(1:gMCS_order, 
                                   function(x) gMCSs.ENSEMBL.mat[which(rowSums(gMCSs.ENSEMBL.mat) == x),])
    gMCSs.ENSEMBL.mat    <- do.call(rbind, gMCSs.ENSEMBL.mat)
  }
  genes.gMCSs.ENSEMBL  <- gMCS.info$genes.gMCSs.ENSEMBL # List of genes (ENSEMBL)
  gMCSs.ENSEMBL.txt    <- gMCS.info$gMCSs.ENSEMBL.txt
  table.genes.HumanGEM <- gMCS.info$table.genes.HumanGEM
  
  # Calculate localT2 thresholds  ####
  # All samples are considered to be from the same cohort. If not, please update in different files
  
  print("Perform the localT2 thesholding")
  
  localT2 <- list()
  
  # withProgress(message = 'Perform the localT2 thesholding',{
  for (cm in levels(sample.cohort)){
    
    print(paste(cm))
    if (localT2_mode=="all_genes_gMCSs"){
      gene.cohort.TPM <- gene.exp[genes.gMCSs.ENSEMBL,, sample.cohort==cm]
    } else {
      gene.cohort.TPM <- gene.exp[table.genes.HumanGEM$ENSEMBL, sample.cohort==cm]
    }
    
    localT2[[cm]]         <- list()
    localT2[[cm]]$OFF     <- quantile(gene.cohort.TPM, 0.25, na.rm = T)
    localT2[[cm]]$ON      <- quantile(gene.cohort.TPM, 0.75, na.rm = T)
    localT2[[cm]]$MAYBEON <- apply(gene.cohort.TPM,1,mean)
    
    print(localT2[[cm]][c(1,2)])
    print(localT2[[cm]]$MAYBEON[1:5])
  }
  # })
  
  
  # Define genes in cohort as ON or OFF ####
  
  if (localT2_mode=="all_genes_gMCSs"){
    gene.ON.OFF <- gene.exp[genes.gMCSs.ENSEMBL,]*0
  } else {
    gene.ON.OFF <- gene.exp[table.genes.HumanGEM$ENSEMBL,]*0
  }
  
  # withProgress(message = 'Define genes in cohort as ON or OFF',{
  for (cm in levels(sample.cohort))
  {
    # n <- 1
    print(paste(cm))
    
    if (localT2_mode=="all_genes_gMCSs"){
      gene.aux.TPM <- gene.exp[genes.gMCSs.ENSEMBL, sample.cohort==cm]
    } else {
      gene.aux.TPM <- gene.exp[table.genes.HumanGEM$ENSEMBL, sample.cohort==cm]
    }
    gene.aux.ON.OFF <-  gene.aux.TPM*0
    
    gene.aux.TH.MAYBEON <- gene.aux.TPM*0
    for (i in 1:dim(gene.aux.TH.MAYBEON)[2])
    { gene.aux.TH.MAYBEON[,i] <- localT2[[cm]]$MAYBEON[rownames(gene.aux.ON.OFF)] }
    
    gene.aux.ON.OFF <- (gene.aux.TPM >= gene.aux.TH.MAYBEON) 
    gene.aux.ON.OFF[gene.aux.TPM > localT2[[cm]]$ON] <- TRUE
    gene.aux.ON.OFF[gene.aux.TPM < localT2[[cm]]$OFF] <- FALSE
    
    gene.ON.OFF[, sample.cohort==cm] <-  gene.aux.ON.OFF*1
    rm(list = c("gene.aux.TPM", "gene.aux.ON.OFF"))
  }
  # })
  
  # table(gene.ON.OFF)
  dim(gene.ON.OFF)
  sum(is.na(gene.ON.OFF))
  
  
  
  # Define first and second gene as ON or OFF ####
  
  
  calculateNumberONgenes <- function(i){
    # i <- 3
    require(Matrix)
    # set all OFF genes as 0 and ON genes as 1
    gMCSs.ENSEMBL.exp <- gMCSs.ENSEMBL.mat
    
    aux <- rownames(gene.ON.OFF)[gene.ON.OFF[,i]>0.5]
    gMCSs.ENSEMBL.exp <- gMCSs.ENSEMBL.exp[,colnames(gMCSs.ENSEMBL.exp) %in% aux]
    
    gene.num.ON.aux <- rowSums(gMCSs.ENSEMBL.exp)
    
    gene.ENSEMBL.ON.aux <- rep("", nrow(gMCSs.ENSEMBL.exp))
    
    if (sum(gMCSs.ENSEMBL.exp)>0){ 
      aux <- as.data.frame(which(gMCSs.ENSEMBL.exp>0, arr.ind=T))
      aux$col <- colnames(gMCSs.ENSEMBL.exp)[aux$col]
      aux <- aggregate(aux$col, by = list(aux$row), paste ,collapse = ",")
      
      gene.ENSEMBL.ON.aux[aux$Group.1] <- aux$x
    }
    
    return(list(gene.num.ON.aux, gene.ENSEMBL.ON.aux))
  }
  
  nn <- ncol(gene.exp)
  
  print(paste("Define genes in cohort as ON or OFF"))
  
  results <- list() 
  for(i in 1:nn) { 
    print(paste('localT2:',i,'/',ncol(gene.exp)))
    results[[i]] <- calculateNumberONgenes(i) 
  }
  
  
  gene.num.ON <- lapply(results, function(x){x[[1]]})
  gene.num.ON <- matrix(unlist(gene.num.ON), ncol = nn)
  
  gene.ENSEMBL.ON <- lapply(results, function(x){x[[2]]})
  gene.ENSEMBL.ON <- matrix(unlist(gene.ENSEMBL.ON), ncol = nn)
  
  colnames(gene.num.ON) = colnames(gene.ENSEMBL.ON) = 
    colnames(gene.exp)
  
  # debug
  # set.seed(28)
  idx.aux <- sample(ncol(gene.num.ON),5)
  print(head(gene.num.ON[,idx.aux]))
  print(head(gene.ENSEMBL.ON[,idx.aux]))
  
  
  
  # Calculate essential genes by sample ####
  
  genes.ENSEMBL.essential <- gMCSs.ENSEMBL.txt[!grepl('--', gMCSs.ENSEMBL.txt)]
  
  
  # single essential genes ##
  
  essential.single.gMCS <- gene.num.ON==1
  essential.single.gMCS.table <- as.data.frame(t(rowsum(as.matrix(t(essential.single.gMCS)*1), group =  t(sample.class))))
  gene.ENSEMBL.ON.only1gen = gene.ENSEMBL.ON
  gene.ENSEMBL.ON.only1gen[ gene.num.ON!=1 ] <- ""
  
  num.essential.gene <- as.data.frame(matrix(0,nrow = length(genes.gMCSs.ENSEMBL),ncol = length(levels(sample.class))+2))
  colnames(num.essential.gene) <- c("gen","num.gMCS",levels(sample.class))
  rownames(num.essential.gene) <- genes.gMCSs.ENSEMBL
  num.essential.gene$gen <- genes.gMCSs.ENSEMBL
  
  mat.essential.gene <- matrix(0,nrow = length(genes.gMCSs.ENSEMBL),ncol = dim(gene.exp)[2])
  colnames(mat.essential.gene) <- colnames(gene.exp)
  rownames(mat.essential.gene) <- genes.gMCSs.ENSEMBL
  
  
  sum(essential.single.gMCS)
  
  # define function to obtain essential genes ##
  
  obtain_essential_single_genes_gMCS <- function(i){
    gen = genes.gMCSs.ENSEMBL[i]
    return(colSums(gene.ENSEMBL.ON.only1gen==gen & essential.single.gMCS))
  }
  
  obtain_essential_single_gMCS_simple_class <- function(i)
  {
    gen = genes.gMCSs.ENSEMBL[i]
    # gen = "ENSG00000136810"
    aux1 <- gene.ENSEMBL.ON.only1gen==gen & essential.single.gMCS
    
    aux2 <- apply(aux1*1, 1, sum)
    aux3 <- aux2/ncol(aux1)
    aux4 <- order(aux3, decreasing = T)
    aux3 <- matrix(aux3, ncol = 1)
    aux5 <- data.frame(gen = gen, task = "combined", gMCS = aux4,aux3[aux4])
    colnames(aux5)[-c(1:3)] <- levels(sample.class)
    aux5 <- aux5[aux5[,4]>0,]
    
    return(aux5)
  }
  
  obtain_essential_single_gMCS_multiple_class <- function(i)
  {
    gen = genes.gMCSs.ENSEMBL[i]
    aux1 <- gene.ENSEMBL.ON==gen & essential.single.gMCS
    
    if (sum(aux1)==0){
      return(NA)
    } else {
      
      aux2 <- t(rowsum(t((aux1)*1), t(sample.class)))
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
  
  # run essential genes
  print(paste("Calculating essential genes for",ncol(gene.exp), "samples:"))
  results <- foreach(i=1:nn) %dopar% obtain_essential_single_genes_gMCS(i) 
  
  #run essential gMCS for each gen
  print(paste("Calculating essential gMCS for each gene for",ncol(gene.exp), "samples:"))
  if (length(levels(sample.class))>1){
    list.gMCS.essential <- foreach(i=1:nn) %dopar% obtain_essential_single_gMCS_multiple_class(i)
  } else {
    list.gMCS.essential <- foreach(i=1:nn) %dopar% obtain_essential_single_gMCS_simple_class(i) 
  }
  
  
  
  mat.essential.gene <- do.call(rbind,results)
  rownames(mat.essential.gene) <- genes.gMCSs.ENSEMBL
  
  
  list.gMCS.essential <- do.call(rbind, list.gMCS.essential)
  list.gMCS.essential <- na.omit(list.gMCS.essential)
  
  
  num.essential.gene$gen <- genes.gMCSs.ENSEMBL
  num.essential.gene$num.gMCS <- apply(gMCSs.ENSEMBL.mat[,genes.gMCSs.ENSEMBL],2,sum)
  if (length(levels(sample.class))>1) {
    num.essential.gene[,levels(sample.class)] <- t(apply(mat.essential.gene, 1, function(x){rowsum((x>0)*1, t(sample.class))}))
  } else {
    num.essential.gene[,levels(sample.class)] <- apply(mat.essential.gene>0, 1, sum)
  }
  
  
  ratio.essential.gene <- num.essential.gene
  for (col in levels(sample.class))
  {
    ratio.essential.gene[,col] <- ratio.essential.gene[,col]/sum(sample.class==col)
  }
  
  ## Set essential genes (gMCS of only 1 gene)
  for (col in levels(sample.class))
  {
    num.essential.gene[genes.ENSEMBL.essential, col] <- sum(sample.class==col)
    ratio.essential.gene[genes.ENSEMBL.essential, col] <- 1
  }
  mat.essential.gene[genes.ENSEMBL.essential, ] <- 1
  
  list.gMCS.essential.mat <- function(g){Matrix::Matrix(gene.ENSEMBL.ON==g & essential.single.gMCS, sparse = T)}
  list.gMCS.essential.mat <- foreach(g = genes.gMCSs.ENSEMBL) %dopar% list.gMCS.essential.mat(g)
  names(list.gMCS.essential.mat) <- genes.gMCSs.ENSEMBL
  for(EssentialGene in genes.ENSEMBL.essential){
    list.gMCS.essential.mat[[EssentialGene]][which(unlist(gMCS.info$gMCSs.ENSEMBL.list) == EssentialGene)[1],] <- rep(TRUE,5)
  }
  
  # close cluster
  stopCluster(cl = cl)
  
  
  # Add the gene information (symbol, entrez, isEssential)
  colnames(num.essential.gene)[colnames(num.essential.gene)=="gen"] <- "ENSEMBL"
  colnames(ratio.essential.gene)[colnames(ratio.essential.gene)=="gen"] <- "ENSEMBL"
  colnames(list.gMCS.essential)[colnames(list.gMCS.essential)=="gen"] <- "ENSEMBL"
  
  num.essential.gene <- merge(table.genes.HumanGEM, num.essential.gene)
  ratio.essential.gene <- merge(table.genes.HumanGEM, ratio.essential.gene)
  list.gMCS.essential <- merge(table.genes.HumanGEM, list.gMCS.essential)
  
  # Export the results 
  return(list(localT2 = localT2,
              #gene.num.ON = gene.num.ON,
              #gene.ENSEMBL.ON = gene.ENSEMBL.ON,
              gene.ON.OFF = gene.ON.OFF,
              num.essential.gene = num.essential.gene,
              ratio.essential.gene = ratio.essential.gene,
              mat.essential.gene = mat.essential.gene,
              list.gMCS.essential = list.gMCS.essential,
              list.gMCS.essential.mat = list.gMCS.essential.mat))
  
}