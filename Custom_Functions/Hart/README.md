# Custom Functions - Hart ##
In this folder the custom functions used to compute Hart's essentiality meta-analysis are kept.

 ## fun_CalculateHartEssentialGenes_gmcsTH_DOM:
The main code to compute the essential genes based on a given threshold. This function is the first step of the pipeline. It needs the gene expression data, the gMCS information, a data.frame with the information regarding the sample classes, the threshold percentage, the number of workers to compute in parallel, the maximum order (length) considered and whether order 1 genes are considered (TRUE by default).

 ## fun_LocalT2CalculateEssential:
 Same code structure as before, but tailored to localT2 thresholding technique.

 ## fun_TweakedHartgMCS_Results:
 This function is used to get the final essentiality results of the pipeline.