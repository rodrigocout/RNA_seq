library("ggfortify")
library("WGCNA")

options(stringsAsFactors = FALSE)
memory.limit()
enableWGCNAThreads()

# --------------------------------------------------------
# Dataset is split based on tissue, using 'OA' for lesioned in sample names
# and 'P' for preserved
# --------------------------------------------------------
splitTissue <- function(datExpr){
  datExprPreserved <- datExpr[grep("P", rownames(datExpr)),]
  datExprOsteo <- datExpr[grep("OA", rownames(datExpr)),]
  print(dim(datExprOsteo))
  tissues <- list("Preserved" = datExprPreserved, "Lesioned" = datExprOsteo)
  return(tissues)
}

# --------------------------------------------------------
# Reads of an average of >= 5 are kept, others are disregarded
# --------------------------------------------------------
filterReads <- function(datExpr){
  vector <- c()
  for (x in 1:(ncol(datExpr)-1)){
    meanReads <- (mean(datExpr[,x]))
    if (meanReads >= 5){
      vector <- c(vector,x)
    }
  }
  datExpr <- data.frame(datExpr[,vector])
  datExpr <- t(datExpr)
  return(datExpr)
}

# --------------------------------------------------------
# Final network containing nodes, edges and weight are exported
# --------------------------------------------------------
exportNetwork <- function(datExpr){
  adj <- adjacency(datExpr, power = 4, type='signed') # Change power according to < 0.8 but closest to it
  TOM <- TOMsimilarity(adj, TOMType = 'signed') # Signed is preferred
  #save(TOM, file="TOM.RData")
  dimnames(TOM) <- list(names(data.frame(datExpr)), names(data.frame(datExpr)))
  vis <- exportNetworkToVisANT(TOM,file = paste("Network.txt", sep=""), weighted = TRUE,threshold = 0)
}

# --------------------------------------------------------
# Co-expression network of lncRNA and mRNA is created 
# --------------------------------------------------------
main <- function(){
  VstNormNoBatch <- data.frame(VstNormNoBatch)
  remove <- which(rownames(VstNormNoBatch) %in% c('__no_feature', '__not_aligned', '__ambiguous', '__too_low_aQual', '__alignment_not_unique')) 
  VstNormNoBatch <- VstNormNoBatch[-c(remove),]
  dataset <- t(VstNormNoBatch)
  
  dataset <- filterReads(dataset)

  dataset2 <- dataset[, which(colnames(dataset) %in% colnames(datExprLnc))]
  datExprLnc2 <- datExprLnc[, which(colnames(datExprLnc) %in% colnames(dataset2))]
  datExpr <- rbind(data.frame(dataset2), datExprLnc2)


  tissues <- splitTissue(t(datExpr))
  dataset <- tissues$Lesioned

  exportNetwork(dataset)
}

main()
