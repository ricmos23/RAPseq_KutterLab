#' Collection of functions to assist in analyses
#' @author carlga ~ Carlos Gallardo
#'


#' Scales data row-wise using selected method
#' @param x A dataframe or matrix
#' @param method Method used for scaling
#' @return Data frame with scaled data
#' 
#' @author carlga
#'
scaleData <- function(x, method = NULL) {
  
  if(is.null(method)) stop('Please, indicate a valid method\n')
  
  if(method == "zscore") {
    row.mean <- apply(x, 1, mean)
    row.sd <- apply(x, 1, sd)
    res <- (x - row.mean) / row.sd
  }
  else if(method == "minmax") {
    row.min <- apply(x, 1, min)
    row.max <- apply(x, 1, max)
    res <- (x - row.min) / (row.max - row.min)
  }
  else if(method == "max") {
    row.max <- apply(x, 1, max)
    res <- x / row.max
  }
  else stop('This method is not available\n')
  
  return(as.data.frame(res))
}

#' Transforms data according to column-wise grouping using a function (e.g. mean, median...)
#' @param x A dataframe or matrix
#' @param group.lbls Vector with grouping of columns
#' @param FUN Function to be used for data transformation (e.g. function(x) apply(x,1,mean))
#' @return Data frame with transformed data
#' 
#' @author carlga
#'
groupTransform <- function(x, group.lbls, FUN) {
  
  group.lbls.uniq <- unique(group.lbls)
  group.lbls.uniq <- split(group.lbls.uniq, 1:length(group.lbls.uniq))
  
  res <- lapply(group.lbls.uniq, function(lbl) FUN(x[, group.lbls==lbl]))
  res <- dplyr::bind_cols(res)
  res <- as.data.frame(res)
  row.names(res) <- row.names(x)
  colnames(res) <- unlist(group.lbls.uniq)
  
  return(res)
}

#' Performs gene set enrichment analysis
#' @param gene.rnk A named vector of scores ranking genes
#' @param gene.sets List with gene sets to test enrichment
#' @param min.size Minimum size of gene sets to test
#' @param max.size Maximum size of gene sets to test
#' @param eps Sets boundary for calculating P value
#' @param nproc Sets BPPARAM to use nproc workers 
#' @return List with GSEA results
#' 
#' @author carlga
#'
runGSEA <- function(gene.rnk, gene.sets, min.size, max.size, eps = 0.0, nproc = 0) {
  require(fgsea)
  
  res <- list()
  res$all <- fgsea::fgsea(stats = gene.rnk,
                          pathways = gene.sets,
                          minSize = min.size,
                          maxSize = max.size,
                          eps = eps,
                          nproc = nproc)
  
  # res$collapsed <- fgsea::collapsePathways(fgseaRes = res$all[order(pval)],
  #                                          pathways = gene.sets,
  #                                          stats = gene.rnk)
  # res$collapsed <- res$all[pathway %in% res$collapsed$mainPathways]
  
  return(res)
}

