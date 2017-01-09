###Compare genes across different stages of tumors
###This module contains code that generates the output of result(DESeq2) comparing stages with each other
library(DESeq2)
##This is the preliminary analysis for the initial over sampling genes found using varselRF

comp.res <- function(dds_comp, contrast, stages.cont, stages.tum)
{
  stages.comp.cont <- lapply(stages.tum, function(x)
  {
    res = results(dds_comp, contrast = c(contrast, x, stages.cont))
    indexes = is.na(res[,5]) | is.na(res[,6])
    res = res[!indexes,]
  })
  names(stages.comp.cont) = stages.tum
  return(stages.comp.cont)
}


get.genes.cont <- function(stage.comp, logfc, adj.pval = NULL, pval = NULL)
{
  genes.list = list()
  if(is.null(adj.pval))
    adj.pval = 1
  if(is.null(pval))
    pval = 1
  for(i in seq_along(stage.comp))
  {
    genes.list[[names(stage.comp)[i]]] = lapply(stage.comp[[i]], function(x)
      {
        get.genes(x,logfc,adj.pval,pval)
    })
    names(genes.list[[names(stage.comp)[i]]]) = names(stage.comp[[i]])
  }
  return(genes.list)
}
get.genes <- function(res, logfc, adj.pval, pval)
{
  return(rownames(res)[abs(res$log2FoldChange) > logfc & res$padj < adj.pval & res$pvalue < pval])
}

get.unique.genes <- function(genes.list)
{
  genes_unique <- sapply(seq(length(genes.list)), function(x)
    {
    g <- Reduce(union, genes.list[-x])
    setdiff(unlist(genes.list[[x]]), g)
  })
  names(genes_unique) = names(genes.list)
  return(genes_unique)
}
