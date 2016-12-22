###Compare genes across different stages of tumors
###This module contains code that generates the output of result(DESeq2) comparing stages with each other
library(DESeq2)
##This is the preliminary analysis for the initial over sampling genes found using varselRF
length(over.sel.genes)

intersect(remove.dots(rownames(res)), over.sel.genes)
intersect(diff.genes[[2]], over.sel.genes)

intersect(remove.dots(rownames(res)[abs(res$log2FoldChange) > 0.5 & res$padj < 0.01 ]), over.sel.genes) #63 genes
intersect(remove.dots(rownames(res)[abs(res$log2FoldChange) > 0.2 & res$padj < 0.01 ]), over.sel.genes) #66 genes
intersect(remove.dots(rownames(res)[abs(res$log2FoldChange) > 0.2 & res$padj < 0.02 ]), over.sel.genes) #67 genes
intersect(remove.dots(rownames(res)[abs(res$log2FoldChange) > 0 & res$padj < 0.5 ]), over.sel.genes) #68 genes


comp.res <- function(dds_comp, contrast, stages.cont, stages.tum)
{
  stages.comp.cont <- lapply(stages.tum, function(x)
  {
    res = results(dds_comp, contrast = c(contrast, x, stages.cont))
    indexes = is.na(res[,5]) & is.na(res[,6])
    res = res[,!indexes]
  })
  names(stages.comp.cont) = stages.tum
  return(stages.comp.cont)
}


stage.comp = list()
stage.comp[['stage i']] <- comp.res(dds_tumor_reported,'stage', 'stage i',c('stage ii', 'stage iii', 'stage iv'))
stage.comp[['stage ii']] <- comp.res(dds_tumor_reported,'stage', 'stage ii',c('stage iii', 'stage iv'))
stage.comp[['stage iii']] <- comp.res(dds_tumor_reported,'stage', 'stage iii',c('stage iv'))

stages.comp.diff <- list()
stages.comp.diff[['stage i']] <- comp.res(dds_tumor_reported[diff.genes[[1]],], 'stage', 'stage i', c('stage ii', 'stage iii', 'stage iv'))
stages.comp.diff[['stage ii']] <- comp.res(dds_tumor_reported[diff.genes[[1]],],'stage', 'stage ii',c('stage iii', 'stage iv'))
stages.comp.diff[['stage iii']] <- comp.res(dds_tumor_reported[diff.genes[[1]],],'stage', 'stage iii',c('stage iv'))

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
genes.list <- get.genes.cont(stage.comp, 2, adj.pval = 0.01)
Reduce(union,Reduce(union,genes.list))

genes.list <- list()

genes.list[['stage i']] =list()
genes.list[['stage i']][['stage ii']] = get.genes(stage.comp$`stage i`$`stage ii`, 2, 0.01, 1)
genes.list[['stage i']][['stage iii']] = get.genes(stage.comp$`stage i`$`stage iii`,2 , 0.01, 1)
genes.list[['stage i']][['stage iv']] = get.genes(stage.comp$`stage i`$`stage iv`, 2, 0.01, 1)
