source('main/updated/initialisation.R')
#net.features.updated <- list()
get.feature.shrunken <- function(data, stages, train.ind.list)
{
  shrunken.features <- list()
  shrunken.features.ob <- get.shrunken.features(train.ind.list, data, stages)
  shrunken.features[['genes.object']] <- shrunken.features.ob$genes
  shrunken.features[['genes']] <- get.genes.shrunken(shrunken.features$genes.object)
  shrunken.features[['atleast_1']] <- get.genes.common(shrunken.features$genes, 1)
  shrunken.features[['atleast_2']] <- get.genes.common(shrunken.features$genes, 2)
  shrunken.features[['atleast_3']] <- get.genes.common(shrunken.features$genes, 3)
  shrunken.features[['atleast_4']] <- get.genes.common(shrunken.features$genes, 4)
  shrunken.features[['genes_stage']] <- get.shrunken.group.stage(shrunken.features$genes.object)
  shrunken.features[['atleast_dfs']] <- get.shrunken.common.stage(shrunken.features$genes_stage, c(1,2,3,4))
  print('Shrunken Completed')
  return(shrunken.features)
}

get.feature.varSelRf <- function(data, stages, train.ind.list)
{
  varselRF.features <- list()
  varselRF.features[['genes.object']] <- do.varselRF(train.ind.list, data, stages)
  varSelRF.features[['genes.list']] <- get.min.oob.varselRf(varselRF.features$genes.object)
  varselRF.features[['atleast_1']] <- get.genes.common(varselRF.features$genes.list,1)
  varselRF.features[['atleast_2']] <- get.genes.common(varselRF.features$genes.list,2)
  varselRF.features[['atleast_3']] <- get.genes.common(varselRF.features$genes.list,3)
  varselRF.features[['atleast_4']] <- get.genes.common(varselRF.features$genes.list,4)
  print('VarSelRF Completed')
  return(varselRF.features)
}

get.feature.deseq2 <- function(dds.tum, stages, train.ind.list)
{
  sample.df <- colData(dds.tum)
  sample.df[,3] <- sapply(as.character(sample.df[,3]), function(stage)
  {
    if(stage == 'stage ii')
      stage = 'stage i'
    else if(stage == 'stage iii')
      stage = 'stage iv'
    else
      stage = stage
  })
  sample.df$stage <- as.factor(sample.df$stage)
  
  dds_obj <- create.Deseq2(train.ind.list, counts(dds.tum),
                           colData = sample.df)
  res.train.dds_obj <- do.Deseq2(dds_obj)
  
  features.deseq2 <- list()
  features.deseq2[['genes.object']] <- res.train.dds_obj
  features.deseq2[['genes.list']] <- get.deseq2.genes(res.train.dds_obj, c(1,1.5,2))
  
  features.deseq2[['atleast_1']] <- sapply(features.deseq2$genes.list, function(x)
  {
    get.genes.common(x, 1)
  })
  features.deseq2[['atleast_2']] <- sapply(features.deseq2$genes.list, function(x)
  {
    get.genes.common(x, 2)
  })
  
  features.deseq2[['atleast_3']] <- sapply(features.deseq2$genes.list, function(x)
  {
    get.genes.common(x, 3)
  })
  features.deseq2[['atleast_4']] <- sapply(features.deseq2$genes.list, function(x)
  {
    get.genes.common(x, 4)
  })
  return(features.deseq2)
  print('Finished DESeq2')
  
}

get.features.sam <- function(data, stages, train.ind.list)
{
  res.sam <- get.sam.features(train.ind.list, data, stages)
  sam.genes <- lapply(res.sam, function(sam.obj)
  {
    a <- rbind(sam.obj$siggenes.table$genes.up, sam.obj$siggenes.table$genes.lo)
    data.frame(GeneID = as.character(a[,1]), Gene.Name = as.character(a[,2]), Score = as.numeric(a[,3]), 
               Fold.Change = as.numeric(a[,4]), q.val = as.numeric(a[,5]), stringsAsFactors = F)
    
  })
  
  features.sam <- list()
  features.sam[['sam.obj']] <- res.sam
  features.sam[['genes.list']] <- get.sam.genes(sam.genes, c(1, 1.5, 2))
  
  features.sam[['atleast_1']] <- sapply(features.sam$genes.list, function(genes.list)
  {
    get.genes.common(genes.list, 1)
  })
  features.sam[['atleast_2']] <- sapply(features.sam$genes.list, function(genes.list)
  {
    get.genes.common(genes.list, 2)
  })
  features.sam[['atleast_3']] <- sapply(features.sam$genes.list, function(genes.list)
  {
    get.genes.common(genes.list, 3)
  })
  features.sam[['atleast_4']] <- sapply(features.sam$genes.list, function(genes.list)
  {
    get.genes.common(genes.list, 4)
  })
  return(features.sam)
}
get.features <- function(data, stages, train.ind.list)
{
  net.features <- list()
  ######Shrunken####
  net.features[['shrunken']] <- get.feature.shrunken(data, stages, train.ind.list)
  
  #####VarSelRF####
  net.features[['varSelRF']] <- get.feature.varSelRf(data, stages, train.ind.list)

  ######DESeq2
  load('environment/accuracy_feature/updated/dds_tumor_reported.RData')
  net.features[['deseq2']] <- get.feature.deseq2(dds_tumor_reported, stages, train.ind.list)
  
  ####SAMSeq
  net.features[['sam']] <- get.features.sam(assay(dds_tumor_reported), stages, train.ind.list)
  
  return(net.features)
}
# net.features.updated <- get.features(vst_tumor_tum, stages.levels.comb, gr.updated.train)
# net.features.trial <- get.features(vst_tumor_tum, stages.levels.comb, gr.trial.train)


get.filter.fea <- function(features.list, folds.list)
{
  filter.fea <- list()
  for(i in seq_along(folds.list))
  {
    filter.fea[[folds.list[i]]] <- lapply(features.list, function(genes.list)
    {
      genes.list[[folds.list[i]]]  
    })
  }
  return(filter.fea)
}
get.class.fea <- function(net.fea)
{
  fea.list <- list()
  for(i in seq_along(net.fea))
  {
    fea.name <- names(net.fea)[i]
    if(fea.name %in% c('shrunken','varSelRF'))
      fea.list[[fea.name]] <- net.fea[[fea.name]][c(3,4,5,6)]
    else
    {
      temp.fea <- get.filter.fea(net.fea[[fea.name]][c(3,4,5,6)], paste0(c(1,1.5,2), ' fold'))
      for(j in seq_along(temp.fea))
        fea.list[[paste0(fea.name, names(temp.fea)[j])]] <- temp.fea[[names(temp.fea)[j]]]
    }
  }
  return(fea.list)
}
