source('main/updated/initialisation.R')
#net.features.updated <- list()
get.features <- function(data, stages, train.ind, file_save)
{
  net.features <- list()
  ######Shrunken####
  net.features[['shrunken']] <- list()
  shrunken.features.ob <- get.shrunken.features(train.ind, data, stages)
  net.features$shrunken[['genes.object']] <- shrunken.features.ob$genes
  net.features$shrunken[['genes']] <- get.genes.shrunken(net.features$shrunken$genes.object)
  net.features$shrunken[['atleast_1']] <- get.genes.common(net.features$shrunken$genes, 1)
  net.features$shrunken[['atleast_2']] <- get.genes.common(net.features$shrunken$genes, 2)
  net.features$shrunken[['atleast_3']] <- get.genes.common(net.features$shrunken$genes, 3)
  net.features$shrunken[['atleast_4']] <- get.genes.common(net.features$shrunken$genes, 4)
  net.features$shrunken[['genes_stage']] <- get.shrunken.group.stage(net.features$shrunken$genes.object)
  net.features$shrunken[['atleast_dfs']] <- get.shrunken.common.stage(net.features$shrunken$genes_stage, c(1,2,3,4))
  save(net.features.updated,file = file_save)
  print('Shrunken Completed')

  #####VarSelRF####
  net.features[['varSelRF']] <- list()
  net.features[['varSelRF']][['genes.object']] <- do.varselRF(gr.trial.train, data, stages)
  net.features$varSelRF[['genes.list']] <- get.min.oob.varselRf(net.features$varSelRF$genes.object)
  net.features$varSelRF[['atleast_1']] <- get.genes.common(net.features$varSelRF$genes.list,1)
  net.features$varSelRF[['atleast_2']] <- get.genes.common(net.features$varSelRF$genes.list,2)
  net.features$varSelRF[['atleast_3']] <- get.genes.common(net.features$varSelRF$genes.list,3)
  net.features$varSelRF[['atleast_4']] <- get.genes.common(net.features$varSelRF$genes.list,4)
  save(net.features,file = file_save)
  print('VarSelRF Completed')

  ######DESeq2
  load('environment/accuracy_feature/updated/dds_tumor_reported.RData')
  sample.df <- colData(dds_tumor_reported)
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

  dds_obj <- create.Deseq2(gr.updated.train, counts(dds_tumor_reported),
                           colData = sample.df)
  res.train.dds_obj <- do.Deseq2(dds_obj)

  net.features[['deseq2']] <- list()
  net.features$deseq2[['genes.object']] <- res.train.dds_obj
  net.features$deseq2[['genes.list']] <- get.deseq2.genes(res.train.dds_obj, c(1,1.5,2))

  net.features$deseq2[['atleast_1']] <- sapply(net.features$deseq2$genes.list, function(x)
  {
    get.genes.common(x, 1)
  })
  net.features$deseq2[['atleast_2']] <- sapply(net.features$deseq2$genes.list, function(x)
  {
    get.genes.common(x, 2)
  })

  net.features$deseq2[['atleast_3']] <- sapply(net.features$deseq2$genes.list, function(x)
  {
    get.genes.common(x, 3)
  })
  net.features$deseq2[['atleast_4']] <- sapply(net.features$deseq2$genes.list, function(x)
  {
    get.genes.common(x, 4)
  })
  print('Finished DESeq2')
#  save(net.features, file = file_save)
  return(net.features)
}
net.features.updated <- get.features(vst_tumor_tum, stages.levels.comb, gr.updated.train, 
                                     'environment/accuracy_feature/updated/net_features_updated.RData')
net.features.trial <- get.features(vst_tumor_tum, stages.levels.comb, gr.trial.train, 
                                     'environment/accuracy_feature/updated/net_features_trial.RData')



