source('main/updated/initialisation.R')

net.features.updated <- list()
######Shrunken####
net.features.updated[['shrunken']] <- list()
shrunken.features.ob <- get.shrunken.features(gr.updated.train, vst_tumor_tum, stages.levels.comb)
net.features.updated$shrunken[['genes.object']] <- shrunken.features.ob$genes
net.features.updated$shrunken[['genes']] <- get.genes.shrunken(net.features.updated$shrunken$genes.object)
net.features.updated$shrunken[['atleast_1']] <- get.genes.common(net.features.updated$shrunken$genes, 1)
net.features.updated$shrunken[['atleast_2']] <- get.genes.common(net.features.updated$shrunken$genes, 2)
net.features.updated$shrunken[['atleast_3']] <- get.genes.common(net.features.updated$shrunken$genes, 3)
net.features.updated$shrunken[['atleast_4']] <- get.genes.common(net.features.updated$shrunken$genes, 4)
net.features.updated$shrunken[['genes_stage']] <- get.shrunken.group.stage(net.features.updated$shrunken$genes.object)
net.features.updated$shrunken[['atleast_dfs']] <- get.shrunken.common.stage(net.features.updated$shrunken$genes_stage, c(1,2,3,4))
save(net.features.updated,file = 'environment/accuracy_feature/updated/net_features_updated.RData')

######VarSelRF####
net.features.updated[['varSelRF']] <- list()
net.features.updated[['varSelRF']][['genes.object']] <- do.varselRF(gr.updated.train, vst_tumor_tum, stages.levels.comb)
net.features.updated$varSelRF[['genes.list']] <- get.min.oob.varselRf(net.features.updated$varSelRF$genes.object)
net.features.updated$varSelRF[['atleast_1']] <- get.genes.common(net.features.updated$varSelRF$genes.list,1)
net.features.updated$varSelRF[['atleast_2']] <- get.genes.common(net.features.updated$varSelRF$genes.list,2)
net.features.updated$varSelRF[['atleast_3']] <- get.genes.common(net.features.updated$varSelRF$genes.list,3)
net.features.updated$varSelRF[['atleast_4']] <- get.genes.common(net.features.updated$varSelRF$genes.list,4)

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
res.train.dds_obj <- do.Deseq2(dds_obj_updated)

net.features.updated[['deseq2']] <- list()
net.features.updated$deseq2[['genes.object']] <- res.train.dds_obj
net.features.updated$deseq2[['genes.list']] <- get.deseq2.genes(res.train.dds_obj, c(1,1.5,2))
net.features.updated$deseq2$genes.list$'3 fold' <- net.features.updated$deseq2$genes.list$`3 fold`$`3 fold`
net.features.updated$deseq2$genes.list$'2.5 fold' <- 
  get.deseq2.genes(res.train.dds_obj, c(2.5))
net.features.updated$deseq2$genes.list$'2.5 fold' <- net.features.updated$deseq2$genes.list$`2.5 fold`$`2.5 fold`

net.features.updated$deseq2[['atleast_1']] <- sapply(net.features.updated$deseq2$genes.list, function(x)
{
  get.genes.common(x, 1)
})
net.features.updated$deseq2[['atleast_2']] <- sapply(net.features.updated$deseq2$genes.list, function(x)
{
  get.genes.common(x, 2)
})

net.features.updated$deseq2[['atleast_3']] <- sapply(net.features.updated$deseq2$genes.list, function(x)
{
  get.genes.common(x, 3)
})
net.features.updated$deseq2[['atleast_4']] <- sapply(net.features.updated$deseq2$genes.list, function(x)
{
  get.genes.common(x, 4)
})

