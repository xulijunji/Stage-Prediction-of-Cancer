source('decomp.R')
source('main/updated/helper_func.R')
source('main/final.R')
#library(pamr)

load('environment/stages.level.comb.RData')
load('environment/accuracy_feature/dds_nor_tum_comb.RData')
load('environment/accuracy_feature/vs_nor_comb.RData')
load('environment/accuracy_feature/tumor_ind_vs.RData')
load('environment/first_trial_shrunken_classifier.RData')
load('environment/accuracy_feature/updated/gr.RData')
load('environment/accuracy_feature/updated/net_features_updated.RData')
load('environment/accuracy_feature/net_features.RData')

vst_tumor_tum <- vs_normal_comb_reported[tumor.ind.vs, ]

gr.updated <- build.groups(260,5)
stage.dist <- get.stage.distribution(gr.updated, stages.levels.comb)
gr.updated.train <- gr.updated[-1]

test.group <- gr.updated[[1]]
train.group <- unlist(gr.updated[-1])

net.features.updated <- list()
######Shrunken####
net.features.updated[['shrunken']] <- list()
shrunken.features.ob <- get.shrunken.features(gr.train, vst_tumor_tum, stages.levels.comb)
net.features.updated$shrunken[['genes.object']] <- shrunken.features.ob$genes
net.features.updated$shrunken[['genes']] <- get.genes.shrunken(net.features.updated$shrunken$genes.object)
net.features.updated$shrunken[['atleast_1']] <- get.genes.common(net.features.updated$shrunken$genes, 1)
net.features.updated$shrunken[['atleast_2']] <- get.genes.common(net.features.updated$shrunken$genes, 2)
net.features.updated$shrunken[['atleast_4']] <- get.genes.common(net.features.updated$shrunken$genes, 4)
net.features.updated$shrunken[['genes_stage']] <- get.shrunken.group.stage(net.features.updated$shrunken$genes.object)
net.features.updated$shrunken[['atleast_dfs']] <- get.shrunken.common.stage(net.features.updated$shrunken$genes_stage, c(1,2,4))
save(net.features.updated,file = 'environment/accuracy_feature/updated/net_features_updated.RData')

######VarSelRF####
net.features.updated[['varSelRF']] <- list()
net.features.updated[['varSelRF']][['genes.object']] <- do.varselRF(gr.train, vst_tumor_tum, stages.levels.comb)
net.features.updated$varSelRF[['genes.list']] <- get.min.oob.varselRf(net.features.updated$varSelRF$genes.object)
net.features.updated$varSelRF[['atleast_1']] <- get.genes.common(net.features.updated$varSelRF$genes.list,1)
net.features.updated$varSelRF[['atleast_2']] <- get.genes.common(net.features.updated$varSelRF$genes.list,2)
net.features.updated$varSelRF[['atleast_3']] <- get.genes.common(net.features.updated$varSelRF$genes.list,3)
net.features.updated$varSelRF[['atleast_4']] <- get.genes.common(net.features.updated$varSelRF$genes.list,4)

######DESeq2
load('environment/accuracy_feature/updated/dds_tumor_reported.RData')
net.features.updated[['deseq2']] <- list()
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
