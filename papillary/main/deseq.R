library(DESeq2)
load('environment/stages.level.comb.RData')
load('environment/dds_tumor_reported_normal_stage.RData')
load('environment/dds_tumor_reported.RData')
load('environment/first_trial_shrunken_classifier.RData')
load('environment/dds_object.RData')
load('environment/accuracy_feature/classifer_list.RData')
load('environment/accuracy_feature/net_features.RData')
source('across_tumors.R')
source('main/final.R')

sample.df <- colData(dds_tumor_reported_normal_stage)
sample.df[,4] <- sapply(as.character(sample.df[,4]), function(stage)
  {
  if(stage == 'stage ii')
    stage = 'stage i'
  else if(stage == 'stage iii')
    stage = 'stage iv'
  else
    stage = stage
})
dds_nor_tum_comb <- 
  DESeqDataSetFromMatrix(counts(dds_tumor_reported_normal_stage),
                         colData = sample.df, design = ~stage.type)
dds_nor_tum_comb <- dds_nor_tum_comb[rowSums(assay(dds_nor_tum_comb)) > 2]
dds_nor_tum_comb <- DESeq(dds_nor_tum_comb)
save(dds_nor_tum_comb, file = 'environment/dds_nor_tum_comb.RData')

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
gr = first.trial$gr
dds_obj <- create.Deseq2(gr, counts(dds_tumor_reported), 
                         colData = sample.df)
res.train.dds_obj <- do.Deseq2(dds_obj)

deseq.genes.list <- lapply(c(1,1.5,2), function(fold)
  {
  lapply(res.train.dds_obj, function(x)
  {
    get.genes(x, fold, 0.05, 0.05)
  })
})
names(deseq.genes.list) <- c("1fold", "1.5fold", "2fold")


net.features[['deseq']] <- list()
net.features$deseq[['genes.object']] <- res.train.dds_obj
net.features$deseq[['genes.list']] <- deseq.genes.list
net.features$deseq[['atleast_1']] <- sapply(net.features$deseq$genes.list, function(x)
  {
  get.genes.common(x, 1)
})
net.features$deseq[['atleast_3']] <- sapply(net.features$deseq$genes.list, function(x)
{
  get.genes.common(x, 3)
})
net.features$deseq[['atleast_5']] <- sapply(net.features$deseq$genes.list, function(x)
{
  get.genes.common(x, 5)
})

length(intersect(net.features$deseq$atleast_1, g2))
sapply(net.features$deseq$genes.list$`1.5fold`, length)

classifier.list[['deseq']] <- list()
classifier.list$deseq[['atleast_1']] <- do.rf(first.trial$gr, 
                                              vs_normal_comb_reported[tumor.ind.vs,], 
                                              net.features$deseq$atleast_1,
                                              stages.levels.comb, list(), 
                                              list())
classifier.list$deseq[['atleast_1']] <- do.knn(first.trial$gr, 
                                                  vs_normal_comb_reported[tumor.ind.vs,], 
                                                  net.features$deseq$atleast_1,
                                                  stages.levels.comb, classifier.list$deseq$atleast_1[[1]], 
                                                  classifier.list$deseq$atleast_1[[2]])
classifier.list$deseq[['atleast_1']] <- do.svm(first.trial$gr, 
                                               vs_normal_comb_reported[tumor.ind.vs,], 
                                               net.features$deseq$atleast_1,
                                               stages.levels.comb, classifier.list$deseq$atleast_1[[1]], 
                                               classifier.list$deseq$atleast_1[[2]])
classifier.list$deseq[['atleast_1']] <- do.naive(first.trial$gr, 
                                               vs_normal_comb_reported[tumor.ind.vs,], 
                                               net.features$deseq$atleast_1,
                                               stages.levels.comb, classifier.list$deseq$atleast_1[[1]], 
                                               classifier.list$deseq$atleast_1[[2]])

classifier.list$deseq[['atleast_3']] <- do.rf(first.trial$gr, 
                                              vs_normal_comb_reported[tumor.ind.vs,], 
                                              net.features$deseq$atleast_3,
                                              stages.levels.comb, list(), 
                                              list())
classifier.list$deseq[['atleast_3']] <- do.naive(first.trial$gr, 
                                                 vs_normal_comb_reported[tumor.ind.vs,], 
                                                 net.features$deseq$atleast_3,
                                                 stages.levels.comb, classifier.list$deseq$atleast_3[[1]], 
                                                 classifier.list$deseq$atleast_3[[2]])
classifier.list$deseq[['atleast_3']] <- do.knn(first.trial$gr, 
                                                 vs_normal_comb_reported[tumor.ind.vs,], 
                                                 net.features$deseq$atleast_3,
                                                 stages.levels.comb, classifier.list$deseq$atleast_3[[1]], 
                                                 classifier.list$deseq$atleast_3[[2]])
classifier.list$deseq[['atleast_3']] <- do.svm(first.trial$gr, 
                                                 vs_normal_comb_reported[tumor.ind.vs,], 
                                                 net.features$deseq$atleast_3,
                                                 stages.levels.comb, classifier.list$deseq$atleast_3[[1]], 
                                                 classifier.list$deseq$atleast_3[[2]])

classifier.list$deseq[['atleast_5']] <- do.rf(first.trial$gr, 
                                              vs_normal_comb_reported[tumor.ind.vs,], 
                                              net.features$deseq$atleast_5,
                                              stages.levels.comb, list(), 
                                              list())
classifier.list$deseq[['atleast_5']] <- do.svm(first.trial$gr, 
                                               vs_normal_comb_reported[tumor.ind.vs,], 
                                               net.features$deseq$atleast_5,
                                               stages.levels.comb, classifier.list$deseq$atleast_5[[1]], 
                                               classifier.list$deseq$atleast_5[[2]])
classifier.list$deseq[['atleast_5']] <- do.knn(first.trial$gr, 
                                               vs_normal_comb_reported[tumor.ind.vs,], 
                                               net.features$deseq$atleast_5,
                                               stages.levels.comb, classifier.list$deseq$atleast_5[[1]], 
                                               classifier.list$deseq$atleast_5[[2]])
classifier.list$deseq[['atleast_5']] <- do.naive(first.trial$gr, 
                                               vs_normal_comb_reported[tumor.ind.vs,], 
                                               net.features$deseq$atleast_5,
                                               stages.levels.comb, classifier.list$deseq$atleast_5[[1]], 
                                               classifier.list$deseq$atleast_5[[2]])
save(net.features, file = 'environment/accuracy_feature/net_features.RData')
save(classifier.list, file = 'environment/accuracy_feature/classifer_list.RData')
