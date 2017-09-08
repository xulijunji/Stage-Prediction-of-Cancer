source('main/updated/initialisation.R')

trial.res <- list()

trial.res[['shrunken']] <- lapply(net.features.trial$shrunken[c(3,4,5,6)], function(genes)
  {
  final.res(vst_tumor_tum, train.trial.ind, test.trial.ind, stages.levels.comb, genes, 10)
})
trial.res[['varselRF']] <- lapply(net.features.trial$varSelRF[c(3,4,5,6)], function(genes)
{
  final.res(vst_tumor_tum, train.trial.ind, test.trial.ind, stages.levels.comb, genes, 10)
})
trial.res[['deseq2']] <- list()


get.deseq.features <- function(features.list)
{
  deseq.features <- list()  
  for(i in seq_along(features.list$deseq2$atleast_1))
  {
    fold = names(features.list$deseq2$atleast_1)[i]
    deseq.features[[fold]] <- sapply(c(3,4,5,6), function(at_ind)
    {
      unlist(features.list$deseq2[[at_ind]][i])
    }, simplify = T)
    names(deseq.features[[fold]]) <- names(features.list$deseq2)[c(3,4,5,6)]
  }
  return(deseq.features)
}

deseq.or.features <- get.deseq.features(net.features.updated)
deseq.features <- get.deseq.features(net.features.trial)

trial.res[['deseq2']][['1 fold']] <- lapply(deseq.features$`1 fold`, function(genes)
                                {final.res(vst_tumor_tum, train.trial.ind, test.trial.ind, 
                                               stages.levels.comb, genes, 10)})
trial.res[['deseq2']][['1 fold']] <- lapply(deseq.features$`1 fold`, function(genes)
{final.res(vst_tumor_tum, train.trial.ind, test.trial.ind, 
           stages.levels.comb, genes, 10)})
trial.res[['deseq2']][['1.5 fold']] <- lapply(deseq.features$`1.5 fold`, function(genes)
{final.res(vst_tumor_tum, train.trial.ind, test.trial.ind, 
           stages.levels.comb, genes, 10)})
trial.res[['deseq2']][['2 fold']] <- lapply(deseq.features$`2 fold`, function(genes)
{final.res(vst_tumor_tum, train.trial.ind, test.trial.ind, 
           stages.levels.comb, genes, 10)})

trial.res.or.fe <- lapply(deseq.features$`2 fold`, function(genes)
{final.res(vst_tumor_tum, train.indexes, test.indexes, 
           stages.levels.comb, genes, 10)})
trial.res.new.feat <- lapply(deseq.or.features$`2 fold`, function(genes)
{final.res(vst_tumor_tum, train.trial.ind, test.trial.ind, 
           stages.levels.comb, genes, 10)})

 <- create.net.df(test.results)
b.cv <- create.net.df(cv.results)

library(ggplot2)
ggplot(b, aes(x=Classifier, y=AUC, fill=Feature_Selection)) +
  geom_boxplot()

micr.class <- final.res(t(merged_final[,1:34]), c(1:34), c(1:34), sample_micro_info$stage, 
                        net.features.trial$varSelRF$atleast_4, 10)

get.intersect <- function(genes.list)
{
  return(Reduce(intersect, genes.list))
}
g1 <- get.intersect(list(fea.trial.list$`sam1.5 fold`$atleast_3, fea.trial.list$`deseq21.5 fold`$atleast_3,
                         fea.trial.list$varSelRF$atleast_3, fea.trial.list$shrunken$atleast_3))
length(g1)
t <- final.res(vst_tumor_tum, train.trial.ind, test.trial.ind, stages.levels.comb, g1, 10)
get.aucs(t[[2]])
1get.aucs(res.wcgna[[2]])
length(intersect(wcgna.genes.ent.ens, g1))
intersect(wcgna.genes.ent.ens, diff.genes$`2`)

res.wcgna.mine <- final.res(vst_tumor_tum, train.trial.ind, test.trial.ind, stages.levels.comb,
                            intersect(g1, wcgna.genes.ent.ens), 10)
get.aucs(res.wcgna.mine[[1]])                  
length(intersect(fea.trial.list$varSelRF$atleast_4, wcgna.genes.ent.ens))
