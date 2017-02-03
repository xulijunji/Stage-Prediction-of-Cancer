source('decomp.R')
source('shrunken/pamr.listgenes.R')
source('after_class.R')
source('main/final.R')
library(DESeq2)

load('environment/req_dfs.RData')
load('environment/stages.level.comb.RData')
load('environment/stage.index.RData')
load('environment/first_trial_shrunken_classifier.RData')
load('environment/accuracy_feature/shrunken_list_vs_normal_reported.RData')

###Saving the results of the first trial into a single list
first.trial <- list()
first.trial[['gr']] <- gr
first.trial[['train']] <- pamr.train.comb
first.trial[['cv']] <- pamr.cv.comb
first.trial[['genes']] <- pamr.genes.comb
first.trial[['auc']] <- pamr.aucs.comb
first.trial[['conf.mat']] <- confusion.mat
first.trial[['eval.mat']] <- eval.mat
save(first.trial, file = 'environment/first_trial_shrunken_classifier.RData')
sapply(first.trial$genes, length)
load('environment/first_trial_shrunken_classifier.RData')

confusion.mat <- list()
eval.mat <- list()

out.list <- do.shrunken(gr, req.dfs$vs, stages.levels.comb, confusion.mat, eval.mat)
genes <- Reduce(intersect, out.list[[3]])
out.list <- do.knn(gr, req.dfs$vs, genes, stages.levels.comb, out.list[[1]], out.list[[2]])
out.list <- do.naive(gr, req.dfs$vs, genes, stages.levels.comb, out.list[[1]], out.list[[2]])
out.list <- do.svm(gr, req.dfs$vs, genes, stages.levels.comb, out.list[[1]], out.list[[2]])
out.list <- do.rf(gr, req.dfs$vs, genes.1, stages.levels.comb, confusion.mat, eval.mat)

genes.group <- pamr.listgene(first.trial$train[[5]], 
                             data = list(x=as.matrix(t(req.dfs$vs[sort(unlist(first.trial$gr[-5])),])),
                                         y=stages.levels.comb[sort(unlist(first.trial$gr[-5]))]), 
                             threshold = first.trial$cv[[5]]$threshold[21],
                             fitcv = first.trial$cv[[5]], genenames = T)
genes.2 <- genes.group[,2][as.numeric(genes.group[,4]) > 0]
genes.1 <- genes.group[,2][as.numeric(genes.group[,3]) > 0]

#########Above is without using normal

##vst with normal
load('environment/dds_tumor_reported_normal_stage.RData')
load('environment/dds_nor_tum_comb.RData')
vst_normal_reported <- t(assay(vst(dds_tumor_reported_normal_stage)))

vs_normal_comb_reported <- t(assay(vst(dds_nor_tum_comb)))

View(vst_normal_reported)
View(req.dfs$vs)
View(colData(dds_tumor_reported_normal_stage))
View(df.stage.tumor.rep)


tumor.ind.vs <- which(colData(dds_tumor_reported_normal_stage)[,2] == 'T')
match.ind <- match(colData(dds_tumor_reported_normal_stage)[,1][tumor.ind.vs],
                   df.stage.tumor.rep$sample.id)

sum(droplevels(colData(dds_tumor_reported_normal_stage)[tumor.ind.vs,4]) ==
  stages.levels)

#########Trial
# o1 <- do.shrunken(first.trial$gr, vst_normal_reported[tumor.ind.vs,], stages.levels.comb, confusion.mat, eval.mat)
# o2 <- do.shrunken(first.trial$gr, vst_normal_reported[tumor.ind.vs,], stages.levels.comb, confusion.mat, eval.mat)
# o3 <- do.shrunken(first.trial$gr, req.dfs$vs, stages.levels.comb, confusion.mat, eval.mat)
# View(o1[[3]][[1]])
# g <- get.genes.shrunken(o1[[3]])
# gd <- get.genes.shrunken(o2[[3]])
# go <- get.genes.shrunken(o3[[3]])
# sapply(g,length)
# length(intersect(g2,
#                  get.intersect.genes(g, c(1,2,3,4,5))))

acc.list <- do.shrunken(first.trial$gr, 
                        vs_normal_comb_reported[tumor.ind.vs,], 
                        stages.levels.comb, list(), list())
features.trials.list <- get.genes.shrunken(acc.list[[3]])
sapply(features.trials.list, length)
length(intersect(g2, features.trials.list[[2]]))

features.normal.rep <- list()
features.normal.rep[[1]] <- get.genes.common(features.trials.list, 1)
features.normal.rep[[3]] <- get.genes.common(features.trials.list, 3)
features.normal.rep[[5]] <- get.genes.common(features.trials.list, 5)

accuracy.list <- list()
accuracy.list[['all']] <- acc.list
accuracy.list[['at_3']] <- acc.list
accuracy.list[['at_1']] <- acc.list
remove(acc.list)

out.list <- do.knn(gr, req.dfs$vs, genes, stages.levels.comb, out.list[[1]], out.list[[2]])
accuracy.list[['all']] <- do.knn(first.trial$gr, 
                      vst_normal_reported[tumor.ind.vs,], 
                      features.normal.rep[[5]],
                      stages.levels.comb, accuracy.list$all[[1]], 
                      accuracy.list$all[[2]])
accuracy.list[['all']] <- do.rf(first.trial$gr, 
                                 vst_normal_reported[tumor.ind.vs,], 
                                 features.normal.rep[[5]],
                             stages.levels.comb, accuracy.list$all[[1]], 
                                 accuracy.list$all[[2]])
accuracy.list[['all']] <- do.naive(first.trial$gr, 
                                 vst_normal_reported[tumor.ind.vs,], 
                                 features.normal.rep[[5]],
                                 stages.levels.comb, accuracy.list$all[[1]], 
                                 accuracy.list$all[[2]])
accuracy.list[['all']] <- do.svm(first.trial$gr, 
                                 vst_normal_reported[tumor.ind.vs,], 
                                 features.normal.rep[[5]],
                                 stages.levels.comb, accuracy.list$all[[1]], 
                                 accuracy.list$all[[2]])


accuracy.list[['at_3']] <-do.knn(first.trial$gr, 
                          vst_normal_reported[tumor.ind.vs,], 
                          features.normal.rep[[3]],
                          stages.levels.comb, accuracy.list$at_3[[1]], 
                          accuracy.list$at_3[[2]])
accuracy.list[['at_3']] <-do.rf(first.trial$gr, 
                                 vst_normal_reported[tumor.ind.vs,], 
                                 features.normal.rep[[3]],
                                 stages.levels.comb, accuracy.list$at_3[[1]], 
                                 accuracy.list$at_3[[2]])
accuracy.list[['at_3']] <-do.naive(first.trial$gr, 
                                 vst_normal_reported[tumor.ind.vs,], 
                                 features.normal.rep[[3]],
                                 stages.levels.comb, accuracy.list$at_3[[1]], 
                                 accuracy.list$at_3[[2]])
accuracy.list[['at_3']] <-do.svm(first.trial$gr, 
                                 vst_normal_reported[tumor.ind.vs,], 
                                 features.normal.rep[[3]],
                                 stages.levels.comb, accuracy.list$at_3[[1]], 
                                 accuracy.list$at_3[[2]])

accuracy.list[['at_1']] <-do.knn(first.trial$gr, 
                                 vst_normal_reported[tumor.ind.vs,], 
                                 features.normal.rep[[1]],
                                 stages.levels.comb, accuracy.list$at_1[[1]], 
                                 accuracy.list$at_1[[2]])
accuracy.list[['at_1']] <-do.rf(first.trial$gr, 
                                 vst_normal_reported[tumor.ind.vs,], 
                                 features.normal.rep[[1]],
                           stages.levels.comb, accuracy.list$at_1[[1]], 
                                 accuracy.list$at_1[[2]])
accuracy.list[['at_1']] <-do.naive(first.trial$gr, 
                                 vst_normal_reported[tumor.ind.vs,], 
                                 features.normal.rep[[1]],
                                 stages.levels.comb, accuracy.list$at_1[[1]], 
                                 accuracy.list$at_1[[2]])
accuracy.list[['at_1']] <-do.svm(first.trial$gr, 
                                 vst_normal_reported[tumor.ind.vs,], 
                                 features.normal.rep[[1]],
                                 stages.levels.comb, accuracy.list$at_1[[1]], 
                                 accuracy.list$at_1[[2]])
save(accuracy.list, file = 'environment/accuracy_list_vs_normal_reported.RData')
load('environment/accuracy_feature/shrunken_list_vs_normal_reported.RData')