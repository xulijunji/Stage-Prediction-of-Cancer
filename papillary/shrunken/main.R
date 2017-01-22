source('decomp.R')
source('shrunken/pamr.listgenes.R')
source('after_class.R')
source('shrunken/final.R')
library(DESeq2)

load('environment/req_dfs.RData')
load('environment/stages.level.comb.RData')
load('environment/stage.index.RData')
load('environment/first_trial_shrunken_classifier.RData')

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
out.list <- do.knn(gr, req.dfs$vs, genes., stages.levels.comb, out.list[[1]], out.list[[2]])
out.list <- do.naive(gr, req.dfs$vs, genes, stages.levels.comb, out.list[[1]], out.list[[2]])
out.list <- do.svm(gr, req.dfs$vs, genes, stages.levels.comb, out.list[[1]], out.list[[2]])
out.list <- do.rf(gr, req.dfs$vs, genes.1, stages.levels.comb, confusion.mat, eval.mat)


##vst with normal
load('environment/dds_tumor_reported_normal_stage.RData')
vst_normal_reported <- t(assay(vst(dds_tumor_reported_normal_stage)))
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

