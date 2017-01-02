load('environment/only_tumor_reported.RData')
load('environment/stages_levels.RData')
load('environment/stages.level.comb.RData')
load('environment/diff_genes.RData')

source('shrunken/pamr.listgenes.R')
library(DESeq2)
library(pROC)
library(pamr)

diff.genes$0.5
########All Classes##############
####All Genes#####
#Fpqm
exp.train.pamr.fpqm = pamr.train(list(x=as.matrix(t(req.dfs$fpqm)), 
                                     y=stages.levels))
exp.cv.pamr.fpqm = pamr.cv(exp.train.pamr.fpqm, list(x=as.matrix(t(req.dfs$fpqm)),
                                                   y=stages.levels))
pamr.plotcv(exp.cv.pamr.fpqm)
pamr.confusion(exp.cv.pamr.fpqm, exp.cv.pamr.fpqm$threshold[13])
multiclass.roc(stages.levels, ordered(exp.cv.pamr.fpqm$yhat[,13]))

#VST
exp.train.pamr.vst = pamr.train(list(x=t(as.matrix(req.dfs$vs)), 
                                     y=stages.levels))
exp.cv.pamr.vst = pamr.cv(exp.train.pamr.vst, list(x=as.matrix(t(req.dfs$vs)),
                                                    y=stages.levels))
pamr.plotcv(exp.cv.pamr.vst)
pamr.confusion(exp.cv.pamr.vst, exp.cv.pamr.vst$threshold[19])
multiclass.roc(stages.levels, ordered(exp.cv.pamr.vst$yhat[,19]))

###Diff Genes 1######
#Fpqm
exp.train.pamr.fpqm.diff1 = pamr.train(list(x=as.matrix(t(req.dfs$fpqm[,diff.genes[[1]]])), 
                                      y=stages.levels))
exp.cv.pamr.fpqm.diff1 = pamr.cv(exp.train.pamr.fpqm.diff1, list(x=as.matrix(t(req.dfs$fpqm[,diff.genes[[1]]])),
                                                     y=stages.levels))
exp.cv.pamr.fpqm.diff1$threshold
pamr.plotcv(exp.cv.pamr.fpqm.diff1)
pamr.confusion(exp.cv.pamr.fpqm.diff1, exp.cv.pamr.fpqm.diff1$threshold[11])
plot.auc(exp.cv.pamr.fpqm.diff1,stages.levels)

#Vst
exp.train.pamr.vst.diff1 = pamr.train(list(x=t(as.matrix(req.dfs$vs[,diff.genes$`1`])), 
                                     y=stages.levels))
exp.cv.pamr.vst.diff1 = pamr.cv(exp.train.pamr.vst.diff1, list(x=t(as.matrix(req.dfs$vs[,diff.genes$`1`])),
                                                   y=stages.levels))
pamr.plotcv(exp.cv.pamr.vst.diff1)
pamr.confusion(exp.cv.pamr.vst.diff1, exp.cv.pamr.vst.diff1$threshold[13])
multiclass.roc(stages.levels, ordered(exp.cv.pamr.vst.diff1$yhat[,19]))


diff.genes[["0.5"]] = remove.dots(rownames(res)[abs(res$log2FoldChange) > 0.5 &  res$padj < 0.01])
intersect(unlist(genes.list$`stage i`), diff.genes$`0.5`)


########2 stages combined###############
