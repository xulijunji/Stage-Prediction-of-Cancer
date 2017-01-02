load('environment/only_tumor_reported.RData')
load('environment/stages_levels.RData')
load('environment/stages.level.comb.RData')
load('environment/diff_genes.RData')

source('shrunken/pamr.listgenes.R')
library(DESeq2)
library(pROC)

install.packages('~/pamr_1.55.tar.gz', repos = NULL, type = 'source')
library(pamr)
exp.train.pamr = pamr.train(list(x=as.matrix(exp_prof_tumor_reported), y=stages.levels), threshold.scale = new.scales)

fpqm.train.pamr = pamr.train(list(x=as.matrix(only.tumor.reported$dfs$fpqm), y=stages.levels))
fpqm.cv.pamr = pamr.cv(fpqm.train.pamr, data = list(x=as.matrix(only.tumor.reported$dfs$fpqm), y=stages.levels), nfold = 10)
pamr.plotcv(fpqm.cv.pamr)
pamr.plotcen(fpqm.train.pamr, data = list(x=as.matrix(only.tumor.reported$dfs$fpqm), y=stages.levels), threshold = 3.21)
pamr.confusion(fpqm.train.pamr, threshold = 4.3)
pamr.plotcvprob(fpqm.cv.pamr, data = list(x=as.matrix(exp_fpqm_tumor_reported), y=stages.levels), threshold = 5.2)
pamr.geneplot(fpqm.train.pamr, data = list(x=as.matrix(exp_fpqm_tumor_reported), y=stages.levels), threshold = 5.2)
shrunken.genes <- pamr.listgene(fpqm.train.pamr, data = list(x=as.matrix(only.tumor.reported$dfs$fpqm), y=stages.levels), threshold = 3)


#######vst
exp.train.pamr.vst = pamr.train(list(x=as.matrix(t(req.dfs$fpqm)), 
                                     y=stages.levels))
exp.cv.pamr.vst = pamr.cv(exp.train.pamr.vst, list(x=as.matrix(t(req.dfs$fpqm)),
                                                   y=stages.levels))
pamr.plotcv(exp.cv.pamr.vst)
pamr.confusion(exp.cv.pamr.vst, threshold = 3.8)
table(stages.levels.comb, exp.cv.pamr.vst$yhat[,20])
typeof(pamr.confusion(exp.cv.pamr.vst,3.71))

sensitivity(table(stages.levels.comb, stages.levels.comb))
a = confusionMatrix(exp.cv.pamr.vst$yhat[,20],stages.levels) 
multiclass.roc(stages.levels, ordered(exp.cv.pamr.vst$yhat[,17]))
##Note that Caret we write the transpose of confusion matrix compared to that of others

# fpqm.log.train.pamr = pamr.train(list(x=as.matrix(exp_fpqm_tumor_log_reported), y=stages.levels))
# fpqm.nt.train.pamr = pamr.train(list(x=as.matrix(assay(only.tumor.reported$dfs$nt)), y=stages.levels))
# 
# 
# pamr.menu(list(x=as.matrix(exp_fpqm_tumor_reported), y=stages.levels))
# pamr.menu(list(x=as.matrix(assay(only.tumor.reported$dfs$nt)), y=stages.levels))
# 
# View(exp.train.pamr$centroids)
# View(exp.train.pamr$yhat)
# exp.train.pamr$threshold
# new.scales <- pamr.adaptthresh(exp.train.pamr, ntries = 10, reduction.factor = 3)
# myresults2 <- pamr.cv(exp.train.pamr, list(x=as.matrix(exp_prof_tumor_reported), y=stages.levels))
# 
# which.min(fpqm.train.pamr$errors)
# fpqm.train.pamr$threshold[62]
# pamr.confusion(fpqm.train.pamr, 4.43)
# 
# fpqm.train.pamr.over = pamr.train(list(x=as.matrix(exp_fpqm_tumor_reported[over.sel.genes, ]), y=stages.levels))
# fpqm.cv.pamr.over = pamr.cv(fpqm.train.pamr.over, list(x=as.matrix(exp_fpqm_tumor_reported[over.sel.genes, ]), y=stages.levels))
# pamr.plotcv(fpqm.cv.pamr.over)
# which.min(fpqm.train.pamr.over$errors)
# fpqm.train.pamr.over$threshold[15]
# pamr.confusion(fpqm.train.pamr.over, 3.41)
# fpqm.train.pamr.over$centroid


#####On combined data set##############
stages.level.3way <- stages.levels
stages.level.3way[stages.level.3way == 'stage ii'] = 'stage i'
stages.level.3way <- droplevels(stages.level.3way)

fpqm.train.pamr.comb <- pamr.train(list(x = as.matrix(only.tumor.reported$dfs$fpqm), y = stages.levels.comb))
fpqm.train.pamr.comb$threshold.scale
fpqm.cv.pamr.comb = pamr.cv(fpqm.train.pamr.comb, data = list(x=as.matrix(only.tumor.reported$dfs$fpqm), y=stages.levels.comb), nfold = 10)
pamr.plotcv(fpqm.cv.pamr.comb)
pamr.plotcen(fpqm.train.pamr.comb, data = list(x=as.matrix(only.tumor.reported$dfs$fpqm), y=stages.levels.comb), threshold = 4)
pamr.confusion(fpqm.train.pamr.comb, threshold = 4.3)
roc(stages.levels.comb, ordered(fpqm.cv.pamr.comb$yhat[,15]))
confusionMatrix(fpqm.cv.pamr.comb$yhat[,15], stages.levels.comb)
#pamr.plotcvprob(fpqm.cv.pamr, data = list(x=as.matrix(exp_fpqm_tumor_reported), y=stages.levels), threshold = 5.2)
pamr.geneplot(fpqm.train.pamr.comb, data = list(x=as.matrix(only.tumor.reported$dfs$fpqm), y=stages.levels.comb), threshold = 3.8)
shrunken.genes <- pamr.listgene(fpqm.train.pamr.comb, data = list(x=as.matrix(only.tumor.reported$dfs$fpqm), y=stages.levels.comb), 
                                threshold = 4.3, fitcv = fpqm.cv.pamr.comb, genenames = T)


vt.train.pamr.comb <- pamr.train(list(x = as.matrix(assay(only.tumor.reported$dfs$vs)), 
                                      y = stages.levels.comb))
vt.cv.pamr.comb = pamr.cv(vt.train.pamr.comb,
                          data = list(x=as.matrix(assay(only.tumor.reported$dfs$vs)),
                                      y=stages.levels.comb), nfold = 10)
pamr.plotcv(vt.cv.pamr.comb)
pamr.confusion(vt.train.pamr.comb, threshold = 4.3)
###Using only 1 fold genes######
fpqm.train.pamr.1fold <- pamr.train(list(x = as.matrix(only.tumor.reported$dfs$fpqm), y = stages.levels),
                                    gene.subset = match(diff.genes[[1]], rownames(only.tumor.reported$dfs$fpqm)))
fpqm.cv.pamr.1fold <- pamr.cv(fpqm.train.pamr.1fold, list(x = as.matrix(only.tumor.reported$dfs$fpqm), y = stages.levels)
                                 )
pamr.plotcv(fpqm.cv.pamr.1fold)
pamr.confusion(fpqm.train.pamr.1fold, threshold = 3)

fpqm.train.pamr.5fold <- pamr.train(list(x = as.matrix(only.tumor.reported$dfs$fpqm[diff.genes[[5]],]), y = stages.levels)
                                 )
fpqm.cv.pamr.5fold <- pamr.cv(fpqm.train.pamr.5fold, list(x = as.matrix(only.tumor.reported$dfs$fpqm[diff.genes[[5]],]),
                                                          y = stages.levels))
pamr.plotcv(fpqm.cv.pamr.5fold)

fpqm.train.pamr.shrunken <-  pamr.train(list(x = as.matrix(only.tumor.reported$dfs$fpqm[as.numeric(shrunken.genes[,1]),]), y = stages.levels))
fpqm.train.cv.shrunken <- pamr.cv(fpqm.train.pamr.shrunken, list(x = as.matrix(only.tumor.reported$dfs$fpqm[as.numeric(shrunken.genes[,1]),]),
                                                                 y = stages.levels))
pamr.plotcv(fpqm.train.cv.shrunken)
pamr.confusion(fpqm.train.pamr.shrunken, threshold = 2.5)

fpqm.diff.train.comb.1fold <- pamr.train(list(x = as.matrix(only.tumor.reported$dfs$fpqm[diff.genes[[1]],]), y = stages.levels.comb)
                                         )
fpqm.cv.pamr.1fold.comb <- pamr.cv(fpqm.diff.train.comb.1fold, list(x = as.matrix(only.tumor.reported$dfs$fpqm), y = stages.levels.comb))
pamr.plotcv(fpqm.cv.pamr.1fold.comb)                                   
pamr.confusion(fpqm.diff.train.comb.1fold, threshold = 2.5)
shrunken.genes.1fold <- pamr.listgene(fpqm.diff.train.comb.1fold, data = list(x=as.matrix(only.tumor.reported$dfs$fpqm[diff.genes[[1]],]),
                                                                              y=stages.levels.comb),genenames = T, threshold = 2.5  )
fpqm.tain.pamr.3way <- pamr.train(list(x = as.matrix(only.tumor.reported$dfs$fpqm), y = stages.level.3way))
fpqm.train.cv.3way <- pamr.cv(fpqm.tain.pamr.3way, list(x = as.matrix(only.tumor.reported$dfs$fpqm),
                                                                 y = stages.level.3way))
pamr.plotcv(fpqm.train.cv.3way)
pamr.confusion(fpqm.tain.pamr.3way, threshold = 2)
