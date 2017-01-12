load('environment/only_tumor_reported.RData')
load('environment/stages_levels.RData')
load('environment/stages.level.comb.RData')
load('environment/diff_genes.RData')

source('shrunken/pamr.listgenes.R')
source('after_class.R')
library(DESeq2)
library(pROC)
library(caret)
install.packages('~/pamr_1.55.tar.gz', repos = NULL, type = 'source')
library(pamr)


####Trial purposes
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

load('environment/only_tumor_reported.RData')
load('environment/stages_levels.RData')
load('environment/stages.level.comb.RData')
load('environment/diff_genes.RData')

source('shrunken/pamr.listgenes.R')



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

########2 stages combined###############

###Fpqm Combined
fpqm.train.pamr.comb <- pamr.train(list(x = as.matrix(only.tumor.reported$dfs$fpqm), y = stages.levels.comb))
fpqm.cv.pamr.comb = pamr.cv(fpqm.train.pamr.comb, data = list(x=as.matrix(only.tumor.reported$dfs$fpqm), 
                                                              y=stages.levels.comb), nfold = 10)
pamr.plotcv(fpqm.cv.pamr.comb)
pamr.plotcen(fpqm.train.pamr.comb, data = list(x=as.matrix(only.tumor.reported$dfs$fpqm), y=stages.levels.comb), threshold = 4)
plot.auc(fpqm.cv.pamr.comb, stages.levels.comb)
pamr.confusion(fpqm.train.pamr.comb, threshold = fpqm.cv.pamr.comb$threshold[16])
roc(stages.levels.comb, ordered(fpqm.cv.pamr.comb$yhat[,16]))
confusionMatrix(fpqm.cv.pamr.comb$yhat[,16], stages.levels.comb)
#pamr.plotcvprob(fpqm.cv.pamr, data = list(x=as.matrix(exp_fpqm_tumor_reported), y=stages.levels), threshold = 5.2)
pamr.geneplot(fpqm.train.pamr.comb, data = list(x=as.matrix(only.tumor.reported$dfs$fpqm), y=stages.levels.comb), threshold = 3.8)
shrunken.genes.fpqm.comb <- pamr.listgene(fpqm.train.pamr.comb, data = list(x=as.matrix(only.tumor.reported$dfs$fpqm), y=stages.levels.comb), 
                                threshold = fpqm.cv.pamr.comb$threshold[16],
                                fitcv = fpqm.cv.pamr.comb, genenames = T)

###Vst combined
vs.train.pamr.comb <- pamr.train(list(x = as.matrix(assay(only.tumor.reported$dfs$vs)), 
                                      y = stages.levels.comb))
vs.cv.pamr.comb = pamr.cv(vs.train.pamr.comb,
                          data = list(x=as.matrix(assay(only.tumor.reported$dfs$vs)),
                                      y=stages.levels.comb), nfold = 10)
plot.auc(vs.cv.pamr.comb, stages.levels.comb)
pamr.confusion(vs.cv.pamr.comb, threshold = vs.cv.pamr.comb$threshold[21])
roc(stages.levels.comb, ordered(vs.cv.pamr.comb$yhat[,16]))
confusionMatrix(vs.cv.pamr.comb$yhat[,16], stages.levels.comb)
shrunken.genes.vs.comb <- pamr.listgene(vs.train.pamr.comb, data = list(x=as.matrix(t(req.dfs$vs)),
                                                                        y=stages.levels.comb), 
                                          threshold = vs.cv.pamr.comb$threshold[16],
                                          fitcv = vs.cv.pamr.comb, genenames = T)

###Using only 1 fold genes######
##Fpqm
fpqm.train.pamr.comb.1fold <- pamr.train(list(x = as.matrix(only.tumor.reported$dfs$fpqm[diff.genes$`1`,]), 
                                              y = stages.levels.comb)
                                    )
fpqm.cv.pamr.comb.1fold <- pamr.cv(fpqm.train.pamr.comb.1fold, 
                                   list(x = as.matrix(only.tumor.reported$dfs$fpqm[diff.genes$`1`,]), 
                                           y = stages.levels.comb))
pamr.plotcv(fpqm.cv.pamr.comb.1fold)
plot.auc(fpqm.cv.pamr.comb.1fold, stages.levels.comb)
pamr.confusion(fpqm.cv.pamr.comb.1fold, threshold = fpqm.cv.pamr.comb.1fold$threshold[14])
roc(stages.levels.comb, ordered(fpqm.cv.pamr.comb.1fold$yhat[,14]))
confusionMatrix(fpqm.cv.pamr.comb.1fold$yhat[,14], stages.levels.comb)
shrunken.genes.fpqm.comb.1fold <- pamr.listgene(fpqm.train.pamr.comb.1fold, 
                                          data = list(x=as.matrix(t(req.dfs$fpqm)),
                                                y=stages.levels.comb), 
                                        threshold = fpqm.cv.pamr.comb.1fold$threshold[14],
                                        fitcv = fpqm.cv.pamr.comb.1fold, genenames = T)

##VS
vs.train.pamr.comb.1fold <- pamr.train(list(x = as.matrix(assay(only.tumor.reported$dfs$vs[diff.genes[[1]],])), 
                                      y = stages.levels.comb))
vs.cv.pamr.comb.1fold = pamr.cv(vs.train.pamr.comb.1fold,
                        data = list(x=as.matrix(assay(only.tumor.reported$dfs$vs[diff.genes$`1`,])),
                                      y=stages.levels.comb), nfold = 10)
plot.auc(vs.cv.pamr.comb.1fold, stages.levels.comb)
pamr.confusion(vs.cv.pamr.comb.1fold, threshold = vs.cv.pamr.comb.1fold$threshold[16])
roc(stages.levels.comb, ordered(vs.cv.pamr.comb.1fold$yhat[,16]))
confusionMatrix(vs.cv.pamr.comb.1fold$yhat[,16], stages.levels.comb)
shrunken.genes.vs.comb.1fold <- pamr.listgene(vs.train.pamr.comb.1fold, 
                                              data = list(x=as.matrix(t(req.dfs$vs[,diff.genes$`1`])),
                                                                        y=stages.levels.comb), 
                                        threshold = vs.cv.pamr.comb.1fold$threshold[16],
                                        fitcv = vs.cv.pamr.comb.1fold, genenames = T)


###Looking at the genes
length(intersect(shrunken.genes.vs.comb[,2], shrunken.genes.vs.comb.1fold[,2]))
sum(shrunken.genes.fpqm.comb.1fold[,2] %in% shrunken.genes.fpqm.comb[,2])
length(intersect(shrunken.genes.fpqm.comb[,2], shrunken.genes.fpqm.comb.1fold[,2]))
length(intersect(shrunken.genes.fpqm.comb[,2], shrunken.genes.vs.comb[,2]))
