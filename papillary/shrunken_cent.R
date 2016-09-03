install.packages('~/pamr_1.55.tar.gz', repos = NULL, type = 'source')
library(pamr)
exp.train.pamr = pamr.train(list(x=as.matrix(exp_prof_tumor_reported), y=stages.levels), threshold.scale = new.scales)

fpqm.train.pamr = pamr.train(list(x=as.matrix(exp_fpqm_tumor_reported), y=stages.levels))
fpqm.cv.pamr = pamr.cv(fpqm.train.pamr, data = list(x=as.matrix(exp_fpqm_tumor_reported), y=stages.levels), nfold = 10)
pamr.plotcv(fpqm.cv.pamr)
pamr.plotcen(fpqm.train.pamr, data = list(x=as.matrix(exp_fpqm_tumor_reported), y=stages.levels), threshold = 3.21)
pamr.confusion(fpqm.train.pamr, threshold = 1)
pamr.plotcvprob(fpqm.cv.pamr, data = list(x=as.matrix(exp_fpqm_tumor_reported), y=stages.levels), threshold = 5.2)
pamr.geneplot(fpqm.train.pamr, data = list(x=as.matrix(exp_fpqm_tumor_reported), y=stages.levels), threshold = 5.2)
listgenes(fpqm.train.pamr, data = list(x=as.matrix(exp_fpqm_tumor_reported), y=stages.levels), threshold = 5.2)

fpqm.log.train.pamr = pamr.train(list(x=as.matrix(exp_fpqm_tumor_log_reported), y=stages.levels))
fpqm.nt.train.pamr = pamr.train(list(x=as.matrix(assay(only.tumor.reported$dfs$nt)), y=stages.levels))


pamr.menu(list(x=as.matrix(exp_fpqm_tumor_reported), y=stages.levels))
pamr.menu(list(x=as.matrix(assay(only.tumor.reported$dfs$nt)), y=stages.levels))

View(exp.train.pamr$centroids)
View(exp.train.pamr$yhat)
exp.train.pamr$threshold
new.scales <- pamr.adaptthresh(exp.train.pamr, ntries = 10, reduction.factor = 3)
myresults2 <- pamr.cv(exp.train.pamr, list(x=as.matrix(exp_prof_tumor_reported), y=stages.levels))

which.min(fpqm.train.pamr$errors)
fpqm.train.pamr$threshold[62]
pamr.confusion(fpqm.train.pamr, 4.43)

fpqm.train.pamr.over = pamr.train(list(x=as.matrix(exp_fpqm_tumor_reported[over.sel.genes, ]), y=stages.levels))
fpqm.cv.pamr.over = pamr.cv(fpqm.train.pamr.over, list(x=as.matrix(exp_fpqm_tumor_reported[over.sel.genes, ]), y=stages.levels))
pamr.plotcv(fpqm.cv.pamr.over)
which.min(fpqm.train.pamr.over$errors)
fpqm.train.pamr.over$threshold[15]
pamr.confusion(fpqm.train.pamr.over, 3.41)
fpqm.train.pamr.over$centroid
