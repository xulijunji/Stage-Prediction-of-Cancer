install.packages('~/pamr_1.55.tar.gz', repos = NULL, type = 'source')
library(pamr)
exp.train.pamr = pamr.train(list(x=as.matrix(exp_prof_tumor_reported), y=stages.levels), threshold.scale = new.scales)
fpqm.train.pamr = pamr.train(list(x=as.matrix(exp_fpqm_tumor_reported), y=stages.levels))
fpqm.log.train.pamr = pamr.train(list(x=as.matrix(exp_fpqm_tumor_log_reported), y=stages.levels))
pamr.menu(list(x=as.matrix(exp_prof_tumor_reported), y=stages.levels))
View(exp.train.pamr$centroids)
View(exp.train.pamr$yhat)
exp.train.pamr$threshold
new.scales <- pamr.adaptthresh(exp.train.pamr, ntries = 10, reduction.factor = 3)
myresults2 <- pamr.cv(exp.train.pamr, list(x=as.matrix(exp_prof_tumor_reported), y=stages.levels))

