load('environment/accuracy_feature/updated/net_features_trial.RData')
load('environment/methylation/net_fea_cpgss.RData')

fea.meth.list <- get.class.fea(net.cps)
fea.meth.list$varSelRF$atleast_4 <- NULL
train.trial.model <- get.train.model(t(assay(meth.tum.rep.data)), train.trial.ind, fea.meth.list, stages.levels.comb)
cv.trial.model <- get.cv.model(t(assay(meth.tum.rep.data)), train.trial.ind, fea.meth.list,
                               stages.levels.comb, tr.model = train.trial.model)
test.pred.trial <- get.test.pred(t(assay(meth.tum.rep.data)), train.trial.ind, test.trial.ind, 
                                 fea.meth.list, stages.levels.comb, train.trial.model, cv.trial.model)
cv.trial.results <- get.cv.res(train.trial.ind, cv.trial.model, stages.levels.comb, names(fea.meth.list))
test.trial.results <- get.test.res(test.pred.trial, test.trial.ind, stages.levels.comb, names(fea.meth.list))

get.aucs(test.trial.results$varSelRF$knn)

sapply(fea.meth.list$shrunken, length)

tr1 <- final.res(t(assay(meth.tum.rep.data)), train.trial.ind, test.trial.ind, stages.levels.comb,
                 net.cps$shrunken$genes[[1]], 10)
get.aucs(tr1[[2]])
