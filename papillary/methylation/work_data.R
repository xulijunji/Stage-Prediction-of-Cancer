source('main/updated/initialisation.R')
source('~/Honours/Stage-Prediction-of-Cancer/papillary/main/updated/feature_sel.R')
load('methylation/net_fea_cpgss.RData')
load('environment/methylation/pheno_meth_tum_rep.RData')

sh.fea <- get.feature.shrunken(t(assay(meth.tum.rep.data)), stages.levels.comb,
                               list(train.trial.ind), 2)


fea.list <- get.class.fea(list(sh.fea))
length(intersect(rownames(meth.tum.rep.data)[rowSds(assay(meth.tum.rep.data)) < 0.05], 
                 sh.fea$genes[[1]]))


comb.data <- cbind(vst_tumor_tum, t(assay(meth.tum.rep.data)))
fea.comb <- c(net.features.trial$shrunken$atleast_3, net.cps$shrunken$atleast_3)
train.trial.model <- get.train.model(comb.data, train.trial.ind, list(sh=list(fea.comb)), stages.levels.comb)
tmp.train <- pamr.train()
cv.trial.model <- get.cv.model(comb.data, train.trial.ind, list(sh=list(fea.comb)),
                               stages.levels.comb, tr.model = train.trial.model)
test.pred.trial <- get.test.pred(comb.data, train.trial.ind, test.trial.ind, list(sh=list(fea.comb)), 
                                 stages.levels.comb, train.trial.model, cv.trial.model)
cv.trial.results <- get.cv.res(train.trial.ind, cv.trial.model, stages.levels.comb, 'sh')
test.trial.results <- get.test.res(test.pred.trial, test.trial.ind, stages.levels.comb, 'sh')

get.aucs(test.trial.results$sh$)

net.cps <- list()
net.cps[['shrunken']] <- get.feature.shrunken(t(assay(meth.tum.rep.data)), stages.levels.comb,
                                              gr.trial.train)
net.cps[['varselrf']] <- get.feature.varSelRf(t(assay(meth.tum.rep.data)), stages.levels.comb,
                                              gr.trial.train)
save(net.cps, file = 'methylation/net_fea_cpgss.')
nrow(meth.tum.rep.data)


pca.beta <- prcomp(assay(meth.tum.rep.data[]))
plotPCA(assay(meth.tum.rep.data[net.cps$varSelRF$atleast_1,]), intgroup = 'stage', title = 'ss', colData = data.frame(stage = stages.levels.comb))
library(ggplot2)
ggplot(pca.beta$rotation, aes = c(1,2))
