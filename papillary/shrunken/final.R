library(pamr)
library(pROC)
library(caret)
library(randomForest)
library(e1071)

source('decomp.R')
source('shrunken/pamr.listgenes.R')
source('after_class.R')

load('environment/req_dfs.RData')
load('environment/stages.level.comb.RData')
load('environment/stage.index.RData')

gr <- build.groups(length(stages.levels.comb), 5)
pamr.train.comb <- list()
pamr.cv.comb <- list()
pamr.genes.comb <- list()
pamr.predicted.comb <- list()
pamr.aucs.comb <- list()
confusion.mat <- list()
eval.mat <- list()

for(i in seq_along(gr))
{
  train.ind <- sort(unlist(gr[-i]))
  test.ind <- sort(unlist(gr[i]))
  pamr.train.comb[[i]] <- pamr.train(list(x = as.matrix(t(req.dfs$vs[train.ind,])), 
                                          y = stages.levels.comb[train.ind]))
  pamr.cv.comb[[i]] <- pamr.cv(pamr.train.comb[[i]], 
          data = list(x = t(as.matrix(req.dfs$vs[train.ind,])),
                  y = stages.levels.comb[train.ind]),nfold = 10)
  
 pamr.aucs.comb[[i]] <- sapply(seq_along(pamr.cv.comb[[i]]$threshold), function(x)
 {
   multiclass.roc(stages.levels.comb[train.ind], ordered(pamr.cv.comb[[i]]$yhat[,x]))$auc
 })
 thr.ind = which.max(pamr.aucs.comb[[i]])
 
 pamr.genes.comb[[i]] <- pamr.listgene(pamr.train.comb[[i]], 
                            data = list(x=as.matrix(t(req.dfs$vs[train.ind,])),
                                           y=stages.levels.comb[train.ind]), 
                              threshold = pamr.cv.comb[[i]]$threshold[thr.ind],
                            fitcv = pamr.cv.comb[[i]], genenames = T)[,2]
 test.pred <- pamr.predict(fit = pamr.train.comb[[i]], 
                              newx = as.matrix(t(req.dfs$vs[test.ind,])),
                        threshold = pamr.cv.comb[[i]]$threshold[thr.ind])
 
 out.list <- get.results(actual.cv = stages.levels.comb[train.ind], 
                         pred.cv = pamr.cv.comb[[i]]$yhat[,thr.ind], 
                         actual.test = stages.levels.comb[test.ind], pred.test = test.pred,
                         conf.mat = confusion.mat, eval.mat = eval.mat, classifier = 'shrunken', i)
 confusion.mat <- out.list[[1]]
 eval.mat <- out.list[[2]]
 
 remove(test.pred, train.ind,test.ind,thr.ind, out.list)
}
g1 <- Reduce(intersect, pamr.genes.comb[c(1,2,4,3)])
g2 <- Reduce(intersect, pamr.genes.comb)
pamr.plotcv(pamr.cv.comb[[1]])
sapply(pamr.aucs.comb, max)
pamr.predicted.comb[[4]]

#######Trial########
g <- pamr.listgene(pamr.train.comb[[1]], 
                   data = list(x=as.matrix(t(req.dfs$vs[unlist(gr[-1]),])),
                               y=stages.levels.comb[unlist(gr[-1])]), 
                   threshold = pamr.cv.comb[[1]]$threshold[15],
                   fitcv = pamr.cv.comb[[i]], genenames = T)[,2]
stages <- pamr.predict(pamr.train.comb[[2]], 
                       as.matrix(t(req.dfs$vs[gr[[2]],])),
                       pamr.cv.comb[[2]]$threshold[17])
table(stages)
table(stages.levels.comb[gr[[2]]], stages)
pamr.confusion(pamr.cv.comb[[1]], threshold = pamr.cv.comb[[1]]$threshold[17])
plot(pamr.aucs.comb[[1]])
##########

####Random Forests
for(i in seq_along(gr))
{
  train.ind <- sort(unlist(gr[-i]))
  test.ind <- sort(unlist(gr[i]))
  
  rf.train <- randomForest(req.dfs$vs[train.ind,g2], stages.levels.comb[train.ind])
  cv.pred <- cv.rf(req.dfs$vs[train.ind,g2], 10, stages.levels.comb[train.ind])
  test.pred <- predict(rf.train, req.dfs$vs[test.ind,g2])
  out.list <- get.results(actual.cv = stages.levels.comb[train.ind], pred.cv = cv.pred, 
                          actual.test = stages.levels.comb[test.ind], pred.test = test.pred,
                          conf.mat = confusion.mat, eval.mat = eval.mat, classifier = 'rf', i)
  confusion.mat <- out.list[[1]]
  eval.mat <- out.list[[2]]
  remove(cv.pred,test.pred, rf.train, test.ind, train.ind, out.list)
}

####SVMs
for(i in seq_along(gr))
{
  train.ind <- sort(unlist(gr[-i]))
  test.ind <- sort(unlist(gr[i]))
  
  svm.model <- svm(req.dfs$vs[train.ind,g2], stages.levels.comb[train.ind])
  cv.pred <- cv.svm(req.dfs$vs[train.ind,g2], 10, stages.levels.comb[train.ind])
  test.pred <- predict(svm.model, req.dfs$vs[test.ind,g2])
  out.list <- get.results(actual.cv = stages.levels.comb[train.ind], pred.cv = cv.pred, 
                          actual.test = stages.levels.comb[test.ind], pred.test = test.pred,
                          conf.mat = confusion.mat, eval.mat = eval.mat, classifier = 'svm', i = i)
  confusion.mat <- out.list[[1]]
  eval.mat <- out.list[[2]]
  remove(svm.model, test.ind, train.ind, test.pred, cv.pred, out.list)
}

####KNN
knn



