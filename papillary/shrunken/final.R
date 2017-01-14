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

confusion.mat[['shrunken']][['cv']] <- list()
confusion.mat[['shrunken']][['test']] <- list()
eval.mat[['shrunken']][['cv']] <- list()
eval.mat[['shrunken']][['test']] <- list()
eval.mat$shrunken$cv[['auc']] <- list()
eval.mat$shrunken$cv[['other']] <- list()
eval.mat$shrunken$test[['auc']] <- list()
eval.mat$shrunken$test[['other']] <- list()


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
 stages <- pamr.predict(fit = pamr.train.comb[[i]], 
                              newx = as.matrix(t(req.dfs$vs[test.ind,])),
                        threshold = pamr.cv.comb[[i]]$threshold[thr.ind])
 #print(table(stages))
 #pamr.predicted.comb[[i]] <- create.mat.error(table(stages.levels.comb[test.ind], stages))
 confusion.mat$shrunken$cv[[i]] <- create.mat.error(table(stages.levels.comb[train.ind], pamr.cv.comb[[i]]$yhat[,thr.ind]))
 eval.mat$shrunken$cv[['auc']][[i]] <- pamr.aucs.comb[[i]][thr.ind]
 eval.mat$shrunken$cv[['other']][[i]] <- confusionMatrix(ordered(pamr.cv.comb[[i]]$yhat[,thr.ind]),
                                                    stages.levels.comb[train.ind])
 confusion.mat$shrunken$test[[i]] <- create.mat.error(confusion.mat$shrunken$cv[[i]])
 eval.mat$shrunken$test[['auc']][[i]] <- multiclass.roc(stages.levels.comb[test.ind], ordered(stages))$auc
 eval.mat$shrunken$test[['other']][[i]] <- confusionMatrix(stages, stages.levels.comb[test.ind])
 remove(stages)
 remove(train.ind)
 remove(test.ind)
 remove(thr.ind)
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
confusion.mat$rf = list()
confusion.mat$rf$cv = list()
confusion.mat$rf$test = list()
eval.mat$rf = list()
eval.mat$rf$cv = list()
eval.mat$rf$test = list()
eval.mat$rf$cv[['auc']] <- list()
eval.mat$rf$cv[['other']] <- list()
eval.mat$rf$test[['auc']] <- list()
eval.mat$rf$test[['other']] <- list()
for(i in seq_along(gr))
{
  train.ind <- sort(unlist(gr[-i]))
  test.ind <- sort(unlist(gr[i]))
  
  rf <- randomForest(req.dfs$vs[train.ind,g2], stages.levels.comb[train.ind])
  confusion.mat$rf$cv[[i]] <- rf$confusion
  eval.mat$rf$cv$auc[[i]] <- multiclass.roc(stages.levels.comb[train.ind],ordered(rf$predicted))$auc
  eval.mat$rf$cv$other[[i]] <- confusionMatrix(rf$predicted, stages.levels.comb[train.ind])
  
  pr <- predict(rf, req.dfs$vs[test.ind,g2])
  confusion.mat$rf$test[[i]] <- create.mat.error(table(stages.levels.comb[test.ind], pr))
  eval.mat$rf$test$auc[[i]] <- multiclass.roc(stages.levels.comb[test.ind],ordered(pr))$auc
  eval.mat$rf$test$other[[i]] <- confusionMatrix(pr, stages.levels.comb[test.ind])
  
  remove(pr, rf, test.ind, train.ind)
}

####SVMs
confusion.mat$svm = list()
confusion.mat$svm$cv = list()
confusion.mat$svm$test = list()
eval.mat$svm = list()
eval.mat$svm$cv = list()
eval.mat$svm$test = list()
eval.mat$svm$cv[['auc']] <- list()
eval.mat$svm$cv[['other']] <- list()
eval.mat$svm$test[['auc']] <- list()
eval.mat$svm$test[['other']] <- list()
for(i in seq_along(gr))
{
  train.ind <- sort(unlist(gr[-i]))
  test.ind <- sort(unlist(gr[i]))
  
  svm.model <- svm(req.dfs$vs[train.ind,g2], stages.levels.comb[train.ind])
  #print(svm.model$levels)
  cv.pred <- cv.svm(req.dfs$vs[train.ind,g2], 10, stages.levels.comb[train.ind])
  #print(cv.pred)
  confusion.mat$svm$cv[[i]] <- create.mat.error(table(stages.levels.comb[train.ind], cv.pred))
  eval.mat$svm$cv$auc[[i]] <- multiclass.roc(stages.levels.comb[train.ind], ordered(cv.pred))$auc
  #print(levels(stages.levels.comb[train.ind]))
  eval.mat$svm$cv$other[[i]] <- confusionMatrix(cv.pred, stages.levels.comb[train.ind])
  
  pr <- predict(svm.model, req.dfs$vs[test.ind,g2])
  
  confusion.mat$svm$test[[i]] <- create.mat.error(table(stages.levels.comb[test.ind], pr))
  eval.mat$svm$test$auc[[i]] <- multiclass.roc(stages.levels.comb[test.ind], ordered(pr))$auc
  eval.mat$svm$test$other[[i]] <- confusionMatrix(pr, stages.levels.comb[test.ind])
  remove(svm.model, test.ind, train.ind, pr, cv.pred)
}

svm.model <- svm(req.dfs$vs[sort(unlist(gr[-4])),g1], stages.levels.comb[sort(unlist(gr[-4]))])
pr <- predict(svm.model, req.dfs$vs[sort(unlist(gr[4])),g1])
table(stages.levels.comb[sort(unlist(gr[4]))], pr)
####KNN
knn