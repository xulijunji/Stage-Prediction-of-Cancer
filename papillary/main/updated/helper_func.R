library(pamr)
library(pROC)
library(caret)
library(randomForest)
library(e1071)
library(DESeq2)
library(caret)
library(class)
source('shrunken/pamr.listgenes.R')


###Feature Selection
get.shrunken.features <- function(gr, data, stages.levels)
{
  pamr.train.list <- list()
  pamr.cv.list <- list()
  pamr.aucs.list <- list()
  pamr.genes.list <- list()
  for(i in seq_along(gr))
  {
    train.ind <- sort(unlist(gr[-i]))
    test.ind <- sort(unlist(gr[i]))
    pamr.train.list[[i]] <- pamr.train(list(x = as.matrix(t(data[train.ind,])), 
                                            y = stages.levels[train.ind]))
    pamr.cv.list[[i]] <- pamr.cv(pamr.train.list[[i]], 
                                 data = list(x = t(as.matrix(data[train.ind,])),
                                             y = stages.levels[train.ind]),nfold = 10)
    
    pamr.aucs.list[[i]] <- sapply(seq_along(pamr.cv.list[[i]]$threshold), function(x)
    {
      multiclass.roc(stages.levels[train.ind], ordered(pamr.cv.list[[i]]$yhat[,x]))$auc
    })
    thr.ind = which.max(pamr.aucs.list[[i]])
    print(pamr.aucs.list[[i]][thr.ind])
    
    pamr.genes.list[[i]] <- pamr.listgene(pamr.train.list[[i]], 
                                          data = list(x=as.matrix(t(data[train.ind,])),
                                                      y=stages.levels[train.ind]), 
                                          threshold = pamr.cv.list[[i]]$threshold[thr.ind],
                                          fitcv = pamr.cv.list[[i]], genenames = T)
  }
  req.list <- list()
  req.list[['genes']] <- pamr.genes.list
  req.list[['train']] <- pamr.train.list
  req.list[['cv']] <- pamr.cv.list
  return(req.list)
}


###Classifiers training
build.rf.classifier <- function(data, train.ind, genes.list, stages.levels)
{
  rf.train.list <- lapply(genes.list, function(genes)
  {
    randomForest(data[train.ind, genes], stages.levels[train.ind])  
  })
  names(rf.train.list) = names(genes.list)
  return(rf.train.list)
}

build.shrunken.classifier <- function(data, train.ind, genes.list, stages.levels)
{
  shrunken.train.list <- lapply(genes.list, function(genes)
    {
    pamr.train(data = list(x = as.matrix(t(data[train.indexes, genes])), y = stages.levels[train.ind]))
  })
  names(shrunken.train.list) <- names(genes.list)
  return(shrunken.train.list)
}

build.svm.classifier <- function(data, train.ind, genes.list,
                                  stages.levels, gamma = 0, kernel = 'linear', cost =1,
                                  class.weights =if(length(levels(stages.levels)) == 4) 
                                c('stage i' = 1, 'stage ii' =1, 'stage iii' = 1, 'stage iv' =1)
                                  else c('stage i' = 1, 'stage iv' = 1))
{
  svm.train.list <- lapply(genes.list, function(genes)
  {
    svm.model <- svm(x = data[train.ind, genes], y = stages.levels[train.ind], kernel = kernel, gamma = gamma)
  })
  names(svm.train.list) <- names(genes.list)
  return(svm.train.list)
}

build.nb.classifier <- function(data, train.ind, genes.list, stages.levels)
{
  nb.train.list <- lapply(seq_along(genes.list), function(i){
    naiveBayes(data[train.ind, genes.list[[i]]], stages.levels[train.ind])
  })
  names(nb.train.list) <- names(genes.list)
  return(nb.train.list)
}

###Cross validation Models
cv.shrunken <- function(data, folds, genes.list, train.models.list, train.ind, stages.levels)
{
  cv.models <- lapply(seq_along(genes.list), function(i)
  {
    pamr.cv(train.models.list[[i]], list(x = t(as.matrix(data[train.ind, genes.list[[i]]])),
                              y = stages.levels[train.ind]), folds)
  })
  names(cv.models) <- names(genes.list)
  
  thr <- find.best.shrunken.threshold(cv.models, stages.levels[train.ind])
  
  pred <- lapply(seq_along(cv.models), function(i)
    {
    cv.models[[i]]$yhat[,thr[i]]
  })
  names(pred) <- names(genes.list)
  
  res <- list()
  res[['cv.models']] <- cv.models
  res[['thr']] <- thr
  res[['pred']] <- pred
  return(res)
}

cv.svm.list <- function(data, folds, genes.list, train.ind, stages.levels,
                        gamma = 0, kernel = 'linear', cost =1,
                        class.weights =if(length(levels(stages.levels)) == 4) 
                          c('stage i' = 1, 'stage ii' =1, 'stage iii' = 1, 'stage iv' =1)
                        else c('stage i' = 1, 'stage iv' = 1))
{
  svm.pred <- lapply(genes.list, function(genes)
    {
    cv.svm(data[train.ind, genes], folds, stages.levels[train.ind])
  })
  names(svm.pred) <- names(genes.list)
  return(svm.pred)
}

cv.rf.list <- function(data, folds, genes.list, train.ind, stages.levels,
                       sampsize = if (replace) nrow(data)
                       else ceiling(.632*nrow(data)))
{
  rf.pred <- lapply(genes.list, function(genes)
    {
    cv.rf(data[train.ind, genes], folds, stages.levels[train.ind])
  })
  names(rf.pred) <- names(genes.list)
  return(rf.pred)
}

cv.nb.list <- function(data, folds, genes.list, train.ind, stages.levels)
{
  nb.pred <- lapply(genes.list, function(genes)
  {
    cv.naiveBayes(data[train.ind, genes], folds, stages.levels[train.ind])
  })
  names(nb.pred) <- names(genes.list)
  return(nb.pred)
}

cv.knn.list <- function(data, folds, genes.list, train.ind, stages.levels)
{
  knn.pred <- lapply(genes.list, function(genes)
  {
    find.best.k(data[train.ind, genes], folds, stages.levels[train.ind])
  })
  names(knn.pred) <- names(genes.list)
  return(knn.pred)
}

###Prediction Test
predict.shrunken <- function(train.models, genes.list, test.data, thresholds.list)
{
  if(length(train.models) != length(genes.list))
    print("havoc")
  shrunken.pred <- lapply(seq_along(train.models), function(i)
    {
    pamr.predict(train.models[[i]], t(test.data[,genes.list[[i]]]), 
                 threshold = train.models[[i]]$threshold[thresholds.list[i]])
  })
  names(shrunken.pred) <- names(train.models)
  return(shrunken.pred)
}

predict.knn <- function(cv.models, genes.list, train.data, test.data, stages.train)
{
  if(length(cv.models) != length(genes.list))
    print("havoc")
  pred.knn <- lapply(seq_along(cv.models), function(i)
  {
    knn(train.data[,genes.list[[i]]], test.data[,genes.list[[i]]], stages.train, 
        cv.models[[i]]$best_k)
  })
  names(pred.knn) <- names(cv.models)
  return(pred.knn)
}

predict.model <- function(train.models, genes.list, test.data)
{
  if(length(train.models) != length(genes.list))
    print("havoc")
  pred <- lapply(seq_along(train.models), function(i)
  {
    predict(train.models[[i]], test.data[,genes.list[[i]]])
  })
  names(pred) <- names(train.models)
  return(pred)
}

#####Results
get.eval.list <- function(actual.stages, predict.list)
{
  eval.list <- lapply(predict.list, function(predicted)
    {
    get.eval(actual.stages, predicted)
  })
}
get.eval <- function(actual.stages, pred.stages)
{
  res <- list()
  res[['conf_mat']] <- create.mat.error(table(actual.stages, pred.stages))
  res[['eval']] <- list()
  res[['eval']][['auc']] <- multiclass.roc(actual.stages, ordered(pred.stages))$auc
  res[['eval']][['other']] <- confusionMatrix(pred.stages, actual.stages)
  return(res)
}
# shrunk.pred <- function(data, test.ind, genes.list, stages.levels, train.model.list)
# {
#   pred.list <- mapply(function(train.model, genes)
#   {
#     pamr.predict(train.model, t(data[test.ind, genes]), )
#   }, train.model.list, genes.list)
# }
final.res <- function(data, train.ind, test.ind, stages.levels, genes, folds)
{
  ##Classifiers
  train.list <- list()
  train.list[['shrunken']] <- pamr.train(list(x = t(as.matrix(data[train.ind, genes])), 
                                              y = stages.levels.comb[train.ind]))
  train.list[['svm']] <-  svm(x = data[train.ind, genes], y = stages.levels.comb[train.ind])
  train.list[['nb']] <- naiveBayes(x = data[train.ind, genes], y = stages.levels.comb[train.ind])
  train.list[['rf']] <- randomForest(x = data[train.ind, genes], y = stages.levels.comb[train.ind])

  ###Predicted
  pred.cv.list <- list()
  pred.cv.list[['shrunken']] <- cv.shrunken(data, folds, list(genes), list(train.list$shrunken), 
                                             train.ind, stages.levels)
  pred.cv.list[['knn']] <- find.best.k(data[train.ind, genes], folds, stages.levels[train.ind])
  pred.cv.list[['rf']] <- cv.rf.list(data, folds, list(genes), train.ind, stages.levels)
  pred.cv.list[['svm']] <- cv.svm.list(data, folds, list(genes), train.ind, stages.levels)
  pred.cv.list[['nb']] <- cv.nb.list(data, folds, list(genes), train.ind, stages.levels)
  
  pred.test.class <- list()
  pred.test.class[['shrunken']] <- pamr.predict(train.list$shrunken, t(as.matrix(data[test.ind, genes])),
                                    threshold = train.list$shrunken$threshold[pred.cv.list$shrunken$thr])
  pred.test.class[['svm']] <- predict(train.list$svm, data[test.ind, genes])
  pred.test.class[['rf']] <- predict(train.list$rf, data[test.ind, genes])
  pred.test.class[['nb']] <- predict(train.list$nb, data[test.ind, genes])
  pred.test.class[['knn']] <- knn(data[train.ind, genes], data[test.ind, genes],
                                  stages.levels[train.ind], pred.cv.list$knn$best_k)
  
  #CV results
  #print(pred.cv.list)
  res.cv <- list()
#   print(length(stages.levels[train.ind]))
#   print(length(pred.cv.list$shrunken$pred[[1]]))
#   print(length(pred.cv.list$knn$pred))
#   print(length(pred.cv.list$svm[[1]]))
  res.cv[['shrunken']] <- get.eval(stages.levels[train.ind], pred.cv.list$shrunken$pred[[1]])
  res.cv[['knn']] <- get.eval(stages.levels[train.ind], pred.cv.list$knn$pred)
  res.cv[c('svm','rf', 'nb')] <- lapply(pred.cv.list[c('svm','rf', 'nb')], function(pred.cv){
                                       get.eval(stages.levels[train.ind], pred.cv[[1]])
                                      })
  ##Test Results
  res.test <- lapply(pred.test.class, function(pred.test){
    get.eval(stages.levels[test.ind], pred.test)
 })
 return(list(res.cv, res.test))
}
