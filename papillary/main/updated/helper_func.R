library(pamr)
library(pROC)
library(caret)
library(randomForest)
library(e1071)
library(DESeq2)
library(caret)
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

###Cross validation Models
cv.shrunken <- function(data, folds, genes.list, train.models.list, train.ind, stages.levels)
{
  cv.models <- lapply(seq(length(genes.list)), function(i)
  {
    pamr.cv(train.models.list[[i]], list(x = t(as.matrix(data[train.ind, genes.list[[i]]])),
                              y = stages.levels[train.ind]), folds)
  })
  names(cv.models) <- names(genes.list)
  return(cv.models)
}


#####Results
shrunken.cv.results <- function(cv.models, stages.levels, train.ind)
{
  thr <- sapply(cv.models, function(model)
  {
    pamr.aucs.comb <- sapply(seq_along(model$threshold), function(x)
    {
      multiclass.roc(stages.levels[train.ind], ordered(model$yhat[,x]))$auc
    })
    which.max(pamr.aucs.comb)
  })
  names(thr) <- names(cv.models)
  
  eval <- lapply(seq_along(cv.models), function(i)
    {
    get.eval(stages.levels[train.ind], cv.models[[i]]$yhat[,thr[i]])
  })
  names(eval) <- names(cv.models)
  
  res <- list()
  res[['thr']] <- thr
  res[['metric']] <- eval
  return(res)
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