library(pamr)
library(pROC)
library(caret)
library(randomForest)
library(e1071)
library(DESeq2)
source('shrunken/pamr.listgenes.R')

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
