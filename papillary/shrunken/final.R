library(pamr)
library(pROC)
library(caret)
library(randomForest)
library(e1071)


pamr.train.comb <- list()
pamr.cv.comb <- list()
pamr.genes.comb <- list()
pamr.predicted.comb <- list()
pamr.aucs.comb <- list()
confusion.mat <- list()
eval.mat <- list()

do.shrunken <- function(gr, data, stages.levels, confusion.mat, eval.mat)
{
  for(i in seq_along(gr))
  {
    train.ind <- sort(unlist(gr[-i]))
    test.ind <- sort(unlist(gr[i]))
    pamr.train.comb[[i]] <- pamr.train(list(x = as.matrix(t(data[train.ind,])), 
                                            y = stages.levels[train.ind]))
    pamr.cv.comb[[i]] <- pamr.cv(pamr.train.comb[[i]], 
            data = list(x = t(as.matrix(data[train.ind,])),
                    y = stages.levels.comb[train.ind]),nfold = 10)
    
   pamr.aucs.comb[[i]] <- sapply(seq_along(pamr.cv.comb[[i]]$threshold), function(x)
   {
     multiclass.roc(stages.levels[train.ind], ordered(pamr.cv.comb[[i]]$yhat[,x]))$auc
   })
   thr.ind = which.max(pamr.aucs.comb[[i]])
   
   pamr.genes.comb[[i]] <- pamr.listgene(pamr.train.comb[[i]], 
                              data = list(x=as.matrix(t(req.dfs$vs[train.ind,])),
                                             y=stages.levels[train.ind]), 
                                threshold = pamr.cv.comb[[i]]$threshold[thr.ind],
                              fitcv = pamr.cv.comb[[i]], genenames = T)[,2]
   test.pred <- pamr.predict(fit = pamr.train.comb[[i]], 
                                newx = as.matrix(t(data[test.ind,])),
                          threshold = pamr.cv.comb[[i]]$threshold[thr.ind])
   
   out.list <- get.results(actual.cv = stages.levels[train.ind], 
                           pred.cv = pamr.cv.comb[[i]]$yhat[,thr.ind], 
                           actual.test = stages.levels[test.ind], pred.test = test.pred,
                           conf.mat = confusion.mat, eval.mat = eval.mat, classifier = 'shrunken', i)

  confusion.mat <- out.list[[1]]
  eval.mat <- out.list[[2]]
  }
  return(list(confusion.mat, eval.mat, pamr.genes.comb))
}


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
do.rf <- function(gr, data, genes, stages.levels, confusion.mat, eval.mat)
{
  for(i in seq_along(gr))
  {
    train.ind <- sort(unlist(gr[-i]))
    test.ind <- sort(unlist(gr[i]))
    
    rf.train <- randomForest(data[train.ind,genes], stages.levels[train.ind])
    cv.pred <- cv.rf(data[train.ind,genes], 10, stages.levels[train.ind])
    test.pred <- predict(rf.train, data[test.ind,genes])
    out.list <- get.results(actual.cv = stages.levels[train.ind], pred.cv = cv.pred, 
                            actual.test = stages.levels[test.ind], pred.test = test.pred,
                            conf.mat = confusion.mat, eval.mat = eval.mat, classifier = 'rf', i)
    confusion.mat <- out.list[[1]]
    eval.mat <- out.list[[2]]
  }
  return(list(confusion.mat, eval.mat))
}
####SVMs
do.svm <- function(gr, data, genes, stages.levels, confusion.mat, eval.mat)
{
  for(i in seq_along(gr))
  {
    train.ind <- sort(unlist(gr[-i]))
    test.ind <- sort(unlist(gr[i]))
    
    svm.model <- svm(data[train.ind,genes], stages.levels[train.ind])
    cv.pred <- cv.svm(data[train.ind,genes], 10, stages.levels[train.ind])
    test.pred <- predict(svm.model, data[test.ind,genes])
    out.list <- get.results(actual.cv = stages.levels[train.ind], pred.cv = cv.pred, 
                            actual.test = stages.levels[test.ind], pred.test = test.pred,
                            conf.mat = confusion.mat, eval.mat = eval.mat, classifier = 'svm', i = i)
    confusion.mat <- out.list[[1]]
    eval.mat <- out.list[[2]]
  }
  return(list(confusion.mat, eval.mat))
}
####KNN
do.knn <- function(gr, data, genes, stages.levels, confusion.mat, eval.mat)
{
  for(i in seq_along(gr))
  {
    train.ind <- sort(unlist(gr[-i]))
    test.ind <- sort(unlist(gr[i]))
    ans <- find.best.k(data[train.ind, genes], 10, stages.levels[train.ind])
    k <- ans[[1]]
    print(k)
    cv.pred <- ans[[2]]
    test.pred <- knn(data[train.ind,genes], data[test.ind,genes], cl = stages.levels[train.ind], k)
    out.list <- get.results(actual.cv = stages.levels[train.ind], pred.cv = cv.pred, 
                            actual.test = stages.levels[test.ind], pred.test = test.pred,
                            conf.mat = confusion.mat, eval.mat = eval.mat, classifier = 'knn', i = i)
    confusion.mat <- out.list[[1]]
    eval.mat <- out.list[[2]]
  }
  return(list(confusion.mat, eval.mat))
}

###Naive Bayes
do.naive <- function(gr, data, genes, stages.levels, confusion.mat, eval.mat)
{
  for(i in seq_along(gr))
  {
    train.ind <- sort(unlist(gr[-i]))
    test.ind <- sort(unlist(gr[i]))
    
    nb.model <- naiveBayes(data[train.ind,genes], stages.levels[train.ind])
    cv.pred <- cv.naiveBayes(data[train.ind,genes], 10, stages.levels[train.ind])
    #print(table(stages.levels[train.ind], cv.pred))
    test.pred <- predict(nb.model, data[test.ind,genes])
    out.list <- get.results(actual.cv = stages.levels[train.ind], pred.cv = cv.pred, 
                            actual.test = stages.levels[test.ind], pred.test = test.pred,
                            conf.mat = confusion.mat, eval.mat = eval.mat, classifier = 'nb', i = i)
    confusion.mat <- out.list[[1]]
    eval.mat <- out.list[[2]]
  }
  return(list(confusion.mat, eval.mat))
}
