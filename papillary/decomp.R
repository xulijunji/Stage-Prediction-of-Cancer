library(pamr)
library(pROC)
library(caret)
library(randomForest)
library(e1071)
source('across_tumors.R')
source('after_class.R')

build.groups <- function(total.samples, num.group)
{
  gr=NULL
  num.sel=floor(total.samples/num.group)
  num.add=total.samples%%num.group
  range=1:total.samples
  rr=total.samples
  for(i in 1:num.group)
  {
    vrem=sample(1:rr,size=num.sel)
    sel=range[vrem]
    gr=c(gr,list(sel))
    range=range[-vrem]
    rr=rr-num.sel
  }
  
  if(num.add>0)
  {
    vrem=sample(1:num.group,num.add)
    for(i in 1:num.add)
    {
      gr[[vrem[i]]]=c(gr[[vrem[i]]],range[i])
    }
  }
  gr <- lapply(gr, sort)
  return(gr)
}

get.stage.distribution <- function(gr, stages)
{
  stage.dist <- lapply(gr, function(x)
    {
    table(stages[x])
  })
  return(stage.dist)
}

check.pca.var <- function(pca.var, train.data, test.data)
{
  if(pca.var)
  {
    pca.list <- do.pca(train.data, test.data)
    train.data <- pca.list[[1]]
    test.data <- pca.list[[2]]
  }
  return(list(train.data, test.data))
}

cv.rf <- function(data, folds, stages.levels,pca.var,sampsize = if (replace) nrow(data) else ceiling(.632*nrow(data)))
{
  library(randomForest)
  total.samp <- length(rownames(data))
  gr <- build.groups(total.samp, folds)
  
  if(length(gr) != folds)
    return('kat gaya')
  classes <- levels(as.factor(stages.levels))
  predicted <- factor(rep(c(classes[1]), total.samp), levels = classes)
  for(i in seq(folds))
  {
    train.index = sort(unlist(gr[-i]))
    test.index = sort(unlist(gr[i]))
    data.list <- get.pca.var(pca.var, data[train.index,], data[test.index,])
    rf.model <- randomForest(x = data.list[[1]], y = stages.levels[train.index])
    pr <- predict(rf.model, data.list[[2]])
    predicted[test.index] <- pr
  }
  return(predicted)
}

cv.svm <- function(data, folds, stages.levels, gamma = 0, pca.var = F, kernel = 'linear', cost =1,
                   class.weights = sapply(levels(as.factor(stages.levels)), function(x){
                     c(x=1)
                   } ))
{
  names(class.weights) = levels(as.factor(stages.levels))
  total.samp <- length(rownames(data))
  gr <- build.groups(total.samp, folds)
  if(length(gr) != folds)
    return('kat gaya')
  
  classes <- levels(as.factor(stages.levels))
  predicted <- factor(rep(c(classes[1]), total.samp), levels = classes)
  for(i in seq(folds))
  {
    train.index = sort(unlist(gr[-i]))
    test.index = sort(unlist(gr[i]))
    data.list <- get.pca.var(pca.var, data[train.index, ], data[test.index, ])
    train.data = data.list[[1]]
    test.data = data.list[[2]]
    
    svm.model <- svm(x = train.data, y = stages.levels[train.index], kernel = kernel, gamma = gamma)
    pr <- predict(svm.model, test.data)
    predicted[test.index] <- pr
  }
  return(predicted)
}

cv.knn <- function(data, folds, stages.levels, k, pca.var = F)
{
  library(class)
  total.samp <- length(rownames(data))
  gr <- build.groups(total.samp, folds)
  
  if(length(gr) != folds)
    return('kat gaya')
  classes <- levels(as.factor(stages.levels))
  predicted <- factor(rep(c(classes[1]), total.samp), levels = classes)
  
  for(i in seq(folds))
  {
    train.index = sort(unlist(gr[-i]))
    test.index = sort(unlist(gr[i]))
    data.list <- get.pca.var(pca.var, data[train.index,], data[test.index,])
    predicted[test.index] <- knn(train = data.list[[1]], test = data.list[[2]], k = k, 
                                 cl = stages.levels[train.index])
  }
  return(predicted)
}

cv.naiveBayes <- function(data, folds, stages.levels, pca.var)
{
  library(e1071)
  total.samp <- length(rownames(data))
  gr <- build.groups(total.samp, folds)
  
  if(length(gr) != folds)
    return('kat gaya')
  classes <- levels(as.factor(stages.levels))
  predicted <- factor(rep(c(classes[1]), total.samp), levels = classes)
  for(i in seq(folds))
  {
    train.index = sort(unlist(gr[-i]))
    test.index = sort(unlist(gr[i]))
    data.list <- get.pca.var(pca.var, data[train.index, ], data[test.index, ], 5)
    train.data = data.list[[1]]
    test.data = data.list[[2]]
    nb.model <- naiveBayes(x = train.data, y = as.factor(stages.levels[train.index]))
    #print(predict(object = nb.model, newdata=data[train.index,]))
    predicted[test.index] <- predict(nb.model, test.data)
  }
  return(predicted)
}

get.results <- function(actual.cv, pred.cv, actual.test, pred.test, conf.mat, eval.mat, classifier, i)
{
  library(caret)
  if(i == 1)
  {
    conf.mat[[classifier]] <- list()
    
    conf.mat[[classifier]][['cv']] <- list()
    conf.mat[[classifier]][['test']] <- list()
    
    eval.mat[[classifier]][['cv']] <- list()
    eval.mat[[classifier]][['test']] <- list()
    
    eval.mat[[classifier]][['cv']][['auc']] <- list()
    eval.mat[[classifier]][['cv']][['other']] <- list()
    
    eval.mat[[classifier]][['test']][['auc']] <- list()
    eval.mat[[classifier]][['test']][['other']] <- list()
  }
  conf.mat[[classifier]][['cv']][[i]] <- create.mat.error(table(actual.cv, pred.cv))
  conf.mat[[classifier]][['test']][[i]] <- create.mat.error(table(actual.test, pred.test))
  
  eval.mat[[classifier]][['cv']][['auc']][[i]] <- multiclass.roc(actual.cv, ordered(pred.cv))$auc
  eval.mat[[classifier]][['cv']][['other']][[i]] <- confusionMatrix(pred.cv, actual.cv)
  
  eval.mat[[classifier]][['test']][['auc']][[i]] <- multiclass.roc(actual.test, ordered(pred.test))$auc
  eval.mat[[classifier]][['test']][['other']][[i]] <- confusionMatrix(pred.test, actual.test)
  return(list(conf.mat, eval.mat))
}

find.best.k <- function(data, folds, stages.levels, pca.var = F)
{
  
  library(pROC)
  cv.pred <- lapply(seq(folds), function(k){
      cv.knn(data, folds, stages.levels, k, pca.var)
      })
  
  aucs <- sapply(seq(folds), function(k)
  {
    multiclass.roc(stages.levels, ordered(cv.pred[[k]]))$auc
  })
  k <- which.max(aucs)
  res <- list()
  res[['best_k']] <- k
  res[['pred']] <- cv.pred[[k]]
  return(res)
}

find.best.shrunken.threshold <- function(cv.models, stages.levels)
{
  thr <- sapply(cv.models, function(model)
  {
    pamr.aucs.comb <- sapply(seq_along(model$threshold), function(x)
    {
      multiclass.roc(stages.levels, ordered(model$yhat[,x]))$auc
    })
    which.max(pamr.aucs.comb)
  })
  names(thr) <- names(cv.models)
 return(thr) 
}

get.intersect.genes <- function(genes.list, indexes)
{
  return(Reduce(intersect, genes.list[indexes]))
}

create.heatmap <- function(data, stages, genes, title, col = NULL, labs = NULL,
                           show_colnames = F, show_rownames = F,
                           cluster_cols = F, cluster_rows = F)
{
  #View(data)
  library(pheatmap)
  stage.ind.levels <- list()
  if(length(stages) != length(rownames(data)))
    print('damn')
  for(i in seq_along(levels(stages)))
  {
    stage = levels(stages)[i]
    stage.ind.levels[[stage]] = which(stages == stage)
  }
  row.indexes <- Reduce(union, stage.ind.levels)
  #print(row.indexes)
  annotation_col = data.frame(names = stages[row.indexes])
  rownames(annotation_col) = rownames(data)[row.indexes]
  #print(stages[row.indexes])
  if(!is.null(labs))
    labs = rownames(data)[row.indexes]
  pheatmap(t(data[row.indexes, genes]),
           labels_col = labs,
           cluster_cols = cluster_cols, 
           cluster_rows = cluster_rows,
           show_rownames = show_rownames,
           show_colnames = show_rownames,
           annotation_col = annotation_col,
           color = col, main = title, fontsize = 10)
}
  
get.genes.common <- function(genes.list, max.no.of.models)
{
  total.genes <- Reduce(union, genes.list)
  genes.req <- c()
  for(i in total.genes)
  {
     count = 0
     for(j in seq_along(genes.list))
     {
       if(i %in%  genes.list[[j]])
         count = count + 1
       if(count == max.no.of.models)
       {
         genes.req <- c( genes.req, i)
         break
       }
     }
   }
   return(genes.req)
}

get.mean <- function(class.list, ind, mode)
{
  req.values <- sapply(class.list, function(x)
    {
    x[[mode]][ind]
  })
  return(mean(req.values))
}

get.aucs <- function(acc.list)
{
  aucs <- sapply(acc.list, function(acc)
    {
    acc$eval$auc
  })
  return(aucs)
}
get.conf.mat <- function(res)
{
  conf.mat <- lapply(res, function(x)
    {
    x$conf_mat
  })
  return(conf.mat)
}
  
