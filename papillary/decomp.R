#library(pamr)
#library(pROC)
#library(caret)

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

cv.rf <- function(data, folds, stages.levels,sampsize = if (replace) nrow(data) else ceiling(.632*nrow(data)))
{
  library(randomForest)
  total.samp <- length(rownames(data))
  gr <- build.groups(total.samp, folds)
  
  if(length(gr) != folds)
    return('kat gaya')
  predicted <- factor(rep(c('stage i'), total.samp), levels = c('stage i', 'stage iv'))
  for(i in seq(folds))
  {
    train.index = sort(unlist(gr[-i]))
    test.index = sort(unlist(gr[i]))
    rf.model <- randomForest(x = data[train.index, ], y = stages.levels[train.index])
    pr <- predict(rf.model, data[test.index, ])
    predicted[test.index] <- pr
  }
  return(predicted)
}

cv.svm <- function(data, folds, stages.levels, gamma = 0, kernel = 'linear', cost =1,
                   class.weights =if(length(levels(stages.levels)) == 4) c('stage i' = 1, 'stage ii' =1, 
                                                                           'stage iii' = 1, 'stage iv' =1)
                   else c('stage i' = 1, 'stage iv' = 1))
{
  library(e1071)
  total.samp <- length(rownames(data))
  gr <- build.groups(total.samp, folds)
  if(length(gr) != folds)
    return('kat gaya')
  
  predicted <- factor(rep(c('stage i'), total.samp), levels = c('stage i', 'stage iv'))
  for(i in seq(folds))
  {
    train.index = sort(unlist(gr[-i]))
    test.index = sort(unlist(gr[i]))
    svm.model <- svm(x = data[train.index, ], y = stages.levels[train.index], kernel = kernel, gamma = gamma)
    pr <- predict(svm.model, data[test.index, ])
    predicted[test.index] <- pr
  }
  return(predicted)
}

cv.knn <- function(data, folds, stages.levels, k)
{
  library(class)
  total.samp <- length(rownames(data))
  gr <- build.groups(total.samp, folds)
  
  if(length(gr) != folds)
    return('kat gaya')
  predicted <- factor(rep(c('stage i'), total.samp), levels = c('stage i', 'stage iv'))
  for(i in seq(folds))
  {
    train.index = sort(unlist(gr[-i]))
    test.index = sort(unlist(gr[i]))
    predicted[test.index] <- knn(train = data[train.index, ], test = data[test.index,], k = k, 
                                 cl = stages.levels[train.index])
  }
  return(predicted)
}

cv.naiveBayes <- function(data, folds, stages.levels)
{
  library(e1071)
  total.samp <- length(rownames(data))
  gr <- build.groups(total.samp, folds)
  
  if(length(gr) != folds)
    return('kat gaya')
  predicted <- factor(rep(c('stage i'), total.samp), levels = c('stage i', 'stage iv'))
  for(i in seq(folds))
  {
    train.index = sort(unlist(gr[-i]))
    test.index = sort(unlist(gr[i]))
    nb.model <- naiveBayes(x = data[train.index, ], y = stages.levels[train.index])
    
    predicted[test.index] <- predict(nb.model, data[test.index,])
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

find.best.k <- function(data, folds, stages.levels)
{
  library(pROC)
  cv.pred <- lapply(seq(10), function(k){
      cv.knn(data, folds, stages.levels, k)
      })
  aucs <- sapply(seq(10), function(k)
  {
    multiclass.roc(stages.levels, ordered(cv.pred[[k]]))$auc
  })
  k <- which.max(aucs)
  return(list(k, cv.pred[[k]]))
}


get.intersect.genes <- function(genes.list, indexes)
{
  return(Reduce(intersect, genes.list[indexes]))
}

create.heatmap <- function(data, stages, genes, col = NULL, labs = NULL,
                           title)
{
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
           show_colnames = F,
           cluster_cols = F, 
           annotation_col = annotation_col,
           color = col, main = title, fontsize = 20)
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


  
