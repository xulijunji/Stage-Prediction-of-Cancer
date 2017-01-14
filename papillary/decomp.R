library(pamr)

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