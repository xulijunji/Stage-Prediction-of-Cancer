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

cv.svm <- function(data, folds, stages.levels, gamma = 0, kernel = 'linear', cost =1,
                   class.weights =if(length(levels(stages.levels)) == 4) c('stage i' = 1, 'stage ii' =1, 
                                                                           'stage iii' = 1, 'stage iv' =1)
                   else c('stage i' = 1, 'stage iv' = 1))
{
  total.samp <- length(rownames(data))
  gr <- build.groups(total.samp, folds)
  if(length(gr) != folds)
    return('kat gaya')
  
  predicted <- rep(c('stage i'), total.samp/2)
  for(i in seq(folds))
  {
    train.index = sort(unlist(gr[-i]))
    test.index = sort(unlist(gr[i]))
    svm.model <- svm(x = data[train.index, ], y = stages.levels[train.index], kernel = kernel, gamma = gamma)
    #print(levels(stages.levels[train.index]))
    #print(svm.model$levels)
    pr <- predict(svm.model, data[test.index, ])
    #print(pr)
    predicted[test.index] <- pr
  }
  predicted <- factor(predicted, labels = levels(stages.levels))
  
  return(predicted)
}