compute.class.means <- function(data, labs)
{
  ###data - N*j, where N number of samples and j number of genes
  classes <- levels(labs)
  class.means <- matrix(nrow = length(classes), ncol = dim(data)[2])
  
  for(i in seq(1:dim(data)[2]))
  {
    class.means[,i] = sapply(classes, function(x)
    {
      indexes = which(labs == x)
      mean.gene.class = mean(data[indexes,i])
    })
  }
  return(class.means)
}

compute.var.gene <- function(data, class.means, labs, gene.index)
{
  ###data - N*j, where N number of samples and j number of genes
  classes <- levels(labs)
  vars.class.genes = c()
  vars.classes.gene <- sapply(classes, function(x)
  {
    class.index <- which(classes == x)
    indexes <- which(labs == x)
    #  print(indexes)
    var.class.gene <- (sum((data[indexes, gene.index] - class.means[class.index, gene.index])^2))/(length(indexes)-1)
    #print(var.class.gene)
  })
  return(vars.classes.gene)
}

compute.ftest <- function(data, labs)
{
  ###data - N*j, where N number of samples and j number of genes
  classes <- levels(labs)
  class.means = compute.class.means(data, labs)
  class.lengths <- table(labs)
  means.overall <- apply(class.means, 2, mean)
  total.samples <- sum(class.lengths)
  total.classes <- length(classes)
  genes.vars.classes <- sapply(c(1:dim(data)[2]), function(i)
  {
    compute.var.gene(data, class.means, labs, i)
  })
  vars.genes <- sapply(c(1:dim(data)[2]), function(i)
  {
    gene.vars <- genes.vars.classes[,i]
    net.sigms <- mapply(function(x,y)
    {
      y*(y-1)*x
    }, gene.vars, class.lengths)
    sum(net.sigms)/(total.samples - total.classes)
  })
  
  f.vals <- sapply(c(1:dim(data)[2]), function(i)
  {
    f.val <- sapply(c(1:total.classes), function(x)
    {
      class.lengths[x]*abs(class.means[x,i] - means.overall[i])
    })
    sum(f.val)/((total.classes-1)*vars.genes[i])
  })
  df = data.frame(indexes = order(f.vals, decreasing = T), fvals = f.vals[order(f.vals, decreasing = T)])
  return(df)
}

cv.svm <- function(data, folds, stages.levels, gamma = 0, kernel = 'linear')
{
  folds.indexes = createFolds(1:length(rownames(data)), folds)
  #print(folds.indexes)
  for(i in seq_along(folds.indexes))
  {
    train.indexes = unlist(folds.indexes[c(-i)])
    #print(train.indexes)
    test.indexes = unlist(folds.indexes[c(i)])
    #print(stages.levels[test.indexes])
    length(intersect(train.indexes,test.indexes)) == 0
    svm.whole <- svm(x = data[train.indexes,],
                     y = stages.levels[train.indexes])
    pred.svm = predict(svm.whole, data[test.indexes,])
    #print(pred.svm)
    print(table(stages.levels[test.indexes], pred.svm))
  }
}

cv.svm.leave.one.out <- function(data, stages.levels, gamma = 0, kernel = 'linear')
{
  output = rep(NA, nrow(data))
  for(i in seq_along(stages.levels))
  {
    train = data[-i,]
    test = matrix(data[i,], ncol = ncol(train))
    #print(test)
    svm.whole <- svm(x = train, y = stages.levels[-i], kernel = kernel, gamma = gamma, degree = 3)
    #print(ncol(test))
    #print(ncol(train))
    output[i] = predict(svm.whole, test)
    #print(pred.svm)
  }
  levels(output) = levels(stages.levels)
  print(unique(output))
  print(table(stages.levels, output))
}