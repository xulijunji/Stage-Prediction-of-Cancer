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
#  print(num.group)
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

create.ordered.annotation <- function(stages, names)
{
  stage.ind.levels <- list()
  for(i in seq_along(levels(stages)))
  {
    stage = levels(stages)[i]
    stage.ind.levels[[stage]] = which(stages == stage)
  }
  row.indexes <- Reduce(union, stage.ind.levels)
  annotation = NULL
  annotation = data.frame(stage = stages[row.indexes])
  rownames(annotation) = names[row.indexes]
  return(annotation)
}

create.heatmap <- function(data, stages, genes, title, col = NULL, labs = NULL,
                           show_colnames = F, show_rownames = F,
                           cluster_cols = F, cluster_rows = F,
                           legend_breaks = NULL, breaks = NULL)
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
  annotation_col = NULL
  annotation_col = data.frame(names = stages[row.indexes])
  rownames(annotation_col) = rownames(data)[row.indexes]
  if(!is.null(labs))
    labs = rownames(data)[row.indexes]
  pheatmap(t(data[row.indexes, genes]),
           labels_col = labs,
           cluster_cols = cluster_cols, 
           cluster_rows = cluster_rows,
           show_rownames = show_rownames,
           show_colnames = show_colnames,
           annotation_col = annotation_col,
           color = col, main = title, fontsize = 10,
           legend_breaks = legend_breaks,
           breaks = breaks)
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

get.accuracy <- function(res)
{
  acc <- lapply(res, function(x)
    {
    x$eval$other$overall[[1]]
  })
  return(acc)
}

get.sens <- function(res)
{
  sens <- lapply(res, function(x)
  {
    x$eval$other$byClass[1]
  })
  return(sens)
}

get.spec <- function(res)
{
  spec <-  lapply(res, function(x)
    {
    x$eval$other$byClass[2]
  })  
  return(spec)
}

get.f <- function(res)
{
  f.val <-  lapply(res, function(x)
  {
    x$eval$other$byClass[7]
  })  
  return(f.val)
}

get.sep.vals <- function(res.eval)
{
  vals <- sapply(res.eval, function(res)
    {
    res[[1]]
  })
  other.val <- sapply(res.eval, function(res)
  {
    res[[2]]
  })
  return(list(met = vals, dev = other.val))
}

get.dfs.res <- function(res.param.eval, type)
{
  sep.vals <- lapply(res.param.eval, get.sep.vals)
  auc.ind <- which(names(sep.vals) == 'aucs')
  f.ind <- which(names(sep.vals) == 'f.val')
  j <- 0
  sep.vals <- lapply(sep.vals, function(x)
    {
    j <<- j + 1
    
    if(j != auc.ind &  j != f.ind)
    {
      if(type == 'mean')
      {
        x <- lapply(x, function(y) y * 100)
        x <- lapply(x, function(y) formatC(y, digits = 2, format = 'f'))
      }
      else 
      {
        x[[1]] <- x[[1]]*100
        x[[2]] <- x[[2]]
        x[[1]] <- formatC(x[[1]], digits = 2, format = 'f')
        x[[2]] <- as.character(x[[2]])
        x <- list(met = x[[1]], dev = x[[2]])
      }
    }
    else
    {
      if(type == 'mean')
        x <- lapply(x, function(y) formatC(y, digits = 4, format = 'f'))
      else
      {
        x[[1]] <- formatC(x[[1]], digits = 4, format = 'f')
        x[[2]] <- as.character(x[[2]])
        x <- list(met = x[[1]], dev = x[[2]])
      }
    }
    
  })
  sep.text <- '\u00b1'
  if(type != 'mean')
    sep.text <- ','
  #print(sep.vals)
  dfs <- data.frame(Accuracy = paste(sep.vals$accuracy$met, sep.vals$accuracy$dev, sep = sep.text),
                    Aucs = paste(sep.vals$aucs$met, sep.vals$aucs$dev, sep = sep.text),
                    Specificity = paste(sep.vals$spec$met, sep.vals$spec$dev, sep = sep.text),
                    Sensitivity = paste(sep.vals$sens$met, sep.vals$sens$dev, sep = sep.text),
                    FValue = paste(sep.vals$f.val$met, sep.vals$f.val$dev, sep = sep.text)
                    )
  #print(dfs)
  rownames(dfs) <- names(res.param.eval[[1]])
  return(dfs)
}
  

publish.results <- function(res.list, type = 'mean', indexes = seq(length(res.list[[1]])))
{
  aucs <- lapply(res.list, function(class.res)
  {
    atleast.aucs <- get.aucs(class.res[indexes])
  })
  
  acc <- lapply(res.list, function(class.res)
  {
    atleast.acc <- get.accuracy(class.res[indexes])
  })
  spec <- lapply(res.list, function(class.res)
    {
     atleast.spec <- get.spec(class.res[indexes])
  })
  sens <- lapply(res.list, function(class.res)
  {
    atleast.spec <- get.sens(class.res[indexes])
  })
  f.val <- lapply(res.list, function(class.res)
  {
    atleast.spec <- get.f(class.res[indexes])
  })
  #print(acc)
  req.list <- list(aucs, acc, spec, sens, f.val)
  final <- lapply(req.list, function(param)
   {
    lapply(param, function(x)
    {
      x <- unlist(x)
      if(type == 'mean')
        list(mean(x), sd(x))
      else if(type == 'max')
       list(max(x), which.max(x))
      else
        list(min(x), which.min(x))
    })
  })
  names(final) <- c('aucs', 'accuracy', 'spec', 'sens', 'f.val')
  acc.res <- get.dfs.res(final, type)
  return(acc.res)
}

write.dfs.csvs <- function(folder, lists)
{
  for(i in seq_along(lists))
    write.csv(lists[[i]], paste(folder, names(lists)[i], '.csv'))
}

change.microarray <- function(micr.df, genes.list, micr.genes, stage, folds)
{
  temp <- lapply(genes.list, function(genes)
    {
    genes <- intersect(micr.genes, genes)
    indexes <- seq(length(rownames(micr.df)))
    res_micr <- final.res(micr.df, train.ind = indexes, test.ind = indexes, 
                          stage, genes, folds)
    res_micr <- res_micr[[1]]
  })
  names(temp) <- names(genes.list)
  names.genes <- names(temp)
  names.class <- names(temp[[1]])
  res <- list()
  res <- lapply(seq_along(temp[[1]]), function(i)
    {
    res[[names.class[i]]] <<- list()
    res[[names.class[i]]] <<- lapply(seq_along(temp), function(j){
      res[[names.class[i]]][[names.genes[j]]] <<- temp[[names.genes[j]]][[names.class[i]]]
    })
    #names(res[[names.class[i]]])  = names.genes
    #print(typeof(res[[names.class[[i]]]]))
  })
  for(i in seq_along(res))
    names(res[[i]]) <- names.genes
  names(res) <- names.class
  return(res)
}

create.net.df <- function(test.ac)
{
  req.df <- list()
  req.df$auc <- c()
  req.df$accuracy <- c()
  req.df$sens <- c()
  req.df$spec <- c()
  req.df$f_val <- c()
  req.df$classifier <- c()
  req.df$feature_sel <- c()
  for(j in seq_along(test.ac))
  {
    test.ob <- test.ac[[j]]
    for(i in seq_along(test.ob))
    {
      req.df$auc <- c(req.df$auc, unlist(get.aucs(test.ob[[i]])))
      req.df$accuracy <- c(req.df$accuracy, unlist(get.accuracy(test.ob[[i]])))
      req.df$sens <- c(req.df$sens, unlist(get.sens(test.ob[[i]])))
      req.df$spec <- c(req.df$spec, unlist(get.spec(test.ob[[i]])))
      req.df$f_val <- c(req.df$f_val, unlist(get.f(test.ob[[i]])))
      req.df$classifier <- c(req.df$classifier,rep(names(test.ob)[i], 
                                  length(unlist(get.f(test.ob[[i]])))))
      req.df$feature_sel <- c(req.df$feature_sel,rep(names(test.ac)[j],
                          length(unlist(get.f(test.ob[[i]])))))
    }
  }
  df <- data.frame(AUC <- req.df$auc,  Accuracy <- req.df$accuracy, 
                   Sensitivity <- req.df$sens, Specificity <- req.df$spec,
                   F_Value <- req.df$f_val, Classifier <- req.df$classifier,
                   Feature_Selection <- req.df$feature_sel)
  colnames(df) <- c('AUC', 'Accuracy', 'Sensitivity', 'Specificity', 'F_val', 
                    'Classifier', 'Feature_Selection')
  return(df)
}
summarySE <- function(data=NULL, measurevar, groupvars=NULL, na.rm=FALSE,
                      conf.interval=.95, .drop=TRUE) {
  library(plyr)
  
  # New version of length which can handle NA's: if na.rm==T, don't count them
  length2 <- function (x, na.rm=FALSE) {
    if (na.rm) sum(!is.na(x))
    else       length(x)
  }
  
  # This does the summary. For each group's data frame, return a vector with
  # N, mean, and sd
  datac <- ddply(data, groupvars, .drop=.drop,
                 .fun = function(xx, col) {
                   c(N    = length2(xx[[col]], na.rm=na.rm),
                     mean = mean   (xx[[col]], na.rm=na.rm),
                     sd   = sd     (xx[[col]], na.rm=na.rm)
                   )
                 },
                 measurevar
  )
  
  # Rename the "mean" column    
  datac <- rename(datac, c("mean" = measurevar))
  
  datac$se <- datac$sd / sqrt(datac$N)  # Calculate standard error of the mean
  
  # Confidence interval multiplier for standard error
  # Calculate t-statistic for confidence interval: 
  # e.g., if conf.interval is .95, use .975 (above/below), and use df=N-1
  ciMult <- qt(conf.interval/2 + .5, datac$N-1)
  datac$ci <- datac$se * ciMult
  
  return(datac)
}
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)
  
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots==1) {
    print(plots[[1]])
    
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}

get.clust.score <- function(data, genes, stages, classes)
{
  library(mclust)
  clust <- kmeans(data[,genes], classes)
  return(adjustedRandIndex(clust$cluster, stages))
}

run_hclust_on_a_matrix <- function(x, cor.method = "spearman", hclust.method = "ward.D2")
{
  return(stats::hclust(as.dist(1 - cor(x, method = cor.method)), method = hclust.method))
}

create.list.venn <- function(features.list, fold, group)
{
  fold <- as.character(fold)
  group <- as.character(group)
  g1 <- features.list$shrunken[[paste0('atleast_',group)]]
  g2 <- features.list$varSelRF[[paste0('atleast_',group)]]
  g3 <- features.list[[paste0('deseq2',fold, ' fold')]][[paste0('atleast_',group)]]
  g4 <- features.list[[paste0('sam',fold, ' fold')]][[paste0('atleast_',group)]]
  req.list <- list(g1,g2,g3,g4)
  names(req.list) <- c('Shrunken', 'VarSelRF', paste0('Deseq2_log2FC ',fold), paste0('SamSeq_log2FC ',fold))
  return(req.list)
}

create.boxplots <- function(gene, data, stage)
{
  library(ggplot2)
  boxplots <- list()
  for(i in seq_along(gene))
  {
    box.df <- data.frame(gene = data[,gene[i]], stage = stage)
    g <- ggplot(box.df, aes(x = stage, y = gene))+geom_boxplot()+
      ggtitle(gene[i])
    boxplots[[i]] <- g
  }
  return(boxplots)
}
