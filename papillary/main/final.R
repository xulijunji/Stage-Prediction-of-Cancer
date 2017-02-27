library(pamr)
library(pROC)
library(caret)
library(randomForest)
library(e1071)
#source('across_tumors.R')


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
                    y = stages.levels[train.ind]),nfold = 10)
    
   pamr.aucs.comb[[i]] <- sapply(seq_along(pamr.cv.comb[[i]]$threshold), function(x)
   {
     multiclass.roc(stages.levels[train.ind], ordered(pamr.cv.comb[[i]]$yhat[,x]))$auc
   })
   thr.ind = which.max(pamr.aucs.comb[[i]])
   
   pamr.genes.comb[[i]] <- pamr.listgene(pamr.train.comb[[i]], 
                              data = list(x=as.matrix(t(data[train.ind,])),
                                             y=stages.levels[train.ind]), 
                                threshold = pamr.cv.comb[[i]]$threshold[thr.ind],
                              fitcv = pamr.cv.comb[[i]], genenames = T)
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

do.varselRF <- function(gr, data, stages)
{
  model.vars <- list()
  model.selec_his <- list()
  for(i in seq_along(gr))
  {
    train.indexes = sort(unlist(gr[-i]))
    var.ob <- varSelRF(xdata = data[train.indexes,], 
                       Class = stages[train.indexes])
    model.vars[[i]] = var.ob
    #model.selec_his[[i]] = var.ob$selec.history
  }
  return(model.vars)
}

create.Deseq2 <- function(gr, counts_data, colData)
{
  deseq_list <- lapply(seq(length(gr)), function(x)
    {
    train.ind  <- sort(unlist(gr[-x]))
    dds_obj <- DESeqDataSetFromMatrix(counts_data[,train.ind],
                            colData = colData[train.ind,], design = ~stage)
    dds_obj <- dds_obj[rowSums(assay(dds_obj)) > 2]
    dds_obj <- DESeq(dds_obj)
  })
  return(deseq_list)
}
do.Deseq2 <- function(dds_list)
{
  deseq.res <- lapply(dds_list, function(dds)
    {
      comp.res(dds, 'stage', 'stage i', 'stage iv')[['stage iv']]
  })
  return(deseq.res)
}
#######Trial########
# g <- pamr.listgene(pamr.train.comb[[1]], 
#                    data = list(x=as.matrix(t(req.dfs$vs[unlist(gr[-1]),])),
#                                y=stages.levels.comb[unlist(gr[-1])]), 
#                    threshold = pamr.cv.comb[[1]]$threshold[15],
#                    fitcv = pamr.cv.comb[[i]], genenames = T)[,2]
# stages <- pamr.predict(pamr.train.comb[[2]], 
#                        as.matrix(t(req.dfs$vs[gr[[2]],])),
#                        pamr.cv.comb[[2]]$threshold[17])
# table(stages)
# table(stages.levels.comb[gr[[2]]], stages)
# pamr.confusion(pamr.cv.comb[[1]], threshold = pamr.cv.comb[[1]]$threshold[17])
# plot(pamr.aucs.comb[[1]])
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

###Fast Filter
do.fast.filter <- function(gr, data, stages)
{
  features.ff <- lapply(gr, function(x)
    {
    train.indexes <- sort(setdiff(unlist(gr), x))
    req.data <- cbind(data[train.indexes,], stages[train.indexes])
    select.fast.filter(req.data,disc.method = 'MDL', threshold = 0.05)
  })
  return(features.ff)
}

###
do.forward.corr <- function(gr, data, stages)
{
  features.for.corr <- lapply(gr, function(x)
    {
    train.indexes <- sort(setdiff(unlist(gr), x))
    req.data <- cbind(data[train.indexes,], stages[train.indexes])
    select.forward.Corr(req.data, disc.method = 'MDL')
  })
  return(features.for.corr)
}

##Get the genes which crosses a certain threshold for a shrunken centroid
get.shrunken.stage.wise.genes <- function(shrunken.gene.df, stage.ind, threshold.ind, threshold)
{
  genes <- shrunken.gene.df[,2][which(as.numeric(shrunken.gene.df[,stage.ind]) > 0 &
                                        as.numeric(shrunken.gene.df[, threshold.ind]) > threshold)]
  return(genes)
}

##Given a list of pamr.list.gene df extract from the genes
get.genes.shrunken <- function(shrunken.genes.df.list)
{
  genes = list()
  for(i in seq_along(shrunken.genes.df.list))
  {
    genes[[i]] = shrunken.genes.df.list[[i]][,2]
  }
  return(genes)
}

##Gives the genes from list of varSelRf object
get.varSelRf.genes <- function(var.ob.list, type = 1, indexes = NULL)
{
  ##type = 1 for minimum OOB error
  ##type = 2 for indexes of OOB in selec history
  
  if(type == 1)
  {
    genes.list <- sapply(var.ob.list, function(x)
    {
      x[[3]]
    })
    return(genes.list)
  }
  else
  {
    genes.list <- mapply(function(x,y)
    {
      x[[y]]
    }, var.ob.list, indexes)
    return(genes.list)
  }
}

###Get the minimum error genes varSelRf
get.min.oob.varselRf <- function(varselRf.ob.list)
{
  sel.genes.list <- sapply(varselRf.ob.list, function(obj)
    {
    obj = obj$selec.history
    min.oob.error <- min(obj[,3])
    best.pos <-
      which(obj[,3] == min.oob.error)[which.min(obj[,1][which(obj[,3] == min.oob.error)])]
    genes <- obj[best.pos, 2]
    genes <- as.character(genes)
    genes <- strsplit(genes, split = ' + ', fixed = T)
  })
  return(sel.genes.list)
}

###Shrunken has 2 stages and for each group finds the gene corresponding to each stage
get.shrunken.group.stage <- function(shrunken.gene.object.list)
{
  genes.1.list <- list()
  genes.2.list <- list()
  for(i in seq_along(shrunken.gene.object.list))
  {
    a <- sapply(c(3,4), function(x)
    {
      get.shrunken.stage.wise.genes(shrunken.gene.object.list[[i]], x, 6, 0)
    }
    )
    genes.1.list[[i]] <- a[[1]]
    genes.2.list[[i]] <- a[[2]]
  }
  for(i in seq_along(genes.1.list))
  {
    for(j in seq_along(genes.2.list))
    {
      if(length(intersect(genes.1.list[[i]], genes.2.list[[j]])) != 0)
        return(-1)
    }
  }
  return(list(genes.1.list, genes.2.list))
}

###Returns a df of the genes common in certain number of groups along with their stages
get.shrunken.common.stage <- function(shrunken.stage.object.list, no.groups.list)
{
  genes.list.1 <- sapply(no.groups.list, function(x)
    {
    get.genes.common(shrunken.stage.object.list[[1]],x)}
    )
  genes.list.2 <- sapply(no.groups.list, function(x)
  {
    get.genes.common(shrunken.stage.object.list[[2]],x)}
  )
  dfs.list <- list()
  for(i in seq_along(genes.list.1))
  {
    d1 <- data.frame(genes = c(genes.list.1[[i]]), stage = rep('1', length(genes.list.1[[i]])))
    d2 <- data.frame(genes = c(genes.list.2[[i]]), stage = rep('2', length(genes.list.2[[i]])))
    dfs.list[[i]] <- rbind(d1,d2)
  }
  names(dfs.list) = paste(rep(c('atleast_'), length(no.groups.list)), no.groups.list)
  return(dfs.list)
}