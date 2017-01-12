library(pamr)
library(pROC)
library(caret)
library(randomForest)
library(e1071)

source('decomp.R')
source('shrunken/pamr.listgenes.R')
source('after_class.R')

load('environment/req_dfs.RData')
load('environment/stages.level.comb.RData')
load('environment/stage.index.RData')

gr <- build.groups(length(stages.levels.comb), 5)
pamr.train.comb <- list()
pamr.cv.comb <- list()
pamr.genes.comb <- list()
pamr.predicted.comb <- list()
pamr.aucs.comb <- list()

for(i in seq_along(gr))
{
  train.ind <- sort(unlist(gr[-i]))
  test.ind <- sort(unlist(gr[i]))
  pamr.train.comb[[i]] <- pamr.train(list(x = as.matrix(t(req.dfs$vs[train.ind,])), 
                                          y = stages.levels.comb[train.ind]))
  pamr.cv.comb[[i]] <- pamr.cv(pamr.train.comb[[i]], 
          data = list(x = t(as.matrix(req.dfs$vs[train.ind,])),
                  y = stages.levels.comb[train.ind]),nfold = 10)
  
 pamr.aucs.comb[[i]] <- sapply(seq_along(pamr.cv.comb[[i]]$threshold), function(x)
 {
   multiclass.roc(stages.levels.comb[train.ind], ordered(pamr.cv.comb[[i]]$yhat[,x]))$auc
 })
 thr.ind = which.max(pamr.aucs.comb[[i]])
 print(thr.ind)
 pamr.genes.comb[[i]] <- pamr.listgene(pamr.train.comb[[i]], 
                            data = list(x=as.matrix(t(req.dfs$vs[train.ind,])),
                                           y=stages.levels.comb[train.ind]), 
                              threshold = pamr.cv.comb[[i]]$threshold[thr.ind],
                            fitcv = pamr.cv.comb[[i]], genenames = T)[,2]
 stages <- pamr.predict(fit = pamr.train.comb[[i]], 
                              newx = as.matrix(t(req.dfs$vs[test.ind,])),
                        threshold = pamr.cv.comb[[i]]$threshold[thr.ind])
 print(table(stages))
 pamr.predicted.comb[[i]] <- table(stages.levels.comb[test.ind], stages)
 remove(stages)
 remove(train.ind)
 remove(test.ind)
 remove(thr.ind)
}
g1 <- Reduce(intersect, pamr.genes.comb[c(1,2,4,3)])
g2 <- Reduce(intersect, pamr.genes.comb)
pamr.plotcv(pamr.cv.comb[[1]])
sapply(pamr.aucs.comb, max)
pamr.predicted.comb[[4]]

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
pamr.confusion(pamr.cv.comb[[2]], threshold = pamr.cv.comb[[2]]$threshold[18])
plot(pamr.aucs.comb[[2]])
##########

####Random Forests
rf.train.conf <- list()
rf.test.conf <- list()
for(i in seq_along(gr))
{
  train.ind <- sort(unlist(gr[-i]))
  test.ind <- sort(unlist(gr[i]))
  rf <- randomForest(req.dfs$vs[train.ind,g2], stages.levels.comb[train.ind])
  rf.train.conf[[i]] <- rf$confusion
  pr <- predict(rf, req.dfs$vs[test.ind,g2])
  rf.test.conf[[i]] <- table(stages.levels.comb[test.ind], pr)
  remove(pr, rf, test.ind, train.ind)
}

####SVMs
svms.test.conf <- list()


####KNN