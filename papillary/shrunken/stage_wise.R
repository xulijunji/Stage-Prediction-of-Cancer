###Contains stage wise separation
library(pamr)
source('shrunken/pamr.listgenes.R')
source('after_class.R')
load('environment/stages_levels.RData')

stage.index.san4 <- sort(Reduce(union, stage.ind[-4]))
mapply(function(x,y)
  {
  sum(stage.ind[[y]] %in% x) == length(x)
}, stage.ind, c(1:3))

vs.train.san4 <- pamr.train(list(x = as.matrix(t(req.dfs$vs[stage.index.san4,])), 
                                 y = droplevels(stages.levels[stage.index.san4])))
vs.cv.san4 <- pamr.cv(vs.train.san4, list(x = as.matrix(t(req.dfs$vs[stage.index.san4,])), 
                                          y = droplevels(stages.levels[stage.index.san4])))
pamr.plotcv(vs.cv.san4)
plot.auc(vs.cv.san4, droplevels(stages.levels[stage.index.san4]))
pamr.confusion(vs.cv.san4, threshold = vs.cv.san4$threshold[21])

stage.index.1.2 <- sort(Reduce(union, stage.ind[c(1,2)]))
vs.train.1.2 <- pamr.train(list(x = as.matrix(t(req.dfs$vs[stage.index.1.2,])), 
                                y = droplevels(stages.levels[stage.index.1.2])))
vs.cv.1.2 <- pamr.cv(vs.train.1.2, list(x = as.matrix(t(req.dfs$vs[stage.index.1.2,])), 
                          y = droplevels(stages.levels[stage.index.1.2])))
pamr.plotcv(vs.cv.1.2)
plot.auc(vs.cv.1.2, droplevels(stages.levels[stage.index.1.2]))
pamr.confusion(vs.cv.1.2, threshold = vs.cv.1.2$threshold[10])
multiclass.roc(droplevels(stages.levels[stage.index.1.2]), ordered(vs.cv.1.2$yhat[,10]))$auc
shrunken.genes.vs.1.2 <- pamr.listgene(vs.train.1.2, data = list(x=as.matrix(t(req.dfs$vs)),
                                                y=droplevels(stages.levels[stage.index.1.2])), 
                                       threshold = vs.cv.1.2$threshold[10],
                                       fitcv = vs.cv.1.2, genenames = T)
View(shrunken.genes.vs.1.2)
shrunken.genes.vs.1.2.st1 <- get.shrunken.stage.wise.genes(shrunken.genes.vs.1.2, 3, 6, 0.8)
shrunken.genes.vs.1.2.st2 <- get.shrunken.stage.wise.genes(shrunken.genes.vs.1.2, 4, 6, 0.8)
