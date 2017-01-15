source('decomp.R')
source('shrunken/pamr.listgenes.R')
source('after_class.R')
source('shrunken/final.R')

load('environment/req_dfs.RData')
load('environment/stages.level.comb.RData')
load('environment/stage.index.RData')


###Saving the results of the first trial into a single list
first.trial <- list()
first.trial[['gr']] <- gr
first.trial[['train']] <- pamr.train.comb
first.trial[['cv']] <- pamr.cv.comb
first.trial[['genes']] <- pamr.genes.comb
first.trial[['auc']] <- pamr.aucs.comb
first.trial[['conf.mat']] <- confusion.mat
first.trial[['eval.mat']] <- eval.mat
save(first.trial, file = 'environment/first_trial_shrunken_classifier.RData')
sapply(first.trial$genes, length)

gr <- build.groups(length(stages.levels.comb), 5)
g1 <- Reduce(intersect, pamr.genes.comb[c(1,2,4,3)])
g2 <- Reduce(intersect, pamr.genes.comb)
confusion.mat <- list()
eval.mat <- list()