source('decomp.R')
source('shrunken/pamr.listgenes.R')
source('after_class.R')
source('shrunken/final.R')

load('environment/req_dfs.RData')
load('environment/stages.level.comb.RData')
load('environment/stage.index.RData')

confusion.mat <- list()
eval.mat <- list()
genes.comb <- list()
gr <- build.groups(length(stages.levels.comb), 5)

genes.group <- pamr.listgene(first.trial$train[[5]], 
                             data = list(x=as.matrix(t(req.dfs$vs[sort(unlist(first.trial$gr[-5])),])),
                                                    y=stages.levels.comb[sort(unlist(first.trial$gr[-5]))]), 
                             threshold = first.trial$cv[[5]]$threshold[21],
                             fitcv = first.trial$cv[[5]], genenames = T)
genes.2 <- genes.group[,2][as.numeric(genes.group[,4]) > 0]
genes.1 <- genes.group[,2][as.numeric(genes.group[,3]) > 0]

out.list <- do.shrunken(gr, req.dfs$vs, stages.levels.comb, confusion.mat, eval.mat)
genes <- Reduce(intersect, out.list[[3]])
out.list <- do.knn(gr, req.dfs$vs, genes, stages.levels.comb, out.list[[1]], out.list[[2]])
out.list <- do.naive(gr, req.dfs$vs, genes, stages.levels.comb, out.list[[1]], out.list[[2]])
out.list <- do.svm(gr, req.dfs$vs, genes, stages.levels.comb, out.list[[1]], out.list[[2]])
out.list <- do.rf(gr, req.dfs$vs, genes, stages.levels.comb, out.list[[1]], out.list[[2]])

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
load('environment/first_trial_shrunken_classifier.RData')


library(RColorBrewer)
display.brewer.all()
col <- brewer.pal(9, 'YlOrRd')

library(pheatmap)
names <- sapply(stages.levels.comb, function(x)
  {
  if(x == 'stage i')
    x='E'
  else
    x = 'L'
})
row.indexes <- Reduce(union, stage.ind[c(1,2,3,4)])

pheatmap(t(req.dfs$vs[row.indexes,]), labels_col = names[row.indexes],  
            fontsize_col = 5, cluster_cols = F, color = col,
         annotation_col = annotation_col )
annotation_col = data.frame(names = names[row.indexes])
rownames(annotation_col) = rownames(req.dfs$vs)[row.indexes]


