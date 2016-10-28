library(randomForest)

find.indexes.samples <- function(df, genes, samples, stages.levels, n, stage, 
                                 sampsize = table(droplevels(stages.levels[samples])))
{
  set.seed(1)
  stage.indexes <- lapply(seq(1:n), function(x)
    {
    z = randomForest(x = t(df[genes,samples]), y = droplevels(stages.levels[samples]),
                     strata = droplevels(stages.levels[samples]), sampsize = sampsize)
    new.indexes = rep('none', 260)
    new.indexes[samples] = as.character(z$predicted)
    #print(new.indexes[1:20])
    which(new.indexes == stage & stages.levels == stage)
  })
  return(stage.indexes)
}

rf.fpqm.1.3.samples <- find.indexes.samples(exp_fpqm_tumor_reported, tumor.fpqm.varSelRF.1.3$selected.vars,
                                            Reduce(union,stage.ind[c(1,3)]), stages.levels, 10, 'stage iii')

rf.fpqm.1.3.samples.strat.1 <- find.indexes.samples(exp_fpqm_tumor_reported,
                                                    tumor.fpqm.varSelRF.1.3$selected.vars,
                                                    Reduce(union,stage.ind[c(1,3)]), sampsize = c(41,51),
                                                    stages.levels, 10, 'stage iii')

rf.fpqm.all.genes.samples <- find.indexes.samples(exp_fpqm_tumor_reported,all.genes,
                                                  c(1:260), stages.levels, 10, 'stage iii')

  
rf.fpqm.all.genes.samples.1.3 <- find.indexes.samples(exp_fpqm_tumor_reported, all.genes,
                                                      Reduce(union,stage.ind[c(1,3)]), stages.levels,
                                                      10, 'stage iii')


rf.fpqm <- find.indexes.samples(exp_fpqm_tumor_reported, remove.dots(tumor.fpqm.varSelRF$selected.vars),
                                c(1:260), stages.levels, 10, 'stage iii')


length(intersect(rf.fpqm.1.3.samples[[1]], rf.fpqm.1.3.samples.strat.1[[1]]))
length(intersect(rf.fpqm.1.3.samples[[1]], rf.fpqm[[1]]))
length(intersect(rf.fpqm.all.genes.samples[[1]], rf.fpqm.1.3.samples[[1]]))
length(intersect(rf.fpqm.all.genes.samples[[1]], rf.fpqm[[1]]))
length(intersect(rf.fpqm.all.genes.samples[[1]], rf.fpqm.all.genes.samples.1.3[[1]]))
length(intersect(rf.fpqm.1.3.samples.strat.1[[1]], rf.fpqm.all.genes.samples.1.3[[1]]))

##Conclusion from above so far all same types getting separated

###PCA analysis
library(FactoMineR)
library(factoextra)

df.pca.fpqm <- req.dfs$fpqm
rownames(df.pca.fpqm) = df.stage.tumor.rep$short.id
pca.3stage <- PCA(df.pca.fpqm[stage.ind[[3]], all.genes], graph = F)
clus <- HCPC(pca.3stage, nb.clust = 2)
plot(clus, choice = 'tree', axes = c(1,2))

fviz_pca_ind(pca.3stage, label = 'ind', axes = c(1,2))

pca.1stage <- PCA(df.pca.fpqm[stage.ind[[1]], all.genes], graph = F)
fviz_pca_ind(pca.1stage, label = 'ind', axes = c(1,2))
clus <- HCPC(pca.1stage, nb.clust = 2)
plot(clus,axes = c(1,2), choice = 'tree')


