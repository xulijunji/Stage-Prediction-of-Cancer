library(DESeq2)
library(VennDiagram)

library(FactoMineR)


library(factoextra)
load('environment/dds_tumor_reported_normal_stage.RData')
load('environment/req_dfs.RData')
source('across_tumors.R')
results.stage.reported.normal <- list()
results.stage.reported.normal[['normal']] <- comp.res(dds_tumor_reported_normal_stage,'stage.type',
                                                      'N', c('stage i', 'stage ii', 'stage iii', 'stage iv'))
results.stage.reported.normal[['stage i']] <- comp.res(dds_tumor_reported_normal_stage,'stage.type',
                                                      'stage i', c( 'stage ii', 'stage iii', 'stage iv'))
genes.list.normal.reported <-list()
genes.list.normal.reported[['N']] <- list()
genes.list.normal.reported$N[['stage i']] <- get.genes(results.stage.reported.normal$normal$`stage i`,
                                                               2, 0.05, 0.05)
genes.list.normal.reported$N[['stage ii']] <- get.genes(results.stage.reported.normal$normal$`stage ii`,
                                                                2, 0.05, 0.05)
genes.list.normal.reported$N[['stage iii']] <- get.genes(results.stage.reported.normal$normal$`stage iii`,
                                                                 2, 0.05, 0.05)
genes.list.normal.reported$N[['stage iv']] <- get.genes(results.stage.reported.normal$normal$`stage iv`,
                                                                2, 0.05, 0.05)
genes.list.normal.reported[['stagei']] <- list()
genes.list.normal.reported$stagei[['stage ii']] <- get.genes(results.stage.reported.normal$`stage i`$`stage ii`,
                                                       2, 0.05, 0.05)
genes.list.normal.reported$stagei[['stage iii']] <- get.genes(results.stage.reported.normal$`stage i`$`stage iii`,
                                                             2, 0.05, 0.05)
genes.list.normal.reported$stagei[['stage iv']] <- get.genes(results.stage.reported.normal$`stage i`$`stage iv`,
                                                             2, 0.05, 0.05)


genes.list.normal.reported.unique = list()
genes.list.normal.reported.unique$'N'  <- get.unique.genes(genes.list.normal.reported$N)
genes.list.normal.reported.unique$'stagei'  <- get.unique.genes(genes.list.normal.reported$stagei)

load('environment/genes_list_DESeq.RData')
##1 vs 2
rf.1.2.rep.unique <- randomForest(req.dfs$vs[stage.ind.1.2, genes.list.normal.reported.unique$N$`stage ii`],
                                  y = droplevels(stages.levels[stage.ind.1.2]))
rf.1.2.unique.stagei <- randomForest(req.dfs$vs[stage.ind.1.2, genes.list.normal.reported$stagei$`stage ii`],
                                     y = droplevels(stages.levels[stage.ind.1.2]))
###PCA plots
# stage.index.1.2 <- unlist(stage.ind[c(1,2)])
# stage.index.1.3 <- unlist(stage.ind[c(1,3)])
# pca.stage <- PCA(req.dfs$vs[stage.index.1.2,unlist(genes.list.normal.reported.unique$N[c(1,2)])])
# fviz_pca_ind(pca.stage, habillage = droplevels(stages.levels[stage.index.1.2]), geom = c('point'))
# 
# 
# genes.list.unique <- list()
# genes.list.unique[['stage i']]  <- get.unique.genes(genes.list$`stage i`)
# pca.stage <- PCA(req.dfs$vs[stage.index.1.2,unlist(genes.list.unique$`stage i`$`stage iv`)])
# fviz_pca_ind(pca.stage, habillage = droplevels(stages.levels[stage.index.1.2]), geom = c('point'))

####Using feature selection and other algorithm on it




