library(varSelRF)
library(randomForest)

load('environment/vs_normal_tumor_repored.RData')
load('environment/sample_info_tumor_rep_normal.RData')
load('environment/req_dfs.RData')
source('../function.R')

###Using gene.list.normal.reported from deseq_class

##Using all genes
stage.ind.san4 <- unlist(stage.ind[-4])
stage.ind.1.2 <- unlist(stage.ind[c(1,2)])

length(unique(unlist(genes.list.normal.reported$N)))
rf.all <- randomForest(req.dfs$vs[,unique(unlist(genes.list.normal.reported$N))],
                       y = stages.levels)
rf.comb <- randomForest(req.dfs$vs[,unique(unlist(genes.list.normal.reported$N))],
                       y = stages.levels.comb)

rf.all.unique <- randomForest(req.dfs$vs[,unlist(genes.list.normal.reported.unique$N)],
                              y = stages.levels)
rf.comb.unique <- randomForest(req.dfs$vs[,unlist(genes.list.normal.reported.unique$N)],
                              y = stages.levels.comb)

varSelRF.rf.all <- varSelRF(req.dfs$vs[,unique(unlist(genes.list.normal.reported$N))],
                            Class  = stages.levels, vars.drop.frac = 0.2)
varSelRF.rf.comb <- varSelRF(req.dfs$vs[,unique(unlist(genes.list.normal.reported$N))],
                            Class  = stages.levels.comb, vars.drop.frac = 0.2)

varSelRF.rf.1.2.3 <- varSelRF(req.dfs$vs[stage.ind.san4,unique(unlist(genes.list.normal.reported$N))],
                              Class  = droplevels(stages.levels[stage.ind.san4]), vars.drop.frac = 0.2)
varSelRF.rf.1.2 <- varSelRF(req.dfs$vs[stage.ind.1.2,unique(unlist(genes.list.normal.reported$N))],
                              Class  = droplevels(stages.levels[stage.ind.1.2]), vars.drop.frac = 0.2)
g <- unlist(strsplit(as.character(varSelRF.rf.1.2$selec.history$Vars.in.Forest[27]), " + ", fixed = T))
rf.1.2.varSelRf.1.2 <- randomForest(req.dfs$vs[stage.ind.1.2, g],
                                    droplevels(stages.levels[stage.ind.1.2]))

###Using unique genes
varSelRF.rf.all.unique <- varSelRF(req.dfs$vs[,unique(unlist(genes.list.normal.reported.unique$N))],
                                   Class  = stages.levels, vars.drop.frac = 0.2)
varSelRF.rf.comb <- varSelRF(req.dfs$vs[,unique(unlist(genes.list.normal.reported.unique$N))],
                             Class  = stages.levels.comb, vars.drop.frac = 0.2)
g <- unlist(strsplit(as.character(varSelRF.rf.comb$selec.history$Vars.in.Forest[16]), " + ", fixed = T))
rf.varSelRf.comb.unique <- randomForest(req.dfs$vs[, g],
                                        stages.levels.comb)
