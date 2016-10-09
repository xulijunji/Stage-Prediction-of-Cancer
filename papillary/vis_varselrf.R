load('environment/tumor_fpqm_comb_boot.RData')
load('environment/tumor_fpqm_comb_varSelRF.RData')
load('environment/tumor_fpqm_log_comb_boot.RData')
load('environment/tumor_fpqm_log_comb_varSelRF.RData')
load('environment/tumor_fpqm_over_boot.RData')
load('environment/tumor_fpqm_over_varSelRF.RData')

library(varSelRF)
library(randomForest)
intersect(tumor.fpqm.comb.varSelRF$selected.vars,tumor.fpqm.log.comb.varSelRF$selected.vars)

write.dfs.csvs <- function(folder, lists)
{
  for(i in seq_along(lists))
    write.csv(lists[[i]], paste(folder, names(lists)[i], '.csv'))
}

tumor.exp.boot$all.data.randomForest$confusion
tumor.fpqm.boot$all.data.randomForest$confusion
tumor.fpqm.log.boot$all.data.randomForest$confusion
tumor.fpqm.log.comb.boot$all.data.randomForest$confusion
tumor.fpqm.comb.boot$all.data.randomForest$confusion
tumor.fpqm.over.boot$all.data.randomForest$confusion

con.list <- lapply(list(tumor.fpqm.boot, tumor.fpqm.log.boot, tumor.fpqm.comb.boot, 
                        tumor.fpqm.log.comb.boot, tumor.fpqm.over.boot), function(x)
                        {
                          x[[5]][[5]]
                        }
                   )
names(con.list) = c('fpqm', 'fpqm_log', 'fpqm_comb', 'fpqm_comb_log', 'fpqm_over')

sel.list <- lapply(list(tumor.fpqm.varSelRF, tumor.fpqm.log.varSelRF, tumor.fpqm.comb.varSelRF, 
                  tumor.fpqm.log.comb.varSelRF, tumor.fpqm.over.varSelRF), function(x)
                  {
                    x[[1]]
                  }
)
names(sel.list) = c('fpqm', 'fpqm_log', 'fpqm_comb', 'fpqm_comb_log', 'fpqm_over')

write.dfs.csvs(folder = '~/Dropbox/honours/sem 7/RNA_Seq/papillary/results/tumor/varSelRF/confusion/', con.list)
write.dfs.csvs(folder = '~/Dropbox/honours/sem 7/RNA_Seq/papillary/results/tumor/varSelRF/selection_history/', sel.list)

get.split <- function(gene.list)
{
  gene.list = unlist(strsplit(as.character(gene.list), split = '+ ', fixed = T))
  gene.list = unlist(strsplit(gene.list, split = ' ', fixed = T))
  return(remove.dots(gene.list))
}

fpqm.genes.selected.his = get.split(sel.list$fpqm$Vars.in.Forest[24])
fpqm.log.genes.selected.his = get.split(sel.list$fpqm_log$Vars.in.Forest[24])
length(intersect(fpqm.log.genes.selected.his, fpqm.genes.selected.his))
remove(fpqm.genes.selected.his, fpqm.log.genes.selected.his)

rf.fpqm.tum.rep.min = randomForest(x = t(exp_fpqm_tumor_reported[intersect(fpqm.log.genes.selected.his,
                                                                           fpqm.genes.selected.his), ]),
                                   y = stages.levels, ntree = 5000, keep.forest = T)

rf.fpqm.log.tum.rep.min = randomForest(x = t(exp_fpqm_tumor_log_reported[intersect(fpqm.log.genes.selected.his,
                                                                           fpqm.genes.selected.his), ]),
                                   y = stages.levels, ntree = 5000, keep.forest = T)

rf.fpqm.tum.rep.min$confusion
rf.fpqm.log.tum.rep.min$confusion
remove(rf.fpqm.log.tum.rep.min, rf.fpqm.tum.rep.min)

length(colnames(exp_fpqm_tumor_reported_over))
rf.fpqm.tum.rep.over <- randomForest(x = t(exp_fpqm_tumor_reported_over[over.sel.genes,]),
                                     y = stages.oversamp, keep.forest = T)
rf.over.predict <- predict(object = rf.fpqm.tum.rep.over, type = 'response',
                                   newdata = t(exp_fpqm_tumor_reported[over.sel.genes,]))

rownames(exp_fpqm_tumor_reported_over) = remove.dots(rownames(exp_fpqm_tumor_reported_over))


fpqm.comb.tum.rep = get.split(sel.list$fpqm_comb$Vars.in.Forest[24])
fpqm.log.comb.tum.rep = get.split(sel.list$fpqm_comb_log$Vars.in.Forest[24])
length(intersect(rf.fpqm.log.comb.tum.rep, rf.fpqm.comb.tum.rep))
remove(fpqm.log.comb.tum.rep, fpqm.comb.tum.rep)

rf.fpqm.comb.tum.rep = randomForest(x = t(exp_fpqm_tumor_reported[intersect(fpqm.log.comb.tum.rep,
                                                                            fpqm.comb.tum.rep), ]),
                                    y = as.factor(stages.levels.comb), ntree = 5000, keep.forest = T)
rf.fpqm.comb.tum.rep$confusion
remove(rf.fpqm.comb.tum.rep,  rf.fpqm.log.comb.tum.rep)

create.confusion.mat <- function(prediction, act.stage)
{
  types <- levels(prediction)
  stages.levels <- as.character(act.stage)
  prediction <- as.character(prediction)
  counts = c()
  conf <- lapply(types, function(x)
    {
      stage.indexes = which(act.stage == x)
      for(i in types)
        counts = c(counts,sum(prediction[stage.indexes] == i))
      print(counts)
  })
  names(conf)= types
  return(conf)
}
####Using genes collected from oversampling

tumor.fpqm.over.varSelRF$selected.vars
over.sel.genes = remove.dots(tumor.fpqm.over.varSelRF$selected.vars)
rf.fpqm.tum.rep.over = randomForest(x = t(exp_fpqm_tumor_reported[over.sel.genes, ]),
                              y = stages.levels, ntree = 5000, keep.forest = T)

rf.fpqm.tum.rep.over = predict(object = tumor.fpqm.over.boot$all.data.randomForest, 
                               newdata = t(exp_fpqm_tumor_log_reported))

rf.fpqm.log.tum.rep.over = randomForest(x = t(exp_fpqm_tumor_log_reported[over.sel.genes, ]),
                              y = stages.levels, ntree = 5000, keep.forest = T)

over.sel.genes.175 = get.split(sel.list$fpqm_over$Vars.in.Forest[23])
rf.fpqm.tum.rep.over.175 = randomForest(x = t(exp_fpqm_tumor_reported[over.sel.genes.175, ]),
             y = stages.levels, ntree = 5000, keep.forest = T)
write.csv(rf.fpqm.tum.rep.over.175$confusion, 'try.csv')

remove(over.sel.genes.175, rf.fpqm.tum.rep.over.175, rf.fpqm.log.tum.rep.over, rf.fpqm.tum.rep.over)

nt_tumor_reported = randomForest(x = t(assay(only.tumor.reported$dfs$nt[diff.genes[[5]], ])),
                                 y = df.stage.tumor.rep$stage, ntree = 5000, keep.forest = T)
nt_tumor_reported$confusion
                                  
rf.fpqm.tum.rep$confusion
rf.fpqm.tum.rep.min$confusion

###Using genes from SMOTE
load('environment/tumor_fpqm_smote_comb_var.RData')
rf.fpqm.smt.comb.smotedata <- randomForest(x = t.exp.fpqm.reported[,tumor.fpqm.smt.comb.varSelRF$selected.vars], y = t.exp.fpqm.reported$comb,
                                           ntree = 5000, keep.forest = T) 
rf.fpqm.smt.comb.smotedata$confusion


##NT
rf.nt.smt.comb.smotedata. <- randomForest(x = smt.comb.nt[,tumor.nt.smt.comb.varSelRF$selected.vars], y = smt.comb.nt[['comb']],
                                          ntree = 5000, keep.forest = T)
rf.nt.predict <- predict(object = rf.nt.smt.comb.smotedata., 
                         newdata = t(assay(only.tumor.reported$dfs$nt[tumor.nt.smt.comb.varSelRF$selected.vars,])))

rf.nt.smote <- randomForest(x= t(assay(only.tumor.reported$dfs$nt[tumor.nt.smt.comb.varSelRF$selected.vars,])),
                            y = df.stage.tumor.rep$stage, ntree = 5000, keep.forest = T)
rf.nt.smote.comb <- randomForest(x= t(assay(only.tumor.reported$dfs$nt[tumor.nt.smt.comb.varSelRF$selected.vars,])),
                            y = df.stage.tumor.rep$stage, ntree = 5000, keep.forest = T)

con.list <- lapply(list(rf.nt.smt.comb.smotedata., rf.nt.smote, rf.nt.smote.comb) ,
                         function(x)
                        {
                          x$confusion
                        }
)
names(con.list) = c('nt_smt_smotedata', 'nt_smt', 'nt_smt_comb')
write.dfs.csvs(folder = '~/Dropbox/honours/sem 7/RNA_Seq/papillary/results/tumor/varSelRF/confusion/', con.list)
conf.pred <- create.confusion.mat(rf.nt.predict, stages.levels.comb)

rf.nt.across.stage <- randomForest(x= t(assay(only.tumor.reported$dfs$nt[Reduce(union,Reduce(union,genes.list)),])),
                            y = df.stage.tumor.rep$stage, ntree = 5000, keep.forest = T)

rf.nt.under.comb <- randomForest(x= t(assay(only.tumor.reported$dfs$nt[tumor.nt.under.comb.varSelRF$selected.vars,])),
                            y = as.factor(stages.levels.comb), ntree = 5000, keep.forest = T)
rf.nt.under.comb.unddata <- randomForest(x= und.comb.nt$X[,tumor.nt.under.comb.varSelRF$selected.vars],
                                 y = und.comb.nt$Y, ntree = 5000, keep.forest = T)
rf.nt.under.predict <- predict(object = rf.nt.under.comb.unddata, 
                               newdata = t(assay(only.tumor.reported$dfs$nt[tumor.nt.under.comb.varSelRF$selected.vars,])))
conf.pred.under.nt <- create.confusion.mat(rf.nt.under.predict, stages.levels.comb)
con.list <- lapply(list(rf.nt.across.stage, rf.nt.under.comb, rf.nt.under.comb.unddata) ,
                   function(x)
                   {
                     x$confusion
                   }
)
names(con.list) = c('across_stage', 'nt_under_comb', 'nt_under_comb_unddata')
write.dfs.csvs(folder = '~/Dropbox/honours/sem 7/RNA_Seq/papillary/results/tumor/varSelRF/confusion/', con.list)

rf.nt.across.stage.comb <- randomForest(x= t(assay(only.tumor.reported$dfs$nt[Reduce(union,Reduce(union,genes.list)),])),
                                   y = as.factor(stages.levels.comb), ntree = 5000, keep.forest = T)

sum(which(rf.nt.predict == 'stage i') == which(stages.levels.comb == 'stage i'))

rf.fpqm.strat.rep <- randomForest(x = t(exp_fpqm_tumor_reported[tumor.fpqm.varSelRF.samp.rep$selected.vars,]),
                                  y= stages.levels, ntree = 5000, keep.forest = T, nodesize = 5, 
                                  strata = stages.levels,  sampsize = c(15,15,10,15))
rf.fpqm.strat.rep$confusion

rf.fpqm.strat.not.rep <- randomForest(x = t(exp_fpqm_tumor_reported[tumor.fpqm.varSelRF.samp.not.rep$selected.vars,]),
                                  y= stages.levels, ntree = 5000, keep.forest = T, nodesize = 1,
                                  strata = stages.levels,  sampsize = c(50,15,30,10), replace = F)
rf.fpqm.strat.not.rep$confusion

rf.fpqm.wt <- randomForest(x = t(exp_fpqm_tumor_reported[tumor.fpqm.varSelRF.wt$selected.vars,]),
                                  y= stages.levels, ntree = 5000, keep.forest = T, nodesize = 50, 
                           wt = c(1,7.8,3,11))
rf.fpqm.wt$confusion

rf.fpqm.comb.strat <- randomForest(x = t(exp_fpqm_tumor_reported[tumor.fpqm.varSelRF.comb.strat$selected.vars,]),
                                   y = stages.levels.comb, strata = stages.levels.comb, sampsize = c(100,66))
rf.fpqm.comb.strat$confusion

rf.fpqm.comb.wt <- randomForest(x = t(exp_fpqm_tumor_reported[tumor.fpqm.varSelRF.comb.wt$selected.vars,]),
                                   y = stages.levels.comb, classwt = c(1,6))
rf.fpqm.comb.wt$confusion

rf.fpqm.1.4 <- randomForest(x = t(exp_fpqm_tumor_reported[tumor.fpqm.varSelRF.1.4$selected.vars,stages.levels[union(stage.ind$`stage i`, stage.ind$`stage iv`)]]),
                            y = droplevels(stages.levels[union(stage.ind$`stage i`, stage.ind$`stage iv`)]))
rf.fpqm.1.4$confusion

rf.fpqm.2.3.4 <- randomForest(x = t(exp_fpqm_tumor_reported[tumor.fpqm.varSelRF.2.3.4$selected.vars, 
                                                            Reduce(union, stage.ind[-1])]),
                              y = droplevels(stages.levels[Reduce(union, stage.ind[-1])]))
rf.fpqm.2.3.4$confusion 

rf.fpqm.2.3.4 <- randomForest(x = t(exp_fpqm_tumor_reported[tumor.fpqm.varSelRF.2.3.4.strat$selected.vars, 
                                                            Reduce(union, stage.ind[-1])]),
                              y = droplevels(stages.levels[Reduce(union, stage.ind[-1])]),
                              strata = droplevels(stages.levels[Reduce(union, stage.ind[-1])]), sampsize = c(22,25,15))
rf.fpqm.2.3.4$confusion

rf.fpqm.1.2 <- randomForest(x = t(exp_fpqm_tumor_reported[tumor.fpqm.varSelRF.1.2$selected.vars,Reduce(union,stage.ind[c(1,2)])]),
                            y = droplevels(stages.levels[Reduce(union,stage.ind[c(1,2)])]),
                            strata = droplevels(stages.levels[Reduce(union,stage.ind[c(1,2)])]), sampsize = c(34,22))
rf.fpqm.1.2$confusion

rf.fpqm.1.2.strata.1 <- randomForest(x = t(exp_fpqm_tumor_reported[tumor.fpqm.varSelRF.1.2.strat.1$selected.vars,Reduce(union,stage.ind[c(1,2)])]),
                            y = droplevels(stages.levels[Reduce(union,stage.ind[c(1,2)])]),
                            strata = droplevels(stages.levels[Reduce(union,stage.ind[c(1,2)])]), sampsize = c(22,22))
rf.fpqm.1.2.strata.1$confusion

rf.fpqm.1.2.strata.2 <- randomForest(x = t(exp_fpqm_tumor_reported[tumor.fpqm.varSelRF.1.2.strat.2$selected.vars,Reduce(union,stage.ind[c(1,2)])]),
                                     y = droplevels(stages.levels[Reduce(union,stage.ind[c(1,2)])]),
                                     strata = droplevels(stages.levels[Reduce(union,stage.ind[c(1,2)])]), sampsize = c(22,22))
rf.fpqm.1.2.strata.2$confusion

rf.fpqm.1.2.strata.3 <- randomForest(x = t(exp_fpqm_tumor_reported[tumor.fpqm.varSelRF.1.2.strat.3$selected.vars,Reduce(union,stage.ind[c(1,2)])]),
                                     y = droplevels(stages.levels[Reduce(union,stage.ind[c(1,2)])]),
                                     strata = droplevels(stages.levels[Reduce(union,stage.ind[c(1,2)])]), sampsize = c(22,22))
rf.fpqm.1.2.strata.3$confusion

rf.fpqm.3.4 <- randomForest(x = t(exp_fpqm_tumor_reported[tumor.fpqm.varSelRF.3.4$selected.vars,Reduce(union,stage.ind[c(3,4)])]),
                            y = droplevels(stages.levels[Reduce(union,stage.ind[c(3,4)])]),
                            strata = droplevels(stages.levels[Reduce(union,stage.ind[c(3,4)])]), sampsize = c(20,15))
rf.fpqm.3.4$confusion

rf.fpqm.3.4.1 <- randomForest(x = t(exp_fpqm_tumor_reported[tumor.fpqm.varSelRF.3.4.strat.1$selected.vars,Reduce(union,stage.ind[c(3,4)])]),
                            y = droplevels(stages.levels[Reduce(union,stage.ind[c(3,4)])]),
                            strata = droplevels(stages.levels[Reduce(union,stage.ind[c(3,4)])]), sampsize = c(20,15))
rf.fpqm.3.4.1$confusion

rf.fpqm.3.4.2 <- randomForest(x = t(exp_fpqm_tumor_reported[tumor.fpqm.varSelRF.3.4.strat.2$selected.vars,Reduce(union,stage.ind[c(3,4)])]),
                              y = droplevels(stages.levels[Reduce(union,stage.ind[c(3,4)])]),
                              strata = droplevels(stages.levels[Reduce(union,stage.ind[c(3,4)])]), sampsize = c(30,15))#remb all
rf.fpqm.3.4.2$confusion

rf.fpqm.2.3 <-  randomForest(x = t(exp_fpqm_tumor_reported[tumor.fpqm.varSelRF.2.3$selected.vars,Reduce(union,stage.ind[c(2,3)])]),
                             y = droplevels(stages.levels[Reduce(union,stage.ind[c(2,3)])]),
                             strata = droplevels(stages.levels[Reduce(union,stage.ind[c(2,3)])]), sampsize = c(20,20))#remb all
rf.fpqm.2.3$confusion

rf.fpqm.2.3.1 <-  randomForest(x = t(exp_fpqm_tumor_reported[tumor.fpqm.varSelRF.2.3.strat.1$selected.vars,Reduce(union,stage.ind[c(2,3)])]),
                             y = droplevels(stages.levels[Reduce(union,stage.ind[c(2,3)])]),
                             strata = droplevels(stages.levels[Reduce(union,stage.ind[c(2,3)])]))
rf.fpqm.2.3.1$confusion

rf.fpqm.1.3 <- randomForest(x = t(exp_fpqm_tumor_reported[tumor.fpqm.varSelRF.1.3$selected.vars,Reduce(union,stage.ind[c(1,3)])]),
                            y = droplevels(stages.levels[Reduce(union,stage.ind[c(1,3)])]),
                            strata = droplevels(stages.levels[Reduce(union,stage.ind[c(1,3)])]),
                            sampsize = c(102,51))
rf.fpqm.1.3$confusion

rf.fpqm.1.3.strat.1 <- randomForest(x = t(exp_fpqm_tumor_reported[tumor.fpqm.varSelRF.1.3.strat.1$selected.vars,Reduce(union,stage.ind[c(1,3)])]),
                            y = droplevels(stages.levels[Reduce(union,stage.ind[c(1,3)])]),
                            strata = droplevels(stages.levels[Reduce(union,stage.ind[c(1,3)])]),
                            sampsize = c(41,51))
rf.fpqm.1.3.strat.1$confusion

rf.fpqm.1.3.strat.2 <- randomForest(x = t(exp_fpqm_tumor_reported[tumor.fpqm.varSelRF.1.3.strat.2$selected.vars,Reduce(union,stage.ind[c(1,3)])]),
                                    y = droplevels(stages.levels[Reduce(union,stage.ind[c(1,3)])]),
                                    strata = droplevels(stages.levels[Reduce(union,stage.ind[c(1,3)])]),
                                    sampsize = c(60,51))
rf.fpqm.1.3.strat.2$confusion

rf.fpqm.2.4 <- randomForest(x = t(exp_fpqm_tumor_reported[tumor.fpqm.varSelRF.2.4$selected.vars, Reduce(union,stage.ind[c(2,4)])]),
                            y = droplevels(stages.levels[Reduce(union,stage.ind[c(2,4)])]))
rf.fpqm.2.4$confusion

all.genes.2.3.4 <-  unique(c(tumor.fpqm.varSelRF.2.4$selected.vars, tumor.fpqm.varSelRF.2.3$selected.vars,
                          tumor.fpqm.varSelRF.3.4.strat.1$selected.vars))
rf.fpqm.2.3.4.comb <- randomForest(x = t(exp_fpqm_tumor_reported[all.genes.2.3.4, 
                                                                 Reduce(union, stage.ind[-1])]),
                                   y = droplevels(stages.levels[Reduce(union, stage.ind[-1])]),
                                   strata = droplevels(stages.levels[Reduce(union, stage.ind[-1])]),
                                   sampsize = c(22,25,15))
rf.fpqm.2.3.4.comb$confusion


all.genes <- unique(c(tumor.fpqm.varSelRF.1.2.strat.3$selected.vars, tumor.fpqm.varSelRF.2.3$selected.vars, 
               tumor.fpqm.varSelRF.3.4.strat.1$selected.vars, tumor.fpqm.varSelRF.1.4$selected.vars,
               tumor.fpqm.varSelRF.1.3.strat.1$selected.vars, tumor.fpqm.varSelRF.2.4$selected.vars))


replicate.conf <- function(N, data, stages.levels, sampsize = c(172,22,51,15))
{
 a <- function()
 {
   rf.fpqm.all.genes <- randomForest(x = data,
                                     y = stages.levels, strata = stages.levels, sampsize = sampsize)
   return(rf.fpqm.all.genes$confusion[,5])
 }
  means <- replicate(N, a())
  print(apply(means,1,mean))
}
replicate.conf(100, t(exp_fpqm_tumor_reported[all.genes, ]), stages.levels, sampsize = c(15,13,14,15))
replicate.conf(100, t(exp_fpqm_tumor_reported[remove.dots(tumor.fpqm.varSelRF$selected.vars), ]), stages.levels, sampsize = c(15,13,14,15))

pca_over_sel_genes = list()
for(i in seq_along(only.tumor.reported$dfs))
{
  pca_over_sel_genes[[names(only.tumor.reported$dfs)[i]]] = 
    plotPCA(only.tumor.reported$dfs[[i]][over.sel.genes, ], intgroup = 'stage',
            title = names(only.tumor.reported$dfs),
            colData = df.stage.tumor.rep)
}
pca_over_sel_genes$fpqm
lda_over_sel_genes = list()
for(i in seq_along(only.tumor.reported$dfs))
{
  mat = only.tumor.reported$dfs[[i]]
  if(typeof(mat) == 'S4')
    mat = assay(mat)
  
  lda1 = lda(t(mat[over.sel.genes, ]), grouping = stages.levels)
  lda1.p = predict(lda1)
  lda1.p$x = data.frame(lda1.p$x)
  lda1.p$x$group = stages.levels
  lda_over_sel_genes[[names(only.tumor.reported$dfs)[i]]] = 
    create.lda.plots(list(lda1.p, lda1), 
                     title = names(only.tumor.reported$dfs)[i])
  remove(mat)
}
lda_over_sel_genes$nt
lda_over_sel_genes$fpqm_log


intersect(diff.genes[[5]], over.sel.genes)
write(over.sel.genes, '')