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