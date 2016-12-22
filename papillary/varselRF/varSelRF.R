install.packages('~/varSelRF_0.7-5.tar.gz', repos = NULL, type = 'source')
library(varSelRF)
library(snow)


df.stages$sample.id == colnames(exp_fpqm_tumor)

stages.levels <- df.stages$stage[indexes.stages.reported]
stages.levels <- droplevels(stages.levels)


cl <- makeMPIcluster(8)
tumor.fpqm.varSelRF <- varSelRF(xdata = t(exp_fpqm_tumor_reported), Class = stages.levels, 
                                 keep.forest = T)
tumor.fpqm.boot <- varSelRFBoot(t(exp_fpqm_tumor_reported), stages.levels, bootnumber = 10, usingCluster = T,
                                TheCluster = cl, srf= tumor.fpqm.varSelRF)
tumor.fpqm.log.varSelRF <- varSelRF(xdata = t(exp_fpqm_tumor_log_reported), Class = stages.levels, 
                                keep.forest = T)
tumor.fpqm.log.boot <- varSelRFBoot(t(exp_fpqm_tumor_log_reported), stages.levels, bootnumber = 10, usingCluster = T,
                                TheCluster = cl, srf= tumor.fpqm.log.varSelRF)
tumor.exp.varSelRF <- varSelRF(xdata = t(exp_prof_tumor_reported), Class = stages.levels, 
                                keep.forest = T)
tumor.exp.boot <- varSelRFBoot(t(exp_prof_tumor_reported), stages.levels, bootnumber = 10, usingCluster = T,
                                TheCluster = cl, srf= tumor.exp.varSelRF)

tumor.fpqm.comb.varSelRF <- varSelRF(xdata = t(exp_fpqm_tumor_reported), Class = stages.levels.comb, 
                                     keep.forest = T)
tumor.fpqm.comb.boot <- varSelRFBoot(t(exp_fpqm_tumor_reported), stages.levels.comb, bootnumber = 10, usingCluster = T,
                                     TheCluster = cl, srf= tumor.fpqm.comb.varSelRF)
tumor.fpqm.log.comb.varSelRF <- varSelRF(xdata = t(exp_fpqm_log_tumor_reported), Class = stages.levels.comb, 
                                         keep.forest = T)
tumor.fpqm.log.comb.boot <- varSelRFBoot(t(exp_fpqm_log_tumor_reported), stages.levels.comb, bootnumber = 10, usingCluster = T,
                                         TheCluster = cl, srf= tumor.fpqm.log.comb.varSelRF)

tumor.fpqm.over.varSelRF <- varSelRF(xdata = t(exp_fpqm_tumor_reported_over), Class = stages.oversamp,
                                     keep.forest = T)
tumor.fpqm.over.boot <- varSelRFBoot(t(exp_fpqm_tumor_reported_over), Class = stages.oversamp, bootnumber = 10, usingCluster = T,
                                     TheCluster = cl, srf= tumor.fpqm.over.varSelRF)

tumor.fpqm.smt.comb.varSelRF <-  varSelRF(xdata = data.frame(smt.comb[, -which(c('comb') == colnames(smt.comb))]),
                                          Class = smt.comb$comb, keep.forest = T)
tumor.nt.smt.comb.varSelRF <- varSelRF(xdata = data.frame(smt.comb.nt[, -which(c('comb') == colnames(smt.comb.nt))]),
                                            Class = smt.comb.nt$comb, keep.forest = T)
tumor.fpqm.under.comb.varSelRF <- varSelRF(xdata = und.comb.fpqm$X, Class = und.comb.fpqm$Y, keep.forest = T)
tumor.nt.under.comb.varSelRF <- varSelRF(xdata = und.comb.nt$X, Class = und.comb.nt$Y, keep.forest = T)

tumor.fpqm.varSelRF.samp.rep <- varSelRF(xdata = t(exp_fpqm_tumor_reported), Class = stages.levels, 
                                sampsize = c(30,22,25,15), strata = stages.levels, replace = T)
tumor.fpqm.varSelRF.samp.not.rep <- varSelRF(xdata = t(exp_fpqm_tumor_reported), Class = stages.levels, 
                                             strata = stages.levels,
                                sampsize = c(25,15,20,10), replace = F)
tumor.fpqm.varSelRF.wt <- varSelRF(xdata = t(exp_fpqm_tumor_reported), Class = stages.levels, 
                                            classwt = c(1,7.8,3.37,11.47))

tumor.fpqm.varSelRF.comb.strat <-  varSelRF(xdata = t(exp_fpqm_tumor_reported), Class = stages.levels.comb,
                                      strata = stages.levels.comb, sampsize = c(100,66))
tumor.fpqm.varSelRF.comb.wt <-  varSelRF(xdata = t(exp_fpqm_tumor_reported), Class = stages.levels.comb,
                                      classwt = c(1,2.5))

tumor.fpqm.varSelRF.1.4 <-  varSelRF(xdata = t(exp_fpqm_tumor_reported[,union(stage.ind$`stage i`, stage.ind$`stage iv`)]),
                                      Class = droplevels(stages.levels[union(stage.ind$`stage i`, stage.ind$`stage iv`)]))
tumor.fpqm.varSelRF.1.4.strat <-  varSelRF(xdata = t(exp_fpqm_tumor_reported[,union(stage.ind$`stage i`, stage.ind$`stage iv`)]),
                                      Class = droplevels(stages.levels[union(stage.ind$`stage i`, stage.ind$`stage iv`)]))
tumor.fpqm.varSelRF.2.3.4 <- varSelRF(xdata = t(exp_fpqm_tumor_reported[,Reduce(union,stage.ind[-1])]),
                                      Class = droplevels(stages.levels[Reduce(union,stage.ind[-1])]))
tumor.fpqm.varSelRF.2.3.4.strat <- varSelRF(xdata = t(exp_fpqm_tumor_reported[,Reduce(union,stage.ind[-1])]),
                                      Class = droplevels(stages.levels[Reduce(union,stage.ind[-1])]),
                                      strata = droplevels(stages.levels[Reduce(union, stage.ind[-1])]), sampsize = c(22,25,15))
tumor.fpqm.varSelRF.1.2 <- varSelRF(xdata = t(exp_fpqm_tumor_reported[,Reduce(union,stage.ind[c(1,2)])]),
                                    Class = droplevels(stages.levels[Reduce(union,stage.ind[c(1,2)])]))
tumor.fpqm.varSelRF.1.2.strat.1 <- varSelRF(xdata = t(exp_fpqm_tumor_reported[,Reduce(union,stage.ind[c(1,2)])]),
                                    Class = droplevels(stages.levels[Reduce(union,stage.ind[c(1,2)])]),
                                    strata = droplevels(stages.levels[Reduce(union,stage.ind[c(1,2)])]), 
                                    sampsize = c(22,22))
tumor.fpqm.varSelRF.1.2.strat.2 <- varSelRF(xdata = t(exp_fpqm_tumor_reported[,Reduce(union,stage.ind[c(1,2)])]),
                                    Class = droplevels(stages.levels[Reduce(union,stage.ind[c(1,2)])]),
                                    strata = droplevels(stages.levels[Reduce(union,stage.ind[c(1,2)])]), 
                                    sampsize = c(22,22))
tumor.fpqm.varSelRF.1.2.strat.3 <- varSelRF(xdata = t(exp_fpqm_tumor_reported[,Reduce(union,stage.ind[c(1,2)])]),
                                    Class = droplevels(stages.levels[Reduce(union,stage.ind[c(1,2)])]),
                                    strata = droplevels(stages.levels[Reduce(union,stage.ind[c(1,2)])]), 
                                    sampsize = c(22,22))

tumor.fpqm.varSelRF.3.4 <- varSelRF(xdata = t(exp_fpqm_tumor_reported[,Reduce(union,stage.ind[c(3,4)])]),
                                    Class = droplevels(stages.levels[Reduce(union,stage.ind[c(3,4)])]))
                                      
tumor.fpqm.varSelRF.3.4.strat.1 <- varSelRF(xdata = t(exp_fpqm_tumor_reported[,Reduce(union,stage.ind[c(3,4)])]),
                                    Class = droplevels(stages.levels[Reduce(union,stage.ind[c(3,4)])]),
                                    strata = droplevels(stages.levels[Reduce(union,stage.ind[c(3,4)])]), 
                                    sampsize = c(20,15))
tumor.fpqm.varSelRF.3.4.strat.2 <- varSelRF(xdata = t(exp_fpqm_tumor_reported[,Reduce(union,stage.ind[c(3,4)])]),
                                            Class = droplevels(stages.levels[Reduce(union,stage.ind[c(3,4)])]),
                                            strata = droplevels(stages.levels[Reduce(union,stage.ind[c(3,4)])]), 
                                            sampsize = c(30,15))

tumor.fpqm.varSelRF.2.3 <- varSelRF(xdata = t(exp_fpqm_tumor_reported[,Reduce(union,stage.ind[c(2,3)])]),
                                    Class = droplevels(stages.levels[Reduce(union,stage.ind[c(2,3)])]))

tumor.fpqm.varSelRF.2.3.strat.1 <- varSelRF(xdata = t(exp_fpqm_tumor_reported[,Reduce(union,stage.ind[c(2,3)])]),
                                            Class = droplevels(stages.levels[Reduce(union,stage.ind[c(2,3)])]),
                                            strata = droplevels(stages.levels[Reduce(union,stage.ind[c(2,3)])]), 
                                            sampsize = c(22,25))
tumor.fpqm.varSelRF.1.3 <- varSelRF(xdata = t(exp_fpqm_tumor_reported[,Reduce(union,stage.ind[c(1,3)])]),
                                    Class = droplevels(stages.levels[Reduce(union,stage.ind[c(1,3)])]))
tumor.fpqm.varSelRF.1.3.strat.2 <- varSelRF(xdata = t(exp_fpqm_tumor_reported[,Reduce(union,stage.ind[c(1,3)])]),
                                    Class = droplevels(stages.levels[Reduce(union,stage.ind[c(1,3)])]),
                                    strata = droplevels(stages.levels[Reduce(union,stage.ind[c(1,3)])]), sampsize = c(102,51))
tumor.fpqm.varSelRF.1.3.strat.1 <- varSelRF(xdata = t(exp_fpqm_tumor_reported[,Reduce(union,stage.ind[c(1,3)])]),
                                            Class = droplevels(stages.levels[Reduce(union,stage.ind[c(1,3)])]),
                                            strata = droplevels(stages.levels[Reduce(union,stage.ind[c(1,3)])]), sampsize = c(60,51))

tumor.fpqm.varSelRF.2.4 <- varSelRF(xdata = t(exp_fpqm_tumor_reported[,Reduce(union,stage.ind[c(2,4)])]),
                                    Class = droplevels(stages.levels[Reduce(union,stage.ind[c(2,4)])]))
                                    
tumor.fpqm.varSelRF.all <- varSelRF(xdata = t(exp_fpqm_tumor_reported[all.genes,]),
                                    Class = stages.levels)
tumor.fpqm.varSelRF.all <- varSelRF(xdata = t(exp_fpqm_tumor_reported[all.genes,]),
                                    Class = stages.levels, strata = stages.levels,
                                    sampsize = c(30,15,20,15))
tumor.fpqm.varSelRF.diff5 <- varSelRF(xdata = t(exp_fpqm_tumor_reported[diff.genes[[5]],]),
                                      Class = stages.levels, strata = stages.levels, sampsize = c(20,16,18,15)
                                      )
rf.fpqm.all <- randomForest(x = t(exp_fpqm_tumor_reported[tumor.fpqm.varSelRF.all$selected.vars,]),
                            y = stages.levels, strata = stages.levels, sampsize = c(30,15,20,15))

save(tumor.fpqm.smt.comb.varSelRF, file = 'environment/tumor_fpqm_smote_comb_var.RData')
save(tumor.nt.smt.comb.varSelRF, file = 'environment/tumor_nt_smote_comb_var.RData')
save(tumor.fpqm.under.comb.varSelRF, file = 'environment/tumor_fpqm_under_comb_var.RData')
save(tumor.nt.under.comb.varSelRF, file = 'environment/tumor_nt_under_comb_var.RData')
save(tumor.fpqm.varSelRF.samp.rep, file = 'environment/tumor_fpqm_var_strat_rep.RData')
save(tumor.fpqm.varSelRF.samp.not.rep, file = 'environment/tumor_fpqm_var_strat_not_rep.RData')
save(tumor.fpqm.varSelRF.wt, file = 'environment/tumor_fpqm_var_wt.RData')
save(tumor.fpqm.varSelRF.1.4, file = 'environment/tumor_fpqm_var_1_4.RData')
save(tumor.fpqm.varSelRF.1.4.strat, file = 'environment/tumor_fpqm_var_1_4_strat.RData')
save(tumor.fpqm.varSelRF.comb.wt, file = 'environment/tumor_fpqm_var_comb_wt.RData')
save(tumor.fpqm.varSelRF.comb.strat, file = 'environment/tumor_fpqm_var_comb_strat.RData')
save(tumor.fpqm.varSelRF.1.2, file = 'environment/tumor_fpqm_var_1_2.RData')
save(tumor.fpqm.varSelRF.1.2.strat.1, file = 'environment/tumor_fpqm_var_1_2_strat_1.RData')
save(tumor.fpqm.varSelRF.1.2.strat.2, file = 'environment/tumor_fpqm_var_1_2_strat_2.RData')
save(tumor.fpqm.varSelRF.1.2.strat.3, file = 'environment/tumor_fpqm_var_1_2_strat_3.RData')
save(tumor.fpqm.varSelRF.3.4, file = 'environment/tumor_fpqm_var_3_4.RData')
save(tumor.fpqm.varSelRF.2.3.4, file = 'environment/tumor_fpqm_var_2_3_4.RData')
save(tumor.fpqm.varSelRF.2.3.4.strat, file = 'environment/tumor_fpqm_var_2_3_4_strat.RData')
save(tumor.fpqm.varSelRF.3.4, file = 'environment/tumor_fpqm_var_3_4.RData')
save(tumor.fpqm.varSelRF.2.4, file = 'environment/tumor_fpqm_var_2_4.RData')
save(tumor.fpqm.varSelRF.1.3, file = 'environment/tumor_fpqm_var_1_3.RData')
save(tumor.fpqm.varSelRF.1.3.strat.1, file = 'environment/tumor_fpqm_var_1_3_strat_1.RData')
save(tumor.fpqm.varSelRF.1.3.strat.2, file = 'environment/tumor_fpqm_var_1_3_strat_2.RData')
stopCluster(cl)

load('environment/tumor_fpqm_var.RData')
load('environment/tumor_fpqm_var_1_2.RData')
load('environment/tumor_fpqm_var_1_2_strat_1.RData')
load('environment/tumor_fpqm_var_1_2_strat_2.RData')
load('environment/tumor_fpqm_var_1_2_strat_3.RData')
load('environment/tumor_fpqm_var_3_4.RData')
load('environment/tumor_fpqm_var_2_3_4.RData')
load('environment/tumor_fpqm_var_2_3_4_strat.RData')
load('environment/tumor_fpqm_var_1_4.RData')
load('environment/tumor_fpqm_var_1_4_strat.RData')
load('environment/tumor_fpqm_var_1_3.RData')
load('environment/tumor_fpqm_var_1_3_strat_1.RData')
load('environment/tumor_fpqm_var_1_3_strat_2.RData')
load('environment/tumor_fpqm_var_2_4.RData')
load('environment/tumor_fpqm_var_2_4.RData')
load('environment/tumor_fpqm_var_2_3.RData')
load('environment/tumor_fpqm_var_2_3_strat1.RData')
load('environment/tumor_fpqm_')
genes.varselrf = list()
genes.varselrf[['fpqm']] = remove.dots(unlist(strsplit(as.character(tumor.fpqm.varSelRF$selec.history$Vars.in.Forest[25]),
                                    " + ", fixed = T)))
genes.varselrf[['var_1_2']] = unlist(strsplit(as.character(tumor.fpqm.varSelRF.1.2$selec.history$Vars.in.Forest[25]),
                                    " + ", fixed = T))

####Trying Stratified sampling####


###
plot.lists <- function(varSelRF.lists, index)
{
  library(varSelRF)
  library(ggplot2)
  len <- length(varSelRF.lists)
  #print(len)
  par(mfrow=c(1,len))
  for(i in varSelRF.lists)
    plot(i,which=index)

}

plot.lists(list(tumor.fpqm.varSelRF, tumor.fpqm.log.varSelRF), index = 1)
plot.lists(list(tumor.fpqm.varSelRF, tumor.fpqm.log.varSelRF), index = 2)

plot.lists(list(tumor.fpqm.comb.varSelRF, tumor.fpqm.log.comb.varSelRF), index = 1)
plot.lists(list(tumor.fpqm.comb.varSelRF, tumor.fpqm.log.comb.varSelRF), index = 2)

plot.lists(list(tumor.fpqm.over.varSelRF), index = 1)
plot.lists(list(tumor.fpqm.over.varSelRF), index = 2)

plot.lists(list(tumor.fpqm.smt.comb.varSelRF, tumor.nt.smt.comb.varSelRF), index = 1)
plot.lists(list(tumor.fpqm.smt.comb.varSelRF, tumor.nt.smt.comb.varSelRF), index = 2)

plot.lists(list(tumor.fpqm.under.comb.varSelRF, tumor.nt.under.comb.varSelRF), index = 1)
plot.lists(list(tumor.fpqm.under.comb.varSelRF, tumor.nt.under.comb.varSelRF), index = 2)


plot(tumor.exp.boot, subject.names = df.stages$short.id)
selProbPlot(tumor.exp.boot)


dir.copy = paste(getwd(), 'results/tumor/varSelRF', sep = '/')

tumor.exp.varSelRF$selected.vars
tumor.fpqm.varSelRF$selected.vars
tumor.fpqm.log.varSelRF$selected.vars

tumor.exp.varSelRF$initialOrderedImportances[1:30]
tumor.fpqm.varSelRF$initialOrderedImportances[1:30]
tumor.fpqm.log.varSelRF$initialOrderedImportances[1:30]

##BootStrap
tumor.exp.boot$bootstrap.pred.error
tumor.fpqm.boot$bootstrap.pred.error
tumor.fpqm.log.boot$bootstrap.pred.error



tumor.exp.boot$all.data.randomForest$importance
tumor.fpqm.boot$all.data.randomForest$importance
tumor.fpqm.log.boot$all.data.randomForest$importance

View(tumor.exp.boot$class.predictions)

tumor.exp.boot$number.of.vars
tumor.fpqm.boot$number.of.vars
tumor.fpqm.log.boot$number.of.vars

tumor.exp.boot$overlap
tumor.fpqm.boot$overlap
tumor.fpqm.log.boot$overlap

save(tumor.exp.boot, file = 'environment/varSelRF_exp_boot.RData')
save(tumor.exp.varSelRF, file = 'environment/varSelRF_exp.RData')

save(tumor.fpqm.log.boot, file = 'environment/varSelRF_fpqm_log_boot.RData')
save(tumor.fpqm.log.varSelRF, file = 'environment/varSelRF_fpqk_log.RData')

save(tumor.fpqm.boot, file = 'environment/varSelRF_fpqm_boot.RData')
save(tumor.fpqm.varSelRF, file = 'environment/varSelRF_fpqm.RData')

save(tumor.fpqm.log.comb.boot, file = 'environment/tumor_fpqm_log_comb_boot.RData')
save(tumor.fpqm.log.comb.varSelRF, file = 'environment/tumor_fpqm_log_comb_varSelRF.RData')
save(tumor.fpqm.comb.boot, file = 'environment/tumor_fpqm_comb_boot.RData')
save(tumor.fpqm.comb.varSelRF, file = 'environment/tumor_fpqm_comb_varSelRF.RData')


save(tumor.indexes, file = 'environment/tumor_indexes.RData')
save(diff.ids, file = 'environment/diff_ids.RData') ##Used for tumor-control
save(indexes.stages.reported, file = 'environment/indexes_stages_reported.RData') ##Contain indexes for 
#which stage has been reported once the tumor.indexes are known

save(exp_prof, file = 'environment/exp_prof.RData' )
save(exp_fpqm, file = 'environment/exp_fpqm.RData' )


###Combining Stages
stages.levels.comb <- sapply(as.character(stages.levels), function(x)
  {
  if(x == 'stage i' || x == 'stage ii')
    x = 'stage i'
  else
    x = 'stage iv'
})
save(stages.levels.comb, file = 'environment/stages.level.comb.RData' )
