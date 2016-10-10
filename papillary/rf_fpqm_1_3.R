rf.fpqm.1.3.samples <- lapply(seq(1:10), function(z)
{
  z = randomForest(x = t(exp_fpqm_tumor_reported[tumor.fpqm.varSelRF.1.3$selected.vars,
                                                                  Reduce(union,stage.ind[c(1,3)])]),
  y=droplevels(stages.levels[Reduce(union,stage.ind[c(1,3)])]),sampsize = c(102,51))

which(z$predicted == 'stage iii' & stages.levels[Reduce(union,stage.ind[c(1,3)])]
                                                      == 'stage iii')
})