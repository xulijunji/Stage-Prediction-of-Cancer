source('main/updated/initialisation.R')
load('environment/accuracy_feature/updated/net_features_trial.RData')
load('environment/accuracy_feature/updated/new_data/fea.RData')
get.stage.distribution(gr.trial.train, stages.levels.comb)
get.stage.distribution(list(train.trial.ind), stages.levels.comb)
get.stage.distribution(list(test.trial.ind), stages.levels.comb)

create.list.venn <- function(features.list, group, fold)
{
  group <- as.character(group)
  fold <- as.character(fold)
  g1 <- features.list$shrunken[[paste0('atleast_',group)]]
  g2 <- features.list$varSelRF[[paste0('atleast_',group)]]
  g3 <- features.list[[paste0('deseq2',fold,' fold')]][[paste0('atleast_',group)]]
  g4 <- features.list[[paste0('sam',fold,' fold')]][[paste0('atleast_',group)]]
  req.list <- list(g1,g2,g3,g4)
  names(req.list) <- c('Shrunken', 'VarSelRF', paste0('Deseq2 log2FC',fold), paste0('SamSeq log2FC ',fold))
  return(req.list)
}


library(VennDiagram)
venn.diagram(create.list.venn(fea.trial.list, 4, 1), filename = '~/h_short/papillary/paper/images/int_atleast4.png',imagetype = 'png')
venn.diagram(create.list.venn(fea.trial.list, 3, 1), filename = '~/h_short/papillary/paper/images/int_atleast3.png', imagetype = 'png')
Reduce(intersect, create.list.venn(fea.trial.list, 4, 1.5))
