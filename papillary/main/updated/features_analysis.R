library(VennDiagram)
sapply(net.features.updated$shrunken[c(3,4,5,6)], length)
sapply(net.features.updated$varSelRF[c(3,4,5,6)], length)
lapply(net.features.updated$deseq2[c(3,4,5,6)], function(x)
  {
  sapply(x, length) 
})

create.list.venn <- function(net.features, fold, group)
{
  fold <- as.character(fold)
  group <- as.character(group)
  g1 <- net.features$shrunken[[paste0('atleast_',group)]]
  g2 <- net.features$varSelRF[[paste0('atleast_',group)]]
  g3 <- net.features$deseq2[[paste0('atleast_',group)]][[paste0(fold, ' fold')]]
  req.list <- list(g1,g2,g3)
  names(req.list) <- c('Shrunken', 'VarSelRF', paste0('Deseq2_',fold))
  return(req.list)
}

venn.diagram(create.list.venn(net.features.updated, 1.5, 1), filename = 'images/tumor/feature_updated/1_1.5fold')
