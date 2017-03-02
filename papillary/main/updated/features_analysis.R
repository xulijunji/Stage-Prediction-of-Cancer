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
#1.5 fold
venn.diagram(create.list.venn(net.features.updated, 1.5, 1), filename = 'images/tumor/feature_updated/1_1.5fold.tiff')
venn.diagram(create.list.venn(net.features.updated, 1.5, 2), filename = 'images/tumor/feature_updated/2_1.5fold.tiff')
venn.diagram(create.list.venn(net.features.updated, 1.5, 3), filename = 'images/tumor/feature_updated/3_1.5fold.tiff')
venn.diagram(create.list.venn(net.features.updated, 1.5, 4), filename = 'images/tumor/feature_updated/4_1.5fold.tiff')

#1 fold
venn.diagram(create.list.venn(net.features.updated, 1, 1), filename = 'images/tumor/feature_updated/1_1fold.tiff')
venn.diagram(create.list.venn(net.features.updated, 1, 2), filename = 'images/tumor/feature_updated/2_1fold.tiff')
venn.diagram(create.list.venn(net.features.updated, 1, 3), filename = 'images/tumor/feature_updated/3_1fold.tiff')
venn.diagram(create.list.venn(net.features.updated, 1, 4), filename = 'images/tumor/feature_updated/4_1fold.tiff')

#2 fold
venn.diagram(create.list.venn(net.features.updated, 2, 1), filename = 'images/tumor/feature_updated/1_2fold.tiff')
venn.diagram(create.list.venn(net.features.updated, 2, 2), filename = 'images/tumor/feature_updated/2_2fold.tiff')
venn.diagram(create.list.venn(net.features.updated, 2, 3), filename = 'images/tumor/feature_updated/3_2fold.tiff')
venn.diagram(create.list.venn(net.features.updated, 2, 4), filename = 'images/tumor/feature_updated/4_2fold.tiff')
