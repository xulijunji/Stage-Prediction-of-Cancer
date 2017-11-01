load('environment/accuracy_feature/updated/new_data/fea.RData')
library(VennDiagram)

colors()[grep('red', colors())]

create.list.venn <- function(features.list, fold, group)
{
  fold <- as.character(fold)
  group <- as.character(group)
  g1 <- features.list$shrunken[[paste0('atleast_',group)]]
  g2 <- features.list$varSelRF[[paste0('atleast_',group)]]
  g3 <- features.list[[paste0('deseq2',fold, ' fold')]][[paste0('atleast_',group)]]
  g4 <- features.list[[paste0('sam',fold, ' fold')]][[paste0('atleast_',group)]]
  req.list <- list(g1,g2,g3,g4)
  names(req.list) <- c('Shrunken', 'VarSelRF', paste0('Deseq2_log2FC ',fold), paste0('SamSeq_log2FC ',fold))
  return(req.list)
}

venn.diagram(create.list.venn(fea.trial.list, 1.5, 4), filename = 'paper/images/figureS3_A', 
                  fill = c("purple", "indianred2", "pink", 'green'), imagetype = 'png',
                main = 'Intersection of genes across different feature selection \nalgorithms for log2FC ', main.cex = 1.5)
venn.diagram(create.list.venn(fea.trial.list, 2, 4), filename = 'paper/images/figureS3_B', 
             fill = c("purple", "indianred2", "pink", 'green'), imagetype = 'png',
             main = 'Intersection of genes across different feature selection \nalgorithms for group atleast 4', main.cex = 1.5)

grid.draw(v)
write(Reduce(intersect, create.list.venn(fea.trial.list, 1, 4)), 'g_int_4.txt')
