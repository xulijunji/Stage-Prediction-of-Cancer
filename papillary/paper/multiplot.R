library(ggplot2)
library(gridExtra)
library(grid)

grid_arrange_shared_legend <- function(plots, ncol = length(plots), nrow = 1, position = c("bottom", "right")) 
{
  
  position <- match.arg(position)
  g <- ggplotGrob(plots[[1]] + theme(legend.position = position))$grobs
  legend <- g[[which(sapply(g, function(x) x$name) == "guide-box")]]
  lheight <- sum(legend$height)
  lwidth <- sum(legend$width)
  gl <- lapply(plots, function(x) x + theme(legend.position="none"))
  gl <- c(gl, ncol = ncol, nrow = nrow)
  
  combined <- switch(position,
                     "bottom" = arrangeGrob(do.call(arrangeGrob, gl),
                                            legend,
                                            ncol = 1,
                                            heights = unit.c(unit(1, "npc") - lheight, lheight)),
                     "right" = arrangeGrob(do.call(arrangeGrob, gl),
                                           legend,
                                           ncol = 2,
                                           widths = unit.c(unit(1, "npc") - lwidth, lwidth)))
  
  grid.newpage()
  grid.draw(combined)
  
  # return gtable invisibly
  invisible(combined)
  
}
# names(test.trial.results) <- c('Shrunken', 'VarSelRF', 'DeSeq2 log2FC 1', 'DeSeq2 log2FC 1.5', 
#                     'DeSeq2 log2FC 2', 'SamSeq log2FC 1', 'SamSeq log2FC 1.5', 'SamSeq log2FC 2')
# test.df <- create.net.df(test.trial.results)
# p1 <- ggplot(test.df, aes(x=Classifier, y=AUC, fill=Feature_Selection)) +
#   geom_boxplot()
# p2 <- ggplot(test.df, aes(x=Classifier, y=Accuracy, fill=Feature_Selection)) +
#   geom_boxplot()
# p3 <- ggplot(test.df, aes(x=Classifier, y=Specificity, fill=Feature_Selection)) +
#   geom_boxplot()
# p4 <- ggplot(test.df, aes(x=Classifier, y=Sensitivity, fill=Feature_Selection)) +
#   geom_boxplot()
# p5 <- ggplot(test.df, aes(x=Classifier, y=F_val, fill=Feature_Selection)) +
#   geom_boxplot()
# 
# grid_arrange_shared_legend(list(p1,p2,p3,p4,p5), ncol = 2, nrow = 3)
# 
# 
# names(cv.micr.res) <- c('Shrunken', 'VarSelRF', 'DeSeq2 log2FC 1', 'DeSeq2 log2FC 1.5', 
#                        'DeSeq2 log2FC 2', 'SamSeq log2FC 1', 'SamSeq log2FC 1.5', 'SamSeq log2FC 2')
# micr.net.df <-create.net.df(cv.micr.res)
# p.micr.auc <-ggplot(micr.net.df, aes(x=Classifier, y=AUC, fill=Feature_Selection)) + geom_boxplot()
# grid_arrange_shared_legend(list(p.micr.auc))

load('environment/accuracy_feature/updated/new_data/cv_results.RData')
cv.net.df <- create.net.df(cv.trial.results)
p1 <- ggplot(cv.net.df, aes(x=Classifier, y=AUC, fill=Feature_Selection)) +
  geom_boxplot()
p2 <- ggplot(cv.net.df, aes(x=Classifier, y=Accuracy, fill=Feature_Selection)) +
  geom_boxplot()
p3 <- ggplot(cv.net.df, aes(x=Classifier, y=Specificity, fill=Feature_Selection)) +
  geom_boxplot()
p4 <- ggplot(cv.net.df, aes(x=Classifier, y=Sensitivity, fill=Feature_Selection)) +
  geom_boxplot()
p5 <- ggplot(cv.net.df, aes(x=Classifier, y=F_val, fill=Feature_Selection)) +
  geom_boxplot()
grid_arrange_shared_legend(list(p1, p2, p3, p4, p5), ncol = 2, nrow = 3)
