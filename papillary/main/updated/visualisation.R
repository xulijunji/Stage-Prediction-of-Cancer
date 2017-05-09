source('main/updated/initialisation.R')
source('../function.R')
load('environment/sample_info.RData')
load('environment/dds.RData')
load('environment/accuracy_feature/vs_nor_comb.RData')
library(RColorBrewer)
library(ggplot2)

vs_match_tum <- vst(dds)
rownames(vs_match_tum) = remove.dots(rownames(vs_match_tum))
plotPCA(vs_match_tum[diff.genes$`5`,], intgroup = 'type', title = '5 fold genes', colData = colData(dds))
plotPCA(t(vs_normal_comb_reported[, diff.genes$`5`]), intgroup = 'stage.type', title = '5 fold genes', colData = sample.info.all.rep)
plotPCA(t(vst_tumor_tum[,diff.genes$`5`]), intgroup = 'stage.type', title = '5 fold genes', 
        colData = sample.info.all.rep[tumor.ind.vs,])
plotPCA(t(vst_tumor_tum[,net.features.updated$varSelRF$atleast_2]), intgroup = 'stage.type', title = '5 fold genes', 
        colData = data.frame(stage.type=stages.levels.comb))

col <- brewer.pal(9, 'PuRd')
create.heatmap(vs_normal_comb_reported, sample.info.all.rep$stage.type, diff.genes$`5`, 
               'Heatmap using 5 fold genes', col = col, 
               cluster_rows = T, cluster_cols = T)
create.heatmap(t(counts(dds_tumor_reported)), stages.levels.comb, colnames(vst_tumor_tum), 
               'Heatmap using 1 fold genes', col = col, 
               cluster_rows = T, cluster_cols = F)


a <- create.net.df(cv.results)
aucs <- summarySE(a, measurevar = 'AUC',
                groupvars = c('Classifier', 'Feature_Selection'))
accuracy <- summarySE(a, measurevar = 'Accuracy',
                  groupvars = c('Classifier', 'Feature_Selection'))
specificity <- summarySE(a, measurevar = 'Specificity',
                  groupvars = c('Classifier', 'Feature_Selection'))
sensitivity <- summarySE(a, measurevar = 'Sensitivity',
                      groupvars = c('Classifier', 'Feature_Selection'))
f_val <- summarySE(a, measurevar = 'F_val',
                  groupvars = c('Classifier', 'Feature_Selection'))

pd <- position_dodge(0.3)
pe.auc <- ggplot(aucs, aes(x=Classifier, y=AUC, colour=Feature_Selection)) + 
  geom_errorbar(aes(ymin=AUC-se, ymax=AUC+se), width=.1,  position=pd) +
  geom_line()+
  geom_point(position=pd, size =3)
pe.acc <- ggplot(accuracy, aes(x=Classifier, y=Accuracy, colour=Feature_Selection)) + 
  geom_errorbar(aes(ymin=Accuracy-se, ymax=Accuracy+se), width=.1,  position=pd) +
  geom_line()+
  geom_point(position=pd, size =3)
pe.sens <- ggplot(sensitivity, aes(x=Classifier, y=Sensitivity, colour=Feature_Selection)) + 
  geom_errorbar(aes(ymin=Sensitivity-se, ymax=Sensitivity+se), width=.1,  position=pd) +
  geom_line()+
  geom_point(position=pd, size =3)
pe.spec <- ggplot(specificity, aes(x=Classifier, y=Specificity, colour=Feature_Selection)) + 
  geom_errorbar(aes(ymin=Specificity-se, ymax=Specificity+se), width=.1,  position=pd) +
  geom_line()+
  geom_point(position=pd, size =3)
pe.fval <- ggplot(f_val, aes(x=Classifier, y=F_val, colour=Feature_Selection)) + 
  geom_errorbar(aes(ymin=F_val-se, ymax=F_val+se), width=.1,  position=pd) +
  geom_line()+
  geom_point(position=pd, size =3)
multiplot(pe.auc, pe.acc, pe.sens, pe.spec, pe.fval, cols =2)

p1 <- ggplot(a, aes(x=Classifier, y=AUC, fill=Feature_Selection)) +
  geom_boxplot()
p2 <- ggplot(a, aes(x=Classifier, y=Accuracy, fill=Feature_Selection)) +
  geom_boxplot()
p3 <- ggplot(a, aes(x=Classifier, y=Specificity, fill=Feature_Selection)) +
  geom_boxplot()
p4 <- ggplot(a, aes(x=Classifier, y=Sensitivity, fill=Feature_Selection)) +
  geom_boxplot()
p5 <- ggplot(a, aes(x=Classifier, y=F_val, fill=Feature_Selection)) +
  geom_boxplot()
multiplot(p1,p2,p3,p4,p5, cols =2)

