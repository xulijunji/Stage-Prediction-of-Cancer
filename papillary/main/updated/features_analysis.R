library(VennDiagram)
load('environment/accuracy_feature/updated/net_features_updated.RData')
load('environment/accuracy_feature/updated/cv_model.RData')
load('environment/accuracy_feature/updated/cv_results.RData')
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

####HeatMaps
library(RColorBrewer)
library(pheatmap)
col <- colorRampPalette(rev(brewer.pal(9, 'RdYlBu')))(100)
df.tumor <- create.ordered.annotation(stages.levels.comb[],
                                      rownames(vst_tumor_tum[,]))
breaks = c(seq(5,10, length.out = 70), seq(10.1,24, length.out = 30))

genes = net.features.updated$deseq2$atleast_1$`1 fold`
g.int = intersect(net.features.updated$deseq2$atleast_2$`1 fold`,
      intersect(net.features.updated$shrunken$atleast_2, net.features.updated$varSelRF$atleast_2))
data = vst_tumor_tum[train.indexes,late.stage.genes]
data = vst_tumor_tum[,wcgna.genes.ent.ens]
clus_rows = run_hclust_on_a_matrix(data)
clus_cols = run_hclust_on_a_matrix(t(data))
pheatmap(t(data),
#         annotation_col = df.tumor, 
#         cluster_rows = clus_rows,
#         cluster_cols = clus_cols,
          cluster_cols = T,
          breaks = breaks,
#          main = 'Heatmap using intersection for atleast 2 and 1 fold',
         show_rownames = F, show_colnames = F)

create.heatmap(vst_tumor_tum, stages.levels.comb, g, 't',
               col = col)
create.heatmap(vst_tumor_tum, stages.levels.comb, g, 
               '2 fold', col = col)

##Classifiers
temp.class <- list()
temp.class[['shrunken']] <- pamr.train(list(x = t(as.matrix(vst_tumor_tum[train.indexes, g])), 
                                            y = stages.levels.comb[train.indexes]))
temp.class[['svm']] <-  svm(x = vst_tumor_tum[train.indexes, g], y = stages.levels.comb[train.indexes])
temp.class[['nb']] <- naiveBayes(x = vst_tumor_tum[train.indexes, g], y = stages.levels.comb[train.indexes])
temp.class[['rf']] <- randomForest(x = vst_tumor_tum[train.indexes, g], y = stages.levels.comb[train.indexes])



###Predicted
pred.cv.class <- list()
pred.cv.class[['shrunken']] <- cv.shrunken(vst_tumor_tum, 10, list(g), list(temp.class$shrunken), 
                                           train.indexes, stages.levels.comb)
pred.cv.class[['knn']] <- find.best.k(vst_tumor_tum[train.indexes, g], 10, stages.levels.comb[train.indexes])

pred.test.class <- list()
pred.test.class[['shrunken']] <- pamr.predict(temp.class$shrunken, t(as.matrix(vst_tumor_tum[test.indexes, g])),
                                     threshold = temp.class$shrunken$threshold[pred.cv.class$shrunken$thr])
pred.test.class[['svm']] <- predict(temp.class$svm, vst_tumor_tum[test.indexes, g])
pred.test.class[['rf']] <- predict(temp.class$rf, vst_tumor_tum[test.indexes, g])
pred.test.class[['nb']] <- predict(temp.class$nb, vst_tumor_tum[test.indexes, g])
pred.test.class[['knn']] <- knn(vst_tumor_tum[train.indexes, g], vst_tumor_tum[test.indexes, g],
                            stages.levels.comb[train.indexes], pred.cv.class$knn$best_k)

##
res.test <- lapply(pred.test.class, function(pred.test){
    get.eval(stages.levels.comb[test.indexes], pred.test)
})

g_1.5_3 <- intersect(intersect(net.features.updated$shrunken$atleast_3, net.features.updated$varSelRF$atleast_3),
               net.features.updated$deseq2$atleast_3$`1.5 fold`) 
g_1.5_2 <- intersect(intersect(net.features.updated$shrunken$atleast_2, net.features.updated$varSelRF$atleast_2),
                     net.features.updated$deseq2$atleast_2$`1.5 fold`) 
g_1.5_1 <- intersect(intersect(net.features.updated$shrunken$atleast_1, net.features.updated$varSelRF$atleast_1),
                     net.features.updated$deseq2$atleast_1$`1.5 fold`) 

g_1_3 <- intersect(intersect(net.features.updated$shrunken$atleast_3, net.features.updated$varSelRF$atleast_3),
                     net.features.updated$deseq2$atleast_3$`1 fold`) 
g_1_2 <- intersect(intersect(net.features.updated$shrunken$atleast_2, net.features.updated$varSelRF$atleast_2),
                   net.features.updated$deseq2$atleast_2$`1 fold`) 
g_1_1 <- intersect(intersect(net.features.updated$shrunken$atleast_1, net.features.updated$varSelRF$atleast_1),
                   net.features.updated$deseq2$atleast_1$`1 fold`) 

g_2_1 <- intersect(intersect(net.features.updated$shrunken$atleast_1, net.features.updated$varSelRF$atleast_1),
                   net.features.updated$deseq2$atleast_1$`2 fold`)
g_2_2 <- intersect(intersect(net.features.updated$shrunken$atleast_2, net.features.updated$varSelRF$atleast_2),
                   net.features.updated$deseq2$atleast_2$`2 fold`)
g_2_3 <- intersect(intersect(net.features.updated$shrunken$atleast_3, net.features.updated$varSelRF$atleast_3),
                   net.features.updated$deseq2$atleast_3$`2 fold`)

a <- final.res(vst_tumor_tum, train.indexes, test.indexes, stages.levels.comb, g, 10)
res_1.5_3 <- final.res(vst_tumor_tum, train.indexes, test.indexes, stages.levels.comb, g_1.5_3, 10)
res_1.5_2 <- final.res(vst_tumor_tum, train.indexes, test.indexes, stages.levels.comb, g_1.5_2, 10)
res_1.5_1 <- final.res(vst_tumor_tum, train.indexes, test.indexes, stages.levels.comb, g_1.5_1, 10)

res_1_3 <- final.res(vst_tumor_tum, train.indexes, test.indexes, stages.levels.comb, g_1_3, 10)
res_1_2 <- final.res(vst_tumor_tum, train.indexes, test.indexes, stages.levels.comb, g_1_2, 10)
res_1_1 <- final.res(vst_tumor_tum, train.indexes, test.indexes, stages.levels.comb, g_1_1, 10)

res_2_1 <- final.res(vst_tumor_tum, train.indexes, test.indexes, stages.levels.comb, g_2_1, 10)
res_2_2 <- final.res(vst_tumor_tum, train.indexes, test.indexes, stages.levels.comb, g_2_2, 10)

get.conf.mat(a[[2]])



###First finding the intersection between the degs b/w normal and tumor and b/w stages
length(g_1_1)
length(intersect(diff.genes$`1`, g_1_1))
length(intersect(diff.genes$`2`, g_1_1))

length(g_1_2)
length(intersect(diff.genes$`1`, g_1_2))
length(intersect(diff.genes$`2`, g_1_2))

length(g_1.5_1)
length(intersect(diff.genes$`1`, g_1.5_1))
length(intersect(diff.genes$`2`, g_1.5_1))

early.stage.genes <- intersect(diff.genes$`1`, g_1.5_1)
late.stage.genes <- setdiff(g_1.5_1, intersect(diff.genes$`1`, g_1.5_1))

col <- colorRampPalette(rev(brewer.pal(9, 'RdYlBu')))(100)
df.nor.tumor <- create.ordered.annotation(sample.info.all.rep$comb.stage,
                                      rownames(vs_normal_comb_reported))
sample.info.all.rep$comb.stage <- as.factor(sapply(sample.info.all.rep$stage.type, function(stage)
  {
  if(stage == 'stage i' | stage == 'stage ii')
    'early'
  else if(stage == 'stage iii' | stage == 'stage iv')
    'late'
  else
    'N'
}))
breaks = c(seq(5,10, length.out = 70), seq(10.1,24, length.out = 30))
data = vs_normal_comb_reported[,late.stage.genes]
clus_rows = run_hclust_on_a_matrix(data)
clus_cols = run_hclust_on_a_matrix(t(data))
pheatmap(t(data),
         annotation_col = df.nor.tumor, 
         #cluster_rows = clus_rows,
         #         cluster_cols = clus_cols,
         breaks = breaks,
         #          main = 'Heatmap using intersection for atleast 2 and 1 fold',
         show_rownames = F, show_colnames = F)



length(intersect(intersect(diff.genes$`1`, g_1_2) , intersect(diff.genes$`1`, g_1.5_1)))

r <- final.res(vst_tumor_tum, train.indexes, test.indexes, stages.levels.comb, 
          g_1_2, 10)
get.aucs(r[[2]])

####WGCNA#######
####Hub Genes
wcgna.genes.ent <- c(22974, 6241, 11065, 990, 699, 10112, 9232, 55165, 9212,
                     55355, 259266, 146909, 64151, 4751, 24137, 890, 55388)
wcgna.genes.ent.ens <- 	c('ENSG00000112984', 'ENSG00000175063',	'ENSG00000186185',
                          'ENSG00000088325', 'ENSG00000090889', 'ENSG00000066279', 
                          'ENSG00000117650', 'ENSG00000138180',	'ENSG00000123485',	
                          'ENSG00000065328', 'ENSG00000171848', 'ENSG00000109805', 
                          'ENSG00000169679', 'ENSG00000145386',	'ENSG00000178999',
                          'ENSG00000164611', 'ENSG00000094804')
  
intersect(g_2_1, wcgna.genes.ent.ens)
intersect(net.features.updated$deseq2$atleast_1$`2 fold`, wcgna.genes.ent.ens)
res.wcgna <- final.res(vst_tumor_tum, train.indexes, test.indexes, stages.levels.comb, wcgna.genes.ent.ens, 
                       10)
g = intersect(net.features.updated$shrunken$atleast_1, diff.genes$`2`)
res.g <- final.res(vst_tumor_tum, train.indexes, test.indexes, stages.levels.comb, g, 
                   10)

get.aucs(res.wcgna[[1]])
get.aucs(res.wcgna[[2]])
get.aucs(res_1_2[[1]])
get.aucs(res_1_2[[2]])
get.aucs(res_1_1[[2]])
get.aucs(res_2_2[[2]])
get.aucs(res_2_1[[1]])
get.aucs(res_1.5_1[[1]])
get.aucs(res_1.5_2[[2]])
get.aucs(res_1.5_3[[2]])
get.aucs(res_1_3[[2]])
get.aucs(res.g[[2]])
intersect(g_1_1, wcgna.genes.ent.ens)

library(ggplot2)

create.boxplots <- function(gene, data, stage)
{
  boxplots <- list()
  for(i in seq_along(gene))
  {
    box.df <- data.frame(gene = data[,gene[i]], stage = stage)
    g <- ggplot(box.df, aes(x = stage, y = gene))+geom_boxplot()+
      ggtitle(gene[i])
    boxplots[[i]] <- g
  }
  return(boxplots)
}
boxplots.g_2_2 <- create.boxplots(g_2_2, vs_normal_comb_reported, sample.info.all.rep$stage.type)
multiplot(boxplots.g_2_2, col = 3)

boxplots.wcgna <- create.boxplots(wcgna.genes.ent.ens, vs_normal_comb_reported, sample.info.all.rep$stage.type)
multiplot(boxplots.wcgna, col = 4)

boxplot.g_1.5_1 <- create.boxplots(g_1.5_1[55], vs_normal_comb_reported, 
                                   sample.info.all.rep$stage.type)
multiplot(boxplot.g_1.5_1, col = 1)
