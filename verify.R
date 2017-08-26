load('environment/accuracy_feature/updated/dds_tumor_reported.RData')
load('environment/accuracy_feature/updated/net_features_updated.RData')
load('environment/accuracy_feature/updated/gr.RData')
temp.df <- t(counts(dds_tumor_reported))

col <- brewer.pal(9, 'RdPu')

##2 fold
create.heatmap(vst_tumor_tum[train.indexes,], stages.levels.comb[train.indexes],
      net.features.updated$deseq2$atleast_3$`2 fold`, 
      'Heatmap using 2 fold genes', col = col, 
      cluster_rows = T, cluster_cols = F)
create.heatmap(vst_tumor_tum[train.indexes,], stages.levels.comb[train.indexes],
               net.features.updated$deseq2$atleast_3$`2 fold`, 
               'Heatmap using 2 fold genes', col = col, 
               cluster_rows = T, cluster_cols = T)

create.heatmap(temp.df[train.indexes,], stages.levels.comb[train.indexes],
               net.features.updated$deseq2$atleast_3$`2 fold`,
               'Heatmap using 2 fold genes', col = col,
               cluster_cols = T,
               #breaks = c(1e+02,1e+03,1e+04, 2e+04)
               )
get.clust.score(vst_tumor_tum[train.indexes,], 
                net.features.updated$deseq2$atleast_3$`2 fold`, 
                stages.levels.comb[train.indexes], 2)

##1.5 fold
create.heatmap(vst_tumor_tum[train.indexes,], stages.levels.comb[train.indexes],
               net.features.updated$deseq2$atleast_3$`1.5 fold`, 
               'Heatmap using 1.5 fold genes', col = col, 
               cluster_rows = T, cluster_cols = F)
get.clust.score(vst_tumor_tum[train.indexes,], 
                net.features.updated$deseq2$atleast_3$`1.5 fold`, 
                stages.levels.comb[train.indexes], 2)

# create.heatmap(vst_tumor_tum[test.indexes,], stages.levels.comb[test.indexes],
#                net.features.updated$deseq2$atleast_3$`1.5 fold`, 
#                'Heatmap using 2 fold genes', col = col, 
#                cluster_rows = T, cluster_cols = T)
create.heatmap(temp.df[train.indexes,], stages.levels.comb[train.indexes],
               net.features.updated$deseq2$atleast_3$`1.5 fold`,
               'Heatmap using 1.5 fold genes', col = col,
               cluster_cols = T,
               breaks = c(1e+02,1e+03,1e+04, 2e+04)
)

##1 fold
create.heatmap(vst_tumor_tum[train.indexes,], stages.levels.comb[train.indexes],
               net.features.updated$deseq2$atleast_3$`1 fold`, 
               'Heatmap using 1 fold genes', col = col, 
               cluster_rows = T, 
               cluster_cols = T)
create.heatmap(temp.df[train.indexes,], stages.levels.comb[train.indexes],
               net.features.updated$deseq2$atleast_3$`1 fold`,
               'Heatmap using 1 fold genes', col = col,
               cluster_cols = F, cluster_rows = T,
               breaks = c(1e+02,1e+03,1e+04, 2e+04)
)
get.clust.score(vst_tumor_tum[train.indexes,], 
                net.features.updated$deseq2$atleast_3$`1 fold`, 
                stages.levels.comb[train.indexes], 2)


##3 fold
create.heatmap(vst_tumor_tum[train.indexes,],
               stages.levels.comb[train.indexes],
               net.features.updated$deseq2$atleast_3$`3 fold`, 
               'Heatmap using 3 fold genes atleast 3', col = col, 
               cluster_rows = T, cluster_cols = T)
create.heatmap(vst_tumor_tum[train.indexes,],
               stages.levels.comb[train.indexes],
               net.features.updated$deseq2$atleast_2$`3 fold`, 
               'Heatmap using 3 fold genes atleast 2', col = col, 
               cluster_rows = T, cluster_cols = T)
create.heatmap(vst_tumor_tum[train.indexes,],
               stages.levels.comb[train.indexes],
               net.features.updated$deseq2$atleast_1$`3 fold`, 
               'Heatmap using 3 fold genes atleast 1', col = col, 
               cluster_rows = T, cluster_cols = T)
create.heatmap(vst_tumor_tum[train.indexes,],
               stages.levels.comb[train.indexes],
               net.features.updated$deseq2$atleast_4$`2.5 fold`, 
               'Heatmap using 3 fold genes', col = col, 
               cluster_rows = T, cluster_cols = T)
get.clust.score(vst_tumor_tum[train.indexes,], 
                net.features.updated$deseq2$atleast_3$`3 fold`, 
                stages.levels.comb[train.indexes], 2)

##2.5 fold
create.heatmap(vst_tumor_tum[train.indexes,],
               stages.levels.comb[train.indexes],
               net.features.updated$deseq2$atleast_3$`2.5 fold`, 
               'Heatmap using 2.5 fold genes', col = col, 
               cluster_rows = T, cluster_cols = T)
get.clust.score(vst_tumor_tum[train.indexes,], 
                net.features.updated$deseq2$atleast_3$`2.5 fold`, 
                stages.levels.comb[train.indexes], 2)


##Shrunken
create.heatmap(vst_tumor_tum[train.indexes,], stages.levels.comb[train.indexes],
               net.features.updated$shrunken$atleast_3, 
               'Heatmap using shrunken genes', col = col, 
               cluster_rows = T, cluster_cols = T)
create.heatmap(temp.df, stages.levels.comb,
               net.features.updated$shrunken$atleast_3, 
               'Heatmap using shrunken genes', col = col, 
               cluster_rows = T, cluster_cols = T,
               breaks = c(1e+02,1e+03,1e+04, 2e+04)
)
get.clust.score(vst_tumor_tum[train.indexes,], 
                net.features.updated$shrunken$atleast_3, 
                stages.levels.comb[train.indexes], 2)

##VarSelRF
create.heatmap(vst_tumor_tum[train.indexes,], stages.levels.comb[train.indexes],
               net.features.updated$varSelRF$atleast_1, 
               'Heatmap using varSelRf genes', col = col, 
               cluster_rows = T, cluster_cols = T)
create.heatmap(temp.df[train.indexes,], stages.levels.comb[train.indexes],
               net.features.updated$varSelRF$atleast_1, 
               'Heatmap using varSelRf genes', col = col, 
               cluster_rows = T, cluster_cols = T,
               breaks = c(1e+02,1e+03,1e+04, 2e+04)
)
get.clust.score(vst_tumor_tum[train.indexes,], 
                net.features.updated$varSelRF$atleast_3, 
                stages.levels.comb[train.indexes], 2)


###Atleast 3
##deseq2 1fold
g.at3.1 <- intersect(net.features.updated$shrunken$atleast_3,
            intersect(net.features.updated$varSelRF$atleast_3,
                      net.features.updated$deseq2$atleast_3$`1 fold`))
create.heatmap(vst_tumor_tum[train.indexes,], stages.levels.comb[train.indexes],
               g.at3.1, 
               'Heatmap using atleast 3 genes 1 fold', col = col, 
               cluster_rows = T, cluster_cols = T,show_rownames = F
                )
get.clust.score(vst_tumor_tum[train.indexes,], 
                net.features.updated$deseq2$atleast_3$`2 fold`, 
                stages.levels.comb[train.indexes], 2)
c.3.1 <- final.res(vst_tumor_tum, train.indexes, test.indexes, stages.levels.comb,
                g.at3.1, 10)
get.clust.score(vst_tumor_tum[train.indexes,], 
                g.at3.1, 
                stages.levels.comb[train.indexes], 2)

##deseq2 1.5 fold
g.at3.1.5 <- intersect(net.features.updated$shrunken$atleast_3,
                   intersect(net.features.updated$varSelRF$atleast_3,
                             net.features.updated$deseq2$atleast_3$`1.5 fold`))
create.heatmap(vst_tumor_tum[train.indexes,], stages.levels.comb[train.indexes],
               g.at3.1.5, 
               'Heatmap using atleast 3 genes 1.5 fold', col = col, 
               cluster_rows = T, cluster_cols = T)
c.3.1_5 <- final.res(vst_tumor_tum, train.indexes, test.indexes, stages.levels.comb,
                     g.at3.1.5, 10)
get.clust.score(vst_tumor_tum[train.indexes,], 
                g.at3.1.5, 
                stages.levels.comb[train.indexes], 2)

##Atleast 2
##deseq2 1 fold
g.at2.1 <- intersect(net.features.updated$shrunken$atleast_2,
                     intersect(net.features.updated$varSelRF$atleast_2,
                        net.features.updated$deseq2$atleast_2$`1 fold`))
create.heatmap(vst_tumor_tum[train.indexes,], stages.levels.comb[train.indexes],
               g.at2.1, 
               'Heatmap using atleast 2 genes 1 fold', col = col, 
               cluster_rows = T, cluster_cols = T)
c.2.1 <- final.res(vst_tumor_tum, train.indexes, test.indexes, stages.levels.comb,
                     g.at2.1, 10)
get.clust.score(vst_tumor_tum[train.indexes,], 
                g.at2.1, 
                stages.levels.comb[train.indexes], 2)

g.at2.1.5 <- intersect(net.features.updated$shrunken$atleast_2,
                     intersect(net.features.updated$varSelRF$atleast_2,
                               net.features.updated$deseq2$atleast_2$`1.5 fold`))
create.heatmap(vst_tumor_tum[train.indexes,], stages.levels.comb[train.indexes],
               g.at2.1.5, 
               'Heatmap using atleast 2 genes 1.5 fold', col = col, 
               cluster_rows = T, cluster_cols = T)
c.2.1_5 <- final.res(vst_tumor_tum, train.indexes, test.indexes, stages.levels.comb,
                   g.at2.1.5, 10)
get.clust.score(vst_tumor_tum[train.indexes,], 
                g.at2.1.5, 
                stages.levels.comb[train.indexes], 2)


