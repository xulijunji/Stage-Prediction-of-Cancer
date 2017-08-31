library(pheatmap)
library(RColorBrewer)
load('environment/accuracy_feature/vs_nor_comb.RData')
load('environment/accuracy_feature/tumor_ind_vs.RData')
load('environment/stages.level.comb.RData')
load('environment/sample_info_tumor_rep_normal.RData')
load('environment/accuracy_feature/updated/net_features_trial.RData')
load('environment/accuracy_feature/updated/net_features_updated.RData')

source('main/updated/initialisation.R')

stage.com.norm <- c()  
for(stage in sample.info.all.rep$stage.type)
{
  if(stage == 'N')
    stage.com.norm <- c(stage.com.norm, 'N')
  else if(stage == 'stage i' | stage == 'stage ii')
    stage.com.norm <- c(stage.com.norm, 'early')
  else
    stage.com.norm <- c(stage.com.norm, 'late')
}
sample.info.all.rep$stage.comb.norm <- as.factor(stage.com.norm)
  
col <- colorRampPalette(rev(brewer.pal(9, 'RdYlBu')))(100)
vst_tumor_tum <- vs_normal_comb_reported[tumor.ind.vs, ]

genes.clus <- net.features.trial$deseq2$atleast_4$`2 fold`
data.clus <- vs_normal_comb_reported[ , genes.clus]
breaks = c(seq(5,7, length.out = 40), seq(7.1,24, length.out = 60))

clus_rows = run_hclust_on_a_matrix(data.clus)
clus_cols = run_hclust_on_a_matrix(t(data.clus))
pheatmap(t(data.clus),
         annotation_col = create.ordered.annotation(sample.info.all.rep$stage.comb.norm,
                                                    rownames(data.clus)),
#                 cluster_rows = clus_rows,
         #         cluster_cols = clus_cols,
         breaks = breaks,
         #          main = 'Heatmap using intersection for atleast 2 and 1 fold',
         show_rownames = F, show_colnames = F)

##On originial data
genes.orig <- net.features.updated$deseq2$atleast_4$`2 fold`
data.orig <- vst_tumor_tum[test.trial.ind,genes.orig]
pheatmap(t(data.orig),
         annotation_col = create.ordered.annotation(stages.levels.comb[test.trial.ind],
                                                    rownames(data.orig)),
         #                 cluster_rows = clus_rows,
         #         cluster_cols = clus_cols,
         breaks = breaks,
         #          main = 'Heatmap using intersection for atleast 2 and 1 fold',
         show_rownames = F, show_colnames = F)

