###heatmaps and visualisation based on various features
genes.comb <- list()
gr <- build.groups(length(stages.levels.comb), 5)
load('environment/accuracy_feature/shrunken_list_vs_normal_reported.RData')
load('environment/accuracy_feature/dds_nor_tum_comb.RData')
load('environment/accuracy_feature/vs_nor_comb.RData')
load('environment/accuracy_feature/tumor_ind_vs.RData')
load('environment/accuracy_feature/net_features.RData')
load('environment/accuracy_feature/classifer_list.RData')
source('decomp.R')
source('main/final.R')

library(DESeq2)
library(RColorBrewer)
library(VennDiagram)
display.brewer.all()
col <- brewer.pal(9, 'PuRd')

create.heatmap(vs_normal_comb_reported[tumor.ind.vs, ], 
               stages.levels.comb,
               col = col,
genes = net.features$shrunken$atleast_dfs$atleast_5$genes[
 net.features$shrunken$atleast_dfs$atleast_5$stage=='1'  
]
                )

create.heatmap(vs_normal_comb_reported[tumor.ind.vs, ], 
               stages.levels.comb,
               col = col,
genes = g, title = 'Heatmap of the 21 genes selected by all feature
selection methods'
)

g = intersect(net.features$deseq$atleast_5$`1.5fold`,
              net.features$shrunken$atleast_5)

venn.diagram(list(net.features$shrunken$atleast_5, net.features$varSelRf$atleast_5,
                  net.features$deseq$atleast_5$`1.5fold`), filename = 'images/tumor/features/int.jpg',
             category.names = c('Shrunken', 'VarSelRF', 'DESeq2'), main = 'Venn Diagram of genes present in 
             all 5 groups using different feature selection algorithms',
             fill = c("skyblue", "pink1", "mediumorchid"), main.col = 'blue4', cat.col = 'blue4',
            main.cex = 1.5, cat.cex = 1.2, cex = 1.2         )
      

get.mean(classifier.list$deseq$atleast_5[[2]]$rf$test$other,
         2, 'byClass')


