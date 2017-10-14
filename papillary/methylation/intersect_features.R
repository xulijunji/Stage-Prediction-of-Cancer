load('environment/methylation/net_fea_cpgss.RData')
load('environment/accuracy_feature/updated/net_features_trial.RData')
load('environment/genes_map.RData')
source('main/updated/initialisation.R')

library("AnnotationDbi")
library("org.Hs.eg.db")
genes.map <- select(org.Hs.eg.db,keys=colnames(vst_tumor_tum),column="SYMBOL", keytype="ENSEMBL")
genes.map <- genes.map[-which(duplicated(genes.map$ENSEMBL)),]
sum(colnames(vst_tumor_tum) == genes.map$ENSEMBL) == ncol(vst_tumor_tum)
rownames(genes.map)  <- genes.map$ENSEMBL
save(genes.map, file = 'environment/genes_map.RData')

####27K##########3
cpg.probe <- rownames(dmp.san.high)
genes <- net.features.trial$varSelRF$atleast_3

length(genes)
length(cpg.probe)
##Intersect genes between cpg and degs
genes.int <- intersect(probes.uniqe[cpg.probe], genes.map[genes,2])
length(genes.int)
cpgs.int <- names(unlist(lapply(genes.int, function(gene) which(probes.uniqe == gene))))
length(cpgs.int)                


###450 k
g <- Reduce(intersect, create.list.venn(fea.trial.list, 1, 4))
intersect(probes.450k.uniqe[rownames(dmp.450k)], genes.map[g,2])
high.beta.diff <- rownames(dmp.450k)[abs(dmp.450k[,5]) > 0.1]
length(high.beta.diff)
intersect(probes.450k.uniqe[high.beta.diff], genes.map[g,2])

ind.high <- match(intersect(probes.450k.uniqe[high.beta.diff], genes.map[g,2]),probes.450k.uniqe[high.beta.diff])
dmp.450k[high.beta.diff, 5][ind.high]
ens.g <- genes.map[match(intersect(probes.450k.uniqe[high.beta.diff], genes.map[g,2]), genes.map[,2]),1]

create.boxplots(ens.g[1], vst_tumor_tum, stages.levels.comb) ##feature_anal

rownames(rowData(mat.450k.req)) <- rowData(mat.450k.req)[,1]
View(rowData(mat.450k.req)[match(high.beta.diff[ind.high], rowData(mat.450k.req)[,1]),])
