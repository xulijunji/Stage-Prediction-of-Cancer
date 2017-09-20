load('environment/methylation/net_fea_cpgss.RData')
load('environment/accuracy_feature/updated/net_features_trial.RData')
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
