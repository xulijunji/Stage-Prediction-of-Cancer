source('~/Honours/Stage-Prediction-of-Cancer/papillary/main/updated/feature_sel.R')
load('environment/methylation/meth_450k.data')
load('environment/sample_info_tumor_rep_normal.RData')
load('environment/stages.level.comb.RData')
nrow(data.450k)
ncol(data.450k)
library(dplyr)
library(minfi)

case.tum.rep.rna.ids <- filter(sample.info.all.rep, type == 'T') %>% select(sample.names) %>% unlist
case.tum.rep.rna.ids <- get.case.from.sample(case.tum.rep.rna.ids)
head(case.tum.rep.rna.ids)
length(case.tum.rep.rna.ids) == 260

meth.tum.inds.450k <-  which(colData(data.450k)[,4] == 'TP' | colData(data.450k)[,4] == 'TAP')
case.meth.tum.ids.450k <- get.case.from.sample(colnames(data.450k)[meth.tum.inds.450k])
length(case.meth.tum.ids.450k)
head(case.meth.tum.ids.450k)  

sum(case.tum.rep.rna.ids %in% case.meth.tum.ids.450k)

meth.rep.tum.inds.450k <- match(intersect(case.tum.rep.rna.ids, case.meth.tum.ids.450k), case.meth.tum.ids.450k)
length(meth.rep.tum.inds.450k)
meth.rep.tum.inds.450k <- meth.tum.inds.450k[meth.rep.tum.inds.450k]
table(colData(data.450k)[meth.rep.tum.inds.450k, 10])

table(sapply(strsplit(colData(data.450k)[meth.rep.tum.inds.450k,2], '-', fixed = T), function(x) x[4]))

meth.450k.tum.rep.data <- data.450k[,meth.rep.tum.inds.450k]

head(assay(meth.450k.tum.rep.data))
meth.450k.tum.rep.data <- meth.450k.tum.rep.data[-(which(is.na(rowSums(assay(meth.450k.tum.rep.data))))), ]
densityPlot(assay(meth.450k.tum.rep.data), colData(meth.450k.tum.rep.data)$gender)

chroms <- as.factor(sapply(strsplit(rowData(meth.450k.tum.rep.data)[,6], ':', fixed = T ), function(x) x[2]))
x.chrom.inds <- which(chroms == levels(chroms)[23])
y.chrom.inds <- which(chroms == levels(chroms)[24])
densityPlot(assay(meth.450k.tum.rep.data)[union(x.chrom.inds, y.chrom.inds),], colData(meth.450k.tum.rep.data)$gender)
meth.450k.tum.rep.data <- meth.450k.tum.rep.data[-union(x.chrom.inds, y.chrom.inds),]
densityPlot(assay(meth.450k.tum.rep.data), colData(meth.450k.tum.rep.data)$gender)
nrow(meth.450k.tum.rep.data)

probes.450k <- strsplit(rowData(meth.450k.tum.rep.data)[,2], ';', fixed = T)
probe.450k.uniq.ind <- which(sapply(sapply(probes.450k, unique), length) == 1)
length(probe.450k.uniq.ind)
probes.450k.uniqe <- sapply(probes.450k[probe.450k.uniq.ind], unique)
names(probes.450k.uniqe) = rowData(meth.450k.tum.rep.data)[,1]

meth.450k.tum.rep.data <- meth.450k.tum.rep.data[probe.450k.uniq.ind,]
densityPlot(assay(meth.450k.tum.rep.data), sampGroups = stage)
##Outlier
high.ind <- order(apply(assay(meth.450k.tum.rep.data), 2, function(col) sum(col < 0.2)), decreasing = T)[1]
densityPlot(assay(meth.450k.tum.rep.data)[, -high.ind])
meth.450k.tum.rep.data <- meth.450k.tum.rep.data[,-high.ind]

mat.450k.req  <- makeGenomicRatioSetFromMatrix(mat = assay(meth.450k.tum.rep.data), 
                                               rownames = rownames(meth.450k.tum.rep.data), 
                                          pData = colData(meth.450k.tum.rep.data), what = 'Beta' )
setdiff(rownames(assay(meth.450k.tum.rep.data)), rownames(mat.450k.req))
rowData(mat.450k.req) <- rowData(meth.450k.tum.rep.data)[match(rownames(mat.450k.req),
  rownames(meth.450k.tum.rep.data)),]

stage.ind.450k <- match(get.case.from.sample(colnames(mat.450k.req)),case.tum.rep.rna.ids)
length(stage.ind.450k) == 249
table(stages.levels.comb[stage.ind.450k])
table(colData(mat.450k.req)[,10])

nrow(assay(mat.450k.req))
dmp.450k <- dmpFinder(dat = assay(mat.450k.req), pheno = stages.levels.comb[stage.ind.450k],
                      type = 'categorical', qCutoff = 0.05 )
stage.meth.450 <- stages.levels.comb[stage.ind.450k]
stagei.ind <- which(stage.meth.450 == 'stage i')
stageiv.ind <- which(stage.meth.450 == 'stage iv')
beta.diff <- rowMeans(assay(mat.450k.req)[rownames(dmp.450k), stagei.ind]) -
    rowMeans(assay(mat.450k.req)[rownames(dmp.450k), stageiv.ind])
dmp.450k <- cbind(dmp.450k, diff = beta.diff)
nrow(dmp.450k)

library(doMC)
library(foreach)
registerDoMC(4)

sum(sort(stage.ind.450k) == stage.ind.450k)
ordered.stage.ind <- order(stage.ind.450k)
names(ordered.stage.ind) <- as.character(stage.ind.450k)

###tran
gr.450k <- list()
train.450k.ind <- intersect(train.trial.ind, stage.ind.450k)
test.450k.ind <- intersect(test.trial.ind, stage.ind.450k)

for(i in seq_along(gr.trial))
  gr.450k[[i]] <- intersect(gr.trial[[i]], stage.ind.450k)

gr.new.450 <- list()
for(i in seq_along(gr.trial))
  gr.new.450[[i]] <- ordered.stage.ind[as.character(gr.450k[[i]])]
get.stage.distribution(gr.new.450, stages.levels.comb[stage.ind.450k])
get.stage.distribution(gr.450k, stages.levels.comb)
gr.450k.train <- gr.450k[-5]

net.450.cpgs <- list()
net.450.cpgs[['shrunken']] <- get.feature.shrunken(t(assay(mat.450k.req)[rownames(dmp.450k1),]), stages.levels.comb[stage.ind.450k], gr.new.450[-5])
library(doParallel)
registerDoParallel(5)
vs <-get.feature.varSelRf(t(assay(mat.450k.req)), stages.levels.comb[stage.ind.450k], gr.new.450[-5])
  save(vs, file = 'environment/methylation/vs.RData')
sum(get.case.from.sample(colnames(mat.450k.req)) == get.case.from.sample(sample.info.all.rep$sample.names[tumor.ind.vs][stage.ind.450k]))
