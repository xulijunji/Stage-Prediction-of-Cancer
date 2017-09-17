load('environment/methylation/whole_data.RData')
load('environment/methylation/meth_tum_rep_inds.RData')
load('environment/stages.level.comb.RData')
library(SummarizedExperiment)
library(minfi)
ncol(data)

table(colData(data)[meth.rep.tum.inds, 10])
##Only tumor samples with known stage
meth.tum.rep.data <- data[,meth.rep.tum.inds] ##from get_data.R
remove(data)

##Removing the ones with na
meth.tum.rep.data <- meth.tum.rep.data[-(which(is.na(rowSums(assay(meth.tum.rep.data))))), ]
densityPlot(assay(meth.tum.rep.data), colData(meth.tum.rep.data)$gender)

chroms <- as.factor(sapply(strsplit(rowData(meth.tum.rep.data)[,6], ':', fixed = T ), function(x) x[2]))
x.chrom.inds <- which(chroms == levels(chroms)[23])
y.chrom.inds <- which(chroms == levels(chroms)[24])
densityPlot(assay(meth.tum.rep.data)[union(x.chrom.inds, y.chrom.inds),], colData(meth.tum.rep.data)$gender)
meth.tum.rep.data <- meth.tum.rep.data[-union(x.chrom.inds, y.chrom.inds),]
densityPlot(assay(meth.tum.rep.data), colData(meth.tum.rep.data)$gender)

hist(rowSds(assay(meth.tum.rep.data)))

probes <- strsplit(rowData(meth.tum.rep.data)[,2], ';', fixed = T)
probe.uniq.ind <- which(sapply(sapply(probes, unique), length) == 1)

meth.tum.rep.data <- meth.tum.rep.data[probe.uniq.ind,]
densityPlot(assay(meth.tum.rep.data), sampGroups = stages.levels.comb)

mat.req  <- makeGenomicRatioSetFromMatrix(mat = assay(meth.tum.rep.data), rownames = rownames(meth.tum.rep.data), 
                                          pData = colData(meth.tum.rep.data), what = 'Beta' )

rowData(mat.req) <- rowData(meth.tum.rep.data)
qc <- getQC(mat.req)
  
densityPlot(assay(data)[,-meth.27.ind]) ## A total of 5601
dmp <- dmpFinder(assay(meth.tum.rep.data), pheno = stages.levels.comb, 
                 'categorical', 0.05)
nrow(dmp)

##Looking at outliers
high.ind <- order(apply(assay(meth.tum.rep.data), 2, function(col) sum(col < 0.2)), decreasing = T)[1:2]
hist(apply(assay(meth.tum.rep.data), 2, function(col) sum(col < 0.2)))
colnames(meth.tum.rep.data)[high.ind]
densityPlot(assay(meth.tum.rep.data[,-high.ind]))
dmp.san.high <-dmpFinder(assay(meth.tum.rep.data[,-high.ind]), pheno = stages.levels.comb[-high.ind], 
                     'categorical', 0.05)
nrow(dmp.san.high)  #6103

length(intersect(rownames(dmp), rownames(dmp.san.high)))  ##5585


##Looking at low standard deviation
hist(rowSds(assay(meth.tum.rep.data)))
length(intersect(rownames(meth.tum.rep.data)[rowSds(assay(meth.tum.rep.data)) < 0.05], 
                 rownames(dmp)))


##DMRs
bmp <- bumphunter(object=mat.req, design = model.matrix(~stages.levels.comb),
                  cutoff = 0.1, B = 1000 )
cl = clusterMaker(seqnames(mat.req), pos = start(mat.req), maxGap = 500)
bf <- blockFinder(object = mat.req, design = model.matrix(~stages.levels.comb), 
                  what = 'Beta', cluster = cl, cutoff = 0.1)
library(doParallel)
registerDoParallel(cores = 8)
View(bmp$table)


###
source('main/updated/initialisation.R')
shrunk.meth <- get.shrunken.features(list(c(1:260)), data = t(assay(mat.req)), stages.levels.comb, 2)
shrunk.meth.diff <- get.shrunken.features(list(c(1:260)), data = t(assay(mat.req)[rownames(dmp),]), 
                                          stages.levels.comb, 2)
cpg.shrunken <- get.genes.shrunken(shrunk.meth$genes)
cpg.shrunken.diff <- get.genes.shrunken(shrunk.meth.diff$genes)
length(intersect(cpg.shrunken[[1]], rownames(dmp)))
