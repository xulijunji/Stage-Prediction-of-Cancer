source('decomp.R')
source('main/final.R')
load('environment/accuracy_feature/net_features.RData')
load('environment/accuracy_feature/dds_nor_tum_comb.RData')

sapply(net.features$shrunken$genes.list, length)
sapply(net.features$varSelRf$genes.list, length)
sapply(net.features$deseq$genes.list$`1.5fold`, length)

##Shrunken and VarSelRF
length(intersect(net.features$shrunken$atleast_1, net.features$varSelRf$atleast_1))
length(intersect(net.features$shrunken$atleast_3, net.features$varSelRf$atleast_3))
length(intersect(net.features$shrunken$atleast_5, net.features$varSelRf$atleast_5))


##Shrunken and DeSeq2
#1 fold
length(intersect(net.features$shrunken$atleast_1, net.features$deseq$atleast_1$`1fold`))
length(intersect(net.features$shrunken$atleast_3, net.features$deseq$atleast_3$`1fold`))
length(intersect(net.features$shrunken$atleast_5, net.features$deseq$atleast_5$`1fold`))

#1.5 fold
length(intersect(net.features$shrunken$atleast_1, net.features$deseq$atleast_1$`1.5fold`))
length(intersect(net.features$shrunken$atleast_3, net.features$deseq$atleast_3$`1.5fold`))
length(intersect(net.features$shrunken$atleast_5, net.features$deseq$atleast_5$`1.5fold`))

#2fold
length(intersect(net.features$shrunken$atleast_1, net.features$deseq$atleast_1$`2fold`))
length(intersect(net.features$shrunken$atleast_3, net.features$deseq$atleast_3$`2fold`))
length(intersect(net.features$shrunken$atleast_5, net.features$deseq$atleast_5$`2fold`))


