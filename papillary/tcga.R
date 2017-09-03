library(TCGAbiolinks)
load('environment/sample_info_tumor_rep_normal.RData')
query <- GDCquery(project = 'TCGA-KIRP',
                  data.category = "Transcriptome Profiling",
                  data.type = "Gene Expression Quantification", 
                  workflow.type = "HTSeq - Counts")
mut <- GDCquery_Maf('KIRP', pipelines = 'mutect2')
cli <- GDCquery_clinic('TCGA-KIRP', 'clinical')
subtype <- TCGAquery_subtype(tumor = 'KIRP')

sample.info.all.rep.tum <- sample.info.all.rep[sample.info.all.rep$type == 'T',]
tcga.barcodes <- sapply(strsplit(sample.info.all.rep.tum$sample.names, '-', fixed = T),
       function(x) paste(x[[1]], x[[2]], x[[3]], sep = '-'))

sum(is.na(match(tcga.barcodes, cli$submitter_id)))

common <- intersect(tcga.barcodes, cli.check$submitter_id)
cli.check <- cli[,c(1,6)]
rownames(cli.check) <- cli.check$submitter_id
cli.check <- cli.check[common,]
cli.check$mine <- sample.info.all.rep.tum$stage.type[match(cli.check$submitter_id, tcga.barcodes)]
