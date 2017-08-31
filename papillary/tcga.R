query <- GDCquery(project = 'TCGA-KIRP',
                  data.category = "Transcriptome Profiling",
                  data.type = "Gene Expression Quantification", 
                  workflow.type = "HTSeq - Counts")
mut <- GDCquery_Maf('KIRP', pipelines = 'mutect2')
cli <- GDCquery_clinic('TCGA-KIRP', 'biospecimen')
subtype <- TCGAquery_subtype(tumor = 'KIRP')
