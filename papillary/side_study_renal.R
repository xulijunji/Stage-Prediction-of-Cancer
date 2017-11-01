load('environment/genes_map.RData')
load('environment/diff_genes.RData')
library(dplyr)
genes.rcc <- read.delim('~/genes_rcc.txt', header = F)
genes.rcc <- as.character(genes.rcc$V1)
which(is.na(match(genes.rcc, genes.map$SYMBOL)))
genes.rcc[19]  <- 'PRR35'
genes.rcc[40] <- 'RBBP8NL'
genes.rcc[43] <- 'PRDM16'
genes.rcc[60] <- 'BRINP3'
genes.rcc[102] <- 'MFSD4A'
genes.rcc[104] <- 'TMEM52B'
genes.rcc[110] <- 'PLEKHD1'
genes.rcc[136] <- 'CLC'
genes.rcc[172] <- 'LINC00922'
genes.rcc[175] <- 'PLPPR5'
genes.rcc[183] <- 'FDCSP'

tolower('FDCSP') %in% tolower(genes.map$SYMBOL)
genes.rcc <- intersect(genes.rcc, genes.map$SYMBOL)
ens.genes.rcc <- genes.map$ENSEMBL[match(genes.rcc, genes.map$SYMBOL)]
df.rcc <- data.frame(genes.names = genes.rcc, ens.ids = ens.genes.rcc,
                     type = as.character(rep(c('up'), length(ens.genes.rcc))),
                     stringsAsFactors = F, row.names = 1)
which(df.rcc$genes.names == 'GRIK3')
df.rcc$type[124:178] <- 'low'
rcc.genes.up <- filter(df.rcc, type == 'up') %>% select(ens.ids) %>% unlist
rcc.genes.low <- filter(df.rcc, type == 'low') %>% select(ens.ids) %>% unlist

length(intersect(as.character(ens.genes.rcc), diff.genes$`1`))
intersect(rcc.genes.low, diff.genes.up$`2`)

tum.nor.low.multicell <- read.delim('~/t2.txt', header = F, stringsAsFactors = F)
tum.union.up.multicell <- read.delim('~/t.txt', header = F, stringsAsFactors = F)
intersect(tum.nor.low.multicell$V1, tum.union.up.multicell$V1)

tum.union.up.dev <- read.delim('~/t3.txt', header = F, stringsAsFactors = F)
intersect(tum.nor.low.multicell$V1, tum.union.up.dev$V1)
intersect(tum.union.up.dev$V1, tum.union.up.multicell$V1)

