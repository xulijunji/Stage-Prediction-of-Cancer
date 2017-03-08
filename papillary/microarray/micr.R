source("http://bioconductor.org/biocLite.R")
biocLite('hgu133plus2.db')
library('hgu133plus2.db')
biocLite("GEOquery")
library(GEOquery)
library(matrixStats)

gds <- getGEO(filename = '~/Downloads/data/papillary/microarray/GDS1344.soft')

gds2 <- getGEO(filename = '~/Downloads/data/papillary/microarray/GDS1344_full.soft')
names(Meta(gds))
Meta(gds)$type
colnames(Table(gds))
gpl570 <- getGEO(filename = '~/Downloads/data/papillary/microarray/gpl570.txt')

eset <- GDS2eSet(gds2, GPL = gpl570)
eset <- data.frame(exprs(eset))
eset$probe_id <- rownames(eset)
annot.df <- Table(gpl570)
entrez_df = as.data.frame(hgu133plus2ENTREZID)
entrez_df.tab <- table(entrez_df$probe_id)
sum(entrez_df.tab > 1)
ensembl_df = as.data.frame(hgu133plus2ENSEMBL2PROBE)
ensembl_df.tab = table(ensembl_df$probe_id)
ensembl_df.tab = ensembl_df.tab[ensembl_df.tab == 1]
ensembl_df <- ensembl_df[match(names(ensembl_df.tab), ensembl_df$probe_id),]
sum(table(ensembl_df$probe_id) > 1)

ind <- match(entrez_df$probe_id, rownames(eset))
eset <- eset[ind,]
merged <- merge(eset, entrez_df, by = 'probe_id')

merged$IQR <- rowIQRs(merged[,c(2:35)])
merged <- merged[order(merged$gene_id),]

final.indexes <- c()
i = 1
while(i <= length(merged$gene_id))
{
  ind.to.check <- c(i)
  j = i+1
  
  while(TRUE)
  {
    if(merged$gene_id[j] != merged$gene_id[j-1])
      break
  
    ind.to.check <- c(ind.to.check, j)
    j = j + 1
  }  
  
  ind.iqr <- ind.to.check[which.max(merged$IQR[ind.to.check])]
  final.indexes <- c(final.indexes, ind.iqr)
  i = j
}

merged_final <- merged[final.indexes,]
ens.id.match <- match(merged_final$probe_id, ensembl_df$probe_id)
merged_final$ens_id <- ensembl_df$ensembl_id[ens.id.match]
merged_final <- merged_final[order(merged_final$ens_id), ]

ensembl=useMart("ensembl")
ensembl=useDataset("hsapiens_gene_ensembl",mart=ensembl)
genes.entrez = getBM(attributes = c('ensembl_gene_id', 'entrezgene'), filters = 'entrezgene',
                     values = merged_final$gene_id, mart = ensembl)
genes.entrez <- genes.entrez[order(genes.entrez$ensembl_gene_id), ]

merged_final <- merged_final[-which(is.na(merged_final$ens_id)),]
sum(table(merged_final$probe_id) == 1) == length(merged_final$probe_id)
sum(table(merged_final$ens_id) == 1) == length(merged_final$probe_id)
sum(table(merged_final$gene_id) == 1) == length(merged_final$probe_id)
ens.dup <- which(table(merged_final$ens_id) > 1)

##genes to be removed (the entrez ids)
###ENSG00000011454 - 23637(kept) 2844(removed)
###ENSG00000086205 - 2346(kept)  219595(removed)
###ENSG00000104064 - 2553(kept)  55056(removed)
###ENSG00000111850 - 57150(kept) 63914(removed)
###ENSG00000127603 - 23499(kept) 643314(removed)
###ENSG00000136213 - 55501(kept) 101927181(removed)
###ENSG00000143226 - 2212(kept)  9103(removed)
###ENSG00000156273 - 571(kept)   100379661(removed)
###ENSG00000230124 - 84320(kept) 100527964(removed)
###ENSG00000233098 - 440416(kept) 339260(removed)
###ENSG00000267322 - 103091864(kept) 677769(removed)
ens.dup <- match(names(ens.dup), merged_final$ens_id)
merged_final$gene_id[ens.dup]

ens.dup.rem <- c(ens.dup[1]+1, ens.dup[2], ens.dup[3]+1, ens.dup[4]+1, ens.dup[5]+1, ens.dup[6],
                 ens.dup[7]+1, ens.dup[8], ens.dup[9], ens.dup[10], ens.dup[11]+1)
merged_final$gene_id[ens.dup.rem]
merged_final <- merged_final[-ens.dup.rem,]
save(merged_final,file = 'environment/accuracy_feature/updated/microarray_df.RData')

#######Note the below story is for what was done earlier
#######Now it has been taken care of ensembl_df where probes mapping to more than 1 ens_id removed
##3742  ENSG00000151079
which(ensembl_df$ensembl_id == 'ENSG00000130035')  ##756 and 26837
ensembl_df$probe_id[c(756,26837)]  ####"1553347_s_at" "220929_at"

which(ensembl_df$ensembl_id == 'ENSG00000151079')  ##only 1 756
ensembl_df$probe_id[756]      ###"1553347_s_at"
which(entrez_df$probe_id =="1553347_s_at") ##748
entrez_df$gene_id[748]  ###3742

###So probe id 1553347_s_at has 2 ens ids namely  ENSG00000130035 and ENSG00000151079
which(genes.entrez$ensembl_gene_id == 'ENSG00000130035') ##6004 6005
genes.entrez$entrezgene[c(6004,6005)]  ##3742 and 26290
which(genes.entrez$ensembl_gene_id == 'ENSG00000151079') ##0
which(genes.entrez$entrezgene == '3742') ##6004
##ENSG00000130035 has 2 probe ids

which(merged_final$ens_id == 'ENSG00000151079') ###not present
which(merged$gene_id == '3742') #19205

##for 'ENSG00000151079' I found out the probe from which I found the entrez 3742
##this 3742 could belong to either of the two ens ids of which we are not sure 
##Online says 3742 belongs to 'ENSG00000151079' whereas ensembl mart says 3742 belongs to 'ENSG00000130035'

##The problem I had that for a given ens id I had 2 entrez ids but now 1 of those entrez ids had some
##other ens id associated with it. Further that ens id was present in ens_df but not in merged_final.
##Thus I had fears of mapping wrong.But it turns out that the entrez id which had another ens id was associated
##with a probe with 2 ens ids. The probe though had only 1 entrez id making that entrez id also having
##2 ens ids. While using match in probe ids from entrez to ensembl I pick ENSG00000130035 and not
##ENSG00000151079 since it comes first. Since this has happened it here it means while using match other probe
##ids that could have multiple ens ids mappings pick only 1 ens id.

##Ideally probe with multiple associations should be removed so if I remove those then I am set anyway

a = ensembl_df$probe_id[match(merged_final$probe_id, ensembl_df$probe_id)]
a <- a[!is.na(a)]
gt_2 <- names(table(ensembl_df$probe_id))[which(table(ensembl_df$probe_id) > 1)]


####Sample analysis
merged_final <- merged_final[,c(2:35,37,1,36,38)]
sample_micro_info <- data.frame(accession = sort(colnames(merged_final)[1:34]))
sample_micro_info$class <- c(2,1,1,1,1,1,1,1,1,1,
                             1,1,1,1,1,2,1,2,2,1,
                             1,1,2,2,2,1,1,1,1,2,
                             2,2,2,2)
sample_micro_info$stage <- c(2,1,1,1,1,2,1,1,1,2,
                             1,1,1,2,1,2,2,1,2,1,
                             1,1,2,2,2,1,1,1,1,2,
                             2,2,2,2)
save(sample_micro_info, file = 'environment/accuracy_feature/updated/col_micro.RData')
