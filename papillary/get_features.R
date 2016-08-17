###Using Feature_Extract
##Env variables saved exp_prof, exp_fpqm
##Proj variables - genes.entrez, sample.file.map, ens.ids, ens.ids.all

##mRNA_counts
wd = '~/Downloads/data/papillary'
total_case_files = 321
cases = 289
mRNA_files_path = paste(wd,'mrna',sep = '/')
file_name = 'metadata.cart.mRna.json'
proj_dir = '~/Dropbox/honours/RNA_Seq/papillary/'
exp_prof = do.everything(wd, mRNA_files_path, total_case_files, file_name, proj_dir)
save(exp_prof, file = 'environment/exp_prof.RData')

remove(exp_prof,exp_fpqm)

##mRNA_fpkm
wd = '~/Downloads/data/papillary'
total_case_files = 321
cases = 289
mRNA_files_path = paste(wd,'fpqm-ql',sep = '/')
file_name = 'metadata.cart.fpqm.json'
proj_dir = '~/Dropbox/honours/RNA_Seq/papillary/'
exp_fpqm <- get.count.matrix(t.ens.ids.all, t.indexes, tfiles, t1, mRNA_files_path, proj_dir)
save(exp_prof, file = 'environment/exp_fpqm.RData')

remove(wd,total_case_files,cases,mRNA_files_path, file_name,proj_dir)
##Testing
#t1 = get.json.mapping(wd, file_name)
#tfiles = get.mRNA_files(mRNA_files_path, total_case_files)
#tfiles = unzip.files(tfiles, mRNA_files_path, proj_dir, total_case_files)
#t.ens.ids.all = get.genes.files(tfiles, mRNA_files_path, proj_dir)
#length(t.g)
#t.g[1]
#t.g = remove.dots(t.ens.ids.all)
#t.genes.entrez <- get.entrez(t.g)
#t.ens.ids <- get.ens.ids.with.entrez(t.genes.entrez)
#test(t.ens.ids, t.genes.entrez)
#t.indexes <- get.indexes(t.ens.ids, t.g)


#exp_fpqm <- exp_fpqm[rowSums(exp_fpqm)>1,]
#exp_fpqm1 <- do.everything(wd, mRNA_files_path, total_case_files, file_name, proj_dir)
#write.csv(exp_fpqm, 'fpqm.csv')
#for(i in seq_along(colnames(exp_fpqm)))
#{
#  print(sum(exp_fpqm[,i] == exp_fpqm1[,i]) == length(exp_fpqm[,1]))
#}