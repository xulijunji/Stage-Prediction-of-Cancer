stable_genes = list()
stable_genes_files = list.files('results/tumor/stable_gene_sel/')
stable_genes_files = stable_genes_files[2:17]
names.stable.genes <- c('1','2','cho','f')
names.dfs <- c('fpqm','fpqm_log','nt','vs')

fvals.df <- lapply(req.dfs, function(x)
      {
        compute.ftest(x, stages.levels)
      }
  )
names(fvals.df) = names(req.dfs)
write.dfs.csvs('results/tumor/stable_gene_sel/',fvals.df)

####Testing dummy vars
A = matrix(nrow = 12, ncol = 4)
A[,1] = c(0.65,0.85,0.9,0.9,1.1,1.5,1.3,1.2,1.5,1.7,2.0,1.6)
A[,2] = c(0.2,1.0,1.2, 0.7,1.4,1.8,0.9, 1.0,1.7,2.1,2.0,1.2)
A[,3] = c(0.2,1.4,0.8, 0.6,2.0,1.2,1.0, 0.8,1.6,2.3,1.9,1.4)
A[,4] = c(0.7,1.0,1.3, 0.5,0.7,1.6,1.2, 1.5,0.2,1.8,0.3,1.2)

exm.lab <- factor(c('stage i', 'stage i', 'stage i', 'stage ii','stage ii','stage ii','stage ii', 
                    'stage iii','stage iii','stage iii','stage iii','stage iii'))

means <- compute.class.means(A, exm.lab)
fvals <- compute.ftest(A, exm.lab)
df <- compute.ftest(A, exm.lab)
remove(ftest,means,vars,A,exm.lab)
remove(req.dfs) ##write_input.R
####Testing ends

j = 1
for(i in stable_genes_files)
{
  quo = ceiling(j / 4)
  rem = j %% 4
  if(rem == 0)
  {
    rem = 4
    f = read.csv(paste('results/tumor/stable_gene_sel/',i,sep = ''), header = T)
  }
  else
    f = read.csv(paste('results/tumor/stable_gene_sel/',i,sep = ''), header = F)
  
  stable_genes[[names.dfs[quo]]][[names.stable.genes[rem]]] = f
  j = j + 1
}
remove(f,i,j,names.dfs,names.stable.genes,stable_genes_files)

length(intersect(stable_genes$fpqm$cho$V1[1:100], stable_genes$fpqm$f$indexes[1:100]))

stable_genes_files_proc = list.files('results/tumor/stable_gene_sel/proc/')
stable_genes_proc = list()
j = 1
names.dfs.proc = c('fpqm','vs')
names.stable.genes.proc = c('1','2','cho')
for(i in stable_genes_files_proc)
{
  quo = ceiling(j / 3)
  rem = j %% 3
  if(rem == 0)
    rem = 3
  f = read.csv(paste('results/tumor/stable_gene_sel/proc/',i,sep = ''), header = F)
  
  stable_genes_proc[[names.dfs.proc[quo]]][[names.stable.genes.proc[rem]]] = f
  print(names.dfs.proc[quo])
  print(names.stable.genes.proc[rem])
  j = j + 1
}
remove(f,i,j,names.dfs.proc,names.stable.genes.proc,stable_genes_files_proc)
