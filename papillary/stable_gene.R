stable_genes = list()
stable_genes_files = list.files('results/tumor/stable_gene_sel/')

j = 1

compute.class.means <- function(data, labs)
{
  ###data - N*j, where N number of samples and j number of genes
  classes <- levels(labs)
  class.means <- matrix(nrow = length(classes), ncol = dim(data)[2])
  
  for(i in seq(1:dim(data)[2]))
  {
    class.means[,i] = sapply(classes, function(x)
    {
      indexes = which(labs == x)
      mean.gene.class = mean(data[indexes,i])
    })
  }
  return(class.means)
}

compute.var.gene <- function(data, class.means, labs, gene.index)
{
  ###data - N*j, where N number of samples and j number of genes
  classes <- levels(labs)
  vars.class.genes = c()
  vars.classes.gene <- sapply(classes, function(x)
  {
      class.index <- which(classes == x)
      indexes <- which(labs == x)
    #  print(indexes)
      var.class.gene <- (sum((A[indexes, gene.index] - class.means[class.index, gene.index])^2))/(length(indexes)-1)
      #print(var.class.gene)
    })
  return(vars.classes.gene)
}

compute.ftest <- function(data, labs)
{
  ###data - N*j, where N number of samples and j number of genes
  classes <- levels(labs)
  class.means = compute.class.means(data, labs)
  class.lengths <- table(labs)
  means.overall <- apply(class.means, 2, mean)
  total.samples <- sum(class.lengths)
  total.classes <- length(classes)
  genes.vars.classes <- sapply(c(1:dim(data)[2]), function(i)
    {
    compute.var.gene(data, class.means, labs, i)
  })
  vars.genes <- sapply(c(1:dim(data)[2]), function(i)
    {
      gene.vars <- genes.vars.classes[,i]
      net.sigms <- mapply(function(x,y)
        {
        y*(y-1)*x
      }, gene.vars, class.lengths)
      sum(net.sigms)/(total.samples - total.classes)
  })
  
  f.vals <- sapply(c(1:dim(data)[2]), function(i)
  {
    f.val <- sapply(c(1:total.classes), function(x)
      {
        class.lengths[x]*(class.means[x,i] - means.overall[i])
    })
    sum(f.val)/((total.classes-1)*vars.genes[j])
  })
  return(f-vals)
}

A = matrix(nrow = 12, ncol = 4)
A[,1] = c(0.65,0.85,0.9,0.9,1.1,1.5,1.3,1.2,1.5,1.7,2.0,1.6)
A[,2] = c(0.2,1.0,1.2, 0.7,1.4,1.8,0.9, 1.0,1.7,2.1,2.0,1.2)
A[,3] = c(0.2,1.4,0.8, 0.6,2.0,1.2,1.0, 0.8,1.6,2.3,1.9,1.4)
A[,4] = c(0.7,1.0,1.3, 0.5,0.7,1.6,1.2, 1.5,0.2,1.8,0.3,1.2)

exm.lab <- factor(c('stage i', 'stage i', 'stage i', 'stage ii','stage ii','stage ii','stage ii', 
                    'stage iii','stage iii','stage iii','stage iii','stage iii'))

means <- compute.class.means(A, exm.lab)
ftest <- compute.ftest(A, exm.lab)

for(i in stable_genes_files)
{
  quo = as.integer(j / 3)
  rem = j %% 3
  if(quo == 0)
  {
     
  }
}