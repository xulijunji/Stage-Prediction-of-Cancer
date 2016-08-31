replace.symbol <- function(sym.list, replace.with, to.replace)
{
  sapply(sym.list, function(x)
    {
    temp <- unlist(strsplit(x, to.replace, fixed = T))
    paste(temp, collapse = '-')
  })
}

get.imp.genes <- function(method, assay, number)
{
  if(typeof(assay) == 'S4')
    assay = assay(assay)
  if(method == 1) ###most expression
    return(order(rowMeans(assay), decreasing = T)[1:number])
  else if(method == 2)
  {
    mads = apply(assay, 1, mad)
    return(order(mads, decreasing = T)[1:number])
  }
  else if(method == 3)
  {
    library(matrixStats)
    print(typeof(assay))
    rv <- rowVars(assay)
    return(order(rv, decreasing=TRUE)[1:number])
  }
    
}

create.list.imp.genes <- function(df.list, number)
{
  top.genes = list()
  for(i in seq_along(df.list))
    top.genes[[names(df.list)[i]]] = get.list.imp.genes(df.list[[i]], c(1,2,3), number)
  return(top.genes)
  
}

get.list.imp.genes <- function(df.list, methods, number)
{
  list.imp.genes <- list()
  for(i in methods)
    list.imp.genes[[i]] = get.imp.genes(i, df.list, number)
  names(list.imp.genes) = c('exp', 'mad', 'var')
  return(list.imp.genes)
}

##Tweeked original multiplot and plotPCA functions
multiplot <- function(list.plots, plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)
  
  # Make a list from the ... arguments and plotlist
  
  plots <- c(list.plots, plotlist)
  

  numPlots = length(plots)
  #print(numPlots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots==1) {
    print(plots[[1]])
    
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                     layout.pos.col = matchidx$col))
    }
  }
}

plotPCA = function(object, intgroup, ntop=500, returnData=FALSE, title, colData = NULL)
{
  library(ggplot2)
  library(matrixStats)
  
  if(typeof(object) == 'S4')
  {
    
    # perform a PCA on the data in assay(x) for the selected genes
    pca <- prcomp(t(assay(object)))
    
    # the contribution to the total variance for each component
    percentVar <- pca$sdev^2 / sum( pca$sdev^2 )
    
    if (!all(intgroup %in% names(colData(object)))) {
      stop("the argument 'intgroup' should specify columns of colData(dds)")
    }
    
    intgroup.df <- as.data.frame(colData(object)[, intgroup, drop=FALSE])
    
    # add the intgroup factors together to create a new grouping factor
    group <- if (length(intgroup) > 1) {
      factor(apply( intgroup.df, 1, paste, collapse=" : "))
    } else {
      colData(object)[[intgroup]]
    }
  }
  else
  {
    
    # perform a PCA on the data in assay(x) for the selected genes
    
    object = na.omit(object)
    pca <- prcomp(t(object))
    
    # the contribution to the total variance for each component
    percentVar <- pca$sdev^2 / sum( pca$sdev^2 )
    
    if (!all(intgroup %in% names(colData))) {
      stop("the argument 'intgroup' should specify columns of colData(dds)")
    }
    
    intgroup.df <- as.data.frame(colData[, intgroup, drop=FALSE])
    
    # add the intgroup factors together to create a new grouping factor
    group <- if (length(intgroup) > 1) {
      factor(apply( intgroup.df, 1, paste, collapse=" : "))
    } else {
      colData[[intgroup]]
    }
  }
  
  # assembly the data for the plot
  d <- data.frame(PC1=pca$x[,1], PC2=pca$x[,2], group=group, intgroup.df, name=colnames(object))
  
  if (returnData) {
    attr(d, "percentVar") <- percentVar[1:2]
    return(d)
  }
  
  ggplot(data=d, aes_string(x="PC1", y="PC2", color="group")) + geom_point(size=3) + 
    ggtitle(title)+
    xlab(paste0("PC1: ",round(percentVar[1] * 100),"% variance")) +
    ylab(paste0("PC2: ",round(percentVar[2] * 100),"% variance")) +
    coord_fixed() 
}

get.list.plotPCA <- function(assay.list, gene.list, colData = NULL, intgroup = NULL)
{
  ##arg - the indexes of list for differential genes and total number for others
  ##assay.list - the list of assays/dfs for which each individual plot shall be made
  ##gene.list - the vector containing the required genes for which we want to plot
  ##colData - col_data of DESeq object to be used for non assay type objects
  ##intrgroup - specifying the type of column in coldData to be used as basis of differentiation
  
  
  ##gets the list of plotPCA based on different categories for a list of dfs for a list of genes
  ##using get.list.plotPCA.df.gene
  list.df.gene.pca <- list()
  for(i in seq_along(gene.list[[1]]))
  {
    print(i)
    list.df.gene.pca[[names(gene.list[[1]])[i]]] = list()
    for(j in seq_along(assay.list))
    {
      list.df.gene.pca[[names(gene.list[[1]])[i]]][[names(assay.list)[j]]]  =
        get.list.plotPCA.df.gene(c(100,500,1000,2000,3000,4000,5000),assay.list[[j]], gene.list[[j]][[i]],
                                 names(gene.list[[1]])[i], 2, colData, intgroup)  
    }
    
  }
  return(list.df.gene.pca)
}

get.list.plotPCA.gene <- function(arg, assay.list, gene.vec, main_string, type = 1, colData = NULL, intgroup = NULL)
{
  ##arg - the indexes of list for differential genes and total number for others
  ##assay.list - the list of assays/dfs for which each individual plot shall be made
  ##gene.vec - the vector containing the required genes for which we want to plot
  ##main_string - Title to be shown on plot
  ##type - 1 if we are using differential genes else 2 for all other genes
  ##colData - col_data of DESeq object to be used for non assay type objects
  ##intrgroup - specifying the type of column in coldData to be used as basis of differentiation
  
  ##gets the list of plotPCA based on different categories for a list of dfs for a given type of genes
  ##using get.list.plotPCA.df.gene
  
  list.dfs.pca = lapply(assay.list, function(x)
  {
    get.list.plotPCA.df.gene(arg, x , gene.vec, main_string, type, colData, intgroup)
  })
  names(list.dfs.pca) = names(assay.list)
  return(list.dfs.pca)
}

get.list.plotPCA.df.gene <- function(arg, assay, gene.list, main_string, type = 1, colData = NULL, intgroup = NULL)
{
  ##arg - the indexes of list for differential genes and total number for others
  ##assay - the df for which plot shall be made
  ##gene.list - the list containing the required genes for which we want to plot
  ##main_string - Title to be shown on plot
  ##type - 1 if we are using differential genes else 2 for all other genes
  ##colData - col_data of DESeq object to be used for non assay type objects
  ##intrgroup - specifying the type of column in coldData to be used as basis of differentiation
  
  ##gets the list of plotPCA based on different categories for a given df for a given type of genes
    if(type == 1)
    {
      li <- lapply(arg, function(x)
        {
        #print(typeof(assay[gene.list[[x]],]))
        #print(x)
        a = plotPCA(assay[gene.list[[x]],], intgroup = intgroup, 
                title = paste(main_string,as.character(x), sep = '_'), colData = colData)
      
        })
     
    }
    else if(type == 2)
    {
      li <- lapply(arg, function(x)
      {
        plotPCA(assay[gene.list[1:x],], intgroup = intgroup, 
                title = paste(main_string,as.character(x), sep = '_'), colData = colData)})
    }
  names(li) <- arg 
  return(li)
}

save.plots.main.wd <- function(image.direct.main, image.direct.type, list.plots, cols)
{
  ##image.direct.main - main working directory
  ##image.direct.type - subdirectory for directory specific to the list working on
  ##list.plots - list of lists that will be fed to multiplot
  ##col - number of columns in multiplot
  
  ###saves the images of multiplots created on list of lists where the list of a list is fed to multiplot
  for(i in seq_along(list.plots))
  {
    file = paste(image.direct.main, image.direct.type, 
                 paste(names(list.plots)[i],'jpg', sep = '.'), sep = '/')
    jpeg(file, width = 1000, height = 1000)
    multiplot(list.plots[[i]], col = 2)
    dev.off()
    
  }
}

save.plots.under.main.wd <- function(image.direct.main, names.image.direct.type, list.plots, names.cols)
{
  ##image.direct.main - main working directory
  ##names.image.direct.type - subdirectory for directory specific to the list working on which will be 
  ##created first and then the image will be saved to
  ##list.plots - list for each criteria containing a list of lists that will be fed to multiplot
  ##names.col - specifying the number of columns for each criteria
  
  ###loops through a list feeding it save.plots.main.wd function creating its directory
  for(i in seq_along(list.plots))
  {
    dir.create(paste(image.direct.main, names.image.direct.type[i], sep = '/'))
    save.plots.main.wd(image.direct.main, names.image.direct.type[i], list.plots[[names.image.direct.type[i]]], names.cols[i])
  }
}

create.lda.plots <- function(lda.list, title)
{
  library(ggplot2)
  lda.ob = lda.list[[2]]
  lda.pred = lda.list[[1]]
  traces = (lda.ob$svd)^2/sum((lda.ob$svd)^2)
  ggplot(data = lda.pred$x, aes_string(x='LD1', y='LD2', color ='group')) +
  coord_fixed() + geom_point(size=2.5) + ggtitle(title) + 
  xlab(paste0(round(traces[1]*100), '%')) +
  ylab(paste0(round(traces[2]*100), '%'))
}
get.list.lda <- function(assay.list, gene.list, groups)
{
  ##arg - the indexes of list for differential genes and total number for others
  ##assay.list - the list of assays/dfs for which each individual plot shall be made
  ##gene.list - the vector containing the required genes for which we want to plot
  ##colData - col_data of DESeq object to be used for non assay type objects
  ##intrgroup - specifying the type of column in coldData to be used as basis of differentiation
  
  
  ##gets the list of plotPCA based on different categories for a list of dfs for a list of genes
  ##using get.list.plotPCA.df.gene
  list.df.gene.lda <- list()
  for(i in seq_along(gene.list[[1]]))
  {
    print(i)
    list.df.gene.lda[[names(gene.list[[1]])[i]]] = list()
    for(j in seq_along(assay.list))
    {
      list.df.gene.lda[[names(gene.list[[1]])[i]]][[names(assay.list)[j]]]  =
        get.list.lda.df.gene(c(100,500,1000,2000,3000,4000,5000),assay.list[[j]], gene.list[[j]][[i]],
                                 names(gene.list[[1]])[i], 2, groups = groups)  
    }
    
  }
  return(list.df.gene.lda)
}

get.list.lda.df.gene <- function(arg, assay, gene.list, main_string, type = 1, groups)
{
  ##arg - the indexes of list for differential genes and total number for others
  ##assay - the df for which plot shall be made
  ##gene.list - the list containing the required genes for which we want to plot
  ##main_string - Title to be shown on plot
  ##type - 1 if we are using differential genes else 2 for all other genes
  ##colData - col_data of DESeq object to be used for non assay type objects
  ##intrgroup - specifying the type of column in coldData to be used as basis of differentiation
  
  ##gets the list of plotPCA based on different categories for a given df for a given type of genes
  library(MASS)
  library(DESeq2)
  if(typeof(assay) == 'S4')
    assay = assay(assay)
  #print(typeof(gene.list))
  if(type == 1)
  {
    li <- lapply(arg, function(x)
    {
      req.data <- t(assay[gene.list[[x]],])
      #tmp <- cor(req.data)
      #req.dat <- req.data[,!apply(tmp,2,function(x) any(abs(x) > 0.99))]
      #print(head(req.data))
      #print(head(req.dat))
      lda1 = lda(req.data, grouping = groups)
      lda1.p = predict(lda1)
      lda1.p$x = data.frame(lda1.p$x)
      lda1.p$x$group = groups
      lda1.list = list(lda1.p, lda1)
      create.lda.plots(lda1.list, 
                       title = paste(main_string,as.character(x), sep = '_'))
    })
    
  }
  else if(type == 2)
  {
    li <- lapply(arg, function(x)
    {
      req.data <- t(assay[gene.list[1:x],])
      #tmp <- cor(req.data)
      #req.data <- req.data[,!apply(tmp,2,function(x) any(abs(x) > 0.99))]
      lda1 = lda(req.data, grouping = groups)
      lda1.p = predict(lda1)
      lda1.p$x = data.frame(lda1.p$x)
      lda1.p$x$group = groups
      lda1.list = list(lda1.p, lda1)
      create.lda.plots(lda1.list, 
                       title = paste(main_string,as.character(x), sep = '_'))
    })
  }
  
  names(li) <- arg 
  return(li)
}

get.list.lda.gene <- function(arg, assay.list, gene.vec, main_string, type = 1, group)
{
  ##arg - the indexes of list for differential genes and total number for others
  ##assay.list - the list of assays/dfs for which each individual plot shall be made
  ##gene.vec - the vector containing the required genes for which we want to plot
  ##main_string - Title to be shown on plot
  ##type - 1 if we are using differential genes else 2 for all other genes
  ##colData - col_data of DESeq object to be used for non assay type objects
  ##intrgroup - specifying the type of column in coldData to be used as basis of differentiation
  
  ##gets the list of plotPCA based on different categories for a list of dfs for a given type of genes
  ##using get.list.plotPCA.df.gene
  
  list.dfs.lda = lapply(assay.list, function(x)
  {
    get.list.lda.df.gene(arg, x , gene.vec, main_string, type, group)
  })
  names(list.dfs.lda) = names(assay.list)
  return(list.dfs.lda)
}
