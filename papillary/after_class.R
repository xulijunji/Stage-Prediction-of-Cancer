####Deals with after classification with the output that comes###########
compute.error.conf.mat <- function(conf.mat)
{
  total.samps <- apply(conf.mat, 1, sum)
  error <- 1-diag(conf.mat)/total.samps
  return(error)
}
create.mat.error <- function(conf.mat)
{
  return(cbind(conf.mat, compute.error.conf.mat(conf.mat)))
}

plot.auc <- function(cv, stages.levels)
{
  thresh <- cv$threshold
  aucs <- sapply(seq_along(thresh), function(x)
  {
    multiclass.roc(stages.levels, ordered(cv$yhat[,x]))$auc
  })
  #print(aucs)
  x <- seq(thresh)
  plot(x, aucs)
}
