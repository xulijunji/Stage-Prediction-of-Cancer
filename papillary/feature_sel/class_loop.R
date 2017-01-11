classifier.loop <- function(dattable,classifiers=c("svm","lda","rf","nsc"),feature.selection=c("aucs","InformationGain"),disc.method="MDL",threshold=0.3,
                  threshold.consis=0,attrs.nominal=numeric(),no.feat=20,flag.feature=TRUE,method.cross=c("leaveOneOut","sub-sampling","fold-crossval"))
{
  dattable[,ncol(dattable)]=as.factor(dattable[,ncol(dattable)])
  names(dattable)[ncol(dattable)] <- "class.label"

  if(flag.feature)
  {
    feature.subset=no.feat #validate all feature subsets
  }
  else
  {
    feature.subset=1   #validate one whole subset
  }
  times.selected <- matrix(0,ncol(dattable)-1,feature.subset)
  up=seq(feature.subset,1,-1)
  dimnames(times.selected)=c(list(colnames(dattable)[-ncol(dattable)]),list(as.character(up)))

  #for leave one out and cross-validation
  classi <- array(NA,c(nrow(dattable),feature.subset,length(classifiers)))
  misclassified <- array(NA,c(nrow(dattable),feature.subset,length(classifiers)))
  if(length(classifiers)==1)
  {
    dim3=list(classifiers)
  }
  else
  {
    dim3=classifiers
  }
  attr(classi, "dimnames")[[3]]=dim3
  attr(misclassified, "dimnames")[[3]]=dim3

  #for sub-sampling
  class.error<-array(0,c(nrow(dattable),feature.subset,length(classifiers)))
  attr(class.error, "dimnames")[[3]]=dim3

  if(method.cross=="fold-crossval")
  {
    num.group=10
    gr=NULL
    num.sel=floor(nrow(dattable)/num.group)
    num.add=nrow(dattable)%%num.group
    range=1:nrow(dattable)
    rr=nrow(dattable)
    for(i in 1:num.group)
    {
      vrem=sample(1:rr,size=num.sel)
      sel=range[vrem]
      gr=c(gr,list(sel))
      range=range[-vrem]
      rr=rr-num.sel
    }

    if(num.add>0)
    {
    vrem=sample(1:num.group,num.add)
    for(i in 1:num.add)
    {
      gr[[vrem[i]]]=c(gr[[vrem[i]]],range[i])
    }
    }

    #ptm <- proc.time()

    error=array(0,c(num.group,feature.subset,length(classifiers)))
    attr(error, "dimnames")[[3]]=dim3

    for (i in 1:num.group){

      dti.tr <- dattable[-gr[[i]],]
      class.tr=dattable[-gr[[i]],ncol(dattable)]
      dti.test<- dattable[gr[[i]],]
      class.test=dattable[gr[[i]],ncol(dattable)]

      selfea <- select.process(dti.tr,method=feature.selection,disc.method=disc.method,threshold=threshold,threshold.consis=threshold.consis,attrs.nominal=attrs.nominal,max.no.features=no.feat)
      if(!is.na(selfea[1]))
      {
      if(flag.feature)
      {
        start=1
      }
      else
      {
        start=length(selfea)
      }

      for(q in start:length(selfea))
      {
        times.selected[selfea[1:q],length(selfea)-q+1] <- times.selected[selfea[1:q],length(selfea)-q+1] + 1

        out=general.fun(classifiers,q,selfea,dti.tr,class.tr,dti.test,class.test)
        classi[gr[[i]],length(selfea)-q+1,classifiers] <- out$all.class.sample
        misclassified[gr[[i]],length(selfea)-q+1,classifiers] <- out$all.misclass.sample
        error[i,length(selfea)-q+1,classifiers]=(length(gr[[i]])-colSums(out$all.misclass.sample))/length(gr[[i]])
      }
      }
      cat(paste("Iteration",i))
    }
    #cat(proc.time() - ptm)
    #select on threshold for feature
    if(flag.feature)
    {
      true.classified=NULL
      for(i in 1:length(selfea))
      {
        true.classified=cbind(true.classified,1-error[,i,classifiers])
      }
    }
    else
    {
      true.classified=1-error[,1,classifiers]
      dim(true.classified)=c(num.group,length(classifiers))
      colnames(true.classified)=attr(error, "dimnames")[[3]]
    }

    classscores <- data.frame(dattable[,ncol(dattable)],classi[,1,],misclassified[,1,])


    misclnames <- rep(".correct",length(classifiers))
    for (i in 1:length(classifiers)){
      misclnames[i] <- paste(classifiers[i],misclnames[i],sep="")
    }
    names(classscores) <- c("True class",classifiers,misclnames)
    res <- list(predictions=classscores,no.selected=times.selected,true.classified=true.classified)


  }
  if(method.cross=="sub-sampling")
  {

    num.sample=100
    num.test=1/10

    times.select.inst <- rep(0,nrow(dattable))

    inddat <- 1:nrow(dattable)
    label=levels(dattable[,ncol(dattable)])

    confus=array(0,c(length(label),length(label),length(classifiers)))
    dimnames(confus)=c(list(label),list(label),list(classifiers))

    index.class=NULL
    size.sample=0
    for(i in 1:length(label))
    {
      index <-subset(inddat,dattable[,ncol(dattable)]==levels(dattable[,ncol(dattable)])[i])
      index.class=c(index.class,list(index))
      size.sample=size.sample+ceiling(length(index)*num.test)
    }

    error=array(0,c(num.sample,feature.subset,length(classifiers)))
    attr(error, "dimnames")[[3]]=dim3

    for (ib in 1:num.sample)
    {
      index.test=NULL
      for(i in 1:length(label))
      {
        index=sample(index.class[[i]],ceiling(num.test*length(index.class[[i]]))) #test cases
        index.test=c(index.test,index)
      }

      times.select.inst[index.test] <- times.select.inst[index.test] + 1

      dti.tr <- dattable[-index.test,]
      class.tr=dattable[-index.test,ncol(dattable)]
      dti.test<- dattable[index.test,]
      class.test=dattable[index.test,ncol(dattable)]

      #make the new function
      selfea <- select.process(dti.tr,method=feature.selection,disc.method=disc.method,threshold=threshold,attrs.nominal=attrs.nominal,max.no.features=no.feat)
      if(!is.na(selfea[1]))
      {
      if(flag.feature)
      {
        start=1
      }
      else
      {
        start=length(selfea)
      }

      for(q in start:length(selfea))
      {
       dtisel <- dti.tr[,c(selfea[1:q],ncol(dti.tr))]
       times.selected[selfea[1:q],length(selfea)-q+1] <- times.selected[selfea[1:q],length(selfea)-q+1] + 1

       out=general.fun(classifiers,q,selfea,dti.tr,class.tr,dti.test,class.test)
       class.error[index.test,length(selfea)-q+1,classifiers]=class.error[index.test,length(selfea)-q+1,classifiers]+out$all.misclass.sample
       error[ib,length(selfea)-q+1,classifiers]=length(index.test)-colSums(out$all.misclass.sample)
       if(q==length(selfea))
       {
       for(j in 1:length(classifiers))
       {
        vrem=table(class.test,out$all.class.sample[,j])
        confus[rownames(vrem),colnames(vrem),classifiers[j]]=confus[rownames(vrem),colnames(vrem),classifiers[j]]+vrem
       }
       }
      }
      }
      cat(paste("Iteration",ib))
    }
    #select on threshold for feature
    if(flag.feature)
    {
      true.classified=NULL
      for(i in 1:length(selfea))
      {
        true.classified=cbind(true.classified,1-error[,i,classifiers]/size.sample)
      }
#       if(length(classifiers)==1)
#       {
#         vrem=matrix(1:ncol(true.classified),ncol(true.classified),1)
#       }
#       else
#       {
#         vrem=sapply(classifiers,function(z) which(colnames(true.classified)==z))
#       }

#       windows()
#       par(mfrow=c(length(classifiers),1))
#       for(i in 1:length(classifiers))
#       {
#         values=true.classified[,vrem[,i]]
#         colnames(values)=seq(length(selfea),1,-1)
#         boxplot(values,main = paste("Cross validation",classifiers[i]),ylab = "Classification accuracy",xlab="n of features",col=i+1)
#       }

      classscores <- data.frame(dattable[,ncol(dattable)],times.select.inst,class.error[,1,])
    }
    else
    {
      # par(mfrow=c(1,1))
      true.classified=1-error[,1,classifiers]/size.sample

      dim(true.classified)=c(num.sample,length(classifiers))
      colnames(true.classified)=attr(error, "dimnames")[[3]]
      #boxplot(true.classified,main = "Cross validation",ylab = "Classification accuracy",xlab="Classifiers",col=3)
      classscores <- data.frame(dattable[,ncol(dattable)],times.select.inst,class.error[,1,])
    }


    time.correct <- rep(".time_correct",length(classifiers))
    for (i in 1:length(classifiers)){
      time.correct[i] <- paste(classifiers[i],time.correct[i],sep="")
    }
    names(classscores) <- c("true.label","time.selected",time.correct)
    res <- list(predictions=classscores,no.selected=times.selected,true.classified=true.classified,confus=confus)

    #to plot the ordered features according number of selections (for no.fea features)
    #ordered=order(times.selected[,1],decreasing=T)
    #data.frame(colnames(dattable)[ordered[1:no.feat]],times.selected[ordered[1:no.feat],1])
  }

  if(method.cross=="leaveOneOut")
  {
    #ptm <- proc.time()

  for (i in 1:nrow(dattable)){

    dti.tr <- dattable[-i,]
    class.tr=dattable[-i,ncol(dattable)]
    dti.test<- dattable[i,]
    class.test=dattable[i,ncol(dattable)]

    selfea <- select.process(dti.tr,method=feature.selection,disc.method=disc.method,threshold=threshold,attrs.nominal=attrs.nominal,max.no.features=no.feat)
    if(!is.na(selfea[1]))
    {
    if(flag.feature)
    {
      start=1
    }
    else
    {
      start=length(selfea)
    }

    for(q in start:length(selfea))
    {
      times.selected[selfea[1:q],length(selfea)-q+1] <- times.selected[selfea[1:q],length(selfea)-q+1] + 1

      out=general.fun(classifiers,q,selfea,dti.tr,class.tr,dti.test,class.test)
      classi[i,length(selfea)-q+1,classifiers] <- out$all.class.sample
      misclassified[i,length(selfea)-q+1,classifiers] <- out$all.misclass.sample
    }
    }
    cat(paste("Iteration",i))
  }
    #cat(proc.time() - ptm)
    #select on threshold for feature
    if(flag.feature)
    {
      true.classified=c()
      for(i in 1:length(selfea))
      {
        if(length(classifiers)==1)
        {
          true.classified=cbind(true.classified,colSums(misclassified[,i,,drop=FALSE])/nrow(dattable))
        }
        else
        {
          true.classified=cbind(true.classified,colSums(misclassified[,i,])/nrow(dattable))
        }
      }
#       dim(true.classified)=c(length(classifiers),length(selfea))
#       matplot(1:length(selfea),t(true.classified),xlab="n of features",ylab="Accuracy",type="b", col=1:length(classifiers), lty=1, pch=1:length(classifiers),
#               bty="n", las=1, main="Classification with n of features",xaxt="n")
#       axis(1,1:length(selfea),seq(length(selfea),1,-1))
#       legend("bottomright", col=1:length(classifiers), classifiers, bg="white", lwd=1, pch=1:length(classifiers))
    }
    else
    {
      if(length(classifiers)==1)
      {
        true.classified=colSums(misclassified[,1,,drop=FALSE])/nrow(dattable)
      }
      else
      {
        true.classified=colSums(misclassified[,1,])/nrow(dattable)
      }
#       matplot(1,true.classified,pch=1:2,ylim=c(0,1))
#       legend("bottomright", col=1:length(classifiers), paste(classifiers,format(true.classified,digits=2)),bg="white", lwd=1, pch=1:length(classifiers))
    }
    classscores <- data.frame(dattable[,ncol(dattable)],classi[,1,],misclassified[,1,])


    misclnames <- rep(".correct",length(classifiers))
    for (i in 1:length(classifiers)){
      misclnames[i] <- paste(classifiers[i],misclnames[i],sep="")
    }
    names(classscores) <- c("True class",classifiers,misclnames)
    res <- list(predictions=classscores,no.selected=times.selected,true.classified=true.classified)
  }

  return(res)
}