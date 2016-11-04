###using genes from stable_gene.R
###req.dfs from write_input.R

library(e1071)
library(caret)
library(class)
load('environment/stable_genes.RData')
load('environment/req_dfs.RData')

##KNN
##Based on 1
pred.knn.fpqm.1 <- knn.cv(req.dfs$fpqm[,stable_genes$fpqm$`1`$V1[1:150]], cl = stages.levels, k = 5)
create.mat.error(table(stages.levels, pred.knn.fpqm.1))

pred.knn.fpqm_log.1 <- knn.cv(req.dfs$fpqm[,stable_genes$fpqm_log$`1`$V1[1:150]], cl = stages.levels, k = 4)
table(stages.levels, pred.knn.fpqm_log.1)

pred.knn.nt.1 <- knn.cv(req.dfs$nt[,stable_genes$nt$`1`$V1[1:150]], cl = stages.levels, k = 5)
table(stages.levels, pred.knn.nt.1)

pred.knn.vs.1 <- knn.cv(req.dfs$vs[,stable_genes$vs$`1`$V1[1:150]], cl = stages.levels, k = 5)
table(stages.levels, pred.knn.vs.1)

#Based on 2
pred.knn.fpqm.2 <- knn.cv(req.dfs$fpqm[,stable_genes$fpqm$`2`$V1[1:150]], cl = stages.levels, k = 3)
table(stages.levels, pred.knn.fpqm.2)

pred.knn.fpqm_log.2 <- knn.cv(req.dfs$fpqm_log[,stable_genes$fpqm_log$`2`$V1[1:150]], cl = stages.levels, k = 4)
table(stages.levels, pred.knn.fpqm_log.2)

pred.knn.nt.2 <- knn.cv(req.dfs$nt[,stable_genes$nt$`2`$V1[1:150]], cl = stages.levels, k = 5)
table(stages.levels, pred.knn.nt.2)

pred.knn.vs.2 <- knn.cv(req.dfs$vs[,stable_genes$vs$`2`$V1[1:150]], cl = stages.levels, k = 5)
table(stages.levels, pred.knn.vs.2)

###Based on cho
pred.knn.fpqm.cho <- knn.cv(req.dfs$fpqm[,stable_genes$fpqm$cho$V1[1:150]], cl = stages.levels, k = 4)
table(stages.levels, pred.knn.fpqm.cho)

pred.knn.fpqm_log.cho <- knn.cv(req.dfs$fpqm_log[,stable_genes$fpqm_log$cho$V1[1:150]], cl = stages.levels, k = 4)
table(stages.levels, pred.knn.fpqm_log.cho)

pred.knn.nt.2 <- knn.cv(req.dfs$nt[,stable_genes$nt$`2`$V1[1:150]], cl = stages.levels, k = 5)
table(stages.levels, pred.knn.nt.2)

pred.knn.vs.2 <- knn.cv(req.dfs$vs[,stable_genes$vs$`2`$V1[1:150]], cl = stages.levels, k = 5)
table(stages.levels, pred.knn.vs.2)

##Based on F
pred.knn.fpqm.f <- knn.cv(req.dfs$fpqm[,stable_genes$fpqm$f$V1[1:150]], cl = stages.levels, k = 4)
table(stages.levels, pred.knn.fpqm.f)

pred.knn.fpqm_log.f <- knn.cv(req.dfs$fpqm_log[,stable_genes$fpqm_log$f$V1[1:150]], cl = stages.levels, k = 4)
table(stages.levels, pred.knn.fpqm_log.f)

pred.knn.nt.2 <- knn.cv(req.dfs$nt[,stable_genes$nt$`2`$V1[1:150]], cl = stages.levels, k = 5)
table(stages.levels, pred.knn.nt.2)

pred.knn.vs.2 <- knn.cv(req.dfs$vs[,stable_genes$vs$`2`$V1[1:150]], cl = stages.levels, k = 5)
table(stages.levels, pred.knn.vs.2)

##On Processed data
pred.knn.fpqm.proc.1 <- knn.cv(req.dfs.rem$fpqm[,stable_genes_proc$fpqm$`1`$V1[1:350]], cl = stages.levels, k = 5)
compute.error.conf.mat(table(stages.levels, pred.knn.vs.2))
####SVMs
##Based on 1
cv.svm(req.dfs$fpqm[,stable_genes$fpqm$`1`$V1[1:150]], 5, stages.levels)
cv.svm(req.dfs$vs[,stable_genes$vs$`1`$V1[1:150]], 5, stages.levels)
l.o.1 = list()
#200
l.o.1[['fpqm']] = lapply(c(50,100,150,200,250),function(x){
                         cv.svm.leave.one.out(req.dfs$fpqm[,stable_genes$fpqm_log$`1`$V1[1:x]], stages.levels)})
#200
l.o.1[['fpqm_log']] = lapply(c(50,100,150,200,250),function(x){
      cv.svm.leave.one.out(req.dfs$fpqm_log[,stable_genes$fpqm$`1`$V1[1:x]], stages.levels)})
#450
l.o.1[['nt']] = lapply(c(50,100,150,200,250),function(x){
      cv.svm.leave.one.out(req.dfs$nt[,stable_genes$nt$`1`$V1[1:x]], stages.levels)})
t1 = lapply(c(300,350,400,450),function(x){
  cv.svm.leave.one.out(req.dfs$nt[,stable_genes$nt$`1`$V1[1:x]], stages.levels)})
#400
l.o.1[['vs']] = lapply(c(50,100,150,200,250),function(x){cv.svm.leave.one.out(req.dfs$vs[,stable_genes$vs$`1`$V1[1:x]], stages.levels)})
t = lapply(c(300,350,400,450),function(x){cv.svm.leave.one.out(req.dfs$vs[,stable_genes$vs$`1`$V1[1:x]], stages.levels)})

cv.svm(req.dfs$vs[,stable_genes$vs$`1`$V1[1:150]], 10, stages.levels)
cv.svm(req.dfs$nt[,stable_genes$nt$`1`$V1[1:150]], 10, stages.levels)

lo.2 = list()
#150
lo.2[['fpqm']] = lapply(c(50,100,150,200,250,300,350,400,450),function(x){cv.svm.leave.one.out(req.dfs$fpqm[,stable_genes$fpqm$`2`$V1[1:x]],  stages.levels)})
#300
lo.2[['fpqm_log']] = lapply(c(50,100,150,200,250,300,350,400,450),function(x){cv.svm.leave.one.out(req.dfs$fpqm_log[,stable_genes$fpqm_log$`2`$V1[1:x]],  stages.levels)})
#300
lo.2[['nt']] = lapply(c(50,100,150,200,250,300,350,400,450),function(x){cv.svm.leave.one.out(req.dfs$nt[,stable_genes$nt$`2`$V1[1:x]],  stages.levels)})
#400
lo.2[['vs']] = lapply(c(50,100,150,200,250,300,350,400,450),function(x){cv.svm.leave.one.out(req.dfs$vs[,stable_genes$vs$`2`$V1[1:x]],  stages.levels)})

cv.svm(req.dfs$fpqm[,stable_genes$fpqm$`2`$V1[1:200]], 5, stages.levels)
cv.svm(req.dfs$fpqm_log[,stable_genes$fpqm_log$`2`$V1[1:200]], 10, stages.levels)
cv.svm(req.dfs$vs[,stable_genes$vs$`2`$V1[1:200]], 10, stages.levels)
cv.svm(req.dfs$nt[,stable_genes$nt$`2`$V1[1:200]], 10, stages.levels)

lo.cho = list()
lo.cho[['fpqm']] = cv.svm.leave.one.out(req.dfs$fpqm[,stable_genes$fpqm$cho$V1[1:200]],  stages.levels)
lo.cho[['fpqm_log']] = cv.svm.leave.one.out(req.dfs$fpqm_log[,stable_genes$fpqm_log$cho$V1[1:200]],  stages.levels)
lo.cho[['nt']] = cv.svm.leave.one.out(req.dfs$nt[,stable_genes$nt$cho$V1[1:200]],  stages.levels)
lo.cho[['vs']] = cv.svm.leave.one.out(req.dfs$vs[,stable_genes$vs$cho$V1[1:200]],  stages.levels)

cv.svm(req.dfs$fpqm[,stable_genes$fpqm$cho$V1[1:200]], 5, stages.levels)
cv.svm(req.dfs$fpqm_log[,stable_genes$fpqm_log$cho$V1[1:200]], 5, stages.levels)
cv.svm(req.dfs$vs[,stable_genes$vs$`1`$cho[1:200]], 5, stages.levels)
cv.svm(req.dfs$nt[,stable_genes$nt$`1`$cho[1:200]], 5, stages.levels)

cv.svm(req.dfs$fpqm[,stable_genes$fpqm$f$indexes[1:150]], 5, stages.levels)
cv.svm(req.dfs$fpqm_log[,stable_genes$fpqm_log$f$indexes[1:150]], 10, stages.levels)
cv.svm(req.dfs$vs[,stable_genes$vs$f$indexes[1:150]], 10, stages.levels)
cv.svm(req.dfs$nt[,stable_genes$nt$f$indexes[1:150]], 10, stages.levels)

lo.proc = list()
lo.proc[['1']] = list()
lo.proc[['1']][['fpqm']] = lapply(c(50,100,150,200,250,300,350,400),function(x){
  cv.svm.leave.one.out(req.dfs.rem$fpqm[,stable_genes_proc$fpqm$`1`$V1[1:x]], stages.levels)})
#400
lo.proc[['1']][['vs']] = lapply(c(50,100,150,200,250,300,350,400),function(x){
  cv.svm.leave.one.out(req.dfs.rem$vs[,stable_genes_proc$vs$`1`$V1[1:x]], stages.levels)})

lo.proc[['2']] = list()
lo.proc[['2']][['fpqm']] = lapply(c(50,100,150,200,250,300,350,400),function(x){
  cv.svm.leave.one.out(req.dfs.rem$fpqm[,stable_genes_proc$fpqm$`2`$V1[1:x]], stages.levels)})
#400
lo.proc[['2']][['vs']] = lapply(c(50,100,150,200,250,300,350,400),function(x){
  cv.svm.leave.one.out(req.dfs.rem$vs[,stable_genes_proc$vs$`2`$V1[1:x]], stages.levels)})

lo.proc[['cho']] = list()
lo.proc[['cho']][['fpqm']] = lapply(c(50,100,150,200,250,300,350,400),function(x){
  cv.svm.leave.one.out(req.dfs.rem$fpqm[,stable_genes_proc$fpqm$cho$V1[1:x]], stages.levels)})
#400
lo.proc[['cho']][['vs']] = lapply(c(50,100,150,200,250,300,350,400),function(x){
  cv.svm.leave.one.out(req.dfs.rem$vs[,stable_genes_proc$vs$cho$V1[1:x]], stages.levels)})


##############Fine Tuning and varying weights and trying kernels#########
####SVM

##Trials
svm.fpqm.1 <- svm(req.dfs$fpqm[,stable_genes$fpqm$`1`$V1[1:150]], stages.levels, cross = 10)
table(stages.levels, svm.fpqm.1$fitted) ###Note this method not good as entire data is used for training

##Figuring out the best cost for linear kernel
##5 is the best cost
svm.tune.cost <- replicate(50,tune(svm, train.y = stages.levels, train.x = as.data.frame(req.dfs$vs[,stable_genes$vs$`1`$V1[1:400]]), 
                               ranges = list(cost = c(0.01,0.1,1,5,10,20,100)))$best.parameters)
##Figuring out the best gamma for radial kernel at cost 5
svm.tune.gama <- replicate(50,tune(svm, train.y = stages.levels, train.x = as.data.frame(req.dfs$vs[,stable_genes$vs$`1`$V1[1:400]]),
                                   cost = 5, ranges = list(gamma = c(0.01,0.1,1,10)), kernel = 'radial')$best.parameters)
table(unlist(svm.tune.gama)) ##0.1 is the winner

svm.tune.weights.1 <- replicate(10,tune(svm, train.y = stages.levels, train.x = as.data.frame(req.dfs$vs[,stable_genes$vs$`1`$V1[1:400]]), 
                               class.weights = c('stage i' = 1, 'stage ii' = 7.8, 'stage iii' = 3.37, 'stage iv' = 11.47), 
                               cost = 5)$best.performance)
svm.tune.weights.2 <- replicate(10,tune(svm, train.y = stages.levels, train.x = as.data.frame(req.dfs$vs[,stable_genes$vs$`1`$V1[1:400]]), 
                                        class.weights = c('stage i' = 11, 'stage ii' = 3, 'stage iii' = 7, 'stage iv' = 1), 
                                        cost = 5)$best.performance)

svm.tune <- tune(svm, train.y = stages.levels, train.x = as.data.frame(req.dfs$fpqm_log[,stable_genes$fpqm_log$cho$V1[1:150]]), 
                 ranges = list(cost = c(1,10,5,20)))

svm.tune <- replicate(100,tune(svm, train.y = stages.levels, train.x = as.data.frame(req.dfs$fpqm[,stable_genes$fpqm_log$f$indexes[1:150]]), 
                               class.weights = c('stage i' = 1, 'stage ii' = 7.8, 'stage iii' = 3.37, 'stage iv' = 11.47), 
                               ranges = list(cost = c(1,10,5,100)))$best.performance)
mean(svm.tune)

###Svms based on above
cv.svm.leave.one.out(req.dfs$vs[,stable_genes$vs$`1`$V1[1:400]], stages.levels, cost = 5)
cv.svm.leave.one.out(req.dfs$vs[,stable_genes$vs$`1`$V1[1:400]], stages.levels, cost = 5, kernel = 'radial',
                     gamma = 0.01)
o1.tune1 = cv.svm.leave.one.out(req.dfs$vs[,stable_genes$vs$`1`$V1[1:400]], stages.levels,
                          class.weights = c('stage i' = 1, 'stage ii' = 7.8, 'stage iii' = 3.37, 'stage iv' = 11.47), 
                                cost = 5)
o1.tune2 = cv.svm.leave.one.out(req.dfs$vs[,stable_genes$vs$`1`$V1[1:400]], stages.levels,
                                class.weights = c('stage i' = 11, 'stage ii' = 3, 'stage iii' = 7, 'stage iv' = 1), 
                                cost = 5)
o1.tune3 = cv.svm.leave.one.out(req.dfs$vs[,stable_genes$vs$`1`$V1[1:400]], stages.levels,
                                class.weights = c('stage i' = 1, 'stage ii' = 3, 'stage iii' = 2, 'stage iv' = 5), 
                                cost = 5)
o1.tune4 = cv.svm.leave.one.out(req.dfs$vs[,stable_genes$vs$`1`$V1[1:400]], stages.levels,
                                class.weights = c('stage i' = 1, 'stage ii' = 20, 'stage iii' = 2, 'stage iv' = 25), 
                                cost = 5)
a = cv.svm(req.dfs$vs[,stable_genes$vs$`1`$V1[1:400]],5, stages.levels,
class.weights = c('stage i' = 1, 'stage ii' = 20, 'stage iii' = 2, 'stage iv' = 25), 
cost = 5)