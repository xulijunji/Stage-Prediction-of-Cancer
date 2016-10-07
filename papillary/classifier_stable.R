###using genes from stable_gene.R
###req.dfs from write_input.R

library(e1071)
library(caret)
library(class)

##KNN
##Based on 1
pred.knn.fpqm.1 <- knn.cv(req.dfs$fpqm[,stable_genes$fpqm$`1`$V1[1:150]], cl = stages.levels, k = 5)
table(stages.levels, pred.knn.fpqm.1)

pred.knn.fpqm_log.1 <- knn.cv(req.dfs$fpqm[,stable_genes$fpqm_log$`1`$V1[1:150]], cl = stages.levels, k = 4)
table(stages.levels, pred.knn.fpqm_log.1)

pred.knn.nt.1 <- knn.cv(req.dfs$fpqm[,stable_genes$nt$`1`$V1[1:150]], cl = stages.levels, k = 5)
table(stages.levels, pred.knn.nt.1)

pred.knn.vs.1 <- knn.cv(req.dfs$fpqm[,stable_genes$vs$`1`$V1[1:150]], cl = stages.levels, k = 5)
table(stages.levels, pred.knn.vs.1)

#Based on 2
pred.knn.fpqm.2 <- knn.cv(req.dfs$fpqm[,stable_genes$fpqm$`2`$V1[1:150]], cl = stages.levels, k = 3)
table(stages.levels, pred.knn.fpqm.2)

pred.knn.fpqm_log.2 <- knn.cv(req.dfs$fpqm[,stable_genes$fpqm_log$`2`$V1[1:150]], cl = stages.levels, k = 4)
table(stages.levels, pred.knn.fpqm_log.2)

pred.knn.nt.2 <- knn.cv(req.dfs$fpqm[,stable_genes$nt$`2`$V1[1:150]], cl = stages.levels, k = 5)
table(stages.levels, pred.knn.nt.2)

pred.knn.vs.2 <- knn.cv(req.dfs$fpqm[,stable_genes$vs$`2`$V1[1:150]], cl = stages.levels, k = 5)
table(stages.levels, pred.knn.vs.2)

###Based on cho
pred.knn.fpqm.cho <- knn.cv(req.dfs$fpqm[,stable_genes$fpqm$cho$V1[1:150]], cl = stages.levels, k = 4)
table(stages.levels, pred.knn.fpqm.cho)

pred.knn.fpqm_log.cho <- knn.cv(req.dfs$fpqm[,stable_genes$fpqm_log$cho$V1[1:150]], cl = stages.levels, k = 4)
table(stages.levels, pred.knn.fpqm_log.cho)

pred.knn.nt.2 <- knn.cv(req.dfs$fpqm[,stable_genes$nt$`2`$V1[1:150]], cl = stages.levels, k = 5)
table(stages.levels, pred.knn.nt.2)

pred.knn.vs.2 <- knn.cv(req.dfs$fpqm[,stable_genes$vs$`2`$V1[1:150]], cl = stages.levels, k = 5)
table(stages.levels, pred.knn.vs.2)

##Based on F
pred.knn.fpqm.f <- knn.cv(req.dfs$fpqm[,stable_genes$fpqm$f$V1[1:150]], cl = stages.levels, k = 4)
table(stages.levels, pred.knn.fpqm.f)

pred.knn.fpqm_log.f <- knn.cv(req.dfs$fpqm[,stable_genes$fpqm_log$f$V1[1:150]], cl = stages.levels, k = 4)
table(stages.levels, pred.knn.fpqm_log.f)

pred.knn.nt.2 <- knn.cv(req.dfs$fpqm[,stable_genes$nt$`2`$V1[1:150]], cl = stages.levels, k = 5)
table(stages.levels, pred.knn.nt.2)

pred.knn.vs.2 <- knn.cv(req.dfs$fpqm[,stable_genes$vs$`2`$V1[1:150]], cl = stages.levels, k = 5)
table(stages.levels, pred.knn.vs.2)


####SVM

##Trials
svm.fpqm.1 <- svm(req.dfs$fpqm[,stable_genes$fpqm$`1`$V1[1:150]], stages.levels, cross = 10)
table(stages.levels, svm.fpqm.1$fitted) ###Note this method not good as entire data is used for training


svm.tune <- replicate(100,tune(svm, train.y = stages.levels, train.x = as.data.frame(req.dfs$fpqm_log[,stable_genes$fpqm_log$cho$V1[1:150]]), 
                           class.weights = c('stage i' = 1, 'stage ii' = 3, 'stage iii' = 2.5, 'stage iv' = 3), 
                 ranges = list(cost = c(1,10,5,100)))$best.performance)

svm.tune <- replicate(100,tune(svm, train.y = stages.levels, train.x = as.data.frame(req.dfs$fpqm_log[,stable_genes$fpqm_log$cho$V1[1:150]]), 
                               class.weights = c('stage i' = 1, 'stage ii' = 1.5, 'stage iii' = 1.2, 'stage iv' = 1.5), 
                               ranges = list(cost = c(1,10,5,100)))$best.performance)

svm.tune <- replicate(100,tune(svm, train.y = stages.levels, train.x = as.data.frame(req.dfs$fpqm_log[,stable_genes$fpqm_log$cho$V1[1:150]]), 
                               class.weights = c('stage i' = 1, 'stage ii' = 1.5, 'stage iii' = 1.2, 'stage iv' = 1.5), 
                               ranges = list(cost = c(1,10,5,100)))$best.performance)

svm.tune <- tune(svm, train.y = stages.levels, train.x = as.data.frame(req.dfs$fpqm_log[,stable_genes$fpqm_log$cho$V1[1:150]]), 
                 ranges = list(cost = c(1,10,5,20)))

svm.tune <- replicate(100,tune(svm, train.y = stages.levels, train.x = as.data.frame(req.dfs$fpqm[,stable_genes$fpqm_log$f$indexes[1:150]]), 
                 class.weights = c('stage i' = 1, 'stage ii' = 7.8, 'stage iii' = 3.37, 'stage iv' = 11.47), 
                 ranges = list(cost = c(1,10,5,100)))$best.performance)
mean(svm.tune)

cv.svm(req.dfs$fpqm_log[,stable_genes$fpqm_log$`1`$V1[1:150]], 10, stages.levels)
cv.svm(req.dfs$vs[,stable_genes$vs$`1`$V1[1:150]], 10, stages.levels)
cv.svm(req.dfs$nt[,stable_genes$nt$`1`$V1[1:150]], 10, stages.levels)


cv.svm.leave.one.out(req.dfs$nt[,stable_genes$nt$`1`$V1[1:150]],  stages.levels)

cv.svm(req.dfs$fpqm[,stable_genes$fpqm$`2`$V1[1:150]], 5, stages.levels)
cv.svm(req.dfs$fpqm_log[,stable_genes$fpqm_log$`2`$V1[1:150]], 10, stages.levels)
cv.svm(req.dfs$vs[,stable_genes$vs$`2`$V1[1:150]], 10, stages.levels)
cv.svm(req.dfs$nt[,stable_genes$nt$`2`$V1[1:150]], 10, stages.levels)

cv.svm(req.dfs$fpqm[,stable_genes$fpqm$cho$V1[1:150]], 5, stages.levels)
cv.svm(req.dfs$fpqm_log[,stable_genes$fpqm_log$cho$V1[1:150]], 10, stages.levels)
cv.svm(req.dfs$vs[,stable_genes$vs$`1`$cho[1:150]], 10, stages.levels)
cv.svm(req.dfs$nt[,stable_genes$nt$`1`$cho[1:150]], 10, stages.levels)

cv.svm(req.dfs$fpqm[,stable_genes$fpqm$f$indexes[1:150]], 5, stages.levels)
cv.svm(req.dfs$fpqm_log[,stable_genes$fpqm_log$f$indexes[1:150]], 10, stages.levels)
cv.svm(req.dfs$vs[,stable_genes$vs$f$indexes[1:150]], 10, stages.levels)
cv.svm(req.dfs$nt[,stable_genes$nt$f$indexes[1:150]], 10, stages.levels)
