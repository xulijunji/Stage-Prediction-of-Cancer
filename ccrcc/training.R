load('../environment/kirc_data.RData')
load('../papillary/environment/accuracy_feature/updated/new_data/fea.RData')
source('../papillary/main/updated/classifier.R')
source('../papillary/decomp.R')
source('../papillary/main/updated/helper_func.R')
source('../papillary/main/final.R')
source('../papillary/paper/multiplot.R')

tr.indexes <- sample(seq(536), 429)
test.indexes <- setdiff(seq(536), tr.indexes)

tr.model <- get.train.model(data = norm.data, tr.ind = tr.indexes, 
                            fea.list = fea.trial.list.crc,
                            stages = as.factor(stages.ccrcc.comb))

cv.model <- get.cv.model(data = norm.data, tr.ind = tr.indexes, 
                         fea.list = fea.trial.list.crc, tr.model = tr.model,
                         stages = as.factor(stages.ccrcc.comb))

test.pred <- get.test.pred(data = norm.data, tr.ind = tr.indexes, te.ind = test.indexes,
                           fea.list = fea.trial.list.crc, stages = as.factor(stages.ccrcc.comb),
                           tr.model = tr.model, cv.model = cv.model)

test.results <- get.test.res(test.pred = test.pred, te.ind = test.indexes,
                             fea.names = names(test.pred), stages = stages.ccrcc.comb)
cv.results <- get.cv.res(tr.ind = tr.indexes, cv.model = cv.model, 
                         stages = stages.ccrcc.comb, fea.names = names(cv.model))
get.aucs(test.results$varSelRF$nb)
ccrcc.test.plot <- create.gridplot(test.results)
ggarrange(ccrcc.test.plot[[1]], ccrcc.test.plot[[2]], ccrcc.test.plot[[3]],
          ccrcc.test.plot[[4]], ccrcc.test.plot[[5]], ncol = 2,
           nrow = 3,
           common.legend = T, legend = 'bottom')

ccrcc.cv.plot <- create.gridplot(cv.results)
ggarrange(ccrcc.cv.plot[[1]], ccrcc.cv.plot[[2]], ccrcc.cv.plot[[3]],
          ccrcc.cv.plot[[4]], ccrcc.cv.plot[[5]], ncol = 2,
          nrow = 3,
          common.legend = T, legend = 'bottom')
