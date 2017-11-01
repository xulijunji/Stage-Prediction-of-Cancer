load('environment/accuracy_feature/updated/new_data/test_pred.RData')
load('environment/accuracy_feature/updated/new_data/fea.RData')
sum(test.pred.trial$varSelRF$shrunken$atleast_2 == test.pred.trial$`deseq21 fold`$shrunken$atleast_4)

g_1fold_1up <- Reduce(intersect, create.list.venn(fea.trial.list.up, 1, 1))
g_1fold_1_int_diff <- intersect(g_1fold_1up, diff.genes.up$`2`)
r_1fold_int <- final.res(vst_tumor_tum, train.trial.ind, test.trial.ind, stages.levels.comb,
                 g_1fold_1up, 10)
micr_1fold <- final.res(micr_df, 1:34, 1:34, as.factor(sample_micro_info$stage),
                                    intersect(micr_genes, g_1fold_1_int_diff), 10)
get.aucs(r_1fold_int[[2]])
get.aucs(micr_1fold[[1]])
length(intersect(g_1fold_1up, fea.trial.list$varSelRF$atleast_3))
length(intersect(g_1fold_1_int_diff, fea.trial.list$varSelRF$atleast_3))
r_int_g_1fold_int_diff <- final.res(vst_tumor_tum, train.trial.ind, test.trial.ind, stages.levels.comb,
                                 g_1fold_1_int_diff, 10)
micr_int_g_1fold_diff <- final.res(micr_df, 1:34, 1:34, as.factor(sample_micro_info$stage),
                intersect(micr_genes, g_1fold_1_int_diff), 10)
get.aucs(r_int_g_1fold_int_diff[[2]])
get.aucs(micr_int_g_1fold_1_var[[1]])

g_1fold_3 <- Reduce(intersect, create.list.venn(fea.trial.list.up, 1, 3))
g_1fold_4 <- Reduce(intersect, create.list.venn(fea.trial.list.up, 1, 4))
g_2fold_1 <- Reduce(intersect, create.list.venn(fea.trial.list, 2, 1))


diff_g_int_3 <- intersect(diff.genes.up$`2`, g_1fold_3)
diff_g_int_4 <- intersect(diff.genes.up$`2`, g_1fold_4)

micr_3 <- final.res(micr_df, 1:34, 1:34, as.factor(sample_micro_info$stage),
    intersect(intersect(intersect(micr_genes, g_1fold_3), diff.genes.up$`2`), diff.genes.up$`2`), 10)
get.aucs(micr_3[[1]])

r_2 <- final.res(vst_tumor_tum, train.trial.ind, test.trial.ind, stages.levels.comb,
                 g_2, 10)
r_3 <- final.res(vst_tumor_tum, train.trial.ind, test.trial.ind, stages.levels.comb,
                 g_1fold_3, 10)

r_1fold_3 <- final.res(vst_tumor_tum, train.trial.ind, test.trial.ind, stages.levels.comb,
                    g_1fold_3, 10)
r_1fold_4 <- final.res(vst_tumor_tum, train.trial.ind, test.trial.ind, stages.levels.comb,
                     g_1fold_4, 10)
r_diff_int_3 <- final.res(vst_tumor_tum, train.trial.ind, test.trial.ind, stages.levels.comb,
                          diff_g_int_3, 10)
r_diff_int_4 <- final.res(vst_tumor_tum, train.trial.ind, test.trial.ind, stages.levels.comb,
                          diff_g_int_4, 10)
micr_diff_int_4 <- final.res(micr_df, 1:34, 1:34, as.factor(sample_micro_info$stage),
                             intersect(diff_g_int_4, micr_genes), 10)
micr_diff_int_3 <- final.res(micr_df, 1:34, 1:34, as.factor(sample_micro_info$stage),
                             intersect(diff_g_int_3, micr_genes), 10)
get.aucs(r_diff_int_4[[2]])
get.aucs(micr_diff_int_3[[1]])

r_diff_r_1 <- final.res(vst_tumor_tum, train.trial.ind, test.trial.ind, stages.levels.comb,
                          intersect(diff.genes.up$`2`, g_1), 10)

create.boxplots(intersect(Reduce(intersect, create.list.venn(fea.trial.list, 1, 1)), diff.genes.low$`2`)[1],
                  vs_normal_comb_reported, sample.info.all.rep$stage.type)

get.aucs(r_1[[2]])
get.aucs(r_2[[2]])
get.aucs(r_int_3[[2]])
get.aucs(r_int_4[[2]])
get.aucs(r_diff_int_3[[2]])
get.aucs(r_diff_r_1[[2]])

g.tt <- setdiff(g_1fold_1_int_diff, fea.trial.list.up$varSelRF$atleast_2)
r_gtt <- final.res(vst_tumor_tum, train.trial.ind, test.trial.ind, stages.levels.comb,
                        g.tt, 10)
get.aucs(r_gtt[[2]])

r_var3_up <- final.res(vst_tumor_tum, train.trial.ind, test.trial.ind, stages.levels.comb,
                       fea.trial.list.up$varSelRF$atleast_2, 10)
get.aucs(r_var3_up[[2]])

g_1fold_2up <- Reduce(intersect, create.list.venn(fea.trial.list.up, 1, 2))
r_1fold_2up <- final.res(vst_tumor_tum, train.trial.ind, test.trial.ind, stages.levels.comb,
                       g_1fold_2up, 10)
get.aucs(r_1fold_int[[2]])
g_int_var <- intersect(g_1fold_1_int_diff, fea.trial.list.up$varSelRF$atleast_2)
r_int_g_1fold_1_var <- final.res(vst_tumor_tum, train.trial.ind, test.trial.ind, stages.levels.comb,
                              g_int_var, 10)
get.aucs(r_int_g_1fold_1_var[[2]])

res.wcgna <- final.res(vst_tumor_tum, train.trial.ind, test.trial.ind, stages.levels.comb,
                       wcgna.genes.ent.ens, 10)
micr.wcgna <- final.res(micr_df, 1:34, 1:34, as.factor(sample_micro_info$stage),
                        intersect(wcgna.genes.ent.ens, micr_genes), 10)
get.aucs(res.wcgna[[2]])
get.aucs(micr.wcgna[[1]])

uni <- union(fea.trial.list$varSelRF$atleast_3, g_1fold_1_int_diff)
r <- final.res(vst_tumor_tum, train.trial.ind, test.trial.ind, stages.levels.comb,
              uni, 10)
get.aucs(r[[2]])
r_4 <- final.res(vst_tumor_tum, train.trial.ind, test.trial.ind, stages.levels.comb,
                 g_1fold_4, 10)
get.aucs(r_4[[2]])
