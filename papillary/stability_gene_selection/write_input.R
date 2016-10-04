req.dfs <- lapply(only.tumor.reported$dfs, function(x)
{
  if(typeof(x) == 'S4')
    t.data.frame(as.data.frame(assay(x)))
  else
    t.data.frame(as.data.frame(x))
})
names(req.dfs) <- names(only.tumor.reported$dfs)
write.dfs.csvs('data',req.dfs)
