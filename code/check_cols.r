A <- readRDS('nhanes_data/A_pool.RDS')
cat('RIDRETH1 in cols:', 'RIDRETH1' %in% names(A), '\n')
print(table(A$RIDRETH1, useNA = 'ifany'))
cat('---check dm/smoke/drink/INDFMPIR/DMDEDUC2 missing---\n')
for (v in c('dm','smoke','drink','INDFMPIR','DMDEDUC2','BMXBMI')) {
  cat(v, ': NA =', sum(is.na(A[[v]])), '\n')
}
