library(survey)
library(splines)
A <- readRDS('nhanes_data/A_pool.RDS')
Ad <- A[!is.na(A$stool_grp) & !is.na(A$tot25d) & complete.cases(A[, c('RIDAGEYR','RIAGENDR','BMXBMI','RIDRETH1')]), ]
sb <- ns(Ad$tot25d, df = 4)
Ad$s1 <- sb[, 1]; Ad$s2 <- sb[, 2]; Ad$s3 <- sb[, 3]; Ad$s4 <- sb[, 4]
des <- svydesign(ids = ~SDMVPSU, strata = ~SDMVSTRA, weights = ~wt_mec6, nest = TRUE, data = Ad)
# locate which formula term is flagged missing
f <- diarrhea ~ s1 + s2 + s3 + s4 + RIDAGEYR + RIAGENDR + RIDRETH1 + BMXBMI + DMDEDUC2 + INDFMPIR + dm + smoke + drink
allv <- all.vars(f[[3]])
missing <- setdiff(allv, names(des$variables))
cat('formula vars missing from design:', missing, '\n')
cat('dm class in design:', class(des$variables$dm), ' | NA count:', sum(is.na(des$variables$dm)), '\n')
cat('drink class:', class(des$variables$drink), '| NA count:', sum(is.na(des$variables$drink)), '\n')
