library(survey); library(dplyr)
A <- readRDS('nhanes_data/A_pool.RDS')
Ad <- A %>% filter(RIDAGEYR >= 20, !is.na(stool_grp), !is.na(tot25d),
                   complete.cases(RIDAGEYR, RIAGENDR, BMXBMI, RIDRETH1))
des <- svydesign(ids=~SDMVPSU, strata=~SDMVSTRA, weights=~wt_mec6, nest=TRUE, data=Ad)
# 1) 加权人口总量 (美国 20+ 成人, 2005-2010 合并代表量)
pop <- svytotal(~rep(1, nrow(Ad)), des, na.rm=TRUE)  # 不可靠; 直接用权重和
pop_total <- sum(Ad$wt_mec6)
cat(sprintf('weighted population (20+ US adults): %.0f\n', pop_total))
# 2) 各表型加权患病率 + 95% CI
Ad$y_const <- as.integer(Ad$stool_grp=='constipation')
Ad$y_diarr <- as.integer(Ad$stool_grp=='diarrhea')
des2 <- update(des, y_const = Ad$y_const, y_diarr = Ad$y_diarr)
for (nm in c('y_const','y_diarr')) {
  m <- svymean(as.formula(paste0('~', nm)), des2)
  ci <- confint(m)
  cat(nm, sprintf(': %.2f%% (95%%CI %.2f-%.2f)\n',
      100*as.numeric(m[1]), 100*ci[1,1], 100*ci[1,2]))
}
