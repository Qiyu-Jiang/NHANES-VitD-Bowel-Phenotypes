# 校准后全表: M1/M2/M3 x (四分位/切点) x (腹泻/便秘) + 基线均值
library(survey); library(dplyr)
A <- readRDS('nhanes_data/A_pool.RDS')
Ad <- A %>% filter(RIDAGEYR >= 20, !is.na(stool_grp), !is.na(tot25d),
                   complete.cases(RIDAGEYR, RIAGENDR, BMXBMI, RIDRETH1)) %>%
  mutate(vtd_q = ntile(tot25d, 4),
         vtd_qf = factor(vtd_q, labels=c('Q1','Q2','Q3','Q4')),
         vtd_cut = ifelse(tot25d >= 50, 1, 0),
         constipation = as.integer(stool_grp=='constipation'),
         diarrhea = as.integer(stool_grp=='diarrhea'),
         dm = ifelse(is.na(dm),0,dm), smoke = ifelse(is.na(smoke),0,smoke))
des <- svydesign(ids=~SDMVPSU, strata=~SDMVSTRA, weights=~wt_mec6, nest=TRUE, data=Ad)
m2s <- 'RIDAGEYR + RIAGENDR + RIDRETH1'
m3s <- paste(m2s, '+ BMXBMI + DMDEDUC2 + INDFMPIR + dm + smoke + drink')
cat('---- FULL OR TABLE (calibrated) ----\n')
for (oc in c('diarrhea','constipation')) {
  for (ex in c('vtd_qf','vtd_cut')) {
    for (mi in list(c('M1','1'), c('M2',m2s), c('M3',m3s))) {
      f <- as.formula(paste0(oc, ' ~ ', ex, ' + ', mi[2]))
      fit <- svyglm(f, design=des, family=quasibinomial())
      co <- summary(fit)$coefficients
      ci <- suppressMessages(exp(confint(fit)))
      vars <- grep(ex, rownames(co), value=TRUE, fixed=TRUE)
      for (v in vars) {
        cat(sprintf('%s | %s | %s | %s: OR=%.3f (%.3f-%.3f) p=%.4f\n',
            oc, mi[1], ex, v, exp(co[v,1]), ci[v,1], ci[v,2], co[v,4]))
      }
    }
  }
  ft <- svyglm(as.formula(paste0(oc, ' ~ vtd_q')), design=des, family=quasibinomial())
  cat('trend M1', oc, ':', sprintf('%.4f', summary(ft)$coefficients['vtd_q',4]), '\n')
  ft3 <- svyglm(as.formula(paste0(oc, ' ~ vtd_q + ', m3s)), design=des, family=quasibinomial())
  cat('trend M3', oc, ':', sprintf('%.4f', summary(ft3)$coefficients['vtd_q',4]), '\n')
}
# 基线均值 by group (校准后 25D)
cat('\n---- baseline 25D means by group (calibrated) ----\n')
m <- svyby(~tot25d, ~stool_grp, des, svymean, na.rm=TRUE)
print(m)
mo <- svymean(~tot25d, des, na.rm=TRUE)
cat('overall mean 25D:', sprintf('%.1f', as.numeric(mo)), '\n')
# 加权患病率 (样本未变, 但正式输出)
for (g in c('constipation','normal','diarrhea')) {
  Ad$y <- as.integer(Ad$stool_grp==g)
  des2 <- update(des, y = Ad$y)
  mm <- svymean(~y, des2); ci <- confint(mm)
  cat(g, ': ', sprintf('%.2f%% (95%%CI %.2f-%.2f)\n', 100*as.numeric(mm[1]), 100*ci[1,1], 100*ci[1,2]), '\n')
}
