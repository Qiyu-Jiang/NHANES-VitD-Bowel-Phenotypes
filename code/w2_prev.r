library(survey); library(dplyr)
A <- readRDS('nhanes_data/A_pool.RDS')
Ad <- A %>% filter(RIDAGEYR >= 20, !is.na(stool_grp), !is.na(tot25d),
                   complete.cases(RIDAGEYR, RIAGENDR, BMXBMI, RIDRETH1)) %>%
  mutate(vtd_q = ntile(tot25d, 4),
         constipation = as.integer(stool_grp=='constipation'),
         diarrhea = as.integer(stool_grp=='diarrhea'))
des <- svydesign(ids=~SDMVPSU, strata=~SDMVSTRA, weights=~wt_mec6, nest=TRUE, data=Ad)
for (oc in c('diarrhea','constipation')) {
  ft <- svyglm(as.formula(paste0(oc, ' ~ vtd_q')), design=des, family=quasibinomial())
  cat('trend M1', oc, ':', sprintf('%.4f', summary(ft)$coefficients['vtd_q',4]), '\n')
}
