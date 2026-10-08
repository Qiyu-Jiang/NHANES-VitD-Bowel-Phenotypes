library(survey); library(dplyr)
A <- readRDS('nhanes_data/A_pool.RDS')
Ad <- A %>% filter(RIDAGEYR >= 20, !is.na(stool_grp), !is.na(tot25d),
                   complete.cases(RIDAGEYR, RIAGENDR, BMXBMI, RIDRETH1)) %>%
  mutate(is_normal = as.integer(stool_grp == 'normal'))
des <- svydesign(ids=~SDMVPSU, strata=~SDMVSTRA, weights=~wt_mec6, nest=TRUE, data=Ad)
m <- svymean(~is_normal, des)
ci <- confint(m)
cat(sprintf('normal prevalence: %.2f%% (95%%CI %.2f-%.2f)\n',
    100*as.numeric(m[1]), 100*ci[1,1], 100*ci[1,2]))
