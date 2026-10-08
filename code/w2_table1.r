library(survey)
library(dplyr)
A <- readRDS('nhanes_data/A_pool.RDS')
Ad <- A %>%
  filter(!is.na(stool_grp), !is.na(tot25d)) %>%
  mutate(
    vtd_cut = ifelse(tot25d >= 50, 1, 0),
    dm = ifelse(is.na(dm), 0, dm), smoke = ifelse(is.na(smoke), 0, smoke),
    sex = factor(RIAGENDR, 1:2, c('Male','Female')),
    race = factor(RIDRETH1, 1:5, c('Mexican','Other-Hisp','NH-White','NH-Black','Other')),
    educ = factor(DMDEDUC2, 1:5, c('lt9th','9-11th','HS','SomeCol','College')),
    diabetes = factor(dm, 0:1, c('No','Yes')),
    smoking = factor(smoke, 0:1, c('Never','Ever'))
  ) %>%
  filter(complete.cases(RIDAGEYR, RIAGENDR, BMXBMI, RIDRETH1))
des <- svydesign(ids = ~SDMVPSU, strata = ~SDMVSTRA, weights = ~wt_mec6, nest = TRUE, data = Ad)

rows <- list()
add_row <- function(var, level, g1, g2, g3, ov) {
  rows[[length(rows)+1]] <<- data.frame(var=var, level=level, constipation=g1, diarrhea=g2, normal=g3, overall=ov)
}
# n (unweighted)
n_g <- table(Ad$stool_grp)
add_row('n (unweighted)', '', n_g['constipation'], n_g['diarrhea'], n_g['normal'], nrow(Ad))
# continuous
for (v in c('RIDAGEYR','BMXBMI','INDFMPIR','tot25d')) {
  m <- svyby(as.formula(paste('~', v)), ~stool_grp, des, svymean, na.rm=TRUE)
  mo <- svymean(as.formula(paste('~', v)), des, na.rm=TRUE)
  lbl <- switch(v, RIDAGEYR='Age, years', BMXBMI='BMI, kg/m2',
                INDFMPIR='Family income-to-poverty ratio', tot25d='Serum 25(OH)D, nmol/L')
  add_row(lbl, 'mean',
          m[m$stool_grp=='constipation', v], m[m$stool_grp=='diarrhea', v],
          m[m$stool_grp=='normal', v], as.numeric(mo))
}
# categorical
for (nm in c('sex','race','educ','diabetes','smoking')) {
  m <- svyby(as.formula(paste('~', nm)), ~stool_grp, des, svymean, na.rm=TRUE)
  mo <- svymean(as.formula(paste('~', nm)), des, na.rm=TRUE)
  lvls <- levels(Ad[[nm]])
  for (lv in lvls) {
    cn <- paste0(nm, lv)
    add_row(nm, lv,
            m[m$stool_grp=='constipation', cn], m[m$stool_grp=='diarrhea', cn],
            m[m$stool_grp=='normal', cn], as.numeric(mo[cn]))
  }
}
res <- do.call(rbind, rows)
write.csv(res, 'W2_results/W2_table1_weighted.csv', row.names=FALSE)
print(res, digits=3)
