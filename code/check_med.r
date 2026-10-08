# W3 中介分析最小验证: 找出 boot 重放失败的精确原因
library(dplyr)
library(foreign)
library(mediation)
A <- readRDS('nhanes_data/A_pool.RDS')
fat_rows <- list()
for (s in c("D","E","F")) {
  sub <- c(D="A_2005_2006", E="A_2007_2008", F="A_2009_2010")[[s]]
  df <- as.data.frame(foreign::read.xport(file.path('nhanes_data', sub, paste0("DR1TOT_", s, ".XPT"))))
  keep <- intersect(c("SEQN","DR1TKCAL","DR1TTFAT","DR1DRSTZ"), names(df))
  fat_rows[[s]] <- df[, keep]
}
fat <- bind_rows(fat_rows) %>%
  filter(DR1DRSTZ == 1) %>% group_by(SEQN) %>% slice(1) %>% ungroup() %>%
  mutate(fat_pct = ifelse(DR1TKCAL > 0, DR1TTFAT * 9 / DR1TKCAL * 100, NA))
Ad <- A %>% left_join(fat, by="SEQN") %>%
  filter(!is.na(stool_grp), !is.na(tot25d)) %>%
  mutate(vtd_cut = ifelse(tot25d >= 50, 1, 0),
         diarrhea = as.integer(stool_grp=="diarrhea"),
         dm = ifelse(is.na(dm),0,dm), smoke = ifelse(is.na(smoke),0,smoke)) %>%
  filter(complete.cases(RIDAGEYR, RIAGENDR, BMXBMI, RIDRETH1))
covs <- "RIDAGEYR + RIAGENDR + RIDRETH1 + DMDEDUC2 + INDFMPIR + dm + smoke + drink"
d <- Ad %>% filter(!is.na(BMXBMI))
cat('n =', nrow(d), '\n')
# 构造两模型 (formula 对象)
med_f <- BMXBMI ~ vtd_cut + RIDAGEYR + RIAGENDR + RIDRETH1 + DMDEDUC2 + INDFMPIR + dm + smoke + drink
out_f <- diarrhea ~ vtd_cut + BMXBMI + RIDAGEYR + RIAGENDR + RIDRETH1 + DMDEDUC2 + INDFMPIR + dm + smoke + drink
m0 <- glm(med_f, data = d)
o0 <- glm(out_f, data = d, family = binomial())
cat('single models OK\n')
# mediate 不带 boot
md <- mediate(m0, o0, treat='vtd_cut', mediator='BMXBMI', boot=FALSE, sims=100)
print(summary(md)$d0); cat('no-boot OK\n')
# mediate 带 boot, 但 formula 用字面量直接写在 glm() 里 (重放时环境里有)
set.seed(1)
m1 <- glm(BMXBMI ~ vtd_cut + RIDAGEYR + RIAGENDR + RIDRETH1 + DMDEDUC2 + INDFMPIR + dm + smoke + drink, data = d)
o1 <- glm(diarrhea ~ vtd_cut + BMXBMI + RIDAGEYR + RIAGENDR + RIDRETH1 + DMDEDUC2 + INDFMPIR + dm + smoke + drink, data = d, family = binomial())
md2 <- mediate(m1, o1, treat='vtd_cut', mediator='BMXBMI', boot=TRUE, sims=50)
cat('boot with literal formula OK\n')
