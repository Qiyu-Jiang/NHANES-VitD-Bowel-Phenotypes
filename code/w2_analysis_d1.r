# =====================================================
# W2 主分析脚本  (w2_analysis_d1.R)
# A 部分: 血清 25(OH)D × Bristol 排便表型 (2005-2010 合并, 加权)
# 前置: w1_merge_d1.R 已跑通 (A_pool.RDS 就绪)
# 产出: W2_results/ 下 OR 表 + 剂量-反应图 + 亚组森林图 + 敏感性表
# 运行: Rscript w2_analysis_d1.R  (在 DELIVERY 目录下)
# 生成: 2026-09-22 · AutoClaw
# =====================================================

library(survey)
library(dplyr)
library(splines)
library(ggplot2)

ROOT <- "nhanes_data"
OUT  <- "W2_results"
dir.create(OUT, showWarnings = FALSE)

A <- readRDS(file.path(ROOT, "A_pool.RDS"))

Ad <- A |>
  filter(!is.na(stool_grp), !is.na(tot25d)) |>
  mutate(
    vtd_cut = ifelse(tot25d >= 50, 1, 0),
    vtd10   = tot25d / 10,   # v5: continuous exposure, per 10 nmol/L
    age_grp = ifelse(RIDAGEYR >= 60, "60plus", "lt60"),
    bmi_grp = case_when(BMXBMI < 25 ~ "normal", BMXBMI < 30 ~ "overweight", TRUE ~ "obese"),
    constipation = as.integer(stool_grp == "constipation"),
    diarrhea     = as.integer(stool_grp == "diarrhea"),
    dm    = ifelse(is.na(dm), 0, dm),
    smoke = ifelse(is.na(smoke), 0, smoke)
  ) |>
  filter(complete.cases(RIDAGEYR, RIAGENDR, BMXBMI, RIDRETH1)) |>
  mutate(
    # v2 (2026-10-05): quartiles computed on the final analytical sample
    # (after core-covariate complete-case filtering), n = 13,266
    vtd_q  = ntile(tot25d, 4),
    vtd_qf = factor(vtd_q, labels = c("Q1(low)", "Q2", "Q3", "Q4(high)"))
  )

# 预计算 4 结点自然样条基列 (避免 ns() 写进 survey 公式触发变量检查兼容性问题)
sb <- ns(Ad$tot25d, df = 4)
Ad$s1 <- sb[, 1]; Ad$s2 <- sb[, 2]; Ad$s3 <- sb[, 3]; Ad$s4 <- sb[, 4]

cat("W2-A 分析样本:", nrow(Ad), "\n")

des <- svydesign(ids = ~SDMVPSU, strata = ~SDMVSTRA, weights = ~wt_mec6,
                 nest = TRUE, data = Ad)

cov_m2 <- "RIDAGEYR + RIAGENDR + RIDRETH1"
cov_m3 <- paste(cov_m2, "+ BMXBMI + DMDEDUC2 + INDFMPIR + dm + smoke + drink")

# -----------------------------------------------------
# 1. 主分析
# -----------------------------------------------------
run_or <- function(outcome, exposure, covs) {
  f <- as.formula(paste0(outcome, " ~ ", exposure, " + ", covs))
  fit <- svyglm(f, design = des, family = quasibinomial())
  co <- summary(fit)$coefficients
  ci <- suppressMessages(exp(confint(fit)))
  vars <- grep(exposure, rownames(co), value = TRUE, fixed = TRUE)
  data.frame(outcome = outcome, term = vars,
             OR = exp(co[vars, 1]), LCL = ci[vars, 1], UCL = ci[vars, 2],
             p = co[vars, 4], row.names = NULL)
}

res <- list()
for (oc in c("diarrhea", "constipation")) {
  res[[paste(oc, "M1_q")]]   <- run_or(oc, "vtd_qf", "1")
  res[[paste(oc, "M3_q")]]   <- run_or(oc, "vtd_qf", cov_m3)
  res[[paste(oc, "M1_cut")]] <- run_or(oc, "vtd_cut", "1")
  res[[paste(oc, "M3_cut")]] <- run_or(oc, "vtd_cut", cov_m3)
  # v5: continuous per-10 nmol/L parameterization (M1 / M3)
  res[[paste(oc, "M1_cont")]] <- run_or(oc, "vtd10", "1")
  res[[paste(oc, "M3_cont")]] <- run_or(oc, "vtd10", cov_m3)
}
res_df <- do.call(rbind, res)
res_df$set <- gsub("^\\d+\\.", "", row.names(res_df))
res_df$set <- gsub("\\.term.*$", "", res_df$set)
row.names(res_df) <- NULL
write.csv(res_df, file.path(OUT, "W2_main_OR_table.csv"), row.names = FALSE)
cat("\n---- 主分析 OR 表 (M3 全调整) ----\n")
print(res_df[grepl("M3", res_df$set), c("set", "term", "OR", "LCL", "UCL", "p")], digits = 3)

f_trend1 <- as.formula(paste("diarrhea ~ vtd_q +", cov_m3))
f_trend2 <- as.formula(paste("constipation ~ vtd_q +", cov_m3))
cat("\n腹泻型趋势 P (Q1-Q4 连续):",
    format.pval(summary(svyglm(f_trend1, design = des, family = quasibinomial()))$coefficients["vtd_q", 4], digits = 2), "\n")
cat("便秘型趋势 P (Q1-Q4 连续):",
    format.pval(summary(svyglm(f_trend2, design = des, family = quasibinomial()))$coefficients["vtd_q", 4], digits = 2), "\n")

# -----------------------------------------------------
# 2. 剂量-反应 (系数 + vcov 手工计算, 无 predict 兼容性问题)
# -----------------------------------------------------
f_s <- as.formula(paste("diarrhea ~ s1 + s2 + s3 + s4 +", cov_m3))
fit_s <- svyglm(f_s, design = des, family = quasibinomial())
beta <- coef(fit_s)
V    <- vcov(fit_s)

mode_val <- function(x) as.numeric(names(which.max(table(x, useNA = "no"))))
ref_cov <- c(RIDAGEYR = median(Ad$RIDAGEYR), RIAGENDR = mode_val(Ad$RIAGENDR),
             RIDRETH1 = mode_val(Ad$RIDRETH1), BMXBMI = median(Ad$BMXBMI, na.rm = TRUE),
             DMDEDUC2 = mode_val(Ad$DMDEDUC2), INDFMPIR = median(Ad$INDFMPIR, na.rm = TRUE),
             dm = 0, smoke = 0, drink = mode_val(Ad$drink))

gx <- seq(floor(quantile(Ad$tot25d, .01, na.rm = TRUE)),
          ceiling(quantile(Ad$tot25d, .99, na.rm = TRUE)), by = 1)
sb_grid <- predict(sb, gx)                       # ns 基在网格点的取值
sb_ref  <- predict(sb, median(Ad$tot25d, na.rm = TRUE))

mk_X <- function(basis_row, tot25d_unused = NULL) {
  c(1, basis_row, ref_cov[["RIDAGEYR"]], ref_cov[["RIAGENDR"]], ref_cov[["RIDRETH1"]],
    ref_cov[["BMXBMI"]], ref_cov[["DMDEDUC2"]], ref_cov[["INDFMPIR"]],
    ref_cov[["dm"]], ref_cov[["smoke"]], ref_cov[["drink"]])
}
X_grid <- t(apply(as.matrix(sb_grid), 1, mk_X))
X_ref  <- mk_X(sb_ref)

lp    <- as.numeric(X_grid %*% beta)
lp_se <- sqrt(rowSums((X_grid %*% V) * X_grid))
lp_ref <- as.numeric(X_ref %*% beta)

rd <- data.frame(tot25d = gx,
                 or  = exp(lp - lp_ref),
                 lcl = exp(lp - 1.96 * lp_se - lp_ref),
                 ucl = exp(lp + 1.96 * lp_se - lp_ref))
saveRDS(rd, file.path(OUT, "W2_rcs_points.RDS"))
p_rcs <- ggplot(rd, aes(x = tot25d, y = or)) +
  geom_ribbon(aes(ymin = lcl, ymax = ucl), fill = "#0f766e", alpha = .15) +
  geom_line(color = "#0f766e", linewidth = .8) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "grey55") +
  scale_x_continuous(breaks = seq(0, 200, 25)) +
  labs(x = "Serum 25(OH)D (nmol/L)", y = "OR of diarrhea-type stool (ref: median)",
       title = "Dose-response: serum 25(OH)D and diarrhea-type stool",
       caption = "NHANES 2005-2010, weighted, 4-knot natural spline, fully adjusted") +
  theme_minimal(base_size = 12)
ggsave(file.path(OUT, "Fig2_dose_response.png"), p_rcs, width = 8, height = 5, dpi = 300)
cat("\n剂量-反应图已出: Fig2_dose_response.png\n")

# -----------------------------------------------------
# 3. 亚组分析 + 交互检验 (腹泻型, vtd_cut)
# -----------------------------------------------------
sub_src <- c(RIAGENDR = "RIAGENDR", age_grp = "RIDAGEYR",
             bmi_grp = "BMXBMI", dm = "dm", smoke = "smoke")
drop_cov <- function(covs, token) {
  parts <- trimws(strsplit(covs, "+", fixed = TRUE)[[1]])
  paste(parts[parts != token], collapse = " + ")
}
sub_one <- function(var, lv, label) {
  cond <- if (is.character(lv)) paste0(var, " == '", lv, "'") else paste0(var, " == ", lv)
  d_sub <- eval(parse(text = paste0("subset(des, ", cond, ")")))
  covs_sub <- drop_cov(cov_m3, sub_src[[var]])
  fit <- svyglm(as.formula(paste("diarrhea ~ vtd_cut +", covs_sub)),
                design = d_sub, family = quasibinomial())
  co <- summary(fit)$coefficients; ci <- suppressMessages(exp(confint(fit)))
  data.frame(subgroup = label, level = as.character(lv), n = nrow(d_sub$variables),
             OR = exp(co["vtd_cut", 1]), LCL = ci["vtd_cut", 1], UCL = ci["vtd_cut", 2],
             p = co["vtd_cut", 4], row.names = NULL)
}
int_p <- function(var) {
  covs_int <- drop_cov(cov_m3, sub_src[[var]])
  f <- as.formula(paste("diarrhea ~ vtd_cut *", var, "+", covs_int))
  fit <- svyglm(f, design = des, family = quasibinomial())
  nm <- grep(paste0("vtd_cut[:.]", var), rownames(summary(fit)$coefficients), value = TRUE)[1]
  summary(fit)$coefficients[nm, 4]
}
sg <- rbind(
  do.call(rbind, lapply(c(1, 2), function(lv) sub_one("RIAGENDR", lv, "Sex"))),
  do.call(rbind, lapply(c("lt60", "60plus"), function(lv) sub_one("age_grp", lv, "Age"))),
  do.call(rbind, lapply(c("normal", "overweight", "obese"), function(lv) sub_one("bmi_grp", lv, "BMI"))),
  do.call(rbind, lapply(c(0, 1), function(lv) sub_one("dm", lv, "Diabetes"))),
  do.call(rbind, lapply(c(0, 1), function(lv) sub_one("smoke", lv, "Smoking")))
)
sg$interaction_p <- unlist(lapply(seq_len(nrow(sg)), function(i) {
  grp_name <- sg$subgroup[i]
  src_var  <- switch(grp_name, Sex = "RIAGENDR", Age = "age_grp", BMI = "bmi_grp",
                     Diabetes = "dm", Smoking = "smoke")
  int_p(src_var)
}))
sg$level <- ifelse(sg$subgroup == "Sex" & sg$level == "1", "Male",
            ifelse(sg$subgroup == "Sex" & sg$level == "2", "Female",
            ifelse(sg$subgroup == "Diabetes" & sg$level == "1", "Yes",
            ifelse(sg$subgroup == "Diabetes" & sg$level == "0", "No",
            ifelse(sg$subgroup == "Smoking" & sg$level == "1", "Ever",
            ifelse(sg$subgroup == "Smoking" & sg$level == "0", "Never",
            ifelse(sg$level == "lt60", "<60 years",
            ifelse(sg$level == "60plus", "≥60 years",
            ifelse(sg$level == "normal", "Normal weight",
            ifelse(sg$level == "overweight", "Overweight",
            ifelse(sg$level == "obese", "Obesity", sg$level)))))))))))
write.csv(sg, file.path(OUT, "W2_subgroup_forest.csv"), row.names = FALSE)
cat("\n---- 亚组分析 (腹泻型, 25D<50 vs >=50 nmol/L) ----\n")
print(sg[, c("subgroup", "level", "n", "OR", "LCL", "UCL", "p", "interaction_p")], digits = 3)

sg$level <- factor(sg$level, levels = rev(c("Male", "Female", "<60 years",
  "≥60 years", "Normal weight", "Overweight", "Obesity", "No", "Yes", "Never", "Ever")))
sg$subgroup <- factor(sg$subgroup, levels = c("Sex", "Age", "BMI", "Diabetes", "Smoking"))
p_f <- ggplot(sg, aes(x = OR, y = level)) +
  geom_point(color = "#0f766e", size = 2.2) +
  geom_errorbarh(aes(xmin = LCL, xmax = UCL), color = "#0f766e", height = .25) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "grey55") +
  facet_grid(subgroup ~ ., scales = "free_y", space = "free_y") +
  labs(x = "OR (95% CI)", y = NULL,
       title = "Subgroup analyses: serum 25(OH)D >=50 vs <50 nmol/L",
       caption = "Diarrhea-type stool, >=50 vs <50 nmol/L (ref: <50). NHANES 2005-2010, weighted") +
  theme_minimal(base_size = 12) +
  theme(strip.text.y = element_text(angle = 0, face = "bold"),
        panel.spacing.y = unit(0.6, "lines"))
ggsave(file.path(OUT, "Fig3_subgroup_forest.png"), p_f, width = 8, height = 5.5, dpi = 300)
cat("亚组森林图已出: Fig3_subgroup_forest.png\n")

# -----------------------------------------------------
# 4. 敏感性分析 (腹泻型, vtd_cut)
# -----------------------------------------------------
sens <- list()
add_sens <- function(nm, fit, term = "vtd_cut") {
  co <- summary(fit)$coefficients; ci <- suppressMessages(exp(confint(fit)))
  sens[[nm]] <<- c(OR = exp(co[term, 1]), LCL = ci[term, 1],
                   UCL = ci[term, 2], p = co[term, 4])
}
Ad_i <- Ad
Ad_i$INDFMPIR[is.na(Ad_i$INDFMPIR)] <- median(Ad_i$INDFMPIR, na.rm = TRUE)
Ad_i$DMDEDUC2[is.na(Ad_i$DMDEDUC2)] <- as.numeric(names(sort(table(Ad_i$DMDEDUC2), decreasing = TRUE))[1])
Ad_i$drink[is.na(Ad_i$drink)] <- 0
des_i <- svydesign(ids = ~SDMVPSU, strata = ~SDMVSTRA, weights = ~wt_mec6, nest = TRUE, data = Ad_i)
add_sens("S1_imputed", svyglm(as.formula(paste("diarrhea ~ vtd_cut +", cov_m3)), design = des_i, family = quasibinomial()))
add_sens("S2_no_laxative", svyglm(as.formula(paste("diarrhea ~ vtd_cut +", cov_m3)),
         design = subset(des, is.na(BHQ100) | BHQ100 != 1), family = quasibinomial()))
for (cyc in unique(Ad$cycle)) {
  add_sens(paste("S3_cycle", cyc),
           svyglm(as.formula(paste("diarrhea ~ vtd_cut +", cov_m3)),
                  design = subset(des, cycle == cyc), family = quasibinomial()))
}
add_sens("S4_trim_extremes", svyglm(as.formula(paste("diarrhea ~ vtd_cut +", cov_m3)),
         design = subset(des, tot25d >= 12.5 & tot25d <= 250), family = quasibinomial()))
add_sens("S5_constipation_cut", svyglm(as.formula(paste("constipation ~ vtd_cut +", cov_m3)),
         design = des, family = quasibinomial()))
# v5: continuous per-10 sensitivity battery (mirrors the cut-point battery, symmetric reporting)
for (cyc in unique(Ad$cycle)) {
  add_sens(paste("S6_cycle10", cyc),
           svyglm(as.formula(paste("diarrhea ~ vtd10 +", cov_m3)),
                  design = subset(des, cycle == cyc), family = quasibinomial()), term = "vtd10")
}
add_sens("S7_trim10", svyglm(as.formula(paste("diarrhea ~ vtd10 +", cov_m3)),
         design = subset(des, tot25d >= 12.5 & tot25d <= 250), family = quasibinomial()), term = "vtd10")

sens_df <- do.call(rbind, sens)
write.csv(data.frame(analysis = row.names(sens_df), sens_df), file.path(OUT, "W2_sensitivity.csv"), row.names = FALSE)
cat("\n---- 敏感性分析 ----\n")
print(sens_df, digits = 3)

# -----------------------------------------------------
# 5. v5 非线性检验: ns(4结点) 正交化非线性残差基 (3 df) 联合 Wald
#    口径锁定 2026-10-07 拍板: P = 0.709 (F(3,34) = 0.46)
# -----------------------------------------------------
Bnl  <- ns(Ad$tot25d, df = 4)
Xlin <- cbind(1, Ad$tot25d)
Bres <- Bnl - Xlin %*% solve(crossprod(Xlin), crossprod(Xlin, Bnl))
qrB  <- qr(Bres); rk <- qrB$rank
Qnl  <- qr.Q(qrB)[, seq_len(rk), drop = FALSE]
for (k in seq_len(rk)) Ad[[paste0("z", k)]] <- Qnl[, k]
des_nl <- svydesign(ids = ~SDMVPSU, strata = ~SDMVSTRA, weights = ~wt_mec6, nest = TRUE, data = Ad)
zt <- paste0("z", seq_len(rk), collapse = " + ")
fit_nl <- svyglm(as.formula(paste("diarrhea ~ vtd10 +", zt, "+", cov_m3)),
                 design = des_nl, family = quasibinomial())
tst_nl <- regTermTest(fit_nl, as.formula(paste("~", zt)), method = "Wald")
cat("\n---- 非线性检验 (ns4 正交化非线性残差基, 3 df Wald) ----\n")
print(tst_nl)
write.csv(data.frame(test = "nonlinearity_ns4_orth_residual_wald",
                     F = tst_nl$Ftest, df = rk, p = as.numeric(tst_nl$p)),
          file.path(OUT, "W2_nonlinearity_test.csv"), row.names = FALSE)

cat("\nW2-A 主分析完成, 全部结果写入 ", OUT, "\n", sep = "")
