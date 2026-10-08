# =====================================================
# W2-B 主分析脚本  (w2b_analysis_d1.R)
# B 部分: 胆囊切除史 × 胆汁酸/肝胆机制标志物 (2017-2020 P 合并, 加权)
# 主结局: ALT/AST/GGT 升高 (ULN 上限定义); 剂量: 术后时长
# 前置: w1_merge_d1.R 已跑通 (B_pool.RDS 就绪)
# 产出: W2_results/ 下 B 部分 OR 表 + 剂量-反应数据 + 中介分析框架数据
# 运行: Rscript w2b_analysis_d1.R  (在 DELIVERY 目录下)
# 生成: 2026-09-22 · AutoClaw
# =====================================================

library(survey)
library(dplyr)
library(splines)
library(ggplot2)
library(foreign)

ROOT <- "nhanes_data"
OUT  <- "W2_results"
dir.create(OUT, showWarnings = FALSE)

B <- readRDS(file.path(ROOT, "B_pool.RDS"))

# ---- 合入肝脏生化与血小板 (从官方 XPT 直读, 与 B_pool 按 SEQN 合并) ----
RD <- function(rel) as.data.frame(foreign::read.xport(file.path(ROOT, rel)))
bio <- RD(file.path("B_pmcq_2017_2020", "P_BIOPRO.XPT"))
bio <- bio[, intersect(c("SEQN","LBXSATSI","LBXSASSI","LBXSGTSI","LBXSBTSI"), names(bio))]
cbc <- RD(file.path("B_pmcq_2017_2020", "P_CBC.XPT"))
cbc <- cbc[, intersect(c("SEQN","LBXPLTSI"), names(cbc))]

Bd <- B |>
  left_join(bio, by = "SEQN") |> left_join(cbc, by = "SEQN") |>
  mutate(
    # 肝胆标志物升高 (ULN 口径, 参照试剂常用上限: ALT>45 U/L, AST>40 U/L, GGT>65 U/L 男 / >42 女)
    alt_high = ifelse(is.na(LBXSATSI), NA, as.integer(LBXSATSI > 45)),
    ast_high = ifelse(is.na(LBXSASSI), NA, as.integer(LBXSASSI > 40)),
    ggt_high = ifelse(is.na(LBXSGTSI), NA,
                      as.integer(LBXSGTSI > ifelse(RIAGENDR == 1, 65, 42))),
    ggt_ln = ifelse(LBXSGTSI > 0, log(LBXSGTSI), NA),
    # 连续结局 (右偏, log 变换)
    alt_ln = ifelse(LBXSATSI > 0, log(LBXSATSI), NA),
    ast_ln = ifelse(LBXSASSI > 0, log(LBXSASSI), NA),
    # 剂量: 术后时长 (B_pool 已构造 yrs_since; 无手术史者为 NA)
    yrs_since = yrs_since,
    ever_surg = ch_surg
  ) |>
  filter(!is.na(ever_surg)) |>
  filter(complete.cases(RIDAGEYR, RIAGENDR, BMXBMI, RIDRETH3)) |>
  mutate(
    age_grp = ifelse(RIDAGEYR >= 60, "60plus", "lt60"),
    vtd_note = "P-cycle VID file name unresolved; A-part carries 25D exposure"
  )

cat("W2-B 分析样本 (暴露非缺失+协变量完整):", nrow(Bd), "\n")
cat("暴露组 (胆囊切除):", sum(Bd$ever_surg == 1), "\n")

des <- svydesign(ids = ~SDMVPSU, strata = ~SDMVSTRA, weights = ~wt_mec4,
                 nest = TRUE, data = Bd)

cov_m2 <- "RIDAGEYR + RIAGENDR + RIDRETH3"
cov_m3 <- paste(cov_m2, "+ BMXBMI + DMDEDUC2 + INDFMPIR + dm + smoke + drink")

run_or <- function(outcome, exposure, covs, desobj) {
  f <- as.formula(paste0(outcome, " ~ ", exposure, " + ", covs))
  fit <- svyglm(f, design = desobj, family = quasibinomial())
  co <- summary(fit)$coefficients
  ci <- suppressMessages(exp(confint(fit)))
  vars <- grep(exposure, rownames(co), value = TRUE, fixed = TRUE)
  data.frame(outcome = outcome, term = vars,
             OR = exp(co[vars, 1]), LCL = ci[vars, 1], UCL = ci[vars, 2],
             p = co[vars, 4], row.names = NULL)
}

# -----------------------------------------------------
# 1. 主分析: 胆囊切除 × 三项肝胆标志物升高 (M1/M3)
# -----------------------------------------------------
res <- list()
for (oc in c("alt_high", "ast_high", "ggt_high")) {
  res[[paste(oc, "M1")]] <- run_or(oc, "ever_surg", "1", des)
  res[[paste(oc, "M3")]] <- run_or(oc, "ever_surg", cov_m3, des)
}
res_df <- do.call(rbind, res)
res_df$set <- gsub("^\\d+\\.", "", row.names(res_df))
res_df$set <- gsub("\\.term.*$", "", res_df$set)
row.names(res_df) <- NULL
write.csv(res_df, file.path(OUT, "W2B_main_OR_table.csv"), row.names = FALSE)
cat("\n---- B 主分析 OR 表 (胆囊切除 vs 无) ----\n")
print(res_df[, c("set", "term", "OR", "LCL", "UCL", "p")], digits = 3)

# -----------------------------------------------------
# 2. 连续结局: 胆囊切除 × log-GGT / log-ALT (线性回归, beta 口径)
# -----------------------------------------------------
run_beta <- function(outcome, exposure, covs) {
  f <- as.formula(paste0(outcome, " ~ ", exposure, " + ", covs))
  fit <- svyglm(f, design = des)
  co <- summary(fit)$coefficients
  ci <- suppressMessages(confint(fit))
  data.frame(outcome = outcome,
             beta = co[exposure, 1], LCL = ci[exposure, 1], UCL = ci[exposure, 2],
             p = co[exposure, 4], row.names = NULL)
}
beta_df <- do.call(rbind, list(
  run_beta("ggt_ln", "ever_surg", cov_m3),
  run_beta("alt_ln", "ever_surg", cov_m3),
  run_beta("ast_ln", "ever_surg", cov_m3)
))
write.csv(beta_df, file.path(OUT, "W2B_continuous_beta.csv"), row.names = FALSE)
cat("\n---- 连续结局 (log 变换, beta) ----\n")
print(beta_df, digits = 3)

# -----------------------------------------------------
# 3. 剂量-时程: 术后时长 (仅手术史阳性者, n=1025) × GGT 升高
# -----------------------------------------------------
d_s <- Bd[Bd$ever_surg == 1 & !is.na(Bd$yrs_since), ]
des_s <- svydesign(ids = ~SDMVPSU, strata = ~SDMVSTRA, weights = ~wt_mec4,
                   nest = TRUE, data = d_s)
cat("\n剂量-时程子样本 (手术史阳性且术后时长可算):", nrow(d_s), "\n")
# 三分位分组
d_s$yr_q <- ntile(d_s$yrs_since, 3)
d_s$yr_qf <- factor(dyrq <- d_s$yr_q, labels = c("early tertile", "mid", "late"))
des_s <- update(des_s, yr_qf = d_s$yr_qf)
res_dose <- do.call(rbind, list(
  run_or("ggt_high", "yr_qf", "1", des_s),
  run_or("ggt_high", "yr_qf", cov_m3, des_s)
))
res_dose$set <- gsub("^\\d+\\.", "", row.names(res_dose))
res_dose$set <- gsub("\\.term.*$", "", res_dose$set)
row.names(res_dose) <- NULL
write.csv(res_dose, file.path(OUT, "W2B_dose_tertile.csv"), row.names = FALSE)
cat("\n---- 术后时长三分位 × GGT 升高 ----\n")
print(res_dose[, c("set", "term", "OR", "LCL", "UCL", "p")], digits = 3)

# 连续剂量 (预计算样条列, 参照 = 中位术后时长; 避免 update(design, data=) 兼容性问题)
sb <- ns(d_s$yrs_since, df = 4)
d_s$s1 <- sb[, 1]; d_s$s2 <- sb[, 2]; d_s$s3 <- sb[, 3]; d_s$s4 <- sb[, 4]
des_s2 <- svydesign(ids = ~SDMVPSU, strata = ~SDMVSTRA, weights = ~wt_mec4,
                    nest = TRUE, data = d_s)
fit_s <- svyglm(as.formula(paste("ggt_high ~ s1 + s2 + s3 + s4 +", cov_m3)),
                design = des_s2, family = quasibinomial())
beta_s <- coef(fit_s); V_s <- vcov(fit_s)
mode_val <- function(x) as.numeric(names(which.max(table(x, useNA = "no"))))
ref_cov <- c(RIDAGEYR = median(d_s$RIDAGEYR), RIAGENDR = mode_val(d_s$RIAGENDR),
             RIDRETH3 = mode_val(d_s$RIDRETH3), BMXBMI = median(d_s$BMXBMI, na.rm = TRUE),
             DMDEDUC2 = mode_val(d_s$DMDEDUC2), INDFMPIR = median(d_s$INDFMPIR, na.rm = TRUE),
             dm = 0, smoke = 0, drink = mode_val(d_s$drink))
mk <- function(b) c(1, b, ref_cov[["RIDAGEYR"]], ref_cov[["RIAGENDR"]], ref_cov[["RIDRETH3"]],
                    ref_cov[["BMXBMI"]], ref_cov[["DMDEDUC2"]], ref_cov[["INDFMPIR"]],
                    ref_cov[["dm"]], ref_cov[["smoke"]], ref_cov[["drink"]])
gx <- seq(0, min(60, max(d_s$yrs_since, na.rm = TRUE)), by = 1)
grid_b <- predict(ns(d_s$yrs_since, df = 4), gx)
ref_b  <- predict(ns(d_s$yrs_since, df = 4), median(d_s$yrs_since, na.rm = TRUE))
Xg <- t(apply(as.matrix(grid_b), 1, mk))
Xr <- mk(ref_b)
lp <- as.numeric(Xg %*% beta_s)
se <- sqrt(rowSums((Xg %*% V_s) * Xg))
lpr <- as.numeric(Xr %*% beta_s)
rd_b <- data.frame(yrs = gx, or = exp(lp - lpr),
                   lcl = exp(lp - 1.96 * se - lpr), ucl = exp(lp + 1.96 * se - lpr))
saveRDS(rd_b, file.path(OUT, "W2B_dose_rcs_points.RDS"))
p_b <- ggplot(rd_b, aes(x = yrs, y = or)) +
  geom_ribbon(aes(ymin = lcl, ymax = ucl), fill = "#0f766e", alpha = .15) +
  geom_line(color = "#0f766e", linewidth = .8) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "grey55") +
  labs(x = "Years since cholecystectomy", y = "OR of elevated GGT (ref: median years)",
       title = "Dose-response: time since cholecystectomy and elevated GGT",
       caption = "NHANES 2017-2020 pre-pandemic, weighted, surgery-positive subgroup") +
  theme_minimal(base_size = 12)
ggsave(file.path(OUT, "FigS1_dose_time_GGT.png"), p_b, width = 8, height = 5, dpi = 300)
cat("B 剂量-时程图已出: FigS1_dose_time_GGT.png\n")

# -----------------------------------------------------
# 4. 中介分析框架数据 (BMI 作为 25D-排便 / 手术-GGT 两条路径的中介候选)
#    mediation 包框架代码: 实际跑通放到 W3 (需要吃透 sensitivity 检验)
# -----------------------------------------------------
cat("\n中介分析框架: BMI 与腰围已在 B_pool/Bd 中, mediation 包代码按手册 W3 模板填写即可\n")

cat("\nW2-B 主分析完成, 全部结果写入 ", OUT, "\n", sep = "")
