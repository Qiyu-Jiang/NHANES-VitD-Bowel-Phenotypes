# =====================================================
# W3 中介分析脚本  (w3_mediation_d1.R)
# A 部分: 血清 25(OH)D → 排便表型, 双中介候选
#   M1 = BMI (生理路径: 脂肪螯合/代谢健康)
#   M2 = 膳食脂肪供能比 (膳食路径)
# 框架: Imai 因果中介分析 (mediation 包), bootstrap 1000
# ⚠ 已实测: mediate(boot=TRUE) 的模型公式必须以字面量直接写入 glm() 调用,
#   用预构造 formula 对象会在 bootstrap 重放时报 object not found
# 前置: w1_merge_d1.R 已跑通
# 产出: W2_results/W3_mediation_results.csv + W3_mediation_report.txt
# 运行: Rscript w3_mediation_d1.R  (在 DELIVERY 目录下)
# 生成: 2026-09-22 · AutoClaw
# =====================================================

library(dplyr)
library(foreign)
library(mediation)

ROOT <- "nhanes_data"
OUT  <- "W2_results"

A <- readRDS(file.path(ROOT, "A_pool.RDS"))

# ---- 膳食脂肪变量合入 (三期 DR1TOT, 单日可信回访) ----
fat_rows <- list()
for (s in c("D", "E", "F")) {
  sub <- c(D = "A_2005_2006", E = "A_2007_2008", F = "A_2009_2010")[[s]]
  df <- as.data.frame(foreign::read.xport(
    file.path(ROOT, sub, paste0("DR1TOT_", s, ".XPT"))))
  keep <- intersect(c("SEQN","DR1TKCAL","DR1TTFAT","DR1DRSTZ"), names(df))
  fat_rows[[s]] <- df[, keep]
}
fat <- bind_rows(fat_rows) |>
  filter(DR1DRSTZ == 1) |>
  group_by(SEQN) |> slice(1) |> ungroup() |>
  mutate(fat_pct = ifelse(DR1TKCAL > 0, DR1TTFAT * 9 / DR1TKCAL * 100, NA)) |>
  dplyr::select(SEQN, fat_pct)

# ---- 分析数据集 (与 W2 主分析对齐) ----
Ad <- A |>
  left_join(fat, by = "SEQN") |>
  filter(!is.na(stool_grp), !is.na(tot25d)) |>
  mutate(
    vtd_cut = ifelse(tot25d >= 50, 1, 0),
    constipation = as.integer(stool_grp == "constipation"),
    diarrhea     = as.integer(stool_grp == "diarrhea"),
    dm    = ifelse(is.na(dm), 0, dm),
    smoke = ifelse(is.na(smoke), 0, smoke)
  ) |>
  filter(complete.cases(RIDAGEYR, RIAGENDR, BMXBMI, RIDRETH1))

covs_full <- "RIDAGEYR + RIAGENDR + RIDRETH1 + DMDEDUC2 + INDFMPIR + dm + smoke + drink"
cat("W3 分析样本:", nrow(Ad), "\n")

# -----------------------------------------------------
# 四组中介: {diarrhea, constipation} × {BMI, fat_pct}
# 公式必须字面量写入 (实测验证, 见文件头注释)
# -----------------------------------------------------
run_med <- function(data, outcome, mediator, sims = 1000, seed = 42) {
  req <- c(outcome, mediator, "vtd_cut",
           covs_full_split[!covs_full_split %in% c("BMXBMI", "fat_pct")])
  d <- data[complete.cases(data[, req]), , drop = FALSE]
  d$w <- d$wt_mec6
  set.seed(seed)
  # ⚠ 实测: mediate(boot=TRUE) 的模型公式必须以字面量直接写入 glm() 调用,
  #   用预构造 formula 对象会在 bootstrap 重放时报 object not found
  if (mediator == "BMXBMI") {
    if (outcome == "diarrhea") {
      m_fit <- glm(BMXBMI ~ vtd_cut + RIDAGEYR + RIAGENDR + RIDRETH1 + DMDEDUC2 + INDFMPIR + dm + smoke + drink, data = d, weights = w)
      o_fit <- glm(diarrhea ~ vtd_cut + BMXBMI + RIDAGEYR + RIAGENDR + RIDRETH1 + DMDEDUC2 + INDFMPIR + dm + smoke + drink, data = d, family = binomial(), weights = w)
    } else {
      m_fit <- glm(BMXBMI ~ vtd_cut + RIDAGEYR + RIAGENDR + RIDRETH1 + DMDEDUC2 + INDFMPIR + dm + smoke + drink, data = d, weights = w)
      o_fit <- glm(constipation ~ vtd_cut + BMXBMI + RIDAGEYR + RIAGENDR + RIDRETH1 + DMDEDUC2 + INDFMPIR + dm + smoke + drink, data = d, family = binomial(), weights = w)
    }
  } else {
    if (outcome == "diarrhea") {
      m_fit <- glm(fat_pct ~ vtd_cut + RIDAGEYR + RIAGENDR + RIDRETH1 + DMDEDUC2 + INDFMPIR + dm + smoke + drink, data = d, weights = w)
      o_fit <- glm(diarrhea ~ vtd_cut + fat_pct + RIDAGEYR + RIAGENDR + RIDRETH1 + DMDEDUC2 + INDFMPIR + dm + smoke + drink, data = d, family = binomial(), weights = w)
    } else {
      m_fit <- glm(fat_pct ~ vtd_cut + RIDAGEYR + RIAGENDR + RIDRETH1 + DMDEDUC2 + INDFMPIR + dm + smoke + drink, data = d, weights = w)
      o_fit <- glm(constipation ~ vtd_cut + fat_pct + RIDAGEYR + RIAGENDR + RIDRETH1 + DMDEDUC2 + INDFMPIR + dm + smoke + drink, data = d, family = binomial(), weights = w)
    }
  }
  med <- mediate(m_fit, o_fit, treat = "vtd_cut", mediator = mediator,
                 boot = TRUE, sims = sims)
  s  <- summary(med)
  # mediation 4.5.1: 置信区间在 summary 对象的 *.ci 字段 (2 列矩阵 lower/upper);
  # confint(med) 不可用 (mediate 对象无 vcov)
  g1 <- function(x) if (is.null(x) || length(x) == 0) NA_real_ else as.numeric(x)[1]
  ci2 <- function(m) if (is.null(m) || length(m) == 0) c(NA_real_, NA_real_) else as.numeric(m)[1:2]
  acme_ci <- ci2(s$d0.ci); ade_ci <- ci2(s$z0.ci); prp_ci <- ci2(s$n0.ci)
  data.frame(
    outcome = outcome, mediator = mediator, n = nrow(d),
    ACME = g1(s$d0), ACME_lcl = g1(acme_ci[1]), ACME_ucl = g1(acme_ci[2]), ACME_p = g1(s$d0.p),
    ADE  = g1(s$z0), ADE_lcl  = g1(ade_ci[1]),  ADE_ucl  = g1(ade_ci[2]),  ADE_p  = g1(s$z0.p),
    Total = g1(s$tau.coef), Total_p = g1(s$tau.p),
    Prop_mediated = g1(s$n0), Prop_lcl = g1(prp_ci[1]), Prop_ucl = g1(prp_ci[2])
  )
}

covs_full_split <- c("vtd_cut","RIDAGEYR","RIAGENDR","RIDRETH1","DMDEDUC2",
                     "INDFMPIR","dm","smoke","drink","BMXBMI","fat_pct")

results <- list()
results[["diarrhea_BMI"]]      <- run_med(Ad, "diarrhea",     "BMXBMI",  sims = 1000)
cat("1/4 diarrhea  x BMI    done\n")
results[["diarrhea_FATPCT"]]   <- run_med(Ad, "diarrhea",     "fat_pct", sims = 1000)
cat("2/4 diarrhea  x fat%   done\n")
results[["constipation_BMI"]]  <- run_med(Ad, "constipation", "BMXBMI",  sims = 1000)
cat("3/4 constipation x BMI done\n")
results[["constipation_FATPCT"]] <- run_med(Ad, "constipation", "fat_pct", sims = 1000)
cat("4/4 constipation x fat% done\n")

res_df <- do.call(rbind, results)
write.csv(res_df, file.path(OUT, "W3_mediation_results.csv"), row.names = FALSE)
cat("\n---- W3 中介分析结果汇总 ----\n")
print(res_df[, c("outcome","mediator","n","ACME","ACME_lcl","ACME_ucl","ACME_p",
                 "ADE","ADE_p","Prop_mediated")], digits = 3)

cat("\n解读提示 (Discussion 写法):\n")
cat("- ACME 95%CI 不含 0 => 间接效应成立; Prop_mediated = 中介占比\n")
cat("- BMI 中介成立 => '部分关联经 adiposity 通路' 与 W2 亚组 BMI 梯度互证\n")
cat("- 脂肪供能比中介成立 => '膳食结构是关联的候选载体'\n")
cat("- 横断面 + mediation 的因果措辞在局限段声明 (手册 W3 预案)\n")

cat("\nW3 完成, 结果写入 ", file.path(OUT, "W3_mediation_results.csv"), "\n", sep = "")
