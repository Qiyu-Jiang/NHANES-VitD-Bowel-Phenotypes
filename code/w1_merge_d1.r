# =====================================================
# W1 合并与清洗脚本  (w1_merge_d1.R)
# 项目: D1 双周期设计 (A: 2005-2010 VitD x Bristol; B: 2017-2020 胆囊切除)
# 前置: nhanes_data/ 已由 download_d1_data.R 下载就绪
# 产出: nhanes_data/A_pool.RDS, B_pool.RDS, W1_verify_summary.csv + 终端校验报告
# 运行: Rscript w1_merge_d1.R  (在 DELIVERY 目录下)
# 生成: 2026-09-22 · AutoClaw · 全部编码以官方 codebook + XPT 实测为准
# =====================================================

library(foreign)
library(dplyr)

ROOT <- "nhanes_data"
RD <- function(rel) as.data.frame(foreign::read.xport(file.path(ROOT, rel)))

# -----------------------------------------------------
# 0. 工具: 动态选列 (P 周期与旧周期变量名存在差异, 统一在此处理)
# -----------------------------------------------------
pick1 <- function(df, cands, newname) {
  col <- intersect(cands, names(df))[1]
  if (is.na(col)) {
    df[[newname]] <- NA_real_
  } else {
    df[[newname]] <- df[[col]]
  }
  df[, c("SEQN", newname), drop = FALSE]
}

recode_yn <- function(x) ifelse(x == 1, 1, ifelse(x %in% c(2, 3), 0, NA))  # 1=Yes, 2/3=No/Borderline

# -----------------------------------------------------
# 1. A 部分: 2005-2010 三期合并 (Bristol 结局主体)
# -----------------------------------------------------
a_dirs <- c(D = "A_2005_2006", E = "A_2007_2008", F = "A_2009_2010")
a_lab  <- c(D = "2005-2006",  E = "2007-2008",  F = "2009-2010")

build_a <- function(sfx) {
  sub <- a_dirs[[sfx]]
  demo <- RD(file.path(sub, paste0("DEMO_", sfx, ".XPT")))
  demo <- demo[, intersect(c("SEQN","RIDAGEYR","RIAGENDR","RIDRETH1","DMDEDUC2",
                             "INDFMPIR","WTMEC2YR","SDMVPSU","SDMVSTRA"), names(demo))]
  bhq  <- RD(file.path(sub, paste0("BHQ_",  sfx, ".XPT")))
  bhq  <- bhq[, intersect(c("SEQN","BHQ010","BHQ020","BHQ030","BHQ040",
                            "BHD050","BHQ060","BHQ070","BHQ080","BHQ090","BHQ100"), names(bhq))]
  bmx  <- RD(file.path(sub, paste0("BMX_",  sfx, ".XPT")))
  bmx  <- bmx[, intersect(c("SEQN","BMXBMI","BMXWAIST"), names(bmx))]
  vid  <- RD(file.path(sub, paste0("VID_",  sfx, ".XPT")))
  vid  <- pick1(vid, c("LBXVIDMS","LBDVIDMS"), "tot25d_raw")   # 总 25(OH)D 原始值
  # NCHS 官方校准: 2005-2006 RIA → LC-MS/MS equivalent
  # (Vitamin D Analytical Note, Oct 2015: LC-MS/MS eq = 8.36753 + 0.97012 × RIA original)
  if (sfx == "D") {
    vid$tot25d_raw <- 8.36753 + 0.97012 * vid$tot25d_raw
  }
  vid$tot25d <- vid$tot25d_raw
  vid <- vid[, c("SEQN", "tot25d")]
  dpq  <- RD(file.path(sub, paste0("DPQ_",  sfx, ".XPT")))
  dpq  <- dpq[, intersect(c("SEQN", paste0("DPQ0", 10:90), "DPQ100"), names(dpq))]
  diq  <- pick1(RD(file.path(sub, paste0("DIQ_", sfx, ".XPT"))), "DIQ010", "DIQ010")
  smq  <- pick1(RD(file.path(sub, paste0("SMQ_", sfx, ".XPT"))), "SMQ020", "SMQ020")
  alq  <- RD(file.path(sub, paste0("ALQ_",  sfx, ".XPT")))
  alq  <- pick1(alq, c("ALQ101","ALQ111"), "alc_ever")

  out <- demo |>
    left_join(bhq, by = "SEQN") |> left_join(bmx, by = "SEQN") |>
    left_join(vid, by = "SEQN") |> left_join(dpq, by = "SEQN") |>
    left_join(diq, by = "SEQN") |> left_join(smq, by = "SEQN") |>
    left_join(alq, by = "SEQN")
  out$cycle <- a_lab[[sfx]]
  out
}

A_pool <- bind_rows(lapply(names(a_dirs), build_a))

DPQ_COLS <- intersect(c(paste0("DPQ0", 10:90), "DPQ100"), names(A_pool))
A_pool$DPQ_ok_flag <- rowSums(is.na(A_pool[, DPQ_COLS, drop = FALSE])) <= 2

A_pool <- A_pool |>
  filter(RIDAGEYR >= 20) |>
  mutate(
    # 权重: BHQ 为 MEC 访谈组件, 且分析使用 MEC 检查项 -> MEC 权重, 3 期合并除以 3
    wt_mec6 = WTMEC2YR / 3,
    # Bristol 便型分组 (官方 codebook: 1-7 型, 99=拒答)
    stool_type = ifelse(BHQ060 %in% c(99, 77777, 99999), NA, BHQ060),
    stool_grp  = case_when(
      stool_type %in% 1:2 ~ "constipation",
      stool_type %in% 3:5 ~ "normal",
      stool_type %in% 6:7 ~ "diarrhea",
      TRUE ~ NA_character_),
    # 排便频率: BHD050 = 过去 7 天排便总次数 (实测取值 1-63, 999=拒答)
    bm_pw = ifelse(BHD050 >= 999, NA, BHD050),
    bm_grp = case_when(
      bm_pw < 7  ~ "low(<1/d)",
      bm_pw <= 14 ~ "normal(1-2/d)",
      bm_pw > 14  ~ "high(>2/d)",
      TRUE ~ NA_character_),
    # FISI 四项 (气/黏液/液体/固体, 1-6 级): 全部应答时求和; 单独取固体失禁指示
    fisi_any_na = is.na(BHQ010) | is.na(BHQ020) | is.na(BHQ030) | is.na(BHQ040),
    fisi_sum = ifelse(fisi_any_na, NA_real_,
                      rowSums(cbind(BHQ010, BHQ020, BHQ030, BHQ040), na.rm = TRUE)),
    fi_solid = ifelse(is.na(BHQ040) | BHQ040 %in% c(7, 9), NA,
                      ifelse(BHQ040 >= 1 & BHQ040 <= 5, 1, 0)),
    # 便急 (BHQ070): 1=Yes
    urgency = recode_yn(BHQ070),
    # 协变量
    dm     = recode_yn(DIQ010),
    smoke  = recode_yn(SMQ020),
    drink  = recode_yn(alc_ever),
    educ   = ifelse(DMDEDUC2 %in% 1:5, DMDEDUC2, NA_real_),
    # PHQ-9: 7/9 -> NA, 缺失 <=2 题时按已答求和 (先在 mutate 外算好列名向量)
    dpq_ok = DPQ_ok_flag,
    phq9 = ifelse(dpq_ok, rowSums(across(all_of(DPQ_COLS)), na.rm = TRUE), NA_real_)
  ) |>
  select(-any_of("phq_any_na"))

# 剔除临时标记列
A_pool$DPQ_ok_flag <- NULL

# -----------------------------------------------------
# 2. B 部分: 2017-2020 疫情前合并周期 (暴露核心)
# -----------------------------------------------------
P <- "B_pmcq_2017_2020"
demoP <- RD(file.path(P, "P_DEMO.XPT"))
demoP <- demoP[, intersect(c("SEQN","RIDAGEYR","RIAGENDR","RIDRETH3","DMDEDUC2",
                             "INDFMPIR","WTMECPRP","SDMVPSU","SDMVSTRA"), names(demoP))]
mcqP  <- RD(file.path(P, "P_MCQ.XPT"))
mcqP  <- mcqP[, intersect(c("SEQN","MCQ550","MCQ560","MCQ570"), names(mcqP))]
bmxP  <- RD(file.path(P, "P_BMX.XPT"))
bmxP  <- bmxP[, intersect(c("SEQN","BMXBMI","BMXWAIST"), names(bmxP))]
diqP  <- pick1(RD(file.path(P, "P_DIQ.XPT")), "DIQ010", "DIQ010")
smqP  <- pick1(RD(file.path(P, "P_SMQ.XPT")), "SMQ020", "SMQ020")
alqP  <- RD(file.path(P, "P_ALQ.XPT"))
alqP  <- pick1(alqP, c("ALQ101","ALQ111"), "alc_ever")
dpqP  <- RD(file.path(P, "P_DPQ.XPT"))
dpqP  <- dpqP[, intersect(c("SEQN", paste0("DPQ0", 10:90), "DPQ100"), names(dpqP))]

B_pool <- demoP |>
  left_join(mcqP, by = "SEQN") |> left_join(bmxP, by = "SEQN") |>
  left_join(diqP, by = "SEQN") |> left_join(smqP, by = "SEQN") |>
  left_join(alqP, by = "SEQN") |> left_join(dpqP, by = "SEQN") |>
  filter(RIDAGEYR >= 20) |>
  mutate(
    wt_mec4 = WTMECPRP,          # P 周期 4 年合并权重列名为 WTMECPRP (实测确认)
    ch_surg  = ifelse(MCQ560 == 1, 1, ifelse(MCQ560 == 2, 0, NA)),
    gallstone = ifelse(MCQ550 == 1, 1, ifelse(MCQ550 == 2, 0, NA)),
    surg_age = ifelse(MCQ570 >= 70000, NA, MCQ570),
    yrs_since = RIDAGEYR - surg_age,
    dm = recode_yn(DIQ010),
    smoke = recode_yn(SMQ020),
    drink = recode_yn(alc_ever),
    educ = ifelse(DMDEDUC2 %in% 1:5, DMDEDUC2, NA_real_)
  )

# -----------------------------------------------------
# 3. 校验报告
# -----------------------------------------------------
wmean <- function(df, var, w = "wt") {
  d2 <- df[!is.na(df[[var]]) & !is.na(df[[w]]), ]
  if (nrow(d2) == 0) return(NA_real_)
  v <- as.numeric(d2[[var]]); ww <- d2[[w]]
  sum(ww * v) / sum(ww)
}
cat("\n================ W1 校验报告 ================\n")
cat("\n---- A 部分 (2005-2010 合并) ----\n")
cat("总行数(20+):", nrow(A_pool), "\n")
print(table(A_pool$cycle))
cat("stool_grp 非缺失:", sum(!is.na(A_pool$stool_grp)), "\n")
A_pool$w1 <- A_pool$wt_mec6
tab <- table(A_pool$stool_grp)
print(tab)
for (g in c("constipation", "normal", "diarrhea")) {
  A_pool$y_g <- ifelse(A_pool$stool_grp == g, 1, 0)
  cat(sprintf("加权患病率 %s: %.2f%%\n", g, 100 * wmean(A_pool, "y_g", "w1")))
}
cat("tot25d 非缺失:", sum(!is.na(A_pool$tot25d)), "| 加权均值:",
    sprintf("%.1f nmol/L", wmean(A_pool, "tot25d", "w1")), "\n")
cat("BM 分组分布:\n"); print(table(A_pool$bm_grp, useNA = "ifany"))
cat("FISI 可算:", sum(!is.na(A_pool$fisi_sum)), "| 固体失禁阳性:",
    sum(A_pool$fi_solid == 1, na.rm = TRUE), "\n")

cat("\n---- B 部分 (2017-2020 P 合并) ----\n")
cat("总行数(20+):", nrow(B_pool), "\n")
cat("胆囊手术 MCQ560: yes=", sum(B_pool$ch_surg == 1, na.rm = TRUE),
    " no=", sum(B_pool$ch_surg == 0, na.rm = TRUE),
    " NA=", sum(is.na(B_pool$ch_surg)), "\n")
B_pool$w1 <- B_pool$wt_mec4
B_pool$y_ch <- B_pool$ch_surg
cat(sprintf("加权胆囊手术史比例: %.2f%%\n", 100 * wmean(B_pool, "y_ch", "w1")))
cat("手术年龄 MCQ570 非缺失:", sum(!is.na(B_pool$surg_age)),
    "| 术后时长中位数:", median(B_pool$yrs_since, na.rm = TRUE), "年\n")
print(summary(B_pool$yrs_since))

# -----------------------------------------------------
# 4. 落盘
# -----------------------------------------------------
saveRDS(A_pool, file.path(ROOT, "A_pool.RDS"))
saveRDS(B_pool, file.path(ROOT, "B_pool.RDS"))
summ <- data.frame(
  item = c("A_rows_20plus", "A_stool_nonmiss", "A_wprev_constipation",
           "A_wprev_diarrhea", "A_tot25d_nonmiss", "B_rows_20plus",
           "B_ch_surg_yes", "B_wprev_surgery"),
  value = c(nrow(A_pool), sum(!is.na(A_pool$stool_grp)),
            round(100 * {A_pool$y_g <- ifelse(A_pool$stool_grp=="constipation",1,0); wmean(A_pool,"y_g","w1")}, 2),
            round(100 * {A_pool$y_g <- ifelse(A_pool$stool_grp=="diarrhea",1,0); wmean(A_pool,"y_g","w1")}, 2),
            sum(!is.na(A_pool$tot25d)), nrow(B_pool),
            sum(B_pool$ch_surg == 1, na.rm = TRUE),
            round(100 * wmean(B_pool, "y_ch", "w1"), 2))
)
write.csv(summ, file.path(ROOT, "W1_verify_summary.csv"), row.names = FALSE)
cat("\n已保存: A_pool.RDS / B_pool.RDS / W1_verify_summary.csv\n")
cat("W1 完成: 数据就绪, 进入 W2 主分析\n")
