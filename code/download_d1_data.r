# =====================================================
# D1 数据一键下载与校验脚本  (download_d1_data.R)
# 项目: NHANES 胆囊切除史 x 排便表型 (目标期刊 Scientific Reports)
# 用法: 在 R/RStudio 中 source("download_d1_data.R")
#       工作目录即数据落盘根目录 (自动创建 nhanes_data/)
# 依赖: 仅需 base R 的 foreign 包 (R 自带); 若用 nhanesA 则 install.packages("nhanesA")
# 产出: nhanes_data/ 下按子文件夹存放全部 .XPT + 终端校验报告
# 生成: 2026-09-22 · AutoClaw · 全部 URL 为 CDC 官方现行路径
# =====================================================

ROOT <- "nhanes_data"

# -----------------------------------------------------
# 1. 文件清单 (分组 / 起始年 / 文件名)
#    官方路径模式: https://wwwn.cdc.gov/Nchs/Data/Nhanes/Public/{起始年}/DataFiles/{文件名}.XPT
# -----------------------------------------------------
GROUPS <- list(
  # -- B 部分·主选: 2017-2020 疫情前合并周期 (官方已合并, 含 MCQ560/MCQ570) --
  #   实测命名规律: 医学史为 P_MCQ, 人口学为 P_DEMO (两种风格并存, 均 HEAD 实测可达)
  #   注: P 周期维生素 D 文件名多轮探针未定, 25(OH)D 暴露属 A 部分 (VID_E/F 已实测),
  #       不阻塞; 若需 P 周期 VID, 到组件页人工确认后加入 files
  B_main = list(start = 2017, sub = "B_pmcq_2017_2020",
    files = c("P_DEMO","P_MCQ","P_BMX","P_BIOPRO","P_CBC","P_DPQ","P_SLQ",
              "P_HSQ","P_DIQ","P_BPQ","P_ALQ","P_SMQ","P_DR1TOT")),
  # -- B 部分·备选: 2017-2018 单周期 (J 后缀, 全部文件已实测存在) --
  B_backup = list(start = 2017, sub = "B_backup_2017_2018",
    files = c("DEMO_J","MCQ_J","BMX_J","BIOPRO_J","CBC_J","VID_J",
              "DPQ_J","SLQ_J","HSQ_J","DIQ_J","BPQ_J","ALQ_J","SMQ_J","DR1TOT_J")),
  # -- A 部分·结局主体: 2005-2010 三期 (Bristol 便型/便次) --
  A_2005 = list(start = 2005, sub = "A_2005_2006",
    files = c("DEMO_D","BHQ_D","BMX_D","VID_D","MCQ_D","DPQ_D","SLQ_D","HSQ_D",
              "DIQ_D","BPQ_D","ALQ_D","SMQ_D","DR1TOT_D")),
  A_2007 = list(start = 2007, sub = "A_2007_2008",
    files = c("DEMO_E","BHQ_E","BMX_E","VID_E","MCQ_E","DPQ_E","SLQ_E","HSQ_E",
              "DIQ_E","BPQ_E","ALQ_E","SMQ_E","DR1TOT_E")),
  A_2009 = list(start = 2009, sub = "A_2009_2010",
    files = c("DEMO_F","BHQ_F","BMX_F","VID_F","MCQ_F","DPQ_F","SLQ_F","HSQ_F",
              "DIQ_F","BPQ_F","ALQ_F","SMQ_F","DR1TOT_F"))
)

# -----------------------------------------------------
# 2. 下载
# -----------------------------------------------------
download_one <- function(file, start_year, sub) {
  dir.create(file.path(ROOT, sub), showWarnings = FALSE, recursive = TRUE)
  dest <- file.path(ROOT, sub, paste0(file, ".XPT"))
  if (file.exists(dest) && file.size(dest) > 10000) {
    message("SKIP (已存在): ", dest); return(TRUE)
  }
  url <- sprintf("https://wwwn.cdc.gov/Nchs/Data/Nhanes/Public/%d/DataFiles/%s.XPT",
                 start_year, file)
  ok <- tryCatch({
    download.file(url, dest, mode = "wb", quiet = TRUE)
    file.exists(dest) && file.size(dest) > 10000
  }, error = function(e) FALSE, warning = function(w) FALSE)
  if (!ok && file.exists(dest)) file.remove(dest)
  message(ifelse(ok, "OK  ", "FAIL"), "  ", file, "  <-  ", url)
  ok
}

results <- data.frame(group = character(), file = character(), ok = logical(),
                      stringsAsFactors = FALSE)
for (nm in names(GROUPS)) {
  g <- GROUPS[[nm]]
  for (f in g$files) {
    ok <- download_one(f, g$start, g$sub)
    results <- rbind(results, data.frame(group = nm, file = f, ok = ok))
  }
}

# -----------------------------------------------------
# 3. 校验 (对照 2026-09-22 实测数字)
# -----------------------------------------------------
expect_rows <- c(
  "B_backup_2017_2018/MCQ_J.XPT"      = 8897,
  "B_backup_2017_2018/BIOPRO_J.XPT"   = 6401,
  "B_backup_2017_2018/CBC_J.XPT"      = 8366,
  "A_2005_2006/BHQ_D.XPT"             = 4773,
  "A_2007_2008/BHQ_E.XPT"             = 5707,
  "A_2009_2010/BHQ_F.XPT"             = 6059,
  "A_2009_2010/DR1TOT_F.XPT"          = 10253
)
# 值校验: 文件 -> c(列名, 取值, 期望计数)
expect_vals <- list(
  "B_backup_2017_2018/MCQ_J.XPT" = list(c("MCQ560", 1, 641), c("MCQ550", 1, 623)),
  "A_2009_2010/DR1TOT_F.XPT"     = list(c("DR1DRSTZ", 1, 9623))
)

cat("\n================ 校验报告 ================\n")
all_ok <- TRUE
for (rel in names(expect_rows)) {
  p <- file.path(ROOT, rel)
  if (!file.exists(p)) { cat(sprintf("MISSING  %s\n", rel)); all_ok <- FALSE; next }
  df <- as.data.frame(foreign::read.xport(p))
  flag <- ifelse(nrow(df) == expect_rows[[rel]], "OK  ", "DIFF ")
  if (flag != "OK  ") all_ok <- FALSE
  cat(sprintf("%s %-36s rows=%d (期望 %d)\n", flag, rel, nrow(df), expect_rows[[rel]]))
  for (ck in expect_vals[[rel]]) {
    col <- ck[1]; val <- as.numeric(ck[2]); expn <- as.numeric(ck[3])
    if (col %in% names(df)) {
      n <- sum(df[[col]] == val, na.rm = TRUE)
      f2 <- ifelse(n == expn, "OK  ", "DIFF ")
      if (f2 != "OK  ") all_ok <- FALSE
      cat(sprintf("      %s %s==%s: %d (期望 %d)\n", f2, col, val, n, expn))
    } else {
      all_ok <- FALSE
      cat(sprintf("      DIFF 列 %s 不存在于 %s\n", col, rel))
    }
  }
}
# BHQ_F 特检: 便秘/腹泻问项应仅存在于 F 期
p_f <- file.path(ROOT, "A_2009_2010/BHQ_F.XPT")
if (file.exists(p_f)) {
  bhq_f <- as.data.frame(foreign::read.xport(p_f))
  has <- c("BHQ080","BHQ090","BHQ100") %in% names(bhq_f)
  cat(sprintf("OK   BHQ_F 含 BHQ080/090/100: %s\n", paste(has, collapse="/")))
}

# 维生素 D 变量名速查 (各期命名有差异, 打印含 VD 的列名辅助定位)
cat("\n---- 维生素 D 变量速查 (各期 25(OH)D 列名) ----\n")
vid_files <- c("B_backup_2017_2018/VID_J.XPT",
               "A_2005_2006/VID_D.XPT","A_2007_2008/VID_E.XPT","A_2009_2010/VID_F.XPT")
for (rel in vid_files) {
  p <- file.path(ROOT, rel)
  if (file.exists(p)) {
    df <- as.data.frame(foreign::read.xport(p))
    cat(rel, "->", paste(grep("VD", names(df), value = TRUE), collapse = ", "), "\n")
  }
}

# -----------------------------------------------------
# 4. 失败文件的人工兜底指引
# -----------------------------------------------------
failed <- results[!results$ok, ]
if (nrow(failed) > 0) {
  cat("\n!!!! 有文件下载失败, 请人工核对文件名后手动下载 !!!!\n")
  cat("1) 打开对应周期组件索引页确认准确文件名:\n")
  cat("   2017-2020: https://wwwn.cdc.gov/nchs/nhanes/continuousnhanes/default.aspx?BeginYear=2017\n")
  cat("   2005-2006: https://wwwn.cdc.gov/nchs/nhanes/continuousnhanes/default.aspx?BeginYear=2005\n")
  cat("   (其余周期同理, 修改 BeginYear)\n")
  cat("2) 或用官网变量搜索: https://wwwn.cdc.gov/nchs/nhanes/search/default.aspx\n")
  cat("3) 失败最常见原因: P 后缀文件个别命名不同 (如 P_MCQ 与 XXX_P 两种风格并存)。\n")
  cat("失败清单:\n"); print(failed)
} else {
  cat("\n全部文件下载成功。\n")
}
cat("\n总状态:", ifelse(all_ok, "校验全部通过, 可以开始分析", "存在 DIFF/MISSING, 见上方报告"), "\n")
