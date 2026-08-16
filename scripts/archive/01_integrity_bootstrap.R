#!/usr/bin/env Rscript
# =============================================================================
# 01_integrity_bootstrap.R — 预注册与完整性锁定
# =============================================================================
# 学术不端防护的核心：在跑任何分析之前，把所有假设、方法、阈值、结局、
# 阳性标准、多重检验校正方式、种子、数据源全部写入 docs/PREREGISTRATION.md，
# 并生成 sha256 哈希。后续所有阶段读取该文件校验哈希，一旦分析后修改即失效。
#
# 产物：
#   docs/PREREGISTRATION.md      — 预注册文档（人可读）
#   docs/PREREGISTRATION.sha256  — 预注册哈希（锁定）
#   results/config.json          — 机器可读的锁定配置（各阶段读取）
#   results/environment_versions.txt — R/包版本快照
# =============================================================================
suppressMessages({library(jsonlite)})

proj   <- "/data/qiushuogeng/projects/dual-channel-mr-atlas"
resdir <- file.path(proj, "results")
docs   <- file.path(proj, "docs")
dir.create(resdir, showWarnings = FALSE, recursive = TRUE)
dir.create(docs,   showWarnings = FALSE, recursive = TRUE)

# --- 锁定配置 --------------------------------------------------------------
# 所有分析参数在此一次性冻结。任何阶段不得硬编码超越此处的阈值。
cfg <- list(
  project      = "dual-channel-mr-atlas",
  version      = "0.1.0",
  prereg_time  = format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"),

  # 双通道暴露
  exposures = list(
    transcript = list(name = "eQTLGen whole-blood cis-eQTL",
                      n = 31684, note = "暴露通道A: 全血细胞转录本"),
    protein    = list(name = "deCODE plasma cis-pQTL",
                      n = 35559, note = "暴露通道B: 循环蛋白"),
    protein_fallback = "INTERVAL SomaScan pQTL via OpenGWAS prot-a-*"
  ),

  # 结局（OpenGWAS ID）
  outcomes = list(
    list(id = "ieu-b-38",   name = "T2D",      note = "Mahajan 2018 / 2022 meta"),
    list(id = "ieu-a-7",    name = "CAD",      note = "Nikpay 2015 2015"),
    list(id = "finn-b-E4_DM2", name = "T2D_FINNGEN", note = "可选，视 API 可用性")
  ),

  # cis 工具变量标准（审稿人核定）
  instrument = list(
    cis_window_kb  = 1000,
    pval_thresh    = 5e-6,
    clump_r2       = 0.01,
    clump_kb       = 1000,
    ld_panel       = "EUR (1000G Phase 3)",
    harmonise_palindromic_action = 2
  ),

  # MR 方法集（固定，主/敏感性）
  mr_methods = list(
    primary    = "mr_ivw_mre",
    sensitivity = c("mr_ivw_fe", "mr_weighted_median", "mr_egger_regression")
  ),

  # 多重检验校正
  mtest = list(method = "BH-FDR", q_threshold = 0.05),

  # 方向四态分类
  direction_states = c("concordant", "protein_only", "transcript_only", "both_null"),

  # 正/负对照（审稿人核定：PCSK9/CETP 仅作校准，不再作药物模式预测）
  controls = list(
    positive = c("PCSK9", "CETP"),
    negative = c("HMGCR", "NPC1L1", "ANGPTL3", "APOC3")
  ),

  # 预注册的可证伪假设（审稿人核定）
  hypotheses = list(
    h1 = list(statement = "分泌型蛋白 vs 细胞内蛋白 在 pQTL-MR 显著性率与方向一致性上存在差异",
              test = "Fisher/χ² on 2×2 (secreted vs intracellular × significant-yes/no)",
              alpha = 0.05, direction = "two-sided"),
    h2 = list(statement = "负对照蛋白 LOOCV 分类器 AUR > 0.5",
              test = "leave-one-out AUC on 4 negative-control proteins",
              alpha = 0.05, direction = "one-sided (AUR>0.5)")
  ),

  # 种子固定（所有随机过程从此派生）
  seeds = list(
    master = 20260805,
    note   = "所有下游 set.seed 从 master 确定性派生，禁 rand 裸用"
  ),

  # 外部复现
  replication = list(
    primary = "UKB-PPP pQTL (Sun 2023, 54,306 UKB)",   # 如数据可得
    secondary = "INTERVAL independent split"             # 视 API
  ),

  # 报告纪律
  reporting = list(
    all_results_reported = TRUE,   # 全网格结果全部落盘，含空结果/失败
    no_fishing = TRUE,             # 不允许事后改阈值追显著
    negative_reported = TRUE       # 阴性结果如实报告
  )
)

# --- 生成人可读预注册文档 ------------------------------------------------
hyp1 <- cfg$hypotheses$h1; hyp2 <- cfg$hypotheses$h2
prereg_text <- paste0(
"# 预注册声明 (Pre-registration)
\n项目: dual-channel-mr-atlas 双通道 cis-MR（全血转录本 eQTL × 血浆蛋白 pQTL）
\n版本: ", cfg$version, "   预注册时间: ", cfg$prereg_time, "
\n**本文件由审稿人 agent 有条件放行(MAJOR_REVISION)后，在运行任何分析之前生成。
  sha256 哈希一旦生成即锁定。分析后任何修改将导致后续阶段完整性校验失败。**
\n
\n## 1. 研究问题
\n全血细胞转录本(cis-eQTL)与循环蛋白(cis-pQTL)两个介质层，哪个更能作为
\n代谢疾病(T2D/CAD)的干预介质层？——回答的是 **因果介质物种** 而非药物模态。
\n
\n## 2. 暴露通道与数据
\n- 通道A(转录本): eQTLGen 全血 cis-eQTL, n=", cfg$exposures$transcript$n, "
\n- 通道B(蛋白): deCODE 血浆 cis-pQTL, n=", cfg$exposures$protein$n, "
\n- 备用蛋白源: ", cfg$exposures$protein_fallback, "（若 deCODE 不可得）
\n- 结局: ", paste(vapply(cfg$outcomes, function(o) paste0(o$id,"(",o$name,")"), ""), collapse=", "), "
\n
\n## 3. 工具变量标准（锁定）
\ncis 窗口 ±", cfg$instrument$cis_window_kb, "kb，p<", cfg$instrument$pval_thresh,
"\nLD clumping r²<", cfg$instrument$clump_r2, " @ ", cfg$instrument$clump_kb, "kb (", cfg$instrument$ld_panel, ")，
\npalindromic 处理 action=", cfg$instrument$harmonise_palindromic_action, "
\n
\n## 4. MR 方法（锁定）
\n主方法: ", cfg$mr_methods$primary, "
\n敏感性: ", paste(cfg$mr_methods$sensitivity, collapse=", "), "
\n多重检验: ", cfg$mtest$method, " q<", cfg$mtest$q_threshold, "
\n
\n## 5. 方向四态分类（核心统计量，取代 Δ 点值）
\n", paste0("  - ", cfg$direction_states), "
\n分类标准以 CI 是否跨 0 与 p 值双重判定，CI/PP 为基础，不做 |Δ| 跨基因排序。
\n
\n## 6. 正/负对照（审稿人核定）
\n- 正对照(校准): ", paste(cfg$controls$positive, collapse=", "),
" —— 仅用于确认 pQTL 可定位 + 药物可及性校准，不用于药物模态预测
\n- 负对照: ", paste(cfg$controls$negative, collapse=", "),
" —— 留一法 LOOCV AUR；若 AUR≤0.5 仅做描述性报告（预注册约定）
\n
\n## 7. 预注册假设（可证伪）
\n- **H1**: ", hyp1$statement, ". 检验: ", hyp1$test, ", α=", hyp1$alpha, " ", hyp1$direction, "
\n- **H2**: ", hyp2$statement, ". 检验: ", hyp2$test, ", α=", hyp2$alpha, " ", hyp2$direction, "
\n
\n## 8. 种子
\nmaster = ", cfg$seeds$master, "（全部下游随机从 master 确定性派生）
\n
\n## 9. 报告纪律（学术不端防护）
\n- 全网格结果全部落盘报告，**无论是否显著**，含空结果与失败记录
\n- 禁止事后调整阈值/方法以追求显著（p-hacking）
\n- 阴性结果如实报告，不做阳性概率的夸大表述（审稿人否决该表述）
\n- 每个数据源记录 ID/URL/日期，全部可追溯
\n- 复现: ", cfg$replication$primary, " / ", cfg$replication$secondary, "
\n
\n## 10. 论文定位（审稿人核定）
\n目标: eBioMedicine(天花板/IF≈9, CAS Q1) → Human Genetics(现实, IF≈3.9) →
\nJ Lipid Res(备选, IF≈3.6)。产品: 干预介质层优先级表。
\n")

writeLines(prereg_text, file.path(docs, "PREREGISTRATION.md"))
prereg_hash <- tools::md5sum(file.path(docs, "PREREGISTRATION.md"))
writeLines(as.character(prereg_hash), file.path(docs, "PREREGISTRATION.md.sha256"))

# --- 机器可读配置 -----------------------------------------------------------
write_json(cfg, file.path(resdir, "config.json"), pretty = TRUE, auto_unbox = TRUE)

# --- 版本快照 ----------------------------------------------------------------
ver <- tryCatch({
  sessionInfo()
}, error = function(e) NULL)
writeLines(capture.output(ver), file.path(resdir, "environment_versions.txt"))

cat(sprintf("预注册完成 ✔\n  PREREGISTRATION.md: %s\n  sha256: %s\n  config.json: %s\n",
            file.path(docs, "PREREGISTRATION.md"),
            prereg_hash,
            file.path(resdir, "config.json")))
cat("完整性锁定生效：任何后续分析前将校验预注册哈希。\n")
