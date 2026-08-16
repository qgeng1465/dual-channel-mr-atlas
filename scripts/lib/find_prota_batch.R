#!/usr/bin/env Rscript
# 分批查询 OpenGWAS 全量 prot-a-* 的 id→trait 映射，定位正/负对照蛋白
suppressPackageStartupMessages(library(ieugwasr))
resdir <- "/data/qiushuogeng/projects/dual-channel-mr-atlas/results"
ids_all <- paste0("prot-a-", 1:5000)
targets <- c(PCSK9="Proprotein convertase subtilisin/kexin type 9",
             CETP="Cholesteryl ester transfer protein",
             HMGCR="glutaryl-CoA reductase",
             NPC1L1="NPC1",
             ANGPTL3="Angiopoietin-related protein 3",
             APOC3="Apolipoprotein C-III")
res <- list()
for (start in seq(1, length(ids_all), by=100)) {
  batch <- ids_all[start:min(start+99, length(ids_all))]
  info <- tryCatch(gwasinfo(batch), error=function(e){cat("批次ERR@",start,":",conditionMessage(e),"\n");NULL})
  if (!is.null(info) && nrow(info)>0) {
    info <- info[info$id %in% batch, c("id","trait","sample_size")]
    res[[length(res)+1]] <- info
  }
  if (start %% 500 == 1) cat("进度:", start, "/", length(ids_all), "\n")
}
prota <- do.call(rbind, res)
prota <- prota[!duplicated(prota$id), ]
saveRDS(prota, file.path(resdir,"prota_index.rds"))
cat("缓存", nrow(prota), "条 prot-a 索引\n")
for (g in names(targets)) {
  hits <- prota[grepl(targets[g], prota$trait, ignore.case=TRUE), ]
  if (nrow(hits)>0) cat(g, "=>", paste(apply(hits,1,function(r) paste0(r["id"]," (",r["trait"],", n=",r["sample_size"],")")), collapse=" | "), "\n")
  else cat(g, "=> 未找到\n")
}
