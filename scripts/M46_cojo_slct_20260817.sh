#!/bin/bash
# =============================================================================
# M46_cojo_slct_20260817.sh - GCTA-COJO conditional/joint fine-mapping (v2)
# =============================================================================
# 目的：审稿人任务 1。6 个 top 位点（RBM6/CD101/CNNM2/PLAUR/RIC8A/LAMC1）做
#   GCTA-COJO：
#   * Pass A --cojo-slct (p<5e-8)：全基因组独立的 index signals
#   * Pass B --cojo-slct (p<1e-4)：位点级独立的区域 signals（这些是亚 GWAS 的
#     效应基因位点，区域阈值更有信息量）
#   * Pass C --cojo-joint（Pass B index SNPs ∪ 候选 lead）：候选 lead 仪器在
#     区域内其他信号条件下的联合条件 p (pJ)
# 区域提取用 EUR.bim 的 hg19 坐标做 rsid 交集（GWAS 汇总为 GRCh38 坐标，
#   位置过滤会跨 build 失效；rsid 匹配是 build 无关的）。
# 输入：data/opengwas/full/{t2d,cad}_full.gz + data/ldref/1kg.v3/EUR.{bed,bim,fam}
# 输出：/data/qiushuogeng/tmp/cojo/{tag}.jma.cojo(.jma/.jmi/.cma) → results/m46_cojo_20260817.csv
# =============================================================================
set -e
cd /data/qiushuogeng/projects/dual-channel-mr-atlas
GCTA=/data/gengqiushuo/home/miniconda3/envs/gcta/bin/gcta64
EUR=data/ldref/1kg.v3/EUR
COJO=/data/qiushuogeng/tmp/cojo
mkdir -p $COJO
rm -f $COJO/*.ma $COJO/*.jma.cojo $COJO/*.cma.cojo $COJO/*.jma $COJO/*.jmi 2>/dev/null

# tag|outcome|chr|from(hg19)|to(hg19)|lead
LOCI="rbm6|t2d|3|49457459|50657459|rs10049087
cd101|t2d|1|116961774|118161774|rs10494191
cnnm2|cad|10|104158209|105358209|rs11191447
plaur|cad|19|43562473|44762473|rs4760
ric8a|cad|11|1|811312|rs6598075
lamc1|cad|1|182453661|183653661|rs10458355"

for L in $LOCI; do
  IFS='|' read TAG OUT CHR FROM TO LEAD <<< "$L"
  echo "=== $TAG x $OUT (chr$CHR:$FROM-$TO) lead=$LEAD ==="
  # 1) region rsids from hg19 bim
  awk -v c=$CHR -v f=$FROM -v t=$TO '$1==c && $4>=f && $4<=t {print $2}' $EUR.bim > $COJO/$TAG.rs
  echo "  region rsids (bim): $(wc -l < $COJO/$TAG.rs)"
  # 2) extract summary by rsid (build-agnostic)
  if [ "$OUT" = "t2d" ]; then
    zcat data/opengwas/full/t2d_full.gz | awk -F'\t' '
      NR==FNR{s[$1]=1; next}
      FNR==1 {print "SNP A1 A2 freq b se p N"; next}
      s[$2] && $19!="NA" && $20!="NA" && $21!="NA" && $18!="NA" \
      {print $2"\t"$6"\t"$5"\t"$18"\t"$19"\t"$20"\t"$21"\t"655666}' $COJO/$TAG.rs - > $COJO/$TAG.ma
  else
    zcat data/opengwas/full/cad_full.gz | awk -F'\t' '
      NR==FNR{s[$1]=1; next}
      FNR==1 {print "SNP A1 A2 freq b se p N"; next}
      # CAD file has two allele blocks (hm_* cols 5-6, lowercase cols 14-15); use the
      # self-consistent lowercase block: effect_allele=$14, other_allele=$15, eaf=$16,
      # beta=$20, se=$21, p=$22 (mixing hm_effect_allele=$6 with eaf=$16 mislabels freq).
      s[$2] && $20!="NA" && $21!="NA" && $22!="NA" && $16!="NA" \
      {print $2"\t"$14"\t"$15"\t"$16"\t"$20"\t"$21"\t"$22"\t"296525}' $COJO/$TAG.rs - > $COJO/$TAG.ma
  fi
  echo "  summary SNPs: $(($(wc -l < $COJO/$TAG.ma)-1)) (lead in file: $(grep -c "^$LEAD\b" $COJO/$TAG.ma))"
  # 3) Pass A: genome-wide index signals
  $GCTA --bfile $EUR --cojo-file $COJO/$TAG.ma --cojo-slct --cojo-p 5e-8 \
        --out $COJO/${TAG}_A > $COJO/${TAG}_A.log 2>&1 || echo "  pass A err"
  grep -E "associated SNPs are selected|No SNPs" $COJO/${TAG}_A.log | tail -1
  # 4) Pass B: locus-level index signals
  $GCTA --bfile $EUR --cojo-file $COJO/$TAG.ma --cojo-slct --cojo-p 1e-4 \
        --out $COJO/${TAG}_B > $COJO/${TAG}_B.log 2>&1 || echo "  pass B err"
  grep -E "associated SNPs are selected|No SNPs" $COJO/${TAG}_B.log | tail -1
  # 5) Pass C: cojo-joint on PassB index SNPs + candidate lead
  awk 'NR>1 {print $1}' $COJO/${TAG}_B.jma.cojo 2>/dev/null > $COJO/${TAG}_B.index
  grep -q "^$LEAD$" $COJO/${TAG}_B.index 2>/dev/null || echo "$LEAD" >> $COJO/${TAG}_B.index
  $GCTA --bfile $EUR --cojo-file $COJO/$TAG.ma --cojo-joint $COJO/${TAG}_B.index \
        --out $COJO/${TAG}_J > $COJO/${TAG}_J.log 2>&1 || echo "  pass C err"
  echo "  joint output: $(ls $COJO/${TAG}_J.jma.cojo 2>/dev/null && wc -l < $COJO/${TAG}_J.jma.cojo)"
done
echo "== DONE M46 COJO =="
