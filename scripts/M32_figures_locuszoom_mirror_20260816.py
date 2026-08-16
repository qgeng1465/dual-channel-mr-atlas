#!/usr/bin/env python3
# =============================================================================
# M32_figures_locuszoom_mirror_20260816.py — Fig 5：区域共定位镜像图（locusZoom 型）
# =============================================================================
# AJHG 审稿人核心期待：逐基因座"上 GWAS −log10p / 下 eQTL −log10p / 共享基因轨道 /
#   PP.H4 直标"的镜像图。6 个代表位点 = 5 个新候选（RBM6/CNNM2/PLAUR/CD101/LAMC1）
#   + 1 个灰区 strong（AP3S2，证明"真实位点单工具 MR 功效不足"）。
# 数据：
#   - LD r²：plink --r2 vs lead eQTL SNP（1000G EUR 基因型，hg19）
#   - eQTL：eQTLGen full cis（bychr 明文块，hg19），Zscore → p
#   - GWAS：OpenGWAS full（t2d/cad 为 hg38，按 rsid 匹配回 1kg hg19 坐标统一坐标系）
#   - 基因轨道：cis-EQTL-significant 中落在窗口内的基因（诚实：仅"有显著 eQTL 的基因"）
# 输出：results/figures/20260816_F5_locuszoom_mirror.png
# =============================================================================
import csv, os, subprocess, sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import font_manager
from scipy.stats import norm as _norm

for f in font_manager.fontManager.ttflist:
    if f.name == "Noto Sans CJK JP":
        plt.rcParams["font.sans-serif"] = [f.name]
        break
plt.rcParams["axes.unicode_minus"] = False

PROJ  = "/data/qiushuogeng/projects/dual-channel-mr-atlas"
RES   = f"{PROJ}/results"
EUR   = f"{PROJ}/data/ldref/1kg.v3/EUR"
BYCHR = "/data/qiushuogeng/tmp/eqtlgen_stable/bychr"
SIG   = "/data/qiushuogeng/tmp/eqtlgen_stable/cis-EQTL-significant.txt.gz"
GWASD = f"{PROJ}/data/opengwas/full"
OUT   = f"{RES}/figures"
TMP   = "/data/qiushuogeng/tmp/mirror"
os.makedirs(OUT, exist_ok=True); os.makedirs(TMP, exist_ok=True)
PLINK = f"{PROJ}/tools/plink"

# ---- 6 个代表位点（坐标 = hg19）----
LOCI = [
    dict(sym="RBM6",  ensg="ENSG00000004534", out="t2d", chr=3,  pos=50057459,  anchor="rs10049087", grey=False),
    dict(sym="CNNM2", ensg="ENSG00000148842", out="cad", chr=10, pos=104758209, anchor="rs11191447", grey=False),
    dict(sym="PLAUR", ensg="ENSG00000011422", out="cad", chr=19, pos=44162473,  anchor="rs4760",     grey=False),
    dict(sym="CD101", ensg="ENSG00000134256", out="t2d", chr=1,  pos=117561774, anchor="rs10494191", grey=False),
    dict(sym="LAMC1", ensg="ENSG00000135862", out="cad", chr=1,  pos=183053661, anchor="rs10458355", grey=False),
    dict(sym="AP3S2", ensg="ENSG00000157823", out="t2d", chr=15, pos=90405702,  anchor="rs10852122", grey=True),
]
WIND = 600_000
CONF = dict(t2d="T2D", cad="CAD")

# ---- PP.H4 查表 ----
def pp4_lookup():
    tab = {}
    for out in ("t2d", "cad", "fbg"):
        for r in csv.DictReader(open(f"{RES}/coloc_full_{out}_20260815.csv")):
            if r["ok"] == "TRUE":
                tab[(r["gene"], out)] = float(r["pp4"])
    return tab
PP4 = pp4_lookup()

# ---- plink LD（每个位点一次）----
def calc_ld(loc):
    tag = loc["sym"].lower()
    rc = subprocess.run([PLINK, "--bfile", EUR, "--chr", str(loc["chr"]),
                         "--from-bp", str(loc["pos"]-WIND), "--to-bp", str(loc["pos"]+WIND),
                         "--r2", "--ld-snp", loc["anchor"], "--ld-window-kb", str(2*WIND//1000),
                         "--ld-window", "99999", "--ld-window-r2", "0",
                         "--out", f"{TMP}/{tag}_ld"],
                        capture_output=True, text=True).returncode
    r2 = {}
    ld_path = f"{TMP}/{tag}_ld.ld"
    if rc != 0 or not os.path.exists(ld_path):
        print(f"  [WARN] {loc['sym']}×{loc['out']} plink LD 失败 (rc={rc})，该位点按无 LD 着色")
        return r2
    with open(ld_path) as f:
        next(f)  # 跳表头
        for line in f:
            p = line.split()
            if len(p) >= 7:
                r2[p[5]] = float(p[6])   # SNP_B, R2
    r2[loc["anchor"]] = 1.0
    return r2

# ---- eQTL：bychr 提取该基因全部 cis 行 ----
def calc_eqtl(loc):
    cmd = f"awk -F'\\t' '$7==\"{loc['ensg']}\"' {BYCHR}/chr{loc['chr']}.tsv"
    rows = []
    with os.popen(cmd) as f:
        for line in f:
            p = line.rstrip("\n").split("\t")
            if len(p) < 10: continue
            snp, sp, z = p[0], int(p[2]), float(p[3])
            if abs(sp - loc["pos"]) <= WIND:
                rows.append((snp, sp, z))
    return rows

# ---- GWAS：区域 rsid 匹配 ----
def calc_gwas(loc, region_rs):
    cols = {}
    with os.popen(f"zcat {GWASD}/{loc['out']}_full.gz") as f:
        hdr = f.readline().rstrip("\n").split("\t")
        for i, c in enumerate(hdr): cols[c] = i
    ri = cols["hm_rsid"]; pi = cols["p_value"]
    g = {}
    with os.popen(f"zcat {GWASD}/{loc['out']}_full.gz") as f:
        first = True
        for line in f:
            if first: first = False; continue
            p = line.rstrip("\n").split("\t")
            if p[ri] in region_rs:
                try: v = float(p[pi])
                except ValueError: continue
                if not np.isnan(v) and v > 0: g[p[ri]] = v
    return g

# ---- 基因轨道：窗口内"有显著 cis-eQTL 的基因" ----
def gene_track(loc):
    out = []
    cmd = (f"zcat {SIG} | awk -F'\\t' '$10=={loc['chr']} && $11>={loc['pos']-WIND} "
           f"&& $11<={loc['pos']+WIND} {{print $9\"\\t\"$11}}' | sort -u")
    with os.popen(cmd) as f:
        for line in f:
            s, p = line.rstrip("\n").split("\t")
            out.append((s, int(p)))
    return out

# ---- LD 调色板（locusZoom 经典，略调色盲安全）----
def ld_color(r2):
    if r2 is None or r2 < 0.2:  return "#9e9e9e"
    if r2 < 0.4: return "#8f9c2f"
    if r2 < 0.6: return "#3aa655"
    if r2 < 0.8: return "#1e90ff"
    return "#8b008b"

# ---- 绘制：3 行 × 2 列位点，每位点垂直堆叠 GWAS(上)/eQTL(下) ----
from matplotlib.lines import Line2D

fig = plt.figure(figsize=(10.2, 12.5), dpi=300)
gs = fig.add_gridspec(6, 2, hspace=0.14, wspace=0.22,
                      height_ratios=[1, 0.72, 1, 0.72, 1, 0.72])  # 3 行位点，每行 GWAS+eQTL

for i, loc in enumerate(LOCI):
    col = i % 2
    row = (i // 2) * 2          # 位点行 0,2,4
    axg = fig.add_subplot(gs[row, col])
    axe = fig.add_subplot(gs[row+1, col], sharex=axg)
    tag = loc["sym"].lower()
    # LD
    r2 = calc_ld(loc)
    # eQTL
    eq = calc_eqtl(loc)
    eq = [(s, p, 2*_norm.sf(abs(z))) for s, p, z in eq]
    # region rsid set（bim，hg19 坐标作 x 轴基准）
    region_rs = {s for s, _, _ in eq}
    # GWAS
    gw = calc_gwas(loc, region_rs)
    # x 轴：bim hg19 位置
    bim_pos = {}
    with open(f"{EUR}.bim") as f:
        for line in f:
            p = line.split()
            if p[1] in region_rs: bim_pos[p[1]] = int(p[3])

    # ---- 上板：GWAS -log10p（LD 着色）----
    gx, gy, gc = [], [], []
    for s, p in gw.items():
        pos = bim_pos.get(s)
        if pos is None: continue
        gx.append(pos); gy.append(-np.log10(p)); gc.append(ld_color(r2.get(s)))
    if gx:
        axg.scatter(gx, gy, s=11, c=gc, edgecolor="none", zorder=2)
        ap = bim_pos.get(loc["anchor"])
        ap_p = gw.get(loc["anchor"])
        if ap is not None and ap_p is not None:
            axg.scatter([ap], [-np.log10(ap_p)], marker="D", s=38, c="#e00000",
                        edgecolor="black", linewidth=0.5, zorder=4)
            axg.annotate(loc["anchor"], (ap, -np.log10(ap_p)), textcoords="offset points",
                         xytext=(3, 3), fontsize=5.6, color="#e00000")
    axg.set_ylabel(f"GWAS $-\\log_{{10}}(p)$\n({CONF[loc['out']]})", fontsize=6.5)
    axg.set_title(f"{loc['sym']}  ({loc['ensg']})", fontsize=8, loc="left", fontweight="bold", pad=1)
    axg.set_ylim(0, max(gy) * 1.18 if gy else 1)
    axg.tick_params(labelsize=6)
    axg.set_xticklabels([])
    axg.spines[["top", "right"]].set_visible(False)

    # ---- 下板：eQTL -log10p（LD 着色），y 顶值设 cap ----
    ex, ey, ec = [], [], []
    for s, p, pv in eq:
        if s not in bim_pos: continue
        ex.append(bim_pos[s]); ey.append(-np.log10(max(pv, 1e-300))); ec.append(ld_color(r2.get(s)))
    cap = 50.0
    ey_capped = [min(v, cap) for v in ey]
    axe.scatter(ex, ey_capped, s=11, c=ec, edgecolor="none", zorder=2)
    ap = bim_pos.get(loc["anchor"])
    ap_z = next((z for s, _, z in eq if s == loc["anchor"]), None)
    if ap is not None and ap_z is not None:
        ap_p = 2*_norm.sf(abs(ap_z))
        ap_y = min(-np.log10(max(ap_p, 1e-300)), cap)
        axe.scatter([ap], [ap_y], marker="D", s=38, c="#e00000", edgecolor="black", linewidth=0.5, zorder=4)
    maxZ = max((z for _, _, z in eq), default=0)
    if max(ey) > cap:
        axe.text(0.99, 0.96, f"top eQTL Z={maxZ:.1f}\n($-\\log_{{10}}(p)$>50 capped)",
                 transform=axe.transAxes, fontsize=5.2, ha="right", va="top", color="#444444")
    axe.set_ylabel("eQTL $-\\log_{10}(p)$", fontsize=6.5)
    axe.tick_params(labelsize=6)
    axe.spines[["top", "right"]].set_visible(False)

    # ---- 面板字母（AJHG 规范：大写粗体左上角）----
    axg.text(0.012, 0.985, "ABCDEF"[i], transform=axg.transAxes, fontsize=10, fontweight="bold",
             color="#111111", va="top")
    # ---- PP.H4 标注（上板左上，面板字母下方）----
    pp = PP4.get((loc["ensg"], loc["out"]))
    lab = f"PP.H4 = {pp:.3f}" if pp is not None else "PP.H4 = –"
    if loc["grey"]:
        lab = f"PP.H4 = {pp:.3f}\ngrey-zone (MR p>0.05)"
    axg.text(0.03, 0.72, lab, transform=axg.transAxes, fontsize=6.2, fontweight="bold",
             color="#c00000" if loc["grey"] else "#1a3a6b", va="top", linespacing=1.25,
             bbox=dict(boxstyle="round,pad=0.25", fc="white", ec="#cccccc", lw=0.5))

    # ---- 基因轨道（每位点 eQTL 面板底部）----
    gt = gene_track(loc)
    ymin, yspan = axe.get_ylim()[0], axe.get_ylim()[1] - axe.get_ylim()[0]
    track_y = ymin - 0.30*yspan
    for sym, pos in gt:
        axe.plot([pos, pos], [track_y, track_y + 0.13*yspan], color="#999999", lw=1.1, solid_capstyle="butt")
        axe.text(pos, track_y - 0.07*yspan, sym, fontsize=4.4, ha="center", va="top", color="#555555")
    # 基因本体粗标
    axe.plot([loc["pos"], loc["pos"]], [track_y + 0.26*yspan, track_y + 0.46*yspan],
             color="#c00000", lw=2.4, solid_capstyle="butt")
    axe.text(loc["pos"], track_y + 0.50*yspan, loc["sym"], fontsize=6.2, ha="center",
             color="#c00000", fontweight="bold")

    # x 轴刻度（每行下板）
    axe.set_xlabel(f"chr{loc['chr']} (hg19, Mb)", fontsize=6.5)
    axe.set_xticks(np.linspace(loc["pos"]-WIND, loc["pos"]+WIND, 5))
    axe.set_xticklabels([f"{v/1e6:.2f}" for v in np.linspace(loc["pos"]-WIND, loc["pos"]+WIND, 5)], fontsize=6)
    for axx in (axg, axe):
        axx.axvline(loc["pos"], color="#cccccc", lw=0.6, ls=":", zorder=1)

# LD 图例（顶部）
legend_handles = [
    Line2D([0], [0], marker="o", color="w", markerfacecolor="#9e9e9e", ms=5, label="r²<0.2"),
    Line2D([0], [0], marker="o", color="w", markerfacecolor="#8f9c2f", ms=5, label="0.2–0.4"),
    Line2D([0], [0], marker="o", color="w", markerfacecolor="#3aa655", ms=5, label="0.4–0.6"),
    Line2D([0], [0], marker="o", color="w", markerfacecolor="#1e90ff", ms=5, label="0.6–0.8"),
    Line2D([0], [0], marker="o", color="w", markerfacecolor="#8b008b", ms=5, label="r²≥0.8"),
    Line2D([0], [0], marker="D", color="w", markerfacecolor="#e00000", ms=6, label="lead eQTL SNP"),
]
fig.legend(handles=legend_handles, loc="upper center", ncol=6, fontsize=6.5, frameon=False,
           bbox_to_anchor=(0.5, 0.998))
fig.suptitle("Representative colocalization loci (GWAS / eQTL mirror plots; LD vs lead eQTL SNP, 1000G EUR)",
             fontsize=9, y=0.995)
out_p = f"{OUT}/20260816_F5_locuszoom_mirror.png"
fig.savefig(out_p, bbox_inches="tight")
plt.close(fig)
print("Fig 5 OK ->", out_p)
for loc in LOCI:
    print(f"  {loc['sym']}×{loc['out']} PP.H4={PP4.get((loc['ensg'], loc['out']))}")
