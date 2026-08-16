# 杠杆 B：脂质结局 MR 验证（2026-08-13，探索性）

> 目的：验证 deCODE 蛋白工具是否有效（蛋白通道能否被救活）。
> 若 PCSK9/APOC3/ANGPTL3 等已知脂质因果基因在 LDL/HDL/TG/apoB 结局下能恢复教科书效应，
> 则证明蛋白通道的"工具侧"没问题，"心血管结局未达标"是结局侧信号问题而非工具坏。

## 方法（M12_lipid_validation.R，探索性）

- 结局：Richardson 2020 UKB（欧洲，与 deCODE 冰岛 n=35,559 **无样本重叠**）
  LDL `ieu-b-110` n=440,546 / HDL `ieu-b-109` n=403,943 / TG `ieu-b-111` n=441,016 / apoB `ieu-b-108` n=439,214
- 方法：与主流程（M3）完全一致 —— deCODE cis ±1Mb、p<5e-6、LD clump r²<0.01@1000kb EUR、
  主 IVW-MRE + FE/WM/Egger 敏感性、nsnp=1→Wald。
- 输出：`results/grid/protein_decode_lipid_mr.csv`（11 蛋白 × 4 脂质结局，全网格含空）

## 结果（主方法 IVW-MRE / 单工具 Wald）

| 蛋白 | LDL-C | HDL-C | TG | apoB | 校准解读 |
|---|---|---|---|---|---|
| **PCSK9** | **+0.284, 5.2e-15** | +0.007, 0.30 | -0.001, 0.86 | **+0.281, 2.4e-15** | **教科书 ✓**（PCSK9↑→LDL/apoB↑；TG/HDL 正确 null） |
| **APOC3** | +0.146, 2.8e-4 | **-0.389, 2.6e-8** | — | +0.192, 3.0e-4 | **教科书 ✓**（APOC3↑→HDL↓；LOF 升 HDL） |
| ANGPTL3 | — | — | +0.099, 6.9e-3 | +0.039, 6.3e-3 | 预期方向但**不稳健**（WM null/Egger NS，cis 多效性） |
| PCK1 | -0.009, 0.86 | -0.127, 6.5e-3 | +0.188, 1.1e-4 | +0.013, 0.80 | 单工具 nominal（PCK1 非脂质基因，疑似假阳性） |
| INSR | -0.006, 0.92 | +0.026, 0.80 | -0.021, 0.92 | -0.020, 0.82 | **正确负对照 ✓**（INSR 非脂质基因全 null） |

### 敏感性（PCSK9/APOC3 稳健性核查）

- **PCSK9 × LDL**：MRE 5.2e-15 / FE 0 / WM 2.2e-8 / Egger 8.2e-6 → **全方法稳健**
- **PCSK9 × apoB**：MRE 2.4e-15 / FE 0 / WM 1.8e-10 / Egger 5.2e-6 → **全方法稳健**
- **APOC3 × HDL**：MRE 2.6e-8 / FE 0 / WM 3.9e-4 / Egger 9.4e-3 → **稳健**
- **APOC3 × LDL**：MRE 2.8e-4 / WM 0.011 / Egger 0.072 → 方向一致但仅 nominal（次要）
- **ANGPTL3 × TG**：MRE 0.0069 / **WM 0.56** / Egger 0.13 → 不稳健，**不声称**

## 诚实解读（结论）

1. **deCODE 蛋白工具本身是有效的**：PCSK9→LDL/apoB（p~1e-15）与 APOC3→HDL（p=2.6e-8）
   完美恢复已知因果方向，且全方法稳健。PCSK9 的 TG/HDL 正确为 null，INSR 全 null（负对照正确）。
2. **蛋白通道在心血管结局（T2D/CAD/FBG）未达标的解释**：不是工具坏，而是结局侧信号更弱/
   多信号（如 PCSK9×CAD 的 HEIDI 异质性 p=1.28e-4 是 CAD 区域多信号，非工具失败）。
   这回答了审稿人"蛋白工具是坏的"质疑——校准证据在脂质空间。
3. **不作新发现声称**：PCSK9→LDL、APOC3→HDL 是已知药物靶点的教科书阳性对照
   （PCSK9 抑制剂/Volanesorsen 表型），只作"管线能正确回收已知生物学"的校准展品。
4. **探索性声明**：脂质结局在预注册之外，为杠杆 B 探索性扩展；预注册 v6 修订需记录结局扩展。
5. **已知局限**：ANGPTL3→TG 不稳健（cis 区多效性，强 LOF 变异 + 多 SNP 异质性），
   APOC3→LDL 仅 nominal，PCK1 单工具 nominal 疑似假阳性。只报铁板校准（PCSK9→LDL/apoB、APOC3→HDL）。

## 下一步（与 M13 UKB-PPP 复现联动）

- M13 pQTL 方向一致率：PCSK9 r=0.338 (p<1e-232)、ANGPTL3 r=0.396 (p<1e-287) 共享变异 beta 高度一致；
  GLP1R/GCG 弱信号 r 不显著（诚实：弱 cis 蛋白跨平台差异大）。
- 合并叙事：**pQTL 工具层（deCODE vs UKB-PPP 共享变异一致）+ 脂质结局校准（PCSK9/APOC3 教科书回收）
  + 心血管结局如实探索** —— 蛋白通道从"未校准"升级为"已校准 + 心血管结局信号弱"的诚实三层。
