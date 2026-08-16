# 129 nominal-sig strong 的 r²-LD 聚类（2026-08-16，A8）

> 覆盖 atlas 129 个 nominal-sig strong（106 已知 + 23 新候选）；灰区 2 个（AP3S2×t2d、ZNF19×cad）为独立已知-GWAS 峰案例（见 R3），未入聚类。

> 参考面板：1000G EUR (hg19)，503 个体；匹配 123/129 个 lead；成对 r² 由基因型剂量相关平方给出。未匹配 6 个（rs147526786, rs17716350, rs1123462, rs12332382, rs77100387, rs4338909）——未在 1kg EUR 面板定位，按独立簇保守处理（不含在共簇计数）。

**r²≥0.8 → 106 个独立簇**（多 SNP 簇 14 个；另有 6 个未匹配 lead 按独立计）
**r²≥0.6 → 104 个独立簇**（敏感性）

## 多 SNP 簇（r²≥0.8）
- known:RRN3:t2d (rs2280018) + known:NTAN1:t2d (rs2280018)
- known:CADM4:cad (rs4760) + new:PLAUR:cad (rs4760)
- known:FBXW7:t2d (rs9991574) + known:DKFZP434I0714:t2d (rs9991574)
- known:MED23:t2d (rs2246012) + known:ARG1:t2d (rs2246012)
- known:MLH3:cad (rs175065) + known:EIF2B2:cad (rs175057)
- known:RFT1:t2d (rs2336725) + known:SERBP1P3:t2d (rs4687701)
- known:NUDT5:t2d (rs11257655) + known:CAMK1D:t2d (rs11257655)
- known:BLOC1S2:t2d (rs17882802) + known:PHBP9:t2d (rs11591741) + new:CWF19L1:t2d (rs17668357)
- known:RP11-332O19.3:cad (rs11191472) + new:CNNM2:cad (rs11191447)
- known:LRRC37A15P:t2d (rs223490) + known:KRT8P46:t2d (rs223490) + known:RP11-10L12.2:t2d (rs223490)
- known:RP11-464F9.1:t2d (rs7908825) + known:RP11-464F9.9:t2d (rs11000760)
- known:RP11-51F16.5:t2d (rs7223412) + known:STRADA:t2d (rs11079508)
- known:Metazoa_SRP:cad (rs11204675) + known:Metazoa_SRP:cad (rs10888385)
- new:TAGLN2:cad (rs2789422) + new:CCDC19:cad (rs2789422) + new:VSIG8:cad (rs2789422)

## 23 新候选的 LD 独立性（相对已知 106）
- new:PLAUR:cad → 与 1 个 lead 共簇
- new:SLC12A3:cad → LD 独立簇
- new:C2orf49:t2d → LD 独立簇
- new:CWF19L1:t2d → 与 2 个 lead 共簇
- new:U6atac:t2d → LD 独立簇
- new:CD101:t2d → LD 独立簇
- new:TAGLN2:cad → 与 2 个 lead 共簇
- new:CCDC19:cad → 与 2 个 lead 共簇
- new:RBM6:t2d → LD 独立簇
- new:CNNM2:cad → 与 1 个 lead 共簇
- new:N4BP2L2:cad → lead 未在 1kg EUR 定位，按独立簇保守计
- new:VSIG8:cad → 与 2 个 lead 共簇
- new:LAMC1:cad → LD 独立簇
- new:TPD52:fbg → LD 独立簇
- new:SENP6:cad → LD 独立簇
- new:RIC8A:cad → LD 独立簇
- new:HMGN3:cad → LD 独立簇
- new:MT3:fbg → LD 独立簇
- new:RPL13:t2d → LD 独立簇
- new:ZBTB46:cad → LD 独立簇
- new:PDCD6:t2d → lead 未在 1kg EUR 定位，按独立簇保守计
- new:CLEC3B:t2d → lead 未在 1kg EUR 定位，按独立簇保守计
- new:ZNF100:cad → lead 未在 1kg EUR 定位，按独立簇保守计

> 解读：15 个'已知位点内新效应基因'中，PLAUR（与 CADM4 同 rs4760）、CWF19L1（与 BLOC1S2/PHBP9 共簇）、
> CNNM2（与 RP11-332O19.3 共簇）在 lead 层面与已知 strong 共簇——属'同一已知信号上的新效应基因提名'；
> TAGLN2/CCDC19/VSIG8 共用同一 eQTL lead rs2789422（cis 多基因共享信号）。其余新候选 lead 与已知 lead 无 r²≥0.8 共簇。

## 诚实边界
- 106 已知 top_snp 为 coloc 峰 SNP、23 新候选为 eQTL lead SNP，两口径在 LD 聚类中仅作'是否同一信号'近似判定。
- r² 阈值 0.8 主口径、0.6 敏感性；单链聚类（重叠团合并）；缺失基因型按对剔除。
- '独立位点' = LD 不连通分量，非'独立因果变异'。
