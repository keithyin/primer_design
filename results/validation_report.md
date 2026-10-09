# 引物面板数据校验报告

生成时间：数据流水线 results/primer_design_results.tsv + designed.json

**结果：16 项通过，0 项失败，共 16 项检查**

## X. 交付数据内部一致性
- [PASS] (X) TSV 行数 == designed.json 有效对数 — tsv=24, json=24
- [PASS] (X) TSV 引物序列与 designed.json 一致
- [PASS] (X) product_size == 扩增子坐标跨度（1-based 闭区间）
- [PASS] (X) 重复区长度与 Repeats 特点中拷贝数自洽
- [PASS] (X) 档位标签与拷贝数一致（扩展检测位除外）
- [PASS] (X) >30 ✓ 行的单元数 >= 30
- [PASS] (X) >30 扩展检测行存在且为 FMR1/AFF2 — n=2
## A. 参考基因组模拟 PCR（不经 primer3 的独立复核）
- [PASS] (A) 24 对引物结合位点/产物长度/包含关系全部复核
## B. 重复 tract 独立重分析（正则/相位法 vs 设计时 walking 法）
- [PASS] (B) 重复拷贝数/中断/纯度独立复核（正则+相位法）
## C. 热力学复核（primer3 calcTm + Biopython Tm_NN 双重计算）
- [PASS] (C1) primer3 calcTm 复算与报告 Tm 一致（≤0.3℃） — 最大偏差 0.046℃
- [INFO] Biopython Tm_NN 独立第二意见：与报告 Tm 最大偏差 1.81℃
- [PASS] (C2) Biopython Tm_NN 偏差 ≤2.0℃（跨引擎容差）
- [PASS] (C3) 引物 GC% 复算一致
- [PASS] (C4) homodimer/heterodimer/hairpin ΔG ≥ -45 kcal/mol 复算
## E. 基因注释交叉核对（UCSC ncbiRefSeqCurated vs Ensembl）
- [INFO] #3 chr12:7628695 Ensembl=ENSG00000303381（无符号），UCSC RefSeq 重叠基因：[]
- [INFO] #5 chr7:100694312 GIGYF1：Ensembl 基因体包含该位点，RefSeq 转录本边界距位点 61bp（数据库模型边界差异）
- [INFO] #11 chr13:53201089 PCDH8P1（transcribed_processed_pseudogene）：RefSeq Curated 无对应模型（Ensembl 独有注释），无法交叉核对
- [INFO] #16 chr1:88275022 PKN2-AS1（lncRNA）：RefSeq Curated 无对应模型（Ensembl 独有注释），无法交叉核对
- [INFO] #17 chr12:48971839 WNT10B：Ensembl 基因体包含该位点，RefSeq 转录本边界距位点 103bp（数据库模型边界差异）
- [INFO] #19 chr2:126245765 Ensembl=ENSG00000290113（无符号），UCSC RefSeq 重叠基因：[]
- [INFO] #21 chr10:104587500 Ensembl=ENSG00000237761（无符号），UCSC RefSeq 重叠基因：[]
- [PASS] (E) UCSC RefSeq 与 Ensembl 注释交叉核对
- [INFO] 注释核对：24/24 行与 RefSeq 无矛盾，其中 17/24 行为基因符号直接匹配；其余为 RefSeq 无模型或仅有边界差异（见上列 INFO）
## D. 全基因组特异性独立重扫（numpy 滚动哈希 vs grep -F）
- [PASS] (D1a) 全长脱靶计数：grep 与 numpy 两种实现完全一致（48/48）
- [INFO] 3'-15mer 脱靶计数：7/48 条 grep 因非重叠匹配语义少计（numpy 按位置计数为准），全部满足 numpy ≥ grep，方向与预期一致
- [PASS] (D1b) 3'-15mer 计数关系 numpy ≥ grep（grep -o 非重叠语义）
- [INFO] 扫描范围：25 个 contig，3.09 Gbp（正链编码，模式含正反两向）
