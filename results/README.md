# GC / GCC 重复序列 PCR 引物面板（GRCh38 / hg38）— 已通过多轮独立校验

> **第 3 版说明**：本轮补做了上一版列为"最大方法学缺口"的**错配容忍 primer-pair 特异性验证**
> （`results/mismatch_pcr.*`），并把 dimer/hairpin 从"仅采用 primer3 约束"升级为**独立实现复算**
> （`results/dimer_check.*`）；同时对参考基因组做了身份校验，并据此修复了一处真实的参考不一致。
> 结论分级与边界见 **`results/final_report.md` §6.4、§13、§14**。
> **配对特异性未发现可扩增脱靶，但这不等于"无非特异扩增"——湿实验仍是最终证据。**

## 需求

为 GC 二核苷酸重复（(GC)n，"GC repeats"）与 GC 三核苷酸重复（(GCC)n 家族，"GCC repeats"）
设计侧翼 PCR 引物，拷贝数覆盖 20–25 / 25–30 / >30 三档，扩增子 ≤1 kb（硬上限 2 kb），
给出 ≥20 对引物，并附设计与验证注意事项。**最终交付 24 对引物，全部经两轮独立校验：
16 项数据校验（validation_report.md，方法 A–E）+ 24 行逐行从头复核（data_review.md，P1–P5）。**

## 关键事实：参考基因组中 GC 重复的真实分布

对 GRCh38 主装配（UCSC hg38，455 条 contig 全部扫描）做序列级串联扫描（TRF simpleRepeat
轨道因打分阈值会漏掉短 GC 重复，故直接扫描序列）。因"含中断"重复的边界对种子阈值存在
±1–2 单元的定义敏感性，下表并列三种种子阈值的稳健结果（已做 15bp/12bp 更低种子的敏感性
检查，见 `data/sensitivity_15bp.json`、`data/sensitivity_dinuc.json`，结论均未再变化）：

| 家族 | 完美串联最长 | 含中断（单元计）最长 | ≥20 拷贝位点数 | ≥30 拷贝位点数 |
|---|---|---|---|---|
| (GC)n 二核苷酸 | **14 拷贝**（chr19:47,271,733） | **18**（chrX:2,574,014 区，含中断口径） | **0**（12bp 种子下仍为 0） | **0** |
| (GCC)n 家族（6 种阅读框为同一重复） | **23 拷贝**（chr7 GIGYF1 的 (CCG) 完美段） | 常规种子（≥10 单元纯核）25；密扫种子（≥8 单元）**32**（chr12:7,628,695） | 常规 6；密扫 **10** | 常规 0；密扫 **1** |

原因：参考基因组代表"健康基线"单倍型。(CGG)>55、(GGC)>45 等患者型长等位不存在于参考装配。
因此档位覆盖如下：
- **"GCC 20–25"**：4 个真实位点（FMR1 22、HMGB3 22、AFF2 21、TMEM185A 20）+ 7 个 16–19 拷贝的
  "最接近"位点（表中标注 △）；
- **"GCC 25–30"**：2 个真实位点（CARM1 25、GIGYF1 25；参考基因组中 25 即为该档最高）；
- **"GCC >30"**：密扫发现真实位点 **chr12:7,628,695（32 单元 = 26 精确 + 6 中断，纯度 81%，
  结构 GCT(GCC)₂CCT(GCC)₆CCT(GCC)₈CCT(GCC)₄CCT(GCC)₆GCT）已纳入面板（第 3 行）**；
  另保留 FMR1、AFF2 作扩展检测设计（GRCh38 参考序列中 22/21 拷贝；1 kb 产物理论上限
  ~114 / ~93 拷贝，**只能覆盖前突变区间的下段，不能覆盖完整 55–200**；若要覆盖 200 拷贝，
  产物上限需放宽到 ~1.26 / ~1.32 kb，见表中 Repeats 特点列）；
- **"(GC)n ≥20"**：参考基因组不存在（完美口径最长 14，含中断口径最长 18）。
  如需 ≥20 拷贝 (GC)n 标准品，请按文末第 8 条合成。

## 设计方法

1. **位点发现**：序列级扫描（正反链、全部阅读框），纯 GC 核种子 + 严格外推
   （≤2 个连续非精确单元，运行精确率 ≥60%，总纯度 ≥75%），硬性不变量
   `tract 长度 == 单元数 × 单元长`（已对全部扫描输出断言）。
   常规种子 ≥10 单元（30bp）；对"≥20/>30 位点盘点"另做 ≥8 单元（24bp）密扫与
   ≥5 单元（15bp）敏感性扫描；二核苷酸行按"最大完美串联"定义（两侧不可再延伸）。
2. **引物设计**：primer3（primer3-py 2.3.1）：
   - 长度 20–28（偏好 22），Tm 57–64（偏好 60℃），GC 35–65%；
   - 50 mM 单价盐 / 1.5 mM Mg²⁺ / 50 nM DNA / 0.6 mM dNTP；
   - 二聚体/发夹约束用 primer3 的**热力学参数族**（primer3-py 默认开启 thermodynamic
     oligo alignment，故 `PRIMER_MAX_SELF_ANY` / `_SELF_END` 等非 `_TH` 键**不生效**，
     `PRIMER_MAX_HAIRPIN` 更是根本不存在的键名、会被静默忽略）。实际生效的是
     `PRIMER_MAX_SELF_ANY_TH` / `_SELF_END_TH` / `_HAIRPIN_TH` 与
     `PRIMER_PAIR_MAX_COMPL_ANY_TH` / `_END_TH`，均取 primer3 默认值 47.0。
     换算成物理量（对 `calc_homodimer`/`calc_hairpin` 拟合标定）约相当于
     同源二聚体 ΔG ≥ −12.5 kcal/mol、发夹 ΔG ≥ −0.8 kcal/mol；
     3' 端稳定性 `PRIMER_MAX_END_STABILITY` ≤ 6.0（primer3 默认 100，此项为实收紧）；
   - `SEQUENCE_TARGET` = 重复 tract（强制扩增子夹住重复），
     `SEQUENCE_EXCLUDED_REGION` = tract ± 相邻周期-2 简单重复（自动扩展防止引物落入
     旁邻 (CA)n/(GT)n 微卫星）；
   - 扩增子阶梯 250–900 → 200–1000 → 200–1600 → 200–2000 bp（实际用到 [250,900]×16、
     [200,1000]×6、[200,900]×2；最终 338–990 bp，未动用 2 kb 上限）；
   - 每位移 10 个候选对；坐标约定（实证验证）：`PRIMER_LEFT` start 为左引物 5' 端 0-based、
     `PRIMER_RIGHT` start 为右引物**最右侧碱基（即右引物的 5' 端）**的 0-based 坐标，
     右引物 3' 端位于 `rp - len + 1`；`产物 = rp + 1 - lp`。
3. **特异性择优**：对全部候选引物做全基因组匹配（正反链全长 + 3' 端 15-mer），
   按"全长脱靶 = 0 优先、3' 15-mer 脱靶最少、产物尽量短"择优。
4. **权威特异性计数**：numpy 滚动哈希全基因组重扫（3.21 Gbp，455 contig），按基因组
   位置计数（含重叠出现）。grep -o 只报非重叠匹配，在 (CA)n 区会少计，仅用于排序。
5. **注释**：Ensembl REST（GRCh38）+ UCSC ncbiRefSeqCurated 交叉核对。

## 数据校验

### 第一轮：方法与数据校验（`results/validation_report.md`，16 项全部通过）

| 方法 | 内容 | 结果 |
|---|---|---|
| X | TSV 与 designed.json 内部一致性、产物长=坐标跨度、拷贝数自洽、档位标签（含 >30 行） | PASS |
| A | 参考基因组模拟 PCR（不经 primer3）：引物位点唯一、产物长度重算、重复区完整包含 | PASS |
| B | 重复 tract 独立重分析（正则+相位法，与设计算法不同实现） | PASS |
| C | 热力学：primer3 calcTm（最大偏差 0.046℃）+ Biopython Tm_NN（≤2.0℃ 容差）+ GC% + ΔG 复算 | PASS |
| D | 双实现全基因组重扫：grep -F 与 numpy 滚动哈希全长脱靶计数完全一致；48 条引物全长脱靶 = 0、每条全基因组完整匹配恰好 1 次 | PASS |
| E | UCSC RefSeq vs Ensembl 注释交叉核对 | PASS |

### 第二轮：交付数据从头复核（`results/data_review.md`，24/24 行通过）

复核脚本 `scripts/12_review_dataset.py` **不引用设计流水线任何代码**，直接从参考序列
重算：引物结合位点（精确搜索+唯一性）、产物坐标、重复结构（逐 unit 复算 + 逆推 walk
边界复现）、Tm/GC 复算、与 numpy 重扫交叉核对。5 个位点另做 unit 级深挖
（FMR1、chr12、CARM1、GIGYF1、CCDC9）。

**两轮校验发现并已修正的问题：**
1. tract 外推记账 bug：中断单元计数与边界推进不一致，导致拷贝数与 tract 长度不自洽
   （如 CARM1 旧值 26 单元/72bp 互相矛盾）；已修复并全面加硬性不变量断言。
2. primer3 `PRIMER_RIGHT` 坐标语义误解：产物长度曾虚增 ~21bp；实证后改为 `rp+1-lp`。
3. grep -o 非重叠匹配语义少计 15-mer 命中；权威计数改为 numpy 按位置计数。
4. 校验器自身若干口径问题（ΔG 单位 cal/mol、genePred 字段名、GC 显示精度 96→95.5）。
5. **第二轮复核发现**：辅助文件 `specificity.tsv`（旧面板引物）、`genome_scan.tsv`、
   `dense_scan.tsv`（修正前算法）为陈旧版本 → 已用当前面板与修正算法重建
   （13_repeat_scan_clean.py / 03_check_specificity.py）；README 中 FMR1/AFF2 扩展容量
   数字（41/107 → **114/93**）与"参考基因组最长 25、无 >30"的结论过窄 → 已修正为
   三种子阈值并列，并**新发现 chr12:7,628,695（32 单元）真实 >30 位点并入面板**。

## 文件

- `results/primer_design_results.tsv` — 24 对引物完整设计表（28 列权威设计记录，下游校验的输入）
- **`results/primer_panel_delivery.md` / `.csv` — 交付物本身**：面向样本制备组的精简表
  （基因 / 重复区 / 预期扩增子 / 引物+坐标 / Tm F/R / 重复特点），由 `scripts/21_build_delivery.py`
  从上面的 28 列表生成。`.csv` 为机器可读版（引物序列与坐标分列，可直接导入合成订单模板）；
  `.md` 为同一内容的可读版 + 紧凑脚注
- `results/primer_panel_delivery_notes.md` — 交付表的配套说明：Tm 计算条件、PCR 条件建议、
  覆盖范围限制（`(GC)n` 无 ≥20 拷贝档、FMR1/AFF2 的 1 kb 容量上限）、特异性与校验状态、
  优先人工确认位点清单
- `results/validation_report.md` — 16 项校验报告；`results/data_review.md` — 24 行从头复核（含 unit 级深挖）
- `results/designed.json` / `results/candidates.json` — 设计中间产物
- `results/specificity.tsv` — grep 版命中统计（排序用，非重叠语义；**本包内该文件的 `end15_*`
  两列由带缺陷的旧版 `03_check_specificity.py` 生成，右引物计数偏低，勿引用；重跑修复后的脚本可重建**）；
  `results/method_D_genomescan.json` — numpy 版权威计数
- `results/annotation.tsv` — 位点基因注释
- **第 3 版新增的验证产物**（详见 `results/final_report.md` §6.4、§14）：
  - `results/dimer_check.tsv` / `.json` — dimer / hairpin **双重**复核：① 直查 primer3 `_TH` 族
    （设计真正约束的量，上限 47.0）无越界，全场最大 45.73；② 独立 NN 复算（不依赖 primer3）
    作为物理量参考，与 primer3 跨引擎 r = 0.936。合计 24 对 21 ok / 3 warn / 0 fail；
  - `results/mismatch_pcr.tsv` / `.json` — 错配容忍的**配对** in-silico PCR（三档参数，
    3.088 Gbp / 24 对 / 最宽松档 20,626 结合位点 / **0 可扩增脱靶**）；
  - `results/reference_check.json` — 参考基因组与 UCSC hg38 的身份校验（长度 + 抽样窗口 + 24 位点）
- `data/genome_scan.tsv` — 完美串联扫描（GC ≥10 / GCC ≥15 单元）；`data/interrupted_scan.tsv` — 常规种子（≥10 单元）含中断扫描；`data/dense_scan.tsv` — 密扫（≥8 单元，含 chr12 >30 位点）；
  `data/repeat_extremes.json`、`data/sensitivity_15bp.json`、`data/sensitivity_dinuc.json` — 极值/敏感性汇总
- **参考基因组**：`data/hg38_oneline/{chrom}.fa` — 25 条主染色体 / 3.088 Gbp，
  来源为 UCSC `hg38.p14.chromFa`（patch 以 `_fix`/`_alt` 独立序列附加，主染色体与 hg38 逐位相同）。
  由 `scripts/00_fetch_hg38.py` 获取组装（并行分块 + md5 + 断点续传，含 UCSC REST API 备用通路），
  由 `scripts/20_verify_reference.py` 对 UCSC 校验。
- `scripts/` — 流水线：`00_fetch_hg38.py`（参考获取）→ `13_repeat_scan_clean.py` / `08` / `08b` /
  `14` / `15`（位点扫描）→ `09_build_panel.py` → `02b_design_panel.py` → `03b_pick_specific.py` →
  `04_annotate.py` → `11_validate_genome_scan.py` → `05_build_report.py`（28 列权威表）→
  **`21_build_delivery.py`（交付表，流水线末端）**；校验：`10_validate.py`、`12_review_dataset.py`、
  `16_release_gate.py`（发布门禁，含 F9 对交付表坐标的基因组复算）、`17_dimer_check.py`、
  `18_mismatch_pcr.py` + `19_validate_mismatch.py`、`20_verify_reference.py`、
  `03_check_specificity.py`（grep 版计数，仅排序用）
  （01/02/06/07/08/08b 为早期探索脚本，保留备查）。
  所有脚本的 `BASE` 由 `__file__` 推导（可用环境变量 `PRIMER_DESIGN_BASE` 覆盖）。

## PCR / 验证注意事项（GC 富集重复）

1. **聚合酶与缓冲液**：GC 富集模板建议用高保真 GC-buffer 聚合酶（Q5 GC Enhancer、
   Phusion GC Buffer、KAPA2G Robust 等），尤其 chr12/GIGYF1/CARM1 这类长 GC 重复位点。
2. **添加剂**：DMSO 3–5% 或 betaine 0.5–1.5 M 常可改善解链；可与 7-deaza-dGTP（3:1）合用。
   建议按 0/2/3/5% DMSO 做梯度优化，而非一次性叠加多种添加剂——不同聚合酶/缓冲液差异很大。
3. **循环策略**：98℃ 变性、退火 = Tm − 3℃ 起步做梯度/touchdown（如 65→58℃）；
   延伸按产物长度常规设置（均 <1 kb）。
4. **滑链（stutter）**：GC 富集串联重复在 PCR 与 Sanger 中**易出现** ±1 单元 mixed peak（(GC)n 二核苷酸
   最严重），定长以主峰/峰包中位为准、正反双向测序；精确定长优先荧光毛细管电泳（GeneScan）。
5. **拷贝数定标**：本表"预期扩增片段"按 **GRCh38 参考序列中该位点的拷贝数**计算（已独立复核；
   注意 GRCh38 是单一单倍型，不代表人群典型等位，repeat 位点高度多态）。实际样本
   allele = 预期产物 ± 拷贝数差 × 单元长（GCC 家族 3bp、GC 家族 2bp）。含中断位点的
   单元数存在 ±1–2 的定义敏感性（表中已给出精确/中断/纯度）。
6. **扩展位点（FMR1/AFF2）**：GRCh38 参考序列中为 22/21 拷贝。FMR1 的生物学核心结构为
   (CGG)₁₀-AGG-(CGG)₉（19 个 CGG + 1 个 AGG 中断 = 20 单元）；本表按 walk tract 记账为
   22 单元（19 精确 + 3 中断，两端各多计 1 个边界单元），两套口径**不是同一个计数**。
   患者扩展等位每 +1 拷贝产物 +3 bp；1 kb 产物理论上限 FMR1 ~114 拷贝、AFF2 ~93 拷贝，
   **故 1 kb 只能覆盖前突变 55–200 的下段，AFF2 1 kb 更无法触及 >200 的扩增区间**；
   若要覆盖 200 拷贝需 ~1.26 kb（FMR1）/ ~1.32 kb（AFF2）。若放宽至 2 kb
   理论上限分别 ~447/~427 拷贝，但长 GC 产物扩增效率需实验验证；>200 拷贝全突变按文献
   改用 MP/Southern/长读长或 AGG 分布测定。
7. **验证建议**：先以参考 DNA（如 NA12878）确认单一条带与预期长度，再引入已知拷贝数
   差异的对照；建议扩增后测序验证侧翼-重复连接处。本表引物已做全基因组唯一性筛选
   （双实现，全长脱靶 = 0），仍建议实验确认。
8. **合成阳性对照**：若需要 25–30 / >30 各档或 ≥20 拷贝 (GC)n 的可扩展标准品，建议用本表
   引物侧翼序列 + 指定拷贝数重复区合成 gBlock/质粒（弥补参考基因组缺失档位的标准做法）。
