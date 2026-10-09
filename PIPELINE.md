# primer_design 流水线说明

本文说明本项目（GC / GCC 重复序列侧翼 PCR 引物面板，GRCh38）的**主流程**，以及**每个核心步骤的实现方式**。

- 交付物：`results/primer_design_results.tsv`（24 对引物 / 48 条）+ 三份报表（`validation_report.md`、`data_review.md`、`final_report.md`）
- 权威脚本链：`09 → 02b → 03b → 04 → 11 → 05 → 10`，另有 `12`、`16` 两道独立复核/门禁
- 扫描类脚本：`08`（常规种子，权威）、`08b`（密扫种子）、`13`（修正计数后的干净重扫）、`14/15`（更低种子敏感性）
- 早期探索脚本（保留备查，已被取代）：`01`、`02`、`06`、`07`、`03`

---

## 0. 总览

```
┌─ 阶段 1 位点发现（序列级全基因组扫描）──────────────────────────────┐
│  08_interrupted_scan.py   种子 [GC]{30,}   → data/interrupted_scan.tsv │
│  08b_dense_scan.py        种子 [GC]{24,}   → data/dense_scan.tsv       │
│  13_repeat_scan_clean.py  修正计数后重扫   → data/genome_scan.tsv      │
│                                            data/repeat_extremes.json   │
│  14/15_sensitivity_*.py   15bp / 12bp 种子 → data/sensitivity_*.json   │
└───────────────────────────────────────────────────────────────────────┘
                              ↓
┌─ 阶段 2 面板选定 ─────────────────────────────────────────────────────┐
│  09_build_panel.py     扫描结果 + 硬编码 GC 位点 → data/panel_loci.json│
└───────────────────────────────────────────────────────────────────────┘
                              ↓
┌─ 阶段 3 引物设计 ─────────────────────────────────────────────────────┐
│  02b_design_panel.py   primer3-py，逐位点多档产物长度 → candidates.json│
└───────────────────────────────────────────────────────────────────────┘
                              ↓
┌─ 阶段 4 特异性择优 ───────────────────────────────────────────────────┐
│  03b_pick_specific.py  grep -oFf 全基因组 → results/designed.json      │
└───────────────────────────────────────────────────────────────────────┘
                              ↓
┌─ 阶段 5 注释 ─────────────────────────────────────────────────────────┐
│  04_annotate.py        Ensembl REST → results/annotation.tsv           │
└───────────────────────────────────────────────────────────────────────┘
                              ↓
┌─ 阶段 6 权威特异性计数 ───────────────────────────────────────────────┐
│  11_validate_genome_scan.py  numpy 滚动哈希 → method_D_genomescan.json │
└───────────────────────────────────────────────────────────────────────┘
                              ↓
┌─ 阶段 7 交付表 ───────────────────────────────────────────────────────┐
│  05_build_report.py      → results/primer_design_results.tsv（28 列）   │
│  21_build_delivery.py    → results/primer_panel_delivery.md/.csv ← ★交付│
└───────────────────────────────────────────────────────────────────────┘
                              ↓
┌─ 阶段 8 三重校验 ─────────────────────────────────────────────────────┐
│  10_validate.py      16 项方法校验 → validation_report.md              │
│  12_review_dataset.py P1–P5 从头复核 → data_review.md                  │
│  16_release_gate.py   F1–F8 发布门禁（非零退出即失败）                  │
└───────────────────────────────────────────────────────────────────────┘
```

所有脚本的 `BASE` 均由 `__file__` 推导（`os.path.dirname(os.path.dirname(os.path.abspath(__file__)))`），可用环境变量 `PRIMER_DESIGN_BASE` 覆盖。

**外部依赖与数据**：`primer3-py`、`numpy`、`biopython`（仅校验用）、GNU `grep`、网络（Ensembl REST、UCSC API）。参考基因组 `data/hg38/chromFa/*.fa`（UCSC hg38 主装配 + alt/patch contig）与扁平化的 `data/hg38_oneline/*.fa`。**注意 `data/` 被 `.gitignore` 忽略，不随包交付**，重跑全流程需先自行下载 hg38。

---

## 阶段 1：位点发现（序列级扫描）

### 为什么不用 TRF 注释轨

`01_scan_repeats.py` 最初走的是 UCSC `simpleRepeat`（TRF）轨道：对 24 条主染色体逐条调 `getData/track` REST 接口，过滤 `copy_num >= 20 且 perMatch >= 90`，把 `period==2 且字符集 ⊆ {G,C}` 记为 `GC` 家族、`period==3` 同理记为 `GCC` 家族，重新按 `copies = (end-start)/len(seq)` 复算拷贝数并打档位标签，输出 `data/candidates.tsv`。

**问题**：TRF 用打分阈值筛选，短 GC 重复会被整段漏掉。因此后续改为**直接扫描序列**。

### 06_genome_scan.py — 完美串联扫描

```python
DINUC_PAT  = re.compile(r"(?:GC){20,}")
TRINUC_PAT = "|".join(f"(?:{u}){{20,}}" for u in ["GCC","CCG","CGC","GGC","CGG","GCG"])
```

- **双链扫描**：反向链命中用 `gs, ge = len(s)-m.end(), len(s)-m.start()` 映回正链坐标。
- **合并重叠命中**：`gs <= merged[-1][1]` 即合并（含首尾相接），因为同一串联区会被多个阅读框各命中一次。
- **择优**：对每个合并区同时按二核苷酸、三核苷酸两个家族算 `analyze_tract`（遍历链 × 单元 × `shift`），取纯度高的那个。
- **过滤**：硬编码 `purity < 80` 丢弃。
- 输出 `data/genome_scan.tsv`。

**局限**（已知）：二核苷酸正则只写了 GC 相位，纯 `(CG)n` 要 21 个单元才能被种子命中；扫描包含 alt/patch contig（未过滤 `_`）。

### 08_interrupted_scan.py — 含中断的 walk 扫描（权威算法）

这是整条流水线**最核心的算法**，其余扫描脚本都是它的变体。

**种子与阅读框**：`CORE_PAT = re.compile(r"[GC]{30,}")`，即至少 30 bp 的纯 G/C 核。在核内用 `best_frame_core` 遍历「3 单元/2 单元 × 双链 × 6 种旋转 × 3 个相位」，取完美匹配 bp 数最多者为阅读框，要求 `perfect_bp >= 30`。

**外推（walk）**——定义于 [08_interrupted_scan.py:48](scripts/08_interrupted_scan.py#L48)，其 docstring 明确写有"停止时失败单元**不消费**"与记账不变量：

```python
def walk(seq, i, unit, direction, max_streak=2):
    # 逐单元前进：等于 unit 记 exact，否则记中断
    # 连续中断达到 max_streak(2) → 停
    # 已走 >= 4 单元且 exact/已走 < min_purity(0.6) → 停
    # 关键：触发停止的那个不精确单元「不消费」(ni -= 1, 不推进 edge)
```

从种子两端各走一次（`direction = -1 / +1`），**中断单元计入 `n_units`**（TRF 式拷贝数口径）。

**记账与不变量**：

```python
n_units = n // ulen + le + li + re_ + ri     # n = 核内完美 bp
exact   = n // ulen + le + re_
purity  = exact / n_units                     # 要求 >= 0.75
assert (tract_e - tract_s) == n_units * ulen  # 硬不变量
```

该 `assert`（[08_interrupted_scan.py:138](scripts/08_interrupted_scan.py#L138)）是修正「中断单元计数与边界推进不一致」bug 后加上的——正是它暴露了 CARM1 旧值「26 单元 / 72 bp」互相矛盾的问题。

**其他**：要求 `n_units >= 16`；去重按 `(chrom, start)` 排序，同染色体上 `start < prev_end + 10` 视为重叠，保留 `n_units` 大者。输出 `data/interrupted_scan.tsv`（列：`chrom start end unit ulen n_units exact_units interrupt_units purity tract_len gc bucket`）。

### 08b_dense_scan.py — 密扫种子

与 08 同构，两处差异，用来捕获常规种子漏掉的位点：

| | 08（常规） | 08b（密集） |
|---|---|---|
| 种子 | `[GC]{30,}` | `[GC]{24,}` |
| `walk` 收尾 | 失败单元不消费 | 额外做一次**尾部中断 trim** |
| 不变量断言 | 有 | 无 |

正因种子放宽到 24 bp，才发现了参考基因组中唯一真实的 `>30` 位点 **chr12:7,628,694-7,628,790（32 单元 = 26 精确 + 6 中断，纯度 81%）**。输出 `data/dense_scan.tsv`。

### 13_repeat_scan_clean.py — 修正计数规则后的重建

第二轮复核发现 `genome_scan.tsv` / `dense_scan.tsv` 是修正前算法生成的陈旧版本，故用本脚本重建。四个部分：

- **Part A** 完美串联扫描：`(?:GC){10,}` / `(?:CG){10,}` 与 6 种三核苷酸 `{15,}`，合并后取「正链 × 各旋转 × 双链」中最长的完美串联；输出**最大完美串联**口径（与 walk 口径不同）。
- **Part B** 二核苷酸含中断 walk：种子为 ≥10 个完美拷贝，从种子两个方向 walk，去重按区间重叠保留单元数最多者；产出 `(GC)n ≥20` 的盘点（结论：**不存在**）。
- **Part C** 密集三核苷酸：`CORE = re.compile(r"[GC]{24,}")`，核内按 `shift` 找连续完美段 `n >= 24`，两端 walk 外推，`purity >= 0.75`、`n_units >= 16`，并带 `assert (tract_e - tract_s) == n_units * 3`；按 `start < prev_end + 10` 去重。重建 `data/dense_scan.tsv`。
- **Part D** 对 `data/interrupted_scan.tsv` 逐行验不变量 `(end - start) == n_units * ulen`。

汇总写 `data/repeat_extremes.json`。

### 14 / 15 — 更低种子的敏感性检查

- `14_sensitivity_seed.py`：三核苷酸，种子 `[GC]{15,}`，只扫正链（GCC 家族 6 种旋转在反向互补下封闭），同样 `walk` + `purity >= 0.75` + `n_units >= 16`，另加边界断言；输出 `data/sensitivity_15bp.json`。
- `15_sensitivity_dinuc.py`：二核苷酸，种子 `(?:GC){6,}` / `(?:CG){6,}`（12 bp），`walk` 的 `ulen` 硬编码为 2，过滤 `units >= 12 且 purity >= 0.75`；输出 `data/sensitivity_dinuc.json`。

**种子逐级放宽（30 → 24 → 15 → 12 bp）的目的是验证「参考基因组中不存在 (GC)n ≥ 20」「>30 档只有 chr12 一个位点」这两条结论不会随阈值变化而翻盘**——实测均未翻盘。

---

## 阶段 2：09_build_panel.py — 面板选定

把扫描结果收敛成一份固定的待设计位点清单 `data/panel_loci.json`（24 条）。

1. **读入** `data/interrupted_scan.tsv`。
2. **注入密扫特有行**：从 `data/dense_scan.tsv` 按 `(chrom, start) == ("chr12", 7628694)` 精确匹配，把那条 32 单元的 `>30` 位点手工追加进池子。
3. **扩展位点优先**：`near(t, chrom, lo, hi)` 按**坐标窗口**匹配 FMR1（chrX 147,912,030–147,912,060）与 AFF2（chrX 148,500,610–148,500,650）——用窗口而非精确坐标，是因为不同扫描版本间 tract 边界会漂移几个 bp。这两行打 `expansion_row: True` 与 `bucket_label: ">30（扩展检测）"`。
4. **强制纳入**：HMGB3（chrX 150,983,370–150,983,410）、TMEM185A（chrX 149,631,700–149,631,740）。
5. **填充**：按 `-n_units` 排序，**每条染色体最多 2 个位点**保证多样性，直到 14 个 GCC 位点。
6. **档位标签**：`20 <= n < 25` → `20-25`；`25 <= n < 30` → `25-30`；`n >= 30` → `>30`；`n < 20` → `20-25` 但标 `△ 最接近（n 拷贝）`。
7. **GC 家族 10 个位点硬编码**：参考基因组不存在 ≥20 拷贝的 (GC)n（最长完美串联 14 拷贝，chr19:47,271,732），故直接写入 10 个 12–14 拷贝的最长真实位点，`purity: 100.0`、`bucket_label: "GC家族参考最长档（<20）"`。

---

## 阶段 3：02b_design_panel.py — primer3 引物设计

### 全局参数（`GLOBAL_ARGS`）

| 项 | 值 |
|---|---|
| 长度 | 20–28，最优 22 |
| Tm | 57–64 ℃，最优 60 ℃ |
| GC% | 35–65，最优 50 |
| 盐/浓度 | 单价盐 50 mM、二价 1.5 mM、DNA 50 nM、dNTP 0.6 mM |
| 自我结构（`_TH` 族） | `SELF_ANY_TH` / `SELF_END_TH` / `HAIRPIN_TH` = 47.0 |
| 配对结构（hetero-dimer） | `PAIR_MAX_COMPL_ANY_TH` / `_END_TH` = 47.0 |
| 3' 端稳定性 | `MAX_END_STABILITY` ≤ 6.0（primer3 默认 100） |
| 连续同碱基 | `MAX_POLY_X` ≤ 5（等同默认） |
| 返回数 | `PRIMER_NUM_RETURN = 10` |

> **必须用 `_TH` 族。** primer3-py 默认开启 thermodynamic oligo alignment，此时
> `PRIMER_MAX_SELF_ANY` / `PRIMER_MAX_SELF_END` 这类非 `_TH` 键**不参与判定**；
> `PRIMER_MAX_HAIRPIN` 更是一个**不存在的键名**（只有 `PRIMER_MAX_HAIRPIN_TH`），
> 而 primer3-py 对未知键**静默忽略**。本项目早期版本正是踩了这个坑：
> 三个键全部失效，实际生效的是 primer3 默认值。物理量换算（实测拟合）：
> `SELF_ANY_TH` 47 ≈ 同源二聚体 ΔG −12.5 kcal/mol，`HAIRPIN_TH` 47 ≈ 发夹 ΔG −0.8 kcal/mol。

### 关键实现点

**① 强制夹住重复区**——这是本项目设计的核心约束：

```python
"SEQUENCE_TARGET":          [(tract_offset, tract_len)],  # 扩增子必须包含 tract
"SEQUENCE_EXCLUDED_REGION": [(tract_offset, tract_len)],  # 引物不得落在 tract 内
```

**② 排除区自动加宽** `widen_exclusion(seq, lo, hi)`：向两侧按 2 bp 步长扩张，只要相邻序列仍是 period-2 简单重复（CA/TG/GT 等微卫星），就并入排除区——避免引物落进旁邻微卫星造成非特异。

**③ 模板窗口回退**：依次尝试 `flank = 800 → 500 → 350 → 250`，取 `[start-flank, end+flank)`；过短则跳过。

**④ 产物长度阶梯**（逐档尝试，首个有结果者胜）：

```python
普通位点      : [250,900] → [200,1000] → [200,1600] → [200,2000]
扩展检测行    : [250,400] → [200,550]  → [200,900] → [200,1600] → [200,2000]
```

扩展位点（FMR1/AFF2）从更短的产物起步，是为了给患者长等位留出扩增余量。

**⑤ 产物长度硬上限**：标准路径 `prod > 1000` 丢弃；只有走了 GC 兜底路径才允许到 2000。

**⑥ GC 兜底**：若所有档位都无结果，放宽为 `MAX_GC 75 / MIN_GC 30 / MIN_TM 56 / MAX_POLY_X 6` 重试一遍。

**⑦ 候选收集**：所有档位的所有对都收进来，按 `PRIMER_PAIR_i_PENALTY` 升序，按 `(left_seq, right_seq)` 去重，返回前 10 对。

### 坐标约定（实证确定的坑）

primer3 的 `PRIMER_RIGHT` 返回值语义容易搞错，代码里写明了实测结论：

```python
# LEFT  (start,len): start = 左引物 5' 端，0-based
# RIGHT (start,len): start = 右引物「最右侧碱基」即 5' 端，0-based
#                     右引物 3' 端位于 rp - len + 1
prod          = rp + 1 - lp
amplicon_start = t_start + lp
amplicon_end   = t_start + rp + 1     # 半开
```

早期误用 `rp + len - lp` 导致产物长度虚增约 21 bp，已在复核中修正。

输出 `results/candidates.json`（结构 `{"locus": {...}, "candidates": [...], "status": "ok"}`）。

---

## 阶段 4：03b_pick_specific.py — 全基因组特异性择优

### 扫描方式

把**所有候选引物**（每个位点最多 10 对 × 2 条）的全长序列、反向互补、3' 端 15-mer 及其反向互补全部写进一个 pattern 文件，对 `data/hg38_oneline/*.fa`（每染色体一行）跑**一次** `grep -oFf`，流式解析输出行 `文件名:匹配串`，同时维护两份计数：

- `counts[match]` —— 全基因组命中总数
- `by_chrom[(match, chrom)]` —— 分染色体命中数（用于扣除自身位点所在染色体）

### 打分规则

对每个候选对：

```python
full_off = (正反链全长命中 − 同染色体命中)  # 左右引物求和
end15_off = (3' 15-mer 命中 − 同染色体命中)

prod_pen = product_size / (400 if 扩展行 else 1000)

score = min(full_off, 5) * 400        # 全长脱靶：一票否决级权重
      + min(end15_off, 60) * 8        # 3' 端脱靶：次级权重
      + PRIMER_PAIR_PENALTY           # primer3 自身代价
      + prod_pen * 20                 # 产物越短越好（扩展行权重加倍）
```

`min(..., cap)` 是为防止极个别候选把分值拉到无穷大而掩盖其他因素。

### 定级

`full_off > 0 或 end15_off > 10` → `ok_check`；否则 `ok`。最终交付 24 行中 **18 行 `ok` + 6 行 `ok_check`**，全部 `full_off == 0`。

输出 `results/designed.json`。

> **注意**：这里的 15-mer 计数用 `grep -o`，只报**非重叠**匹配，在 `(CA)n` 这类区域会少计。因此它**只用于排序**，不作为权威计数——权威计数见阶段 6。

---

## 阶段 5：04_annotate.py — 基因注释

- 接口：Ensembl REST `GET /overlap/region/human/{chrom}:{start}-{end}?feature=gene`，`retries=3`、`timeout=30`、间隔 `sleep(0.3)`。
- 坐标：内部 0-based 半开 → Ensembl 1-based 闭区间（`s1, e1 = start + 1, end`）。
- **有重叠基因**：优先 `protein_coding`，同类中取重叠长度最小者（`gkey` 返回 `(biotype 优先级, 重叠长度)`）。
- **无重叠基因**：向两侧扩 100 kb 再查，取距离最近者，标 `intergenic(upstream/downstream Nbp of GENE)`。
- 输出 `results/annotation.tsv`。

---

## 阶段 6：11_validate_genome_scan.py — 权威全基因组计数（方法 D）

完全独立于阶段 4 的 grep 实现（不同工具、不同算法），用于交叉验证。

**算法**：

1. **base-4 编码**：`CODE[256]` 查表把 A/C/G/T 映射到 0–3，其他字符映射到 4。15-mer 编码为 uint32（`4^15 < 2^32`，**单射**，无碰撞）。
2. **查询表**：对 48 条交付引物各取 3' 端 15-mer，正反两向编码，`qmap: code → [(引物索引, 方向)]`；排序成 `codes_arr`。
3. **逐 contig 分块扫描**（`CHUNK = 20_000_000`）：
   - 向量化 K 次乘加得到每个窗口的编码：`code = code * 4 + c[j:j+n]`；
   - 用 `cumsum` 做**有效性掩码**，排除任何跨越非 ACGT 字符的窗口；
   - `np.searchsorted` + 比较定位命中（O(n log m)）；
4. **命中验证与延伸**：每个命中按方向回推全长起始位置，做**逐字符 `==` 比对**确认是否构成全长引物匹配。
5. **按基因组位置计数**（重叠出现也各计一次），这是与 `grep -o` 的语义差异所在。

输出 `results/method_D_genomescan.json`：`{scanned_files, scanned_bp, elapsed_s, rows:[{no, id, chrom, fwd:{genome_full_hits, full_off, end15_hits, end15_off}, rev:{...}}]}`。

**结果**：扫描 455 contig / 3.21 Gbp，48 条引物全长命中各恰好 1 次（自身位点），全长脱靶 0。

---

## 阶段 7：05_build_report.py → 21_build_delivery.py — 表生成

合并三个来源：

| 来源 | 提供 |
|---|---|
| `results/designed.json` | 引物序列、坐标、Tm、GC、产物长度、`ok/ok_check` |
| `results/annotation.tsv` | 基因符号、biotype、位置上下文 |
| `results/method_D_genomescan.json` | **权威**的 `*_full_off` / `*_end15_off` / `genome_full_hits` |

要点：

- **坐标 1-based 化**：`repeat_start = locus.start + 1`、`repeat_end = locus.end`、`amplicon_start = pair.amplicon_start + 1`、`amplicon_end = pair.amplicon_end`。
- **`repeat_feature` 特征串**：拼成中文字符串供人读，同时被后续校验脚本用正则反解——例如 `GCC家族·(CGG)n×22（精确19+中断3，纯度86.4%），66bp，GC 100.0%`。校验脚本用 `×(\d+)`、`精确(\d+)\+中断(\d+)`、`纯度([\d.]+)%`、`），([\d.]+)bp`、`GC ([\d.]+)%` 解析它。
- **扩展位点容量计算**：`head = (1000 - (product_size - tract_len)) // 3`，即 1 kb 产物扣掉侧翼后理论可容纳的拷贝数（FMR1 ~114、AFF2 ~93）。
- **用方法 D 覆盖 grep 计数**，并**重算 flag**（`pair_full_off > 0 或 pair_end15_off > 10` → `ok_check`）。
- 与 method_D 对齐时带 `assert row["no"] == drow["no"]`。

输出 `results/primer_design_results.tsv`（**内部权威设计记录**，28 列）。

### 21_build_delivery.py — 交付表（流水线末端）

把上面 28 列的表压缩成样本制备组使用的 8 列交付格式，输出
`results/primer_panel_delivery.md`（表格 + 紧凑脚注）与 `.csv`（同内容，引物序列与坐标分列，
可直接导入合成订单模板）。

**两项派生量在这里现算**（上游任何产物里都没有）：

```python
# 引物结合坐标（1-based 闭区间，hg38 正链）
fwd_start = amplicon_start
fwd_end   = amplicon_start + len(fwd) - 1
rev_start = amplicon_end - len(rev) + 1
rev_end   = amplicon_end
```

- `Repeats 特点` 由上游 `repeat_feature` 重新压缩：`(CGG)n 重复，hg38 参考等位 22（精确 19 + 中断 3，纯度 86.4%）`。
  **与上游的 `repeat_feature` 是两套格式**——后者被 `10/12/16` 用正则反解，两者分开维护。
- 无基因符号的位点（#3、#19、#21，均为 lncRNA 基因体内）加 `※` 脚注；假基因（#11 `PCDH8P1`）加 `†`。
- 脚本内自带**排序不变量**断言：`amplicon_start ≤ fwd_start ≤ fwd_end < repeat_start ≤ repeat_end < rev_start ≤ rev_end ≤ amplicon_end`，违反即非零退出。
- 交付表的这两列坐标**由 `16_release_gate.py` 的 F9 直接从参考基因组复算核对**，不依赖上游产物——
  即坐标错、几何不自洽、序列与坐标不符都会被门禁拦下（已做负向测试验证 F9 非恒真）。

---

## 阶段 8：三重校验

三层校验**互相独立**，且都不复用设计代码。

### 10_validate.py — 方法校验（16 项）

| 方法 | 实现要点 |
|---|---|
| **X** 内部一致性 | TSV 与 `designed.json` 行数/序列一致；`product_size == amplicon_end − amplicon_start + 1`；重复区长度与特征串拷贝数自洽；档位标签与拷贝数一致；`>30 ✓` 行单元数 ≥ 30；扩展检测行恰为 FMR1 + AFF2 |
| **A** 参考基因组模拟 PCR | **不经 primer3**：直接从 FASTA 取 `[amplicon_start, amplicon_end)`，断言 `amp.startswith(fwd)`、`amp.endswith(rc(rev))`、各自在扩增子内恰好出现 1 次、长度相符、重复区被完整包含、引物与重复区不重叠、引物只含 ACGT |
| **B** 重复 tract 独立重分析 | 用**两种与设计算法不同**的方法复算拷贝数：① 正则找最大连续精确段 `best_run`；② 逐相位对齐计数 `phase_max`。断言 `best_run <= exact <= phase_max`、`exact + interrupts == units`、纯度自洽、GC 行确为完美串联 |
| **C** 热力学 | `primer3.calc_tm`（同引擎，容差 0.3 ℃，实测最大偏差 0.046 ℃）+ `Bio.SeqUtils.MeltingTemp.Tm_NN`（跨引擎，容差 2.0 ℃，实测 1.81 ℃）+ GC% 复算 + homo/hetero/hairpin ΔG 复算（阈值 −45000 cal/mol） |
| **D** 全基因组双实现对照 | 读 `method_D_genomescan.json`，断言 grep 与 numpy 的**全长脱靶完全一致**，且 `numpy_end15 >= grep_end15`（方向性验证非重叠语义） |
| **E** 注释交叉核对 | 调 UCSC `ncbiRefSeqCurated` 轨道（**只查位点 ±50 kb 邻域**，避免整染色体转储），与 Ensembl 基因符号比对；区分「直接匹配 / 转录本边界差异 / Ensembl 独有注释（RefSeq 无模型）」三类，仅真正的矛盾才失败 |

### 12_review_dataset.py — 从零复核（P1–P5）

显式声明"**不引用设计流水线任何代码**"，从 GRCh38 FASTA 重算每个数字：

- **P1** 用 `str.find` 穷举定位引物结合位点（要求唯一），从结合位点**重新推导**产物坐标与 TSV 比对。
- **P2** 按文档化的计数规则**重新实现 `walk`**，并且**从每一个内部最大精确段出发做边界复现**——要求至少有一个种子能走出完全相同的边界（`l_edge == rs0 and r_edge == re1`）与计数（`ne_l + ne_r + L == exact`、`ni_l + ni_r == inter`）。这样才证明交付的 tract 边界不是某次扫描的偶然产物。
- **P3** 长度/GC/Tm 复算 + 配对 ΔTm ≤ 2.0 ℃。
- **P4** 与 `method_D_genomescan.json` 交叉核对，并**重新施加 flag 规则**验证 `ok` / `ok_check` 标签。
- **P5** 5 个位点（FMR1、chr12 >30、CARM1、GIGYF1、CCDC9）做**单元级深挖**：逐单元列出序列/状态/GC 数（tract 外各多列 4 个单元），并给出 6 种阅读框 × 3 个相位的精确单元计数表。

阅读框来源优先取扫描文件（`interrupted_scan.tsv` + `dense_scan.tsv`）的 `unit` 列——因为 chr12 那条只存在于密扫文件，若只看常规扫描会因缺阅读框而把整行每个单元都判为中断。

以 `sys.exit(1 if n_issue else 0)` 返回状态。

### 16_release_gate.py — 发布门禁（F1–F8）

在**干净进程**中跑，触碰每一个交付文件，任一失败即 `exit(1)`。（`results/README.md` 记为"27 项端到端检查"；源码中实际有 25 个 `check()` 调用点，其中 F1 的清单检查按 6 个文件循环，故执行次数取决于运行环境。本机 `data/` 缺失，无法实测该数字。）

- **F1** 文件清单与行数（`designed.json` 24、`candidates.json` 24、`method_D_genomescan.json` 24、`annotation.tsv` 24、`specificity.tsv` 48、`data/panel_loci.json` 24）；文件缺失时**报告而非崩溃**（因为 `data/` 可能不在检出中）。
- **F2** TSV 结构：必需列非空、48 条引物两两不同、长度 20–28 且全为 ACGT、坐标为整数且 `product == 跨度`。
- **F3/F4** 逐行新读基因组验证结合位点/产物跨度/包含关系；GC 行验完美串联**且两侧不可延伸**（最大性）；GCC 行验 `精确+中断 == 单元数`。
- **F5** 五份文件交叉一致，且 `genome_full_hits == 1`。
- **F6** 设计约束：产物 ≤ 1 kb、Tm ∈ [57,64]、ΔTm ≤ 2.0、GC ∈ [35,65]、全长脱靶合计 = 0、flag 分布恰为 18 `ok` + 6 `ok_check`。
- **F7** 覆盖度：GCC 完全符合档位 ≥ 7 行、GC 家族 10 行、扩展检测 2 行、chr12 `>30 ✓` 行存在且为 32 单元。
- **F8** 报表自洽：`validation_report.md` 含"16 项通过，0 项失败"且无 `FAIL`；`data_review.md` 含"24 行中 **24 行**通过"且无 `FAIL`；README 数字同步；8 个权威脚本存在。

---

## 核心算法小结

### `walk` —— 含中断串联的外推规则

```
输入：序列、起点、单元、方向
状态：ne(精确数) ni(中断数) streak(连续中断计数) edge(边界)
循环：
  取下一个单元
    匹配 unit      → ne++, streak = 0
    不匹配         → ni++, streak++
                     streak >= 2                        → 回退该单元，停
                     ne+ni >= 4 且 ne/(ne+ni) < 0.6     → 回退该单元，停
    edge 前进一个单元长度
输出：ne, ni, edge
```

**两个必须一致的口径**（不一致就是 bug）：

1. `edge` 前进的碱基数 **严格等于** `(ne + ni) × ulen`；
2. 触发停止的那个单元**不计入也不推进**（`ni -= 1` 后再 `break`）。

`assert (tract_e − tract_s) == n_units × ulen` 就是这条不变量的机器化检查。

### 两套并存的重复口径

| 口径 | 定义 | 用在哪 |
|---|---|---|
| **最大完美串联** | 纯完美重复，两侧不能再延伸一个单元 | GC 家族（二核苷酸）10 行；`genome_scan.tsv` |
| **含中断 walk** | 允许 ≤2 连续非精确单元，总纯度 ≥75%，中断计入单元数 | GCC 家族（三核苷酸）14 行；`interrupted_scan.tsv` / `dense_scan.tsv` |

同一 tract 在两套口径下计数可以不同（例如 FMR1：文献口径 20 单元、walk 口径 22 单元），文档与报表中已分别标注。

### 三个关键实证结论

1. **(GC)n ≥ 20 在参考基因组中不存在**。完美口径最长 14（chr19:47,271,732），含中断口径最长 18；种子降到 12 bp 仍是 0。
2. **GCC 家族 `>30` 只有一个真实位点**：chr12:7,628,694（32 单元 = 26 精确 + 6 中断，纯度 81%），只有把种子放宽到 24 bp 才扫得到。
3. **参考基因组是"健康基线"单倍型**，患者型长等位（FMR1 >55、AFF2 >200）不在装配中，故 FMR1/AFF2 只能作扩展检测设计。

---

## 脚本总览

| 脚本 | 输入 | 输出 | 角色 |
|---|---|---|---|
| `01_scan_repeats.py` | UCSC REST | `data/candidates.tsv` | 早期：TRF 轨道（已被取代） |
| `02_design_primers.py` | `candidates.tsv` | `designed.json` | 早期：**已损坏/被 02b 取代**，见下 |
| `02b_design_panel.py` | `data/panel_loci.json` | `results/candidates.json` | **引物设计（权威）** |
| `03_check_specificity.py` | `designed.json` | `results/specificity.tsv` | grep 版计数（仅排序用） |
| `03b_pick_specific.py` | `candidates.json` | `results/designed.json` | **特异性择优（权威）** |
| `04_annotate.py` | `designed.json` | `results/annotation.tsv` | **基因注释** |
| `05_build_report.py` | designed + annotation + method_D | `primer_design_results.tsv` | 28 列权威设计记录 |
| `21_build_delivery.py` | `primer_design_results.tsv` | `primer_panel_delivery.{md,csv}` | **交付表（流水线末端）** |
| `06_genome_scan.py` | `chromFa/*.fa` | `data/genome_scan.tsv` | 早期：完美串联扫描 |
| `07_availability.py` | `chromFa/*.fa` | 仅 stdout | 早期：可用性分布 |
| `08_interrupted_scan.py` | `chromFa/*.fa` | `data/interrupted_scan.tsv` | **含中断扫描（权威）** |
| `08b_dense_scan.py` | `chromFa/*.fa` | `data/dense_scan.tsv` | **密扫（24 bp 种子）** |
| `09_build_panel.py` | scan TSV | `data/panel_loci.json` | **面板选定** |
| `10_validate.py` | 交付表 + designed | `validation_report.md` | 方法校验 16 项 |
| `11_validate_genome_scan.py` | designed + 基因组 | `method_D_genomescan.json` | **权威全基因组计数** |
| `12_review_dataset.py` | 交付表 + 基因组 | `data_review.md` | 从头复核 P1–P5 |
| `13_repeat_scan_clean.py` | `chromFa/*.fa` | 重建 genome/dense scan + `repeat_extremes.json` | **修正计数后的重扫** |
| `14_sensitivity_seed.py` | `chromFa/*.fa` | `data/sensitivity_15bp.json` | 15 bp 种子敏感性 |
| `15_sensitivity_dinuc.py` | `chromFa/*.fa` | `data/sensitivity_dinuc.json` | 12 bp 种子敏感性 |
| `16_release_gate.py` | 全部交付文件 | stdout + 退出码 | **发布门禁 F1–F8** |

---

## 已知问题与注意事项

1. **`02_design_primers.py` 已损坏且与下游不兼容**，已被 `02b` + `03b` 取代：
   - [02_design_primers.py:148](scripts/02_design_primers.py#L148) 调用 `analyze_tract(tmpl, tract_offset, tract_len, fam)`，但 `fam` 在 `design_one` 作用域内**从未定义**（它只是 `canonical_motif` 内的局部变量，`main` 里在第 254 行另有一个同名局部）。这会在每个位点抛 `NameError`，被 `main` 的 `try/except Exception` 吞掉，结果是**全部位点 `status = "FAILED: NameError(...)"`**，脚本无法产出任何成功设计。
   - 它写出的是 `{"cand": ..., "pair": ...}` 结构，而下游 `03b` / `04` / `05` / `11` 一律读 `item["locus"]`，schema 不兼容。
   - 其 docstring 声称还会写 `results/loci_designed.tsv`，实际 `main()` 从未写该文件。

2. **`03_check_specificity.py` 的 3' 端 15-mer 键**：必须用 `seq_f[-15:]` 与 `rc(seq_f[-15:])` 组 pattern；早期误用 `seq_r[-15:]`，导致查询的是"5' 端 15-mer 的反向互补"（该串从不存在于 pattern 集合中），**右引物的 3' 端命中数被静默记成 0**。现版本已修复（[03_check_specificity.py:105-110](scripts/03_check_specificity.py#L105-L110) 有说明）。但 `results/specificity.tsv` 是旧版产物，其 `end15_*` 两列**不可引用**，重跑脚本可重建。

3. **`grep -o` 的非重叠语义**：只统计非重叠匹配，在串联重复区会少计。因此 grep 版计数只用于**排序**，一切对外数字以 `method_D_genomescan.json`（numpy 按位置计数）为准。核对关系是 `numpy_end15 >= grep_end15`，方向必须成立。

4. **`data/` 被 `.gitignore` 忽略，不随仓库交付**。本地已重新下载并重建（`data/hg38/chromFa` 与 `data/hg38_oneline` 各 25 条主染色体，3.088 Gbp；扫描产物由 `08`/`13`/`14`/`15` 重建）。但只要 `data/` 不在版本控制内，新克隆的检出仍然无法重跑 `02b` 之后的任何阶段。`16_release_gate.py` 与 `12_review_dataset.py` 对此做了降级处理（报告缺失而非崩溃），受影响的检查项会记为失败。

   > **参考口径差异（重跑时必须知道）**：当前 `data/hg38_oneline` 是 **25 条主染色体 / 3.088 Gbp**，
   > 而历史上 `method_D_genomescan.json` 曾按 **455 contig / 3.209 Gbp**（含 alt/patch/unplaced 支架）
   > 扫描。两者相差约 121 Mb，会使 3' 15-mer 计数出现 ±1 量级的差异（全长精确匹配数不变）。
   > 交付表的历史 `end15_*` 数值来自 455-contig 口径，重跑后的数值来自 25-染色体口径，
   > **两者不可直接比较**；引用时应以同一次运行的口径为准。参考身份已由
   > `20_verify_reference.py` 对 UCSC hg38 校验通过。

5. **`08` 与 `08b` 的 `walk` 实现不完全相同**（08b 缺"失败单元不消费"的记账、另加尾部 trim、且没有不变量断言），因此两者对同一 tract 的中断计数可能不同。权威的 `interrupted_scan.tsv` 来自 08；`dense_scan.tsv` 已由 `13` 用统一规则重建。

6. **`08` 内部有两个不同的精确度阈值**：`walk` 里中止用的 `min_purity = 0.6`，与最终过滤用的 `purity >= 0.75`，二者含义不同，不要混淆。

7. **坐标约定**：所有内部输出为 0-based 半开；对外交付表（TSV）为 1-based 闭区间；`14/15` 仅在打印到 stdout 时用 `start + 1`。

8. **交付表的 `repeat_feature` 是中文特征串**，被多个校验脚本用正则反解。改动其格式会连带打破 `10_validate.py`、`12_review_dataset.py`、`16_release_gate.py` 的解析。

9. **`02b_design_panel.py` 中三个函数是死代码**：`pick_best()`（[:179](scripts/02b_design_panel.py#L179)）无任何调用点，`dinuc_run_penalty()`（[:158](scripts/02b_design_panel.py#L158)）只被 `pick_best` 调用，`gc_fallback()`（[:67](scripts/02b_design_panel.py#L67)）也从未被调用（`design()` 在 [:149](scripts/02b_design_panel.py#L149) 附近把同样的字典内联了一遍）。因此那个"惩罚 3' 端 12 bp 内二核苷酸串联"的自定义逻辑**从未生效**，实际用的是纯 primer3 penalty 排序。保留而不删除是因为删改会改变候选排序行为；如需启用，须作为一次有意的行为变更来对待。

10. **`_TH` 参数族问题（已修复，记录备查）**：`02b` 早期设的是 `PRIMER_MAX_SELF_ANY` / `_SELF_END` / `_HAIRPIN`。在 primer3-py 默认的热力学比对模式下，前两者不参与判定，`PRIMER_MAX_HAIRPIN` 甚至不是合法键名，而 primer3-py **静默忽略未知键**——三个约束全部失效，实际生效的是 primer3 默认值（`*_TH` = 47）。现已改为显式设置 `_TH` 族与 `PRIMER_PAIR_MAX_COMPL_*_TH`。修复前后 `candidates.json` **逐字节一致**，证明这是纯粹的配置如实化、无行为变更。
