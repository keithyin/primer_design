#!/usr/bin/env python3
"""Build the sample-prep delivery table -- the panel's final output.

Reads results/primer_design_results.tsv (the authoritative 28-column table) and
emits the reduced delivery format used by the sample-preparation group:

  # | 基因/区域 | 参考位置(重复区, hg38) | 预期扩增子 |
  Forward Primer (5'->3') | Reverse Primer (5'->3') | Tm F/R (C) | Repeats 特点

Two derived quantities are computed here, because neither is stored upstream:

  * **primer binding coordinates**.  The upstream table stores only the amplicon
    span; the primer intervals follow from it and the primer lengths:
      F: [amplicon_start,                        amplicon_start + len(F) - 1]
      R: [amplicon_end - len(R) + 1,             amplicon_end]
    Both are 1-based closed on the hg38 plus strand, matching the IGV/UCSC
    display convention, and are re-derived from the reference genome by
    16_release_gate.py (F9).

  * **a compact `Repeats 特点` string**, e.g.
      (CGG)n 重复，hg38 参考等位 22（精确 19 + 中断 3，纯度 86.4%）
    This is deliberately NOT the upstream `repeat_feature` string, which the
    validators regex-parse; the two are maintained separately.

Writes results/primer_panel_delivery.md (table + footnotes) and
results/primer_panel_delivery.csv (same content, one field per column).
"""
import csv
import os
import re
import sys

BASE = os.environ.get("PRIMER_DESIGN_BASE") or os.path.dirname(
    os.path.dirname(os.path.abspath(__file__)))

SRC = os.path.join(BASE, "results", "primer_design_results.tsv")
OUT_MD = os.path.join(BASE, "results", "primer_panel_delivery.md")
OUT_CSV = os.path.join(BASE, "results", "primer_panel_delivery.csv")

# Loci annotated only by an Ensembl gene ID (no symbol in either Ensembl or
# UCSC RefSeq) get a footnote marker rather than a bare ENSG id in the table.
NO_SYMBOL_FOOTNOTE = "※"
PSEUDOGENE_FOOTNOTE = "†"
_PSEUDOGENE_BIOTYPES = {"transcribed_processed_pseudogene",
                        "transcribed_unprocessed_pseudogene",
                        "processed_pseudogene", "unprocessed_pseudogene",
                        "pseudogene"}

EXPANSION_NOTE = {
    "FMR1": "脆性X综合征位点；患者等位可扩展至 >200 拷贝，hg38 为健康基线",
    "AFF2": "FRAXE 脆性位点；患者等位可扩展至 >200 拷贝，hg38 为健康基线",
}


def comma(n):
    return f"{int(n):,}"


def region(chrom, start, end):
    """hg38 1-based closed interval, thousands-separated, en dash."""
    return f"{chrom}:{comma(start)}–{comma(end)}"


def repeat_note(row):
    """Compact `Repeats 特点` string built from the upstream repeat_feature."""
    feat = row["repeat_feature"]
    m = re.search(r"\(([GC]{2,3})\)n×(\d+)", feat)
    if not m:
        raise ValueError(f"unparsable repeat_feature for #{row['no']}: {feat!r}")
    unit, units = m.group(1), int(m.group(2))
    ex = re.search(r"精确(\d+)\+中断(\d+)", feat)
    if ex:
        exact, inter = int(ex.group(1)), int(ex.group(2))
        purity = re.search(r"纯度([\d.]+)%", feat)
        tail = (f"（精确 {exact} + 中断 {inter}，"
                f"纯度 {purity.group(1)}%）")
    else:
        tail = "（完美串联，纯度 100%）"
    note = f"({unit})n 重复，hg38 参考等位 {units}{tail}"
    extra = EXPANSION_NOTE.get(row["gene"])
    if extra:
        note += f"。{extra}"
    return note


def gene_label(row):
    gene = row["gene"]
    if gene.startswith("ENSG"):
        return gene + NO_SYMBOL_FOOTNOTE
    if row.get("biotype") in _PSEUDOGENE_BIOTYPES:
        return gene + PSEUDOGENE_FOOTNOTE
    return gene


def build_rows():
    with open(SRC) as f:
        src = list(csv.DictReader(f, delimiter="\t"))
    rows = []
    for r in src:
        fwd, rev = r["fwd"], r["rev"]
        a0, a1 = int(r["amplicon_start"]), int(r["amplicon_end"])
        rows.append({
            "no": r["no"],
            "gene": gene_label(r),
            "gene_raw": r["gene"],
            "family": r["family"],
            "chrom": r["chrom"],
            "repeat_start": int(r["repeat_start"]),
            "repeat_end": int(r["repeat_end"]),
            "product_size": int(r["product_size"]),
            "amplicon_start": a0,
            "amplicon_end": a1,
            "fwd": fwd,
            "fwd_start": a0,
            "fwd_end": a0 + len(fwd) - 1,
            "rev": rev,
            "rev_start": a1 - len(rev) + 1,
            "rev_end": a1,
            "fwd_tm": r["fwd_tm"],
            "rev_tm": r["rev_tm"],
            "repeats": repeat_note(r),
        })
    return rows


def write_markdown(rows):
    out = [
        "# GC / GCC 重复序列侧翼 PCR 引物面板 — 交付表",
        "",
        "**参考基因组**：GRCh38 / hg38（UCSC hg38.p14 主染色体，已对 UCSC 逐位校验）  ",
        f"**交付规模**：{len(rows)} 对引物（"
        f"{sum(1 for r in rows if r['family'] == 'GCC')} 个 GCC 家族位点 + "
        f"{sum(1 for r in rows if r['family'] == 'GC')} 个 GC 家族位点），全部 22 nt  ",
        "**坐标约定**：hg38，1-based 闭区间（与 IGV / UCSC 显示一致）",
        "",
        "| # | 基因 / 区域 | 参考位置（重复区, hg38） | 预期扩增子 | "
        "Forward Primer (5′→3′) | Reverse Primer (5′→3′) | "
        "Tm F/R (°C) | Repeats 特点 |",
        "|---|---|---|---|---|---|---|---|",
    ]
    for r in rows:
        out.append(
            f"| {r['no']} | {r['gene']} "
            f"| {region(r['chrom'], r['repeat_start'], r['repeat_end'])} "
            f"| {r['product_size']} bp "
            f"| `{r['fwd']}` {region(r['chrom'], r['fwd_start'], r['fwd_end'])} "
            f"| `{r['rev']}` {region(r['chrom'], r['rev_start'], r['rev_end'])} "
            f"| {r['fwd_tm']} / {r['rev_tm']} | {r['repeats']} |")
    out += [
        "",
        "**脚注**",
        "",
        f"- **{NO_SYMBOL_FOOTNOTE}** Ensembl 基因 ID；该位点落在 lncRNA 基因体内，"
        "Ensembl 与 UCSC RefSeq 均未给出基因符号。",
        f"- **{PSEUDOGENE_FOOTNOTE}** 已加工假基因（transcribed_processed_pseudogene）。",
        "- 反向引物列为**直接送合成的寡核苷酸序列（5′→3′）**；其右侧坐标是该引物在 "
        "hg38 **正链**上的结合区间（即扩增子末端 `len(R)` 个碱基处），不是引物自身的方向。",
        f"- Tm 由 primer3 2.3.1 按 SantaLucia 最近邻模型计算，条件为 50 mM 单价盐 / "
        f"1.5 mM Mg²⁺ / 50 nM 引物 / 0.6 mM dNTP；Biopython `Tm_NN` 独立复算最大偏差 1.81 °C。"
        "**与其他文献面板的 Tm 不可直接比较**（盐条件不同可差 5–10 °C）。",
        "- 面板的覆盖范围限制（`(GC)n` 档参考基因组不存在 ≥20 拷贝、`>30` 档仅 1 个真实位点、"
        "FMR1/AFF2 在 1 kb 产物下的扩增容量上限）、PCR 条件建议与优先确认位点，"
        "见 `results/primer_panel_delivery_notes.md`。",
        "",
    ]
    with open(OUT_MD, "w") as f:
        f.write("\n".join(out))


CSV_COLS = ["no", "gene", "chrom", "repeat_start", "repeat_end", "product_size",
            "fwd", "fwd_start", "fwd_end", "rev", "rev_start", "rev_end",
            "fwd_tm", "rev_tm", "repeats"]


def write_csv(rows):
    with open(OUT_CSV, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=CSV_COLS, delimiter=",",
                           extrasaction="ignore", lineterminator="\n")
        w.writeheader()
        for r in rows:
            w.writerow({k: (r[k].replace(NO_SYMBOL_FOOTNOTE, "")
                            .replace(PSEUDOGENE_FOOTNOTE, "") if k == "gene"
                            else r[k]) for k in CSV_COLS})


def main():
    rows = build_rows()
    write_markdown(rows)
    write_csv(rows)
    n_gc = sum(1 for r in rows if r["family"] == "GC")
    print(f"wrote results/primer_panel_delivery.md and .csv "
          f"({len(rows)} pairs: {len(rows) - n_gc} GCC + {n_gc} GC)")
    # primer coordinates must sit strictly inside the amplicon and bracket the
    # repeat tract -- a cheap invariant the release gate re-checks against the
    # genome itself.
    bad = [r["no"] for r in rows
           if not (r["amplicon_start"] <= r["fwd_start"] <= r["fwd_end"]
                   < r["repeat_start"] <= r["repeat_end"]
                   < r["rev_start"] <= r["rev_end"] <= r["amplicon_end"])]
    if bad:
        print(f"ERROR: primer/repeat ordering violated on rows {bad}",
              file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
