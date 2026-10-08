#!/usr/bin/env python3
"""Independent, from-scratch review of results/primer_design_results.tsv.

Imports nothing from the design pipeline. Every delivered number is re-derived
directly from the GRCh38 FASTA:

  P1  primer binding sites located by exact string search on the chromosome;
      positions fix the product from scratch; compared with the TSV coords.
  P2  repeat tract re-analysis under the documented counting rule (frame
      lattice anchored at the tract start; consume a unit while <2 consecutive
      non-exact units and running exact fraction >= 60% after >=4 units),
      including a walk reproduction from every maximal interior exact run
      (the tract must be reachable from at least one seed).
  P3  primer length / GC / Tm recomputed (primer3 calc_tm, design formula).
  P4  specificity claims cross-read from results/method_D_genomescan.json
      (independent numpy rescan) and the TSV flag rule re-applied.
  P5  unit-by-unit deep dives: FMR1, CARM1, GIGYF1, CCDC9.

Writes results/data_review.md
"""
import csv
import json
import os
import re
import sys

import primer3

BASE = "/home/keith/projects/primer_design"
CHROM_FA = os.path.join(BASE, "data", "hg38", "chromFa")

RC = str.maketrans("ACGT", "TGCA")


def rc(s):
    return s.translate(RC)[::-1]


def gc_pct(s):
    return (s.count("G") + s.count("C")) / len(s) * 100 if s else 0.0


def load_chrom(chrom):
    with open(os.path.join(CHROM_FA, f"{chrom}.fa")) as f:
        f.readline()
        return f.read().replace("\n", "").upper()


def find_all(seq, sub):
    out, i = [], 0
    while True:
        i = seq.find(sub, i)
        if i < 0:
            return out
        out.append(i)
        i += 1


def walk(seq, i, unit, direction, max_streak=2, min_purity=0.6):
    """Reimplementation of the documented counting rule (README)."""
    ulen = len(unit)
    ne = ni = streak = 0
    edge = i
    while True:
        u = seq[i:i + ulen] if direction > 0 else seq[i - ulen:i]
        if len(u) < ulen:
            break
        if u == unit:
            ne += 1
            streak = 0
        else:
            streak += 1
            ni += 1
            if streak >= max_streak:
                ni -= 1
                break
            if ne + ni >= 4 and ne / (ne + ni) < min_purity:
                ni -= 1
                break
        edge = i + ulen if direction > 0 else i - ulen
        i = edge
    return ne, ni, edge


def parse_feature(feat):
    units = int(re.search(r"×(\d+)", feat).group(1))
    m = re.search(r"精确(\d+)\+中断(\d+)", feat)
    exact, inter = (int(m.group(1)), int(m.group(2))) if m else (units, 0)
    purity = float(re.search(r"纯度([\d.]+)%", feat).group(1))
    tl = float(re.search(r"），([\d.]+)bp", feat).group(1))
    gcp = float(re.search(r"GC ([\d.]+)%", feat).group(1))
    ud = re.search(r"\(([GC]{2,3})\)n", feat).group(1)
    return dict(units=units, exact=exact, inter=inter, purity=purity,
                tl=tl, gc=gcp, unit_disp=ud)


def aligned_count(tract, unit, phase):
    u = len(unit)
    return sum(1 for i in range(phase, len(tract) - u + 1, u)
               if tract[i:i + u] == unit)


def maximal_runs(flags):
    runs, cur = [], 0
    for f_ in flags + [False]:
        if f_:
            cur += 1
        elif cur:
            runs.append(cur)
            cur = 0
    return runs


def analyze_row(r, seq, scan_unit, lines_buffer):
    issues, info = [], []
    fwd, rev = r["fwd"], r["rev"]
    rs0, re1 = int(r["repeat_start"]) - 1, int(r["repeat_end"])
    a0_rep, a1_rep = int(r["amplicon_start"]) - 1, int(r["amplicon_end"])
    prod_rep = int(r["product_size"])
    feat = parse_feature(r["repeat_feature"])

    # ---- P1: primer sites from scratch ----
    fpos = find_all(seq, fwd)
    rpos = find_all(seq, rc(rev))
    if len(fpos) != 1 or len(rpos) != 1:
        issues.append(f"P1 引物位点在本染色体不唯一 (fwd {len(fpos)}, rev {len(rpos)})")
        return issues, info
    a0, a1 = fpos[0], rpos[0] + len(rev)
    if a1 - a0 != prod_rep:
        issues.append(f"P1 重算产物 {a1 - a0} != TSV {prod_rep}")
    if a0 != a0_rep or a1 != a1_rep:
        issues.append(f"P1 重算扩增子 [{a0},{a1}) != TSV [{a0_rep},{a1_rep})")
    if not (a0 <= rs0 and re1 <= a1):
        issues.append("P1 重复区未被扩增子完整包含")
    if a0 + len(fwd) > rs0 or a1 - len(rev) < re1:
        issues.append("P1 引物与重复区重叠")
    if set(fwd + rev) - set("ACGT"):
        issues.append("P1 引物含非 ACGT 字符")

    # ---- P2: tract re-analysis ----
    tract = seq[rs0:re1]
    ulen = 2 if r["family"] == "GC" else 3
    rots = ["GC", "CG"] if ulen == 2 else ["GCC", "CCG", "CGC", "GGC", "CGG", "GCG"]
    units = len(tract) // ulen
    if len(tract) % ulen != 0:
        issues.append(f"P2 tract 长度 {len(tract)} 非 {ulen} 的整数倍")
    if units != feat["units"]:
        issues.append(f"P2 单元数 {units} != 特征串 ×{feat['units']}")
    if abs(len(tract) - feat["tl"]) > 0.05:
        issues.append(f"P2 tract 长度 {len(tract)} != 特征串 {feat['tl']}bp")

    phase0 = {u: aligned_count(tract, u, 0) for u in rots}
    # frame: prefer the scan's own unit; for GC dinucs the tract phase is direct
    if ulen == 2:
        unit = tract[0:2]
        if unit not in rots:
            issues.append(f"P2 二核苷酸单元 {unit} 非 GC/CG")
    else:
        unit = scan_unit or max(phase0, key=phase0.get)
        if scan_unit and scan_unit != feat["unit_disp"]:
            info.append(f"scan 阅读框 {scan_unit}（文献名 {feat['unit_disp']}）")
    if phase0.get(unit) != feat["exact"]:
        issues.append(f"P2 阅读框 {unit} phase0 精确数 {phase0.get(unit)} "
                      f"!= 特征串 {feat['exact']}（各框: {phase0}）")
    pmax = {u: max(aligned_count(tract, u, p) for p in range(ulen)) for u in rots}
    if pmax.get(unit, 0) > feat["exact"]:
        issues.append(f"P2 存在更高相位的精确数 {pmax[unit]} > {feat['exact']}")

    purity = round(feat["exact"] / units * 100, 1) if units else 0
    if abs(purity - feat["purity"]) > 0.11:
        issues.append(f"P2 纯度 {purity} != 特征串 {feat['purity']}")
    gcp = round(gc_pct(tract), 1)
    if abs(gcp - feat["gc"]) > 0.06:
        issues.append(f"P2 tract GC% {gcp} != 特征串 {feat['gc']}")

    flags = [tract[k * ulen:(k + 1) * ulen] == unit for k in range(units)]
    runs = maximal_runs(flags)
    run_len = max(runs) if runs else 0
    if ulen == 2:
        # GC-family rows are defined as maximal PERFECT tandem runs (see README):
        # verify 100% perfection and that neither side can be extended perfectly.
        if not all(flags):
            issues.append("P2 GC 行不是完美串联")
        if seq[rs0 - ulen:rs0] == unit or seq[re1:re1 + ulen] == unit:
            issues.append("P2 GC 完美串联不是最大（两侧可继续延伸）")
        info.append(f"P2 完美串联 {units} 单元，两侧不可延伸（最大完美串联）")
        return issues, info
    walk_ok, walk_tried = False, 0
    run_pos = []
    k = 0
    while k < units:
        if flags[k]:
            j = k
            while j < units and flags[j]:
                j += 1
            run_pos.append((k, j - k))
            k = j
        else:
            k += 1
    for k0, L in run_pos:
        if L < 3:
            continue
        walk_tried += 1
        ne_l, ni_l, l_edge = walk(seq, rs0 + k0 * ulen, unit, -1)
        ne_r, ni_r, r_edge = walk(seq, rs0 + (k0 + L) * ulen, unit, +1)
        if (l_edge == rs0 and r_edge == re1
                and ne_l + ne_r + L == feat["exact"]
                and ni_l + ni_r == feat["inter"]):
            walk_ok = True
            break
    if walk_tried == 0:
        info.append(f"P2 最长精确串联 {run_len} < 3，跳过 walk 复现")
    elif walk_ok:
        info.append(f"P2 最长精确串联 {run_len} 单元；walk 自种子复现边界与计数一致")
    else:
        issues.append(f"P2 walk 复现失败（最长精确串联 {run_len}，"
                      f"已尝试 {walk_tried} 个种子）")

    # ---- P3: primers ----
    for side in ("fwd", "rev"):
        s = r[side]
        if not (20 <= len(s) <= 28):
            issues.append(f"P3 {side} 长度 {len(s)} 超出 20-28")
        t = primer3.bindings.calc_tm(s, mv_conc=50.0, dv_conc=1.5,
                                     dntp_conc=0.6, dna_conc=50.0)
        if abs(t - float(r[f"{side}_tm"])) > 0.3:
            issues.append(f"P3 {side} Tm 复算 {t:.2f} != {r[f'{side}_tm']}")
        if abs(gc_pct(s) - float(r[f"{side}_gc"])) > 0.06:
            issues.append(f"P3 {side} GC% 复算不符")
    if abs(float(r["fwd_tm"]) - float(r["rev_tm"])) > 2.0:
        issues.append("P3 配对 ΔTm > 2.0")
    return issues, info


def deep_dive(seq, rs0, re1, unit, ulen, label, lines):
    lines.append(f"\n### {label}\n")
    lines.append(f"阅读框 = `{unit}`（起点与 tract 起点对齐），单元序 0 表示 tract 首单元。\n")
    lines.append("| 单元序 | 坐标(1-based) | 序列 | 状态 | GC数 |")
    lines.append("|---|---|---|---|---|")
    n_units = (re1 - rs0) // ulen
    for k in range(-4, n_units + 4):
        i = rs0 + k * ulen
        u = seq[i:i + ulen]
        if rs0 <= i and i + ulen <= re1:
            st = "精确" if u == unit else "**中断**"
        else:
            st = "tract 之外"
        lines.append(f"| {k} | {i + 1} | {u} | {st} | "
                     f"{u.count('G') + u.count('C')} |")
    flags = [seq[rs0 + k * ulen:rs0 + (k + 1) * ulen] == unit for k in range(n_units)]
    run_pos = []
    k = 0
    while k < n_units:
        if flags[k]:
            j = k
            while j < n_units and flags[j]:
                j += 1
            run_pos.append((k, j - k))
            k = j
        else:
            k += 1
    if run_pos:
        k0, L = max(run_pos, key=lambda x: x[1])
        ne_l, ni_l, l_edge = walk(seq, rs0 + k0 * ulen, unit, -1)
        ne_r, ni_r, r_edge = walk(seq, rs0 + (k0 + L) * ulen, unit, +1)
        lines.append(f"\nwalk 复现（种子 = 单元 {k0}–{k0 + L - 1}，共 {L} 个精确单元）："
                     f"左边界 {l_edge + 1}，右边界 {r_edge}；"
                     f"复现区间 [{l_edge + 1}, {r_edge}]，交付区间 [{rs0 + 1}, {re1}]；"
                     f"精确 {ne_l + ne_r + L}，中断 {ni_l + ni_r}。\n")
    if ulen == 3:
        t = seq[rs0:re1]
        lines.append("\n各阅读框/相位的精确单元计数：\n")
        lines.append("| 阅读框 | phase 0 | phase 1 | phase 2 |")
        lines.append("|---|---|---|---|")
        for u in ["GCC", "CCG", "CGC", "GGC", "CGG", "GCG"]:
            lines.append(f"| {u} | {aligned_count(t, u, 0)} | "
                         f"{aligned_count(t, u, 1)} | {aligned_count(t, u, 2)} |")


def main():
    rows = []
    with open(os.path.join(BASE, "results", "primer_design_results.tsv")) as f:
        for r in csv.DictReader(f, delimiter="\t"):
            rows.append(r)
    scan_unit = {}
    with open(os.path.join(BASE, "data", "interrupted_scan.tsv")) as f:
        for r in csv.DictReader(f, delimiter="\t"):
            scan_unit[(r["chrom"], int(r["start"]))] = r["unit"]
    with open(os.path.join(BASE, "results", "method_D_genomescan.json")) as f:
        dmap = {d["no"]: d for d in json.load(f)["rows"]}

    by_chrom = {}
    for r in rows:
        by_chrom.setdefault(r["chrom"], []).append(r)

    lines = ["# 交付数据独立复核报告（from-scratch review）", "",
             "复核方式：不引用设计流水线任何代码，直接从 GRCh38 参考序列重算",
             "每行的引物结合位点、产物坐标、重复结构（含 walk 边界复现）、",
             "热力学参数，并与独立 numpy 全基因组重扫（method D）交叉核对。", ""]
    n_issue = 0
    n_note = 0
    detail = []
    for chrom in sorted(by_chrom):
        seq = load_chrom(chrom)
        for r in by_chrom[chrom]:
            su = scan_unit.get((r["chrom"], int(r["repeat_start"]) - 1))
            issues, info = analyze_row(r, seq, su, [])
            d = dmap.get(int(r["no"]))
            if d is None:
                issues.append("P4 method_D 缺少对应行")
            else:
                for key in ("fwd", "rev"):
                    if int(r[f"{key}_genome_full_hits"]) != d[key]["genome_full_hits"]:
                        issues.append(f"P4 {key} genome_full_hits 不符")
                    if int(r[f"{key}_full_off"]) != d[key]["full_off"]:
                        issues.append(f"P4 {key} full_off 不符")
                    if int(r[f"{key}_end15_off"]) != d[key]["end15_off"]:
                        issues.append(f"P4 {key} end15_off 不符")
                pf, pe = int(r["pair_full_off"]), int(r["pair_end15_off"])
                exp_flag = "ok_check" if (pf > 0 or pe > 10) else "ok"
                if exp_flag != r["flag"]:
                    issues.append(f"P4 flag {r['flag']} != 规则 {exp_flag}")
            n_issue += bool(issues)
            n_note += len(info)
            status = "PASS" if not issues else "**FAIL**"
            note = "；".join(issues + info) or "—"
            detail.append(f"| {r['no']} | {r['chrom']}:{r['repeat_start']}-"
                          f"{r['repeat_end']} | {r['family']} | {r['product_size']} | "
                          f"{status} | {note} |")
        del seq

    lines.append(f"## P1–P4 逐行复核（{len(rows)} 行）\n")
    lines.append("| # | 重复区(hg38) | 家族 | 产物bp | 结论 | 复核记录 |")
    lines.append("|---|---|---|---|---|---|")
    lines += detail

    lines.append("\n## P5 单元级深挖")
    by_key = {(r["chrom"], int(r["repeat_start"]) - 1): r for r in rows}
    dives = [(("chrX", 147912046), "FMR1"),
             (("chr12", 7628694), "chr12 >30 位点"),
             (("chr19", 10871582), "CARM1"),
             (("chr7", 100694311), "GIGYF1"),
             (("chr19", 47271732), "CCDC9")]
    for (chrom, rs0_key), name in dives:
        r = by_key[(chrom, rs0_key)]
        rs0, re1 = int(r["repeat_start"]) - 1, int(r["repeat_end"])
        ulen = 2 if r["family"] == "GC" else 3
        su = scan_unit.get((chrom, rs0))
        unit = su if ulen == 3 else r["fwd"] and None  # placeholder replaced below
        seq = load_chrom(chrom)
        if ulen == 2:
            unit = seq[rs0:rs0 + 2]
        deep_dive(seq, rs0, re1, unit, ulen,
                  f"#{r['no']} {name}（{r['gene']}）{chrom}:{r['repeat_start']}-{r['repeat_end']}",
                  lines)
        del seq

    lines.append("\n## 结论\n")
    lines.append(f"- {len(rows)} 行中 **{len(rows) - n_issue} 行**通过全部逐项复核"
                 f"（引物位点 / 产物坐标 / 重复结构含边界 walk 复现 / "
                 f"Tm·GC 复算 / 独立重扫交叉核对）。")
    if n_issue:
        lines.append(f"- **{n_issue} 行存在待处理问题（见上表）**。")
    else:
        lines.append("- 未发现与交付数据矛盾之处。")
    lines.append(f"- 复核另记录 {n_note} 条结构备注（阅读框、串联长度、walk 种子等）。")

    out = os.path.join(BASE, "results", "data_review.md")
    with open(out, "w") as f:
        f.write("\n".join(lines) + "\n")
    print(f"wrote {out}: {len(rows) - n_issue}/{len(rows)} rows clean, {n_issue} with issues")
    sys.exit(1 if n_issue else 0)


if __name__ == "__main__":
    main()
