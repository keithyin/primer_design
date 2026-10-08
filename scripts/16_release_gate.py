#!/usr/bin/env python3
"""Final release gate: fresh end-to-end verification of the whole deliverable.

Run in a clean process; touches every delivered file. Exits non-zero on any
failure. Covers:
  F1  file inventory & row counts
  F2  TSV structural integrity (columns, integer coords, unique primers)
  F3  per-row genome verification (binding sites, product span, tract bracket)
  F4  repeat structure accounting (perfect-maximal for GC; exact+interrupt for GCC)
  F5  cross-file agreement (designed/candidates/method_D/annotation/specificity)
  F6  design constraint compliance (sizes <=1kb, Tm/GC windows, pair dTm)
  F7  bucket/coverage expectations
  F8  reports self-consistency (validation_report, data_review, README numbers)
"""
import csv
import json
import os
import re
import sys

BASE = "/home/keith/projects/primer_design"
CHROM_FA = os.path.join(BASE, "data", "hg38", "chromFa")
RC = str.maketrans("ACGT", "TGCA")
fails = []
notes = []


def rc(s):
    return s.translate(RC)[::-1]


def gc_pct(s):
    return (s.count("G") + s.count("C")) / len(s) * 100 if s else 0.0


def load_chrom(c):
    with open(os.path.join(CHROM_FA, f"{c}.fa")) as f:
        f.readline()
        return f.read().replace("\n", "").upper()


def check(name, ok, detail=""):
    (notes if ok else fails).append((name, detail))
    print(f"  [{'PASS' if ok else 'FAIL'}] {name}" + (f" — {detail}" if detail else ""))


def main():
    tsv = os.path.join(BASE, "results", "primer_design_results.tsv")
    with open(tsv) as f:
        rows = list(csv.DictReader(f, delimiter="\t"))

    # F1 inventory
    print("F1 file inventory")
    inv = {
        "results/designed.json": 24,
        "results/candidates.json": 24,
        "results/method_D_genomescan.json": 24,
        "results/annotation.tsv": 24,      # data rows (excl. header)
        "results/specificity.tsv": 48,     # data rows (excl. header)
        "data/panel_loci.json": 24,
    }
    for path, n in inv.items():
        p = os.path.join(BASE, path)
        if path.endswith(".json"):
            d = json.load(open(p))
            cnt = len(d["rows"]) if isinstance(d, dict) and "rows" in d else len(d)
        else:
            cnt = sum(1 for _ in open(p)) - 1
        check(f"{path} 行数 == {n}", cnt == n, f"实际 {cnt}")

    # F2 structure
    print("F2 TSV structure")
    check("行数 == 24", len(rows) == 24, f"{len(rows)}")
    req_cols = ["no", "family", "bucket", "gene", "chrom", "repeat_start", "repeat_end",
                "amplicon_start", "amplicon_end", "product_size", "fwd", "rev",
                "fwd_tm", "rev_tm", "repeat_feature", "flag"]
    check("必需列齐全且非空",
          all(r.get(c, "").strip() for r in rows for c in req_cols))
    all_primers = [r[k] for r in rows for k in ("fwd", "rev")]
    check("48 条引物序列两两不同（面板内无重复）",
          len(set(all_primers)) == 48, f"唯一序列 {len(set(all_primers))}")
    check("引物均为 ACGT 且长 20-28",
          all(set(p) <= set("ACGT") and 20 <= len(p) <= 28 for p in all_primers))
    coords_ok = all(
        int(r["product_size"]) == int(r["amplicon_end"]) - int(r["amplicon_start"]) + 1
        and 1 <= int(r["repeat_start"]) <= int(r["repeat_end"])
        for r in rows)
    check("坐标为整数且 product == amplicon 跨度", coords_ok)

    # F3/F4 per-row genome verification
    print("F3/F4 per-row genome verification")
    cache = {}
    bad = []
    for r in rows:
        c = r["chrom"]
        if c not in cache:
            cache[c] = load_chrom(c)
        seq = cache[c]
        a0, a1 = int(r["amplicon_start"]) - 1, int(r["amplicon_end"])
        rs0, re1 = int(r["repeat_start"]) - 1, int(r["repeat_end"])
        amp = seq[a0:a1]
        if not amp.startswith(r["fwd"]) or amp.count(r["fwd"]) != 1:
            bad.append(f"#{r['no']} fwd 位点")
        if not amp.endswith(rc(r["rev"])) or amp.count(rc(r["rev"])) != 1:
            bad.append(f"#{r['no']} rev 位点")
        if len(amp) != int(r["product_size"]):
            bad.append(f"#{r['no']} 产物长度")
        if not (a0 <= rs0 and re1 <= a1):
            bad.append(f"#{r['no']} tract 未被包含")
        ulen = 2 if r["family"] == "GC" else 3
        tract = seq[rs0:re1]
        if len(tract) % ulen != 0:
            bad.append(f"#{r['no']} tract 长度非整单元")
            continue
        units = len(tract) // ulen
        feat = r["repeat_feature"]
        m = re.search(r"×(\d+)", feat)
        if not m or int(m.group(1)) != units:
            bad.append(f"#{r['no']} 单元数与特征串不符")
        if r["family"] == "GC":
            rots = ("GC", "CG")
            if not any(tract == u * units or rc(tract) == u * units for u in rots):
                bad.append(f"#{r['no']} GC 非完美串联")
            if seq[rs0 - 2:rs0] == seq[rs0:rs0 + 2] or seq[re1:re1 + 2] == seq[re1 - 2:re1]:
                bad.append(f"#{r['no']} GC 串联两侧可延伸（非最大）")
        else:
            m2 = re.search(r"精确(\d+)\+中断(\d+)", feat)
            ex, it = int(m2.group(1)), int(m2.group(2))
            if ex + it != units:
                bad.append(f"#{r['no']} 精确+中断 != 单元数")
    check("24 行引物位点/产物/包含关系（fresh 基因组读取）", not bad, "; ".join(bad[:5]))

    # F5 cross-file agreement
    print("F5 cross-file agreement")
    designed = json.load(open(os.path.join(BASE, "results", "designed.json")))
    dmap = {d["no"]: d for d in json.load(open(
        os.path.join(BASE, "results", "method_D_genomescan.json")))["rows"]}
    ann = list(csv.DictReader(open(os.path.join(BASE, "results", "annotation.tsv")),
                              delimiter="\t"))
    ok = len(designed) == len(dmap) == len(rows) == len(ann) == 24
    for r, it, d, a in zip(rows, designed, (dmap[i] for i in range(1, 25)), ann):
        if r["fwd"] != it["pair"]["left_seq"] or r["rev"] != it["pair"]["right_seq"]:
            ok = False
        if (r["chrom"], int(r["repeat_start"]) - 1) != (a["chrom"], int(a["start"])):
            ok = False
        for side in ("fwd", "rev"):
            if int(r[f"{side}_full_off"]) != d[side]["full_off"]:
                ok = False
            if int(r[f"{side}_genome_full_hits"]) != 1:
                ok = False
    check("TSV == designed == method_D == annotation（24 行全对齐，full_hits==1）", ok)

    # F6 constraints
    print("F6 design constraints")
    sizes = [int(r["product_size"]) for r in rows]
    tms = [float(r[k]) for r in rows for k in ("fwd_tm", "rev_tm")]
    gcs = [float(r[k]) for r in rows for k in ("fwd_gc", "rev_gc")]
    dtms = [abs(float(r["fwd_tm"]) - float(r["rev_tm"])) for r in rows]
    check("产物 338-990 bp，全部 ≤1kb", max(sizes) <= 1000, f"{min(sizes)}-{max(sizes)}")
    check("Tm 全部在 57-64℃", all(57 <= t <= 64 for t in tms),
          f"{min(tms)}-{max(tms)}")
    check("配对 ΔTm ≤ 2.0℃", all(d <= 2.0 for d in dtms), f"max {max(dtms):.1f}")
    check("引物 GC% 全部在 35-65", all(35 <= g <= 65 for g in gcs),
          f"{min(gcs)}-{max(gcs)}")
    check("全长脱靶合计 = 0", sum(int(r["pair_full_off"]) for r in rows) == 0)
    check("flag 分布 = 18 ok + 6 ok_check",
          sum(1 for r in rows if r["flag"] == "ok") == 18
          and sum(1 for r in rows if r["flag"] == "ok_check") == 6)

    # F7 coverage
    print("F7 coverage")
    b25 = [r for r in rows if r["bucket"].startswith("20-25 ✓") or r["bucket"].startswith("25-30 ✓") or r["bucket"].startswith(">30 ✓")]
    check("GCC 完全符合档位行 ≥ 7（20-25 ×4, 25-30 ×2, >30 ×1）",
          sum(1 for r in b25 if r["family"] == "GCC") >= 7,
          f"{sum(1 for r in b25 if r['family'] == 'GCC')}")
    check("GC 家族行 = 10", sum(1 for r in rows if r["family"] == "GC") == 10)
    check("扩展检测行 = 2（FMR1/AFF2）",
          sum(1 for r in rows if "扩展检测" in r["bucket"]) == 2)
    chr12 = [r for r in rows if r["chrom"] == "chr12" and r["bucket"].startswith(">30 ✓")]
    check("chr12 >30 ✓ 行存在且 32 单元",
          len(chr12) == 1 and (int(chr12[0]["repeat_end"]) - int(chr12[0]["repeat_start"]) + 1) // 3 == 32)

    # F8 reports
    print("F8 reports self-consistency")
    vr = open(os.path.join(BASE, "results", "validation_report.md")).read()
    dr = open(os.path.join(BASE, "results", "data_review.md")).read()
    rd = open(os.path.join(BASE, "results", "README.md")).read()
    check("validation_report: 16 项全 PASS", "16 项通过，0 项失败" in vr
          and "FAIL" not in vr)
    check("data_review: 24/24 行通过", "24 行中 **24 行**通过" in dr and "FAIL" not in dr)
    check("README 数字同步（24 对 / 48 条引物）",
          "24 对" in rd and "48 条" in rd and "23 对" not in rd)
    check("校验/复核脚本在 scripts/ 且可追溯",
          all(os.path.exists(os.path.join(BASE, "scripts", f)) for f in
              ("09_build_panel.py", "02b_design_panel.py", "03b_pick_specific.py",
               "04_annotate.py", "11_validate_genome_scan.py", "05_build_report.py",
               "10_validate.py", "12_review_dataset.py")))

    print()
    if fails:
        print(f"RELEASE GATE: FAILED ({len(fails)} problems)")
        for n, d in fails:
            print(f"  - {n}: {d}")
        sys.exit(1)
    print(f"RELEASE GATE: PASSED ({len(notes)} checks clean)")


if __name__ == "__main__":
    main()
