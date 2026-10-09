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

BASE = os.environ.get("PRIMER_DESIGN_BASE") or os.path.dirname(
    os.path.dirname(os.path.abspath(__file__)))
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
        if not os.path.exists(p):
            # data/ is gitignored and may be absent from a checkout; report it
            # instead of dying with a traceback mid-gate.
            check(f"{path} 行数 == {n}", False, "文件缺失")
            continue
        try:
            if path.endswith(".json"):
                with open(p) as f:
                    d = json.load(f)
                cnt = len(d["rows"]) if isinstance(d, dict) and "rows" in d else len(d)
            else:
                cnt = sum(1 for _ in open(p)) - 1
        except Exception as e:
            check(f"{path} 行数 == {n}", False, f"读取失败: {e}")
            continue
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
    genome_gap = []
    for r in rows:
        c = r["chrom"]
        if c not in cache:
            try:
                cache[c] = load_chrom(c)
            except OSError:
                # no local genome: record it once per chromosome, skip the row
                if c not in genome_gap:
                    genome_gap.append(c)
                cache[c] = None
        seq = cache[c]
        if seq is None:
            continue
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
    if genome_gap:
        check(f"参考基因组可用（{CHROM_FA}）", False,
              f"缺失染色体: {','.join(genome_gap[:5])}"
              f"{' 等' if len(genome_gap) > 5 else ''}；本项及以下基因组相关检查无法执行")
    check("24 行引物位点/产物/包含关系（fresh 基因组读取）",
          not bad and not genome_gap,
          "; ".join(bad[:5]) if bad else
          ("参考基因组不可用，本项未能执行" if genome_gap else ""))

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
               "10_validate.py", "12_review_dataset.py", "21_build_delivery.py")))

    # F9 the delivery table -- the artifact the sample-prep group actually uses.
    # Its two primer-coordinate columns are derived (not stored upstream), and
    # its `Repeats 特点` string is a re-formatted derivative, so both are
    # re-derived here from the reference genome and from the authoritative TSV.
    print("F9 delivery table")
    md_path = os.path.join(BASE, "results", "primer_panel_delivery.md")
    csv_path = os.path.join(BASE, "results", "primer_panel_delivery.csv")
    if not (os.path.exists(md_path) and os.path.exists(csv_path)):
        check("交付表 .md / .csv 存在", False, "文件缺失")
    else:
        check("交付表 .md / .csv 存在", True)
        with open(csv_path, newline="") as f:
            dl = list(csv.DictReader(f))
        check("交付表 CSV 行数 == 24", len(dl) == 24, f"{len(dl)}")
        need = ["no", "gene", "chrom", "repeat_start", "repeat_end", "product_size",
                "fwd", "fwd_start", "fwd_end", "rev", "rev_start", "rev_end",
                "fwd_tm", "rev_tm", "repeats"]
        check("交付表列齐全且非空",
              all(k in dl[0] for k in need)
              and all(str(r.get(k, "")).strip() for r in dl for k in need))

        tsv_by_no = {r["no"]: r for r in rows}
        bad = []
        for r in dl:
            t = tsv_by_no.get(r["no"])
            if t is None:
                bad.append(f"#{r['no']} 在权威 TSV 中不存在")
                continue
            if r["fwd"] != t["fwd"] or r["rev"] != t["rev"]:
                bad.append(f"#{r['no']} 引物序列与权威 TSV 不一致")
            if (r["chrom"], int(r["repeat_start"]), int(r["repeat_end"])) != \
               (t["chrom"], int(t["repeat_start"]), int(t["repeat_end"])):
                bad.append(f"#{r['no']} 重复区坐标与权威 TSV 不一致")
            if abs(float(r["fwd_tm"]) - float(t["fwd_tm"])) > 0.05 or \
               abs(float(r["rev_tm"]) - float(t["rev_tm"])) > 0.05:
                bad.append(f"#{r['no']} Tm 与权威 TSV 不一致")
            # internal geometry of the delivered row
            if int(r["fwd_end"]) - int(r["fwd_start"]) + 1 != len(r["fwd"]):
                bad.append(f"#{r['no']} F 坐标跨度 != 引物长度")
            if int(r["rev_end"]) - int(r["rev_start"]) + 1 != len(r["rev"]):
                bad.append(f"#{r['no']} R 坐标跨度 != 引物长度")
            if int(r["rev_end"]) - int(r["fwd_start"]) + 1 != int(r["product_size"]):
                bad.append(f"#{r['no']} 产物长度 != 引物坐标跨度")
            if not (int(r["fwd_end"]) < int(r["repeat_start"])
                    and int(r["repeat_end"]) < int(r["rev_start"])):
                bad.append(f"#{r['no']} 引物与重复区重叠")
            # fresh genome check of BOTH primer binding intervals
            c = r["chrom"]
            if c not in cache:
                try:
                    cache[c] = load_chrom(c)
                except OSError:
                    cache[c] = None
            seq = cache[c]
            if seq is None:
                continue
            if seq[int(r["fwd_start"]) - 1:int(r["fwd_end"])] != r["fwd"]:
                bad.append(f"#{r['no']} F 坐标处序列 != 正向引物")
            if seq[int(r["rev_start"]) - 1:int(r["rev_end"])] != rc(r["rev"]):
                bad.append(f"#{r['no']} R 坐标处序列 != 反向引物互补")
        check("交付表 24 行坐标/序列/几何（fresh 基因组复算）",
              not bad, "; ".join(bad[:5]))

        with open(md_path) as f:
            md_text = f.read()
        n_md_rows = sum(1 for ln in md_text.splitlines()
                        if ln.startswith("| ") and ln.split("|")[1].strip().isdigit())
        check("交付表 Markdown 含 24 行", n_md_rows == 24, f"{n_md_rows}")
        check("交付表 Markdown 含坐标与 Tm 条件脚注",
              "Tm F/R (°C)" in md_text and "参考基因组" in md_text
              and "primer_panel_delivery_notes.md" in md_text)

    print()
    if fails:
        print(f"RELEASE GATE: FAILED ({len(fails)} problems)")
        for n, d in fails:
            print(f"  - {n}: {d}")
        sys.exit(1)
    print(f"RELEASE GATE: PASSED ({len(notes)} checks clean)")


if __name__ == "__main__":
    main()
