#!/usr/bin/env python3
"""Assemble the final primer-panel locus list.

- GCC family (pure-GC trinucleotide): top real tracts (16-25 units) from
  interrupted_scan.tsv plus the dense-scan chr12 >30-unit tract; FMR1 + AFF2
  flagged as expansion-detecting rows (patient alleles reach >200 copies;
  reference baseline given). HMGB3 / TMEM185A (chrX, >=20 units) force-included.
- GC family (dinucleotide): longest real perfect tracts (12-14 units).

Writes data/panel_loci.json
"""
import csv
import json
import os

BASE = os.environ.get("PRIMER_DESIGN_BASE") or os.path.dirname(
    os.path.dirname(os.path.abspath(__file__)))


def main():
    tracts = []
    with open(os.path.join(BASE, "data", "interrupted_scan.tsv")) as f:
        header = f.readline().rstrip("\n").split("\t")
        for line in f:
            tracts.append(dict(zip(header, line.rstrip("\n").split("\t"))))
    for t in tracts:
        for k in ("start", "end", "ulen", "n_units", "exact_units", "interrupt_units"):
            t[k] = int(t[k])
        for k in ("purity", "tract_len", "gc"):
            t[k] = float(t[k])

    # the dense scan (seed >= 8 units) finds one real >30-unit tract that the
    # regular (seed >= 10 units) scan misses: chr12:7,628,694-7,628,790
    with open(os.path.join(BASE, "data", "dense_scan.tsv")) as f:
        for row in csv.DictReader(f, delimiter="\t"):
            if (row["chrom"], int(row["start"])) == ("chr12", 7628694):
                tracts.append({
                    "chrom": "chr12", "start": 7628694, "end": 7628790,
                    "unit": row["unit"], "ulen": 3,
                    "n_units": int(row["n_units"]),
                    "exact_units": int(row["exact_units"]),
                    "interrupt_units": int(row["interrupt_units"]),
                    "purity": float(row["purity"]),
                    "tract_len": float(row["tract_len"]),
                    "gc": float(row["gc"]), "bucket": ">30",
                })

    # ---- GCC family: top tracts with chromosome diversity ----
    # expansion loci matched by coordinate window (tract bounds shift a few bp
    # between scan versions); then rank by units, keep max 2 per chromosome
    def near(t, chrom, lo, hi):
        return t["chrom"] == chrom and lo <= t["start"] <= hi

    expansion_specs = [
        ("chrX", 147912030, 147912060, "FMR1",
         "fragile-X (CGG) locus; reference 22 units, premutation >55, "
         "full mutation >200 in patients"),
        ("chrX", 148500610, 148500650, "AFF2",
         "FRAXE fragile-site (CCG) locus; reference 21 units, expansion >200 in patients"),
    ]
    force_specs = [("chrX", 150983370, 150983410),   # HMGB3 CCG x22
                   ("chrX", 149631700, 149631740)]   # TMEM185A GCC x20

    pool = sorted(tracts, key=lambda t: -t["n_units"])
    picked = [t for t in pool
              if any(near(t, c, lo, hi) for c, lo, hi, _, _ in expansion_specs)]
    chrom_cnt = {}
    for t in picked:
        chrom_cnt[t["chrom"]] = chrom_cnt.get(t["chrom"], 0) + 1
    for t in pool:
        if t in picked:
            continue
        forced = any(near(t, c, lo, hi) for c, lo, hi in force_specs)
        if forced and len(picked) < 14:
            picked.append(t)
            chrom_cnt[t["chrom"]] = chrom_cnt.get(t["chrom"], 0) + 1
            continue
        c = chrom_cnt.get(t["chrom"], 0)
        if c >= 2:
            continue
        picked.append(t)
        chrom_cnt[t["chrom"]] = c + 1
        if len(picked) == 14:
            break

    loci = []
    for t in picked:
        hit = next(((g, note) for c, lo, hi, g, note in expansion_specs
                    if near(t, c, lo, hi)), None)
        if hit:
            gene, note = hit
            loci.append({
                "id": f"{t['chrom']}:{t['start']}", "chrom": t["chrom"],
                "start": t["start"], "end": t["end"],
                "motif": "GCC", "unit": t["unit"], "ulen": t["ulen"],
                "n_units": t["n_units"], "exact_units": t["exact_units"],
                "interrupt_units": t["interrupt_units"], "purity": t["purity"],
                "tract_len": t["tract_len"], "gc": t["gc"],
                "bucket_label": ">30（扩展检测）",
                "bucket_note": note, "expansion_row": True, "gene_hint": gene,
            })
        else:
            n = t["n_units"]
            if 20 <= n < 25:
                label, flag = "20-25", "✓ 完全符合"
            elif 25 <= n < 30:
                label, flag = "25-30", "✓ 完全符合"
            elif n >= 30:
                label, flag = ">30", "✓ 完全符合（密扫种子≥8单元检出）"
            else:
                label, flag = "20-25", f"△ 最接近（{n}拷贝；参考基因组该档真实最长22）"
            loci.append({
                "id": f"{t['chrom']}:{t['start']}", "chrom": t["chrom"],
                "start": t["start"], "end": t["end"],
                "motif": "GCC", "unit": t["unit"], "ulen": t["ulen"],
                "n_units": t["n_units"], "exact_units": t["exact_units"],
                "interrupt_units": t["interrupt_units"], "purity": t["purity"],
                "tract_len": t["tract_len"], "gc": t["gc"],
                "bucket_label": label, "bucket_note": flag,
            })

    # ---- GC family: longest perfect dinucleotide tracts ----
    dinuc = [
        ("chr19", 47271732, 47271760, 14),
        ("chr1", 88275021, 88275047, 13),
        ("chr12", 48971838, 48971864, 13),
        ("chr4", 61532531, 61532555, 12),
        ("chr2", 126245764, 126245790, 13),
        ("chr1", 83850240, 83850264, 12),
        ("chr10", 104587499, 104587523, 12),
        ("chr17", 47351304, 47351328, 12),
        ("chr5", 139652147, 139652171, 12),
        ("chr16", 73070915, 73070939, 12),
    ]
    for chrom, s, e, n in dinuc:
        loci.append({
            "id": f"{chrom}:{s}", "chrom": chrom, "start": s, "end": e,
            "motif": "GC", "unit": "GC", "ulen": 2,
            "n_units": n, "exact_units": n, "interrupt_units": 0,
            "purity": 100.0, "tract_len": e - s, "gc": 100.0,
            "bucket_label": "GC家族参考最长档（<20）",
            "bucket_note": "参考基因组不存在 ≥20 拷贝的 (GC)n（全基因组最长14拷贝）；"
                           "本行为该家族最长真实位点之一",
        })

    with open(os.path.join(BASE, "data", "panel_loci.json"), "w") as f:
        json.dump(loci, f, indent=1, ensure_ascii=False)
    print(f"panel: {len(loci)} loci "
          f"(GCC {sum(1 for l in loci if l['motif'] == 'GCC')}, "
          f"GC {sum(1 for l in loci if l['motif'] == 'GC')})")


if __name__ == "__main__":
    main()
