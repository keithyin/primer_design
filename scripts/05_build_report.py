#!/usr/bin/env python3
"""Build results/primer_design_results.tsv (the authoritative deliverable table).

Inputs: results/designed.json (chosen pairs), results/annotation.tsv (gene
context), results/method_D_genomescan.json (authoritative specificity counts).
README.md is maintained alongside this pipeline (not generated here).
"""
import json
import os

BASE = "/home/keith/projects/primer_design"

# canonical GC-trinucleotide expansion loci (literature nomenclature),
# matched by coordinate window (tract bounds may shift a few bp between scans)
EXPANSION_SPECS = [
    ("chrX", 147912030, 147912060, "FMR1", "CGG",
     "脆性X综合征（FXTAS/FXS）；前突变55-200，全突变>200"),
    ("chrX", 148500610, 148500650, "AFF2", "CCG",
     "FRAXE脆性位点/智力障碍；扩展>200"),
]


def expansion_of(locus):
    for c, lo, hi, gene, unit, disorder in EXPANSION_SPECS:
        if locus["chrom"] == c and lo <= locus["start"] <= hi:
            return {"gene": gene, "lit_unit": unit, "disorder": disorder}
    return None


def load_tsv(path):
    with open(path) as f:
        header = f.readline().rstrip("\n").split("\t")
        return [dict(zip(header, l.rstrip("\n").split("\t"))) for l in f if l.strip()]


def main():
    with open(os.path.join(BASE, "results", "designed.json")) as f:
        raw = json.load(f)
    ann = {(r["chrom"], int(r["start"]), int(r["end"])): r
           for r in load_tsv(os.path.join(BASE, "results", "annotation.tsv"))}

    rows = []
    for idx0, item in enumerate(raw, 1):
        if item.get("status") not in ("ok", "ok_check"):
            continue
        locus, pair = item["locus"], item["pair"]
        a = ann.get((locus["chrom"], locus["start"], locus["end"]), {})
        exp = expansion_of(locus)
        gene = exp["gene"] if exp else a.get("gene", ".")
        if locus["motif"] == "GCC":
            if exp:
                unit_disp = f"({exp['lit_unit']})n"
            else:
                unit_disp = f"({locus['unit']})n"
            rep = (f"GCC家族·{unit_disp}×{locus['n_units']}"
                   f"（精确{locus['exact_units']}+中断{locus['interrupt_units']}，"
                   f"纯度{locus['purity']:.1f}%），{int(locus['tract_len'])}bp，GC {locus['gc']:.1f}%")
        else:
            rep = (f"GC家族·(GC)n×{locus['n_units']}（完美串联，纯度100%），"
                   f"{int(locus['tract_len'])}bp，GC 100%")
        # bucket
        n = locus["n_units"]
        if exp:
            bucket = "20-25 ✓ ／ >30扩展检测※"
            head = int((1000 - (pair["product_size"] - locus["tract_len"])) // 3)
            rep += f"；扩展检测：1kb内可容纳至~{head}拷贝（{exp['disorder']}）"
        elif locus["motif"] == "GCC" and n >= 30:
            bucket = ">30 ✓"
        elif locus["motif"] == "GCC" and 25 <= n < 30:
            bucket = "25-30 ✓"
        elif locus["motif"] == "GCC" and n >= 20:
            bucket = "20-25 ✓"
        elif locus["motif"] == "GCC":
            bucket = f"GCC 20-25 △最接近（{n}拷贝）"
        else:
            bucket = "GC家族参考最长（<20）△"

        rows.append({
            "no": len(rows) + 1,
            "family": locus["motif"],
            "bucket": bucket,
            "gene": gene,
            "biotype": a.get("biotype", "."),
            "locus_context": a.get("context", "."),
            "chrom": locus["chrom"],
            "repeat_start": locus["start"] + 1,   # 1-based
            "repeat_end": locus["end"],
            "amplicon_start": pair["amplicon_start"] + 1,
            "amplicon_end": pair["amplicon_end"],
            "product_size": pair["product_size"],
            "fwd": pair["left_seq"],
            "rev": pair["right_seq"],
            "fwd_tm": round(pair["left_tm"], 1),
            "rev_tm": round(pair["right_tm"], 1),
            "fwd_gc": round(pair["left_gc"], 1),
            "rev_gc": round(pair["right_gc"], 1),
            "fwd_full_off": pair.get("left_full_off", ""),
            "fwd_end15_off": pair.get("left_end15_off", ""),
            "rev_full_off": pair.get("right_full_off", ""),
            "rev_end15_off": pair.get("right_end15_off", ""),
            "pair_full_off": pair.get("full_off", ""),
            "pair_end15_off": pair.get("end15_off", ""),
            "repeat_feature": rep,
            "flag": item["status"],
        })

    # authoritative specificity metrics: position-based genome rescan (method D,
    # 11_validate_genome_scan.py) -- counts overlapping occurrences, unlike grep -o
    dpath = os.path.join(BASE, "results", "method_D_genomescan.json")
    if os.path.exists(dpath):
        with open(dpath) as f:
            d = json.load(f)
        for row, drow in zip(rows, d["rows"]):
            assert row["no"] == drow["no"]
            for side in ("fwd", "rev"):
                row[f"{side}_full_off"] = drow[side]["full_off"]
                row[f"{side}_end15_off"] = drow[side]["end15_off"]
                row[f"{side}_genome_full_hits"] = drow[side]["genome_full_hits"]
            row["pair_full_off"] = row["fwd_full_off"] + row["rev_full_off"]
            row["pair_end15_off"] = row["fwd_end15_off"] + row["rev_end15_off"]
            if row["pair_full_off"] > 0 or row["pair_end15_off"] > 10:
                row["flag"] = "ok_check"

    cols = list(rows[0].keys())
    with open(os.path.join(BASE, "results", "primer_design_results.tsv"), "w") as f:
        f.write("\t".join(cols) + "\n")
        for r in rows:
            f.write("\t".join(str(r[c]) for c in cols) + "\n")
    n = len(rows)
    fam = sum(1 for r in rows if r["family"] == "GCC")
    print(f"wrote results/primer_design_results.tsv ({n} pairs: {fam} GCC + {n - fam} GC)")


if __name__ == "__main__":
    main()
