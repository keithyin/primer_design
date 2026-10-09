#!/usr/bin/env python3
"""Sequence-level genome-wide scan (hg38, UCSC chromFa) for perfect tandem
GC-dinucleotide and pure-GC trinucleotide repeat tracts with >= 20 units.

TRF's simpleRepeat track uses a score threshold that misses short GC tracts,
so we scan the sequence directly. Both strands are scanned; overlapping frame
matches are merged into single tracts.

Writes data/genome_scan.tsv
"""
import os
import re
import sys

BASE = os.environ.get("PRIMER_DESIGN_BASE") or os.path.dirname(
    os.path.dirname(os.path.abspath(__file__)))
CHROM_FA = os.path.join(BASE, "data", "hg38", "chromFa")
OUT = os.path.join(BASE, "data", "genome_scan.tsv")

DINUC_PAT = re.compile(r"(?:GC){20,}")
TRINUC_ROTATIONS = ["GCC", "CCG", "CGC", "GGC", "CGG", "GCG"]
TRINUC_PAT = re.compile("|".join(f"(?:{u}){{20,}}" for u in TRINUC_ROTATIONS))

RC = str.maketrans("ACGT", "TGCA")


def revcomp(s):
    return s.translate(RC)[::-1]


def analyze_tract(tract, fam):
    """Best rotational frame on either strand; returns perfect copies & purity."""
    units = ("GC", "CG") if fam == "GC" else tuple(TRINUC_ROTATIONS)
    ulen = len(units[0])
    best = (0, units[0], 0)  # perfect_bp, unit, shift
    for strand in ("+", "-"):
        t = tract if strand == "+" else revcomp(tract)
        for unit in units:
            for shift in range(ulen):
                i, n = shift, 0
                while i + ulen <= len(t) and t[i:i + ulen] == unit:
                    n += ulen
                    i += ulen
                if n > best[0]:
                    best = (n, unit, shift)
    perfect_bp, unit, _ = best
    gc = (tract.count("G") + tract.count("C")) / len(tract) * 100
    return {"unit": unit, "perfect_copies": perfect_bp // ulen,
            "perfect_bp": perfect_bp, "purity": round(perfect_bp / len(tract) * 100, 1),
            "tract_len": len(tract), "gc": round(gc, 1),
            "n_units": len(tract) // ulen}


def bucket(n):
    if n < 25:
        return "20-25"
    if n < 30:
        return "25-30"
    return ">30"


def main():
    tracts = {}  # chrom -> list of (s, e)
    for fn in sorted(os.listdir(CHROM_FA)):
        if not fn.endswith(".fa"):
            continue
        chrom = fn[:-3]
        with open(os.path.join(CHROM_FA, fn)) as f:
            f.readline()
            seq = f.read().replace("\n", "").upper()
        hits = []
        for strand, s in (("+", seq), ("-", revcomp(seq))):
            for pat, fam in ((DINUC_PAT, "GC"), (TRINUC_PAT, "GCC")):
                for m in pat.finditer(s):
                    # map reverse-strand match back to forward coordinates
                    if strand == "+":
                        gs, ge = m.start(), m.end()
                    else:
                        gs, ge = len(s) - m.end(), len(s) - m.start()
                    hits.append((gs, ge))
        if not hits:
            continue
        hits.sort()
        merged = []
        for gs, ge in hits:
            if merged and gs <= merged[-1][1]:
                merged[-1][1] = max(merged[-1][1], ge)
            else:
                merged.append([gs, ge])
        tracts[chrom] = merged
        print(f"{chrom}: {len(merged)} tracts", file=sys.stderr)

    rows = []
    for chrom, merged in tracts.items():
        with open(os.path.join(CHROM_FA, f"{chrom}.fa")) as f:
            f.readline()
            seq = f.read().replace("\n", "").upper()
        for gs, ge in merged:
            tract = seq[gs:ge]
            # test both families, keep the one with the higher purity
            a_d = analyze_tract(tract, "GC")
            a_t = analyze_tract(tract, "GCC")
            a = a_d if a_d["purity"] >= a_t["purity"] else a_t
            if a["purity"] < 80:
                continue
            rows.append({
                "chrom": chrom, "start": gs, "end": ge, "motif": "GC" if a["unit"] in ("GC", "CG") else "GCC",
                "unit": a["unit"], "n_units": a["n_units"], "perfect_copies": a["perfect_copies"],
                "purity": a["purity"], "tract_len": a["tract_len"], "gc": a["gc"],
                "bucket": bucket(a["n_units"]),
            })

    rows.sort(key=lambda r: (r["motif"], r["bucket"], r["chrom"], r["start"]))
    with open(OUT, "w") as f:
        f.write("chrom\tstart\tend\tmotif\tunit\tn_units\tperfect_copies\tpurity\ttract_len\tgc\tbucket\n")
        for r in rows:
            f.write("\t".join(str(r[k]) for k in
                    ["chrom", "start", "end", "motif", "unit", "n_units", "perfect_copies",
                     "purity", "tract_len", "gc", "bucket"]) + "\n")

    from collections import Counter
    print(Counter((r["motif"], r["bucket"]) for r in rows), file=sys.stderr)
    print(f"total: {len(rows)} -> {OUT}", file=sys.stderr)


if __name__ == "__main__":
    main()
