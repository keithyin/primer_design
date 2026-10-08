#!/usr/bin/env python3
"""Empirical availability analysis: how many GC-dinuc / GC-trinuc repeat tracts
exist in hg38 at various perfection thresholds, incl. interrupted tracts.
"""
import os
import re
import sys
from collections import Counter

BASE = "/home/keith/projects/primer_design"
CHROM_FA = os.path.join(BASE, "data", "hg38", "chromFa")

RC = str.maketrans("ACGT", "TGCA")
P_GC8 = re.compile(r"(?:GC){8,}|(?:CG){8,}")
P_GCWIN = re.compile(r"(?:[GC]{2}){20,}")   # >=40bp contiguous G/C-only
P_T3WIN = re.compile(r"(?:[GC]{3}){15,}")   # >=45bp contiguous G/C-only


def revcomp(s):
    return s.translate(RC)[::-1]


def best_frame(tract, ulen):
    units = {"GC", "CG"} if ulen == 2 else {"GCC", "CCG", "CGC", "GGC", "CGG", "GCG"}
    best = (0, "")
    for strand in ("+", "-"):
        t = tract if strand == "+" else revcomp(tract)
        for unit in units:
            for shift in range(ulen):
                i, n = shift, 0
                while i + ulen <= len(t) and t[i:i + ulen] == unit:
                    n += ulen
                    i += ulen
                if n > best[0]:
                    best = (n, unit)
    return best  # perfect_bp, unit


def main():
    dist_dinuc = Counter()
    rows = []
    for fn in sorted(os.listdir(CHROM_FA)):
        if not fn.endswith(".fa") or "_" in fn:
            continue  # primary chromosomes only
        chrom = fn[:-3]
        with open(os.path.join(CHROM_FA, fn)) as f:
            f.readline()
            seq = f.read().replace("\n", "").upper()
        seen = []
        for s in (seq, revcomp(seq)):
            for m in P_GC8.finditer(s):
                n = (m.end() - m.start()) // 2
                dist_dinuc[n] += 1
            for m in P_GCWIN.finditer(s):
                gs = m.start() if s is seq else len(s) - m.end()
                ge = m.end() if s is seq else len(s) - m.start()
                if any(abs(gs - s0) < 50 for s0, _ in seen):
                    continue
                seen.append((gs, ge))
                tract = seq[gs:ge]
                pb_d, u_d = best_frame(tract, 2)
                pb_t, u_t = best_frame(tract, 3)
                if pb_d / len(tract) >= pb_t / len(tract):
                    ulen, perfect_bp, unit = 2, pb_d, u_d
                else:
                    ulen, perfect_bp, unit = 3, pb_t, u_t
                n_units = len(tract) // ulen
                rows.append({
                    "chrom": chrom, "start": gs, "end": ge,
                    "unit": unit, "n_units": n_units,
                    "perfect_copies": perfect_bp // ulen,
                    "purity": round(perfect_bp / len(tract) * 100, 1),
                    "gc": round((tract.count("G") + tract.count("C")) / len(tract) * 100, 1),
                })
        print(f"{chrom} done", file=sys.stderr)

    print("\n== perfect dinucleotide (GC/CG)n distribution (n>=8, both strands counted) ==")
    for n in sorted(dist_dinuc):
        if n >= 8:
            print(f"  n={n}: {dist_dinuc[n]}")

    print("\n== GC-only windows >=40bp (contiguous G/C, any phase) ==")
    print("total windows:", len(rows))
    by_unit = Counter(r["unit"] for r in rows)
    print("by unit:", dict(by_unit))
    print("\nwindows with n_units >= 15:")
    for r in sorted(rows, key=lambda r: -r["n_units"])[:40]:
        print(f"  {r['chrom']}:{r['start']}-{r['end']} unit={r['unit']} n_units={r['n_units']} "
              f"perfect={r['perfect_copies']} purity={r['purity']}% gc={r['gc']}%")


if __name__ == "__main__":
    main()
