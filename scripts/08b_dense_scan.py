#!/usr/bin/env python3
"""Find GC-repeat tracts genome-wide (strict periodicity + limited interruptions).

Seed: maximal pure-GC core >= 30 bp. Pick the best periodic frame in the core.
Extend outward unit-by-unit: a unit counts as exact if it equals the frame unit,
otherwise it is an interruption. Stop after 2 consecutive interruptions; require
overall exactness >= 75%. n_units includes interruption units (TRF-like
copy-number convention).

Writes data/dense_scan.tsv
"""
import os
import re
import sys
from collections import Counter

BASE = "/home/keith/projects/primer_design"
CHROM_FA = os.path.join(BASE, "data", "hg38", "chromFa")
OUT = os.path.join(BASE, "data", "dense_scan.tsv")

RC = str.maketrans("ACGT", "TGCA")
TRINUC = ["GCC", "CCG", "CGC", "GGC", "CGG", "GCG"]
DINUC = ["GC", "CG"]


def revcomp(s):
    return s.translate(RC)[::-1]


def best_frame_core(core):
    """Best (perfect_bp, unit, ulen) over strands/frames within a pure-GC core."""
    best = (0, "", 0)
    for ulen, units in ((3, TRINUC), (2, DINUC)):
        for strand in ("+", "-"):
            t = core if strand == "+" else revcomp(core)
            for unit in units:
                for shift in range(ulen):
                    i, n = shift, 0
                    while i + ulen <= len(t) and t[i:i + ulen] == unit:
                        n += ulen
                        i += ulen
                    if n > best[0]:
                        best = (n, unit, ulen)
    return best


def walk(seq, i, unit, direction, max_streak=2):
    """Walk from i (exclusive edge) in unit steps.

    direction=-1: units are seq[i-ulen:i]; +1: seq[i:i+ulen].
    Returns (n_exact, n_interrupt, new_edge). Stops after max_streak consecutive
    interruptions or if exactness drops below 60% after >= 4 units. Trailing
    interruption units are trimmed.
    """
    n_exact = n_int = streak = 0
    edge = i
    while True:
        u = seq[i:i + len(unit)] if direction > 0 else seq[i - len(unit):i]
        if len(u) < len(unit):
            break
        if u == unit:
            n_exact += 1
            streak = 0
        else:
            n_int += 1
            streak += 1
            if streak >= max_streak:
                break
        walked = n_exact + n_int
        if walked >= 4 and n_exact / walked < 0.6:
            break
        edge = i + len(unit) if direction > 0 else i - len(unit)
        i = edge
    # trim trailing interruptions
    while n_int > 0 and n_exact + n_int > 0:
        u = seq[edge - len(unit):edge] if direction > 0 else seq[edge:edge + len(unit)]
        if u == unit:
            break
        n_int -= 1
        edge -= len(unit) if direction > 0 else -len(unit)
    return n_exact, n_int, edge


def bucket(n):
    if n < 25:
        return "20-25"
    if n < 24:
        return "25-30"
    return ">30"


CORE_PAT = re.compile(r"[GC]{24,}")


def main():
    rows = []
    for fn in sorted(os.listdir(CHROM_FA)):
        if not fn.endswith(".fa") or "_" in fn:
            continue
        chrom = fn[:-3]
        with open(os.path.join(CHROM_FA, fn)) as f:
            f.readline()
            seq = f.read().replace("\n", "").upper()
        for m in CORE_PAT.finditer(seq):
            core_s, core_e = m.start(), m.end()
            core = seq[core_s:core_e]
            perfect_bp, unit, ulen = best_frame_core(core)
            if perfect_bp < 30 or not unit:
                continue
            # frame starts inside the core: choose the shift maximizing tract span
            best = None
            for shift in range(ulen):
                i, n = shift, 0
                while i + ulen <= len(core) and core[i:i + ulen] == unit:
                    n += ulen
                    i += ulen
                if n < 24:
                    continue
                start = core_s + shift
                le, li, l_edge = walk(seq, start, unit, -1)
                re_, ri, r_edge = walk(seq, start + n, unit, +1)
                n_units = n // ulen + le + li + re_ + ri
                exact = n // ulen + le + re_
                purity = exact / n_units if n_units else 0
                if purity < 0.75:
                    continue
                score = n_units
                if best is None or score > best[0]:
                    best = (score, start, n, le, li, re_, ri, l_edge, r_edge, purity)
            if best is None:
                continue
            _, start, perfect_n, le, li, re_, ri, l_edge, r_edge, purity = best
            n_units = perfect_n // ulen + le + li + re_ + ri
            if n_units < 16:
                continue
            tract_s, tract_e = min(l_edge, start), max(r_edge, start + perfect_n)
            tract = seq[tract_s:tract_e]
            rows.append({
                "chrom": chrom, "start": tract_s, "end": tract_e,
                "unit": unit, "ulen": ulen, "n_units": n_units,
                "exact_units": perfect_n // ulen + le + re_,
                "interrupt_units": li + ri,
                "purity": round(purity * 100, 1),
                "tract_len": len(tract),
                "gc": round((tract.count("G") + tract.count("C")) / len(tract) * 100, 1),
                "bucket": bucket(n_units),
            })

    rows.sort(key=lambda r: (r["chrom"], r["start"]))
    dedup = []
    for r in rows:
        if dedup and dedup[-1]["chrom"] == r["chrom"] and r["start"] < dedup[-1]["end"] + 10:
            if r["n_units"] > dedup[-1]["n_units"]:
                dedup[-1] = r
            continue
        dedup.append(r)

    dedup.sort(key=lambda r: (r["ulen"], r["bucket"], -r["n_units"]))
    with open(OUT, "w") as f:
        cols = ["chrom", "start", "end", "unit", "ulen", "n_units", "exact_units",
                "interrupt_units", "purity", "tract_len", "gc", "bucket"]
        f.write("\t".join(cols) + "\n")
        for r in dedup:
            f.write("\t".join(str(r[c]) for c in cols) + "\n")

    print(Counter((r["unit"], r["bucket"]) for r in dedup), file=sys.stderr)
    print(f"total: {len(dedup)} -> {OUT}", file=sys.stderr)
    for r in dedup[:60]:
        print(f"  {r['chrom']}:{r['start']}-{r['end']} unit={r['unit']} n_units={r['n_units']} "
              f"exact={r['exact_units']} int={r['interrupt_units']} purity={r['purity']}% "
              f"len={r['tract_len']} gc={r['gc']}% [{r['bucket']}]", file=sys.stderr)


if __name__ == "__main__":
    main()
