#!/usr/bin/env python3
"""Sensitivity check: trinuc GC-family tract scan with a 15bp pure-core seed.

Purpose: establish whether below the 24bp/30bp seeds any longer tract appears
(i.e., how robust the ">30" claim is). Same walk/filters as
13_repeat_scan_clean.py Part C.

Prints loci with n_units >= 20 and the overall maximum.
"""
import os
import re
import sys

BASE = os.environ.get("PRIMER_DESIGN_BASE") or os.path.dirname(
    os.path.dirname(os.path.abspath(__file__)))
CHROM_FA = os.path.join(BASE, "data", "hg38", "chromFa")
TRINUC = ["GCC", "CCG", "CGC", "GGC", "CGG", "GCG"]


def load_chrom(chrom):
    with open(os.path.join(CHROM_FA, f"{chrom}.fa")) as f:
        f.readline()
        return f.read().replace("\n", "").upper()


def walk(seq, i, unit, direction, max_streak=2, min_purity=0.6):
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


def gc_pct(s):
    return (s.count("G") + s.count("C")) / len(s) * 100 if s else 0.0


def main():
    CORE = re.compile(r"[GC]{15,}")
    rows = []
    for fn in sorted(os.listdir(CHROM_FA)):
        if not fn.endswith(".fa") or "_" in fn:
            continue
        chrom = fn[:-3]
        seq = load_chrom(chrom)
        for m in CORE.finditer(seq):
            core_s = m.start()
            core = seq[m.start():m.end()]
            best = None
            for unit in TRINUC:
                for shift in range(3):
                    i, n = shift, 0
                    while i + 3 <= len(core) and core[i:i + 3] == unit:
                        n += 3
                        i += 3
                    if n < 15:
                        continue
                    start = core_s + shift
                    le, li, l_edge = walk(seq, start, unit, -1)
                    re_, ri, r_edge = walk(seq, start + n, unit, +1)
                    n_units = n // 3 + le + li + re_ + ri
                    exact = n // 3 + le + re_
                    purity = exact / n_units if n_units else 0
                    if purity < 0.75:
                        continue
                    if best is None or n_units > best[0]:
                        best = (n_units, start, n, l_edge, r_edge, purity, unit)
            if best is None:
                continue
            n_units, start, perfect_n, l_edge, r_edge, purity, unit = best
            if n_units < 16:
                continue
            tract_s, tract_e = min(l_edge, start), max(r_edge, start + perfect_n)
            assert (tract_e - tract_s) == n_units * 3
            rows.append({"chrom": chrom, "start": tract_s, "end": tract_e,
                         "unit": unit, "n_units": n_units, "purity": round(purity * 100, 1),
                         "tract_len": tract_e - tract_s,
                         "gc": round(gc_pct(seq[tract_s:tract_e]), 1)})
        del seq
        print(f"{chrom}: {len(rows)} cum", file=sys.stderr)
    rows.sort(key=lambda r: (r["chrom"], r["start"]))
    dedup = []
    for r in rows:
        if dedup and dedup[-1]["chrom"] == r["chrom"] and r["start"] < dedup[-1]["end"] + 10:
            if r["n_units"] > dedup[-1]["n_units"]:
                dedup[-1] = r
            continue
        dedup.append(r)
    ge20 = sorted([r for r in dedup if r["n_units"] >= 20], key=lambda r: -r["n_units"])
    print(f"\n15bp-seed scan: {len(dedup)} loci >=16 units")
    print(f"max units: {ge20[0]['n_units'] if ge20 else 0}")
    print(f"loci with >=20 units: {len(ge20)}")
    for r in ge20:
        print(f"  {r['chrom']}:{r['start']+1}-{r['end']} unit={r['unit']} "
              f"n={r['n_units']} purity={r['purity']}% gc={r['gc']}%")
    import json
    with open(os.path.join(BASE, "data", "sensitivity_15bp.json"), "w") as f:
        json.dump({"n_loci_ge16": len(dedup), "max_units": ge20[0]["n_units"] if ge20 else 0,
                   "n_ge20": len(ge20), "ge20_list": ge20}, f, indent=1)


if __name__ == "__main__":
    main()
