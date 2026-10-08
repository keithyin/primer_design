#!/usr/bin/env python3
"""Sensitivity check for GC dinucleotide repeats with a 12bp pure-core seed
(6 units). Confirms the "no (GC)n tract reaches 20 units" claim under a
much lower threshold than the 20bp seed used in 13_repeat_scan_clean.py.
"""
import json
import os
import re
import sys

BASE = "/home/keith/projects/primer_design"
CHROM_FA = os.path.join(BASE, "data", "hg38", "chromFa")


def load_chrom(chrom):
    with open(os.path.join(CHROM_FA, f"{chrom}.fa")) as f:
        f.readline()
        return f.read().replace("\n", "").upper()


def walk(seq, i, unit, direction, max_streak=2, min_purity=0.6):
    ulen = 2
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


def main():
    loci = []
    for fn in sorted(os.listdir(CHROM_FA)):
        if not fn.endswith(".fa") or "_" in fn:
            continue
        chrom = fn[:-3]
        seq = load_chrom(chrom)
        seen = []
        for u in ("GC", "CG"):
            for m in re.finditer(f"(?:{u}){{6,}}", seq):
                s = m.start()
                ne_l, ni_l, le = walk(seq, s, u, -1)
                ne_r, ni_r, re_ = walk(seq, s, u, +1)
                units = ne_l + ni_l + ne_r + ni_r
                exact = ne_l + ne_r
                if units < 12:
                    continue
                purity = exact / units
                if purity < 0.75:
                    continue
                seen.append((le, re_, units, exact, u))
        # dedupe overlapping, keep max units
        seen.sort(key=lambda x: -x[2])
        kept = []
        for le, re_, units, exact, u in seen:
            if any(le < kr and kl < re_ for kl, kr, _, _, _ in kept):
                continue
            kept.append((le, re_, units, exact, u))
        for le, re_, units, exact, u in kept:
            loci.append({"chrom": chrom, "start": le, "end": re_, "unit": u,
                         "total_units": units, "exact_units": exact,
                         "purity": round(exact / units * 100, 1)})
        del seq
        print(f"{chrom}: {len(loci)} cum", file=sys.stderr)
    loci.sort(key=lambda x: -x["total_units"])
    top = loci[:15]
    ge20 = [x for x in loci if x["total_units"] >= 20]
    print(f"\n12bp-seed dinuc scan: {len(loci)} loci >=12 units")
    print(f"max total units: {loci[0]['total_units'] if loci else 0}")
    print(f"loci >=20 units: {len(ge20)}")
    for x in top:
        print(f"  {x['chrom']}:{x['start']+1}-{x['end']} unit={x['unit']} "
              f"units={x['total_units']} exact={x['exact_units']} purity={x['purity']}%")
    with open(os.path.join(BASE, "data", "sensitivity_dinuc.json"), "w") as f:
        json.dump({"n_loci": len(loci), "max_units": loci[0]["total_units"] if loci else 0,
                   "n_ge20": len(ge20), "top15": top}, f, indent=1)


if __name__ == "__main__":
    main()
