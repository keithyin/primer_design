#!/usr/bin/env python3
"""Clean re-scans with the corrected counting rule.

Part A  perfect tandem scan (>=10 dinuc copies / >=15 trinuc copies, all
        rotations, forward strand covers both strands by rotational symmetry)
        -> overwrites data/genome_scan.tsv with invariant-consistent rows.
Part B  interrupted walk scan for GC/CG dinucleotides (seed >= 10 perfect
        copies) -> max total units and the >=20 list (expected: none).
Part C  dense-seed trinuc scan (seed [GC]{24,}) with the corrected walk ->
        overwrites data/dense_scan.tsv; reports >= 20-unit loci.
Part D  invariant check on data/interrupted_scan.tsv (len == units * ulen).

Writes data/repeat_extremes.json with the summary statistics.
"""
import json
import os
import re
import sys

BASE = "/home/keith/projects/primer_design"
CHROM_FA = os.path.join(BASE, "data", "hg38", "chromFa")
TRINUC = ["GCC", "CCG", "CGC", "GGC", "CGG", "GCG"]
DINUC = ["GC", "CG"]


def load_chrom(chrom):
    with open(os.path.join(CHROM_FA, f"{chrom}.fa")) as f:
        f.readline()
        return f.read().replace("\n", "").upper()


def gc_pct(s):
    return (s.count("G") + s.count("C")) / len(s) * 100 if s else 0.0


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


def merge_spans(spans):
    spans = sorted(spans)
    merged = []
    for s, e in spans:
        if merged and s <= merged[-1][1]:
            merged[-1][1] = max(merged[-1][1], e)
        else:
            merged.append([s, e])
    return [tuple(m) for m in merged]


def part_a(chroms):
    rows = []
    for chrom in chroms:
        seq = load_chrom(chrom)
        spans = []
        for u in DINUC:
            for m in re.finditer(f"(?:{u}){{10,}}", seq):
                spans.append((m.start(), m.end()))
        for u in TRINUC:
            for m in re.finditer(f"(?:{u}){{15,}}", seq):
                spans.append((m.start(), m.end()))
        for s, e in merge_spans(spans):
            tract = seq[s:e]
            best = (0, "", "")
            for fam, units in (("GC", DINUC), ("GCC", TRINUC)):
                for u in units:
                    for t in (tract, tract.translate(str.maketrans("ACGT", "TGCA"))[::-1]):
                        for m in re.finditer(f"(?:{u})+", t):
                            n = (m.end() - m.start()) // len(u)
                            if n > best[0]:
                                best = (n, u, fam)
            n, u, fam = best
            if n < (10 if fam == "GC" else 15):
                continue
            rows.append({"chrom": chrom, "start": s, "end": e, "family": fam,
                         "unit": u, "perfect_copies": n, "tract_len": e - s,
                         "gc": round(gc_pct(tract), 1)})
        print(f"A {chrom}: {len(rows)} cumulative rows", file=sys.stderr)
        del seq
    rows.sort(key=lambda r: (r["family"], -r["perfect_copies"], r["chrom"], r["start"]))
    out = os.path.join(BASE, "data", "genome_scan.tsv")
    with open(out, "w") as f:
        f.write("chrom\tstart\tend\tfamily\tunit\tperfect_copies\ttract_len\tgc\n")
        for r in rows:
            f.write("\t".join(str(r[k]) for k in
                    ["chrom", "start", "end", "family", "unit", "perfect_copies",
                     "tract_len", "gc"]) + "\n")
    dinuc = [r for r in rows if r["family"] == "GC"]
    trinuc = [r for r in rows if r["family"] == "GCC"]
    stats = {
        "max_perfect_dinuc": max((r["perfect_copies"] for r in dinuc), default=0),
        "max_perfect_dinuc_loci": [(r["chrom"], r["start"] + 1) for r in dinuc
                                   if r["perfect_copies"] == max((x["perfect_copies"] for x in dinuc), default=0)],
        "dinuc_perfect_ge20": sum(1 for r in dinuc if r["perfect_copies"] >= 20),
        "max_perfect_trinuc": max((r["perfect_copies"] for r in trinuc), default=0),
        "max_perfect_trinuc_loci": [(r["chrom"], r["start"] + 1, r["unit"]) for r in trinuc
                                    if r["perfect_copies"] == max((x["perfect_copies"] for x in trinuc), default=0)],
        "trinuc_perfect_ge20": sum(1 for r in trinuc if r["perfect_copies"] >= 20),
    }
    print(f"A done: {len(rows)} rows -> {out}", file=sys.stderr)
    print(f"A stats: {stats}", file=sys.stderr)
    return stats


def part_b(chroms):
    loci = []
    for chrom in chroms:
        seq = load_chrom(chrom)
        seeds = []
        for u in DINUC:
            for m in re.finditer(f"(?:{u}){{10,}}", seq):
                seeds.append((m.start(), u))
        consumed = []
        for s, u in seeds:
            # walk both directions from the seed start: the right walk covers
            # the seed run itself (every unit consumed, exact or interrupt)
            ne_l, ni_l, l_edge = walk(seq, s, u, -1)
            ne_r, ni_r, r_edge = walk(seq, s, u, +1)
            units = ne_l + ni_l + ne_r + ni_r
            consumed.append((l_edge, r_edge, units, u))
        # dedupe overlapping consumed spans keeping the longest
        consumed.sort(key=lambda x: -x[2])
        kept = []
        for l, r, n, u in consumed:
            if any(l < kr and kl < r for kl, kr, _, _ in kept):
                continue
            kept.append((l, r, n, u))
        for l, r, n, u in kept:
            loci.append({"chrom": chrom, "start": l, "end": r, "total_units": n, "unit": u})
        print(f"B {chrom}: {len(loci)} cumulative loci", file=sys.stderr)
        del seq
    loci.sort(key=lambda x: -x["total_units"])
    ge20 = [x for x in loci if x["total_units"] >= 20]
    stats = {
        "max_dinuc_interrupted_units": loci[0]["total_units"] if loci else 0,
        "max_dinuc_interrupted_loci": [(x["chrom"], x["start"] + 1, x["total_units"])
                                       for x in loci[:5]],
        "dinuc_interrupted_ge20": len(ge20),
        "dinuc_interrupted_ge20_list": [(x["chrom"], x["start"] + 1, x["total_units"])
                                        for x in ge20[:20]],
    }
    print(f"B stats: {stats}", file=sys.stderr)
    return stats


def part_c(chroms):
    rows = []
    CORE = re.compile(r"[GC]{24,}")
    for chrom in chroms:
        seq = load_chrom(chrom)
        for m in CORE.finditer(seq):
            core_s, core_e = m.start(), m.end()
            core = seq[core_s:core_e]
            best = None
            for unit in TRINUC:
                ulen = 3
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
                    n_units = n // 3 + le + li + re_ + ri
                    exact = n // 3 + le + re_
                    purity = exact / n_units if n_units else 0
                    if purity < 0.75:
                        continue
                    if best is None or n_units > best[0]:
                        best = (n_units, start, n, le, li, re_, ri, l_edge, r_edge,
                                purity, unit)
            if best is None:
                continue
            (n_units, start, perfect_n, le, li, re_, ri, l_edge, r_edge,
             purity, unit) = best
            if n_units < 16:
                continue
            tract_s, tract_e = min(l_edge, start), max(r_edge, start + perfect_n)
            assert (tract_e - tract_s) == n_units * 3, \
                f"{chrom}:{tract_s}-{tract_e} invariant violated"
            tract = seq[tract_s:tract_e]
            rows.append({"chrom": chrom, "start": tract_s, "end": tract_e,
                         "unit": unit, "ulen": 3, "n_units": n_units,
                         "exact_units": perfect_n // 3 + le + re_,
                         "interrupt_units": li + ri,
                         "purity": round(purity * 100, 1), "tract_len": len(tract),
                         "gc": round(gc_pct(tract), 1),
                         "bucket": "20-25" if n_units < 25 else
                                   ("25-30" if n_units < 30 else ">30")})
        print(f"C {chrom}: {len(rows)} cumulative rows", file=sys.stderr)
        del seq
    # dedupe overlaps
    rows.sort(key=lambda r: (r["chrom"], r["start"]))
    dedup = []
    for r in rows:
        if dedup and dedup[-1]["chrom"] == r["chrom"] and r["start"] < dedup[-1]["end"] + 10:
            if r["n_units"] > dedup[-1]["n_units"]:
                dedup[-1] = r
            continue
        dedup.append(r)
    out = os.path.join(BASE, "data", "dense_scan.tsv")
    cols = ["chrom", "start", "end", "unit", "ulen", "n_units", "exact_units",
            "interrupt_units", "purity", "tract_len", "gc", "bucket"]
    with open(out, "w") as f:
        f.write("\t".join(cols) + "\n")
        for r in dedup:
            f.write("\t".join(str(r[c]) for c in cols) + "\n")
    ge20 = sorted([r for r in dedup if r["n_units"] >= 20],
                  key=lambda r: -r["n_units"])
    stats = {
        "max_dense_trinuc_units": max((r["n_units"] for r in dedup), default=0),
        "dense_trinuc_ge20": len(ge20),
        "dense_trinuc_ge20_list": [(r["chrom"], r["start"] + 1, r["n_units"],
                                    r["exact_units"], r["purity"]) for r in ge20],
    }
    print(f"C done: {len(dedup)} rows -> {out}", file=sys.stderr)
    print(f"C stats: {stats}", file=sys.stderr)
    return stats


def part_d():
    bad = []
    with open(os.path.join(BASE, "data", "interrupted_scan.tsv")) as f:
        for i, line in enumerate(f):
            if i == 0:
                continue
            c = line.rstrip("\n").split("\t")
            if (int(c[2]) - int(c[1])) != int(c[5]) * int(c[4]):
                bad.append(c[0] + ":" + c[1])
    n_rows = i
    print(f"D invariant check: {n_rows} rows, {len(bad)} violations", file=sys.stderr)
    return {"interrupted_scan_rows": n_rows, "interrupted_scan_invariant_violations": bad}


def main():
    chroms = [f[:-3] for f in sorted(os.listdir(CHROM_FA))
              if f.endswith(".fa") and "_" not in f]
    stats = {}
    stats.update(part_a(chroms))
    stats.update(part_b(chroms))
    stats.update(part_c(chroms))
    stats.update(part_d())
    with open(os.path.join(BASE, "data", "repeat_extremes.json"), "w") as f:
        json.dump(stats, f, indent=1, ensure_ascii=False)
    print("\n=== summary ===")
    print(json.dumps(stats, indent=1, ensure_ascii=False))


if __name__ == "__main__":
    main()
