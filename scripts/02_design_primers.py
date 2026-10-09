#!/usr/bin/env python3
"""Select repeat loci, verify tract copy number from the genome, design flanking
primers with primer3.

Reads   data/candidates.tsv   (from 01_scan_repeats.py)
Uses    data/hg38/chromFa/*.fa (local hg38) when available, else UCSC API
Writes  results/loci_designed.tsv and results/designed.json
"""
import functools
import json
import os
import sys
import urllib.request

import primer3

BASE = os.environ.get("PRIMER_DESIGN_BASE") or os.path.dirname(
    os.path.dirname(os.path.abspath(__file__)))
GENOME_DIR = os.path.join(BASE, "data", "hg38")
CHROM_FA = os.path.join(GENOME_DIR, "chromFa")

# ---------------------------------------------------------------- genome access

_chrom_cache = {}


@functools.lru_cache(maxsize=4)
def load_chrom(chrom):
    """Return (uppercase, gapless) sequence of a chromosome from local hg38 fasta."""
    if chrom in _chrom_cache:
        return _chrom_cache[chrom]
    path = os.path.join(CHROM_FA, f"{chrom}.fa")
    if not os.path.exists(path):
        sys.exit(f"local genome missing: {path}")
    seq = []
    with open(path) as f:
        for line in f:
            if line.startswith(">"):
                continue
            seq.append(line.strip())
    s = "".join(seq).upper()
    _chrom_cache[chrom] = s
    return s


def fetch_flank_ucsc(chrom, start, end):
    url = (f"https://api.genome.ucsc.edu/getData/sequence?genome=hg38;"
           f"chrom={chrom};start={start};end={end}")
    with urllib.request.urlopen(url, timeout=60) as r:
        return json.load(r)["dna"].upper()


def get_template(chrom, start, end):
    """Fetch [start,end) 0-based half-open, preferring local genome."""
    if os.path.isdir(CHROM_FA):
        return load_chrom(chrom)[start:end]
    return fetch_flank_ucsc(chrom, start, end)


# ---------------------------------------------------------------- tract analysis

def revcomp(s):
    return s.translate(str.maketrans("ACGT", "TGCA"))[::-1]


def canonical_motif(unit):
    """Canonical rotational frame: GC family and GCC family."""
    rot = {unit[i:] + unit[:i] for i in range(len(unit))}
    for fam, members in (("GC", {"GC", "CG"}), ("GCC", {"GCC", "CCG", "CGC", "GGC", "CGG", "GCG"})):
        if rot & members:
            return fam, sorted(rot & members)[0]
    return None, unit


def analyze_tract(seq, tract_offset, tract_len, fam):
    """Count perfect copies of the repeat unit inside the tract.

    The genome strand may carry the repeat in any rotational frame; we test all
    rotations of the family's members on both strands and keep the best.
    """
    tract = seq[tract_offset:tract_offset + tract_len]
    units = ({"GC", "CG"} if fam == "GC" else {"GCC", "CCG", "CGC", "GGC", "CGG", "GCG"})
    best = None
    for strand in ("+", "-"):
        t = tract if strand == "+" else revcomp(tract)
        for unit in units:
            # count leading perfect copies from frame 0
            n = 0
            while n + len(unit) <= len(t) and t[n:n + len(unit)] == unit:
                n += len(unit)
            # also try shifting start by 1 or 2 (frame rotation inside the tract)
            for shift in range(len(unit)):
                n2 = 0
                i = shift
                while i + len(unit) <= len(t) and t[i:i + len(unit)] == unit:
                    n2 += len(unit)
                    i += len(unit)
                n_copy = n2 // len(unit)
                if best is None or n_copy > best["perfect"]:
                    best = {"strand": strand, "unit": unit,
                            "perfect": n_copy, "bp_pure": n2,
                            "purity": round(n2 / len(tract) * 100, 1)}
    best["tract_len"] = len(tract)
    best["gc"] = round((t.count("G") + t.count("C")) / len(tract) * 100, 1)
    return best


# ---------------------------------------------------------------- primer design

GLOBAL_ARGS = {
    "PRIMER_OPT_SIZE": 22, "PRIMER_MIN_SIZE": 20, "PRIMER_MAX_SIZE": 28,
    "PRIMER_OPT_TM": 60.0, "PRIMER_MIN_TM": 57.0, "PRIMER_MAX_TM": 64.0,
    "PRIMER_MIN_GC": 35.0, "PRIMER_OPT_GC": 50.0, "PRIMER_MAX_GC": 65.0,
    "PRIMER_SALT_MONOVALENT": 50.0, "PRIMER_SALT_DIVALENT": 1.5,
    "PRIMER_DNA_CONC": 50.0,
    "PRIMER_MAX_SELF_ANY": 45.0, "PRIMER_MAX_SELF_END": 35.0,
    "PRIMER_MAX_HAIRPIN": 45.0,
    "PRIMER_MAX_END_STABILITY": 6.0,
    "PRIMER_MAX_POLY_X": 5,
    "PRIMER_NUM_RETURN": 10,
    "PRIMER_PICK_LEFT_PRIMER": 1, "PRIMER_PICK_RIGHT_PRIMER": 1,
    "PRIMER_PICK_INTERNAL_OLIGO": 0,
}

PRODUCT_RANGES = [[[250, 900]], [[200, 1000]], [[200, 1600]], [[200, 2000]]]


def design_one(locus, cand, tract_info, flank=800):
    """Try successive product-size ranges until primer3 returns a pair."""
    chrom, start, end = locus["chrom"], locus["start"], locus["end"]
    t_start, t_end = start - flank, end + flank
    tmpl = get_template(chrom, t_start, t_end)
    # re-locate tract inside template (tract itself, ±2 bp tolerance)
    off = start - t_start
    tract_offset = max(0, off - 2)
    tract_len = min(len(tmpl), off + (end - start) + 2) - tract_offset
    if off < flank - 20 or len(tmpl) < (end - start) + 200:
        flank = min(flank, off, len(tmpl) - (end - start) - 100)
        if flank < 150:
            return None, "template too short"
        return design_one(locus, cand, tract_info, flank)

    seq_args = {
        "SEQUENCE_ID": f"{chrom}_{start}_{end}",
        "SEQUENCE_TEMPLATE": tmpl,
        "SEQUENCE_EXCLUDED_REGION": [tract_offset, tract_len],
    }
    tract_info = analyze_tract(tmpl, tract_offset, tract_len, fam)
    for rng in PRODUCT_RANGES:
        ga = dict(GLOBAL_ARGS)
        ga["PRIMER_PRODUCT_SIZE_RANGE"] = rng
        res = primer3.bindings.design_primers(seq_args, ga)
        n = res.get("PRIMER_PAIR_NUM_RETURNED", 0)
        if n > 0:
            best = pick_best_pair(res)
            if best is not None:
                # exact amplicon coordinates (0-based half-open on the chromosome)
                lp, ll = best["left_pos"]
                rp, rl = best["right_pos"]
                best["product_size"] = rp + rl - lp
                best["amplicon_start"] = t_start + lp
                best["amplicon_end"] = t_start + rp + rl
                return (best, rng, tract_info), None
    # GC-rich fallback: allow higher primer GC% and lower Tm floor
    for rng in [[[250, 900]], [[200, 1600]]]:
        ga = dict(GLOBAL_ARGS)
        ga["PRIMER_PRODUCT_SIZE_RANGE"] = rng
        ga["PRIMER_MAX_GC"] = 75.0
        ga["PRIMER_MIN_GC"] = 30.0
        ga["PRIMER_MIN_TM"] = 56.0
        ga["PRIMER_MAX_POLY_X"] = 6
        res = primer3.bindings.design_primers(seq_args, ga)
        n = res.get("PRIMER_PAIR_NUM_RETURNED", 0)
        if n > 0:
            best = pick_best_pair(res)
            if best is not None:
                lp, ll = best["left_pos"]
                rp, rl = best["right_pos"]
                best["product_size"] = rp + rl - lp
                best["amplicon_start"] = t_start + lp
                best["amplicon_end"] = t_start + rp + rl
                best["fallback_gc"] = True
                return (best, rng, tract_info), None
    return None, "primer3 found no pair"


def pick_best_pair(res):
    """Prefer pairs whose primers stay clear of the repeat and are close to 60 C."""
    best, best_score = None, None
    for i in range(res["PRIMER_PAIR_NUM_RETURNED"]):
        lp, rp = res[f"PRIMER_LEFT_{i}"], res[f"PRIMER_RIGHT_{i}"]
        lt, rt = res[f"PRIMER_LEFT_{i}_TM"], res[f"PRIMER_RIGHT_{i}_TM"]
        penalty = res[f"PRIMER_PAIR_{i}_PENALTY"]
        tm_pen = abs(lt - 60) + abs(rt - 60)
        score = penalty + tm_pen * 2
        if best_score is None or score < best_score:
            best_score = score
            best = {
                "left_seq": res[f"PRIMER_LEFT_{i}_SEQUENCE"],
                "right_seq": res[f"PRIMER_RIGHT_{i}_SEQUENCE"],
                "left_pos": lp, "right_pos": rp,
                "left_tm": lt, "right_tm": rt,
                "left_gc": res[f"PRIMER_LEFT_{i}_GC_PERCENT"],
                "right_gc": res[f"PRIMER_RIGHT_{i}_GC_PERCENT"],
                "penalty": penalty,
            }
    return best


# ---------------------------------------------------------------- selection

def main():
    # load candidates
    cands = []
    with open(os.path.join(BASE, "data", "candidates.tsv")) as f:
        header = f.readline().rstrip("\n").split("\t")
        for line in f:
            cands.append(dict(zip(header, line.rstrip("\n").split("\t"))))
    for c in cands:
        c["start"], c["end"] = int(c["start"]), int(c["end"])
        c["copies"] = float(c["copies"])
        c["perMatch"] = float(c["perMatch"])
        c["copyNum_trf"] = float(c["copyNum_trf"])

    # plan: 4 loci per (motif, bucket) => 24 pairs
    plan = []
    for motif in ("GC", "GCC"):
        for bucket in ("20-25", "25-30", ">30"):
            pool = [c for c in cands if c["motif"] == motif and c["bucket"] == bucket]
            # prefer perfect (perMatch 100), spread across chromosomes
            pool.sort(key=lambda c: (-c["perMatch"], c["chrom"], abs(c["copies"] - {"20-25": 22.5, "25-30": 27.5, ">30": 40}[bucket])))
            used_chrom = set()
            picked = []
            for c in pool:
                if c["chrom"] in used_chrom and len(picked) < 4:
                    # allow same chrom only if necessary
                    continue
                # enforce min genomic spacing between picked loci on same chrom
                if any(p["chrom"] == c["chrom"] and abs(p["start"] - c["start"]) < 1e6 for p in picked):
                    continue
                picked.append(c)
                used_chrom.add(c["chrom"])
                if len(picked) == 4:
                    break
            for c in picked:
                plan.append(c)
            print(f"plan {motif} {bucket}: {len(picked)} loci (pool {len(pool)})", file=sys.stderr)

    print(f"total planned loci: {len(plan)}", file=sys.stderr)

    rows, raw = [], []
    for i, c in enumerate(plan):
        motif = c["motif"]
        fam, _ = canonical_motif(c["unit"])
        print(f"[{i + 1}/{len(plan)}] {c['chrom']}:{c['start']}-{c['end']} {motif}", file=sys.stderr)
        out = {"cand": c, "status": "not_run"}
        try:
            (pair, rng, tinfo), err = design_one(c, c, None)
        except Exception as e:
            (pair, rng, tinfo), err = None, None, repr(e)
        if pair is None:
            out["status"] = f"FAILED: {err}"
            print(f"   FAILED: {err}", file=sys.stderr)
        else:
            out["status"] = "ok"
            out["pair"] = pair
            out["product_range"] = rng
            out["tract"] = tinfo
            print(f"   ok  L={pair['left_seq']} ({pair['left_tm']:.1f}C)  "
                  f"R={pair['right_seq']} ({pair['right_tm']:.1f}C)  "
                  f"prod={pair['product_size']}bp  tract={tinfo['unit']}x{tinfo['perfect']} "
                  f"(purity {tinfo['purity']}%)  range={rng}", file=sys.stderr)
        raw.append(out)

    with open(os.path.join(BASE, "results", "designed.json"), "w") as f:
        json.dump(raw, f, indent=1)
    print("wrote results/designed.json", file=sys.stderr)


if __name__ == "__main__":
    main()
