#!/usr/bin/env python3
"""Mismatch-tolerant primer-PAIR in-silico PCR, genome-wide (report section 6.4).

This is the check the panel was missing.  Exact-match uniqueness says nothing
about mispriming: a primer does not need a perfect match to be extended.  What
produces an extra band is a primer PAIR -- two sites, each tolerated by its
primer, facing each other at a workable distance.

Method.  For every contig, in both strands:
  1. find every site where a primer's 3'-terminal `seed` bases match exactly
     (seeding is what makes a 3.2 Gbp mismatch search tractable);
  2. extend each seed to the full primer and count mismatches; keep sites with
     <= `mismatches` total, and record the contiguous exact run at the 3' end;
  3. pair F and R sites that face each other with a product length in range,
     in both configurations (plus-strand and minus-strand products).

A single rolling-code pass per contig serves every seed length: for a window of
length S with base-4 code c, the prefix of length s is c // 4**(S-s) and the
suffix is c % 4**s, so shorter seeds need no extra scan.

The on-target product of every pair is a built-in control: it must come out
perfect (0 mismatches, both 3' runs = full primer length, both coordinates
matching the deliverable).  Anything else that pairs is a putative off-target.

Writes results/mismatch_pcr.tsv and results/mismatch_pcr.json.

Usage: python3 scripts/18_mismatch_pcr.py [--configs 15:3,12:4,10:6]
"""
import argparse
import csv
import json
import os
import sys
import time

import numpy as np

BASE = os.environ.get("PRIMER_DESIGN_BASE") or os.path.dirname(
    os.path.dirname(os.path.abspath(__file__)))
GENOME_ONE = os.path.join(BASE, "data", "hg38_oneline")
CHUNK = 20_000_000
MIN_PRODUCT = 100
MAX_PRODUCT = 3000
MAX_HITS_PER_PRIMER = 400_000

RC = str.maketrans("ACGT", "TGCA")
CODE = np.full(256, 4, dtype=np.uint8)
for _i, _ch in enumerate("ACGT"):
    CODE[ord(_ch)] = _i
    CODE[ord(_ch.lower())] = _i


def rc(s):
    return s.translate(RC)[::-1]


def kmer_code(s):
    c = 0
    for ch in s:
        v = int(CODE[ord(ch)])
        if v > 3:
            raise ValueError("non-ACGT in seed")
        c = c * 4 + v
    return c


def load_primers():
    p = os.path.join(BASE, "results", "primer_design_results.tsv")
    with open(p) as f:
        rows = list(csv.DictReader(f, delimiter="\t"))
    return [{"no": int(r["no"]), "gene": r["gene"], "chrom": r["chrom"],
             "fwd": r["fwd"].upper(), "rev": r["rev"].upper(),
             "amp_start": int(r["amplicon_start"]),
             "amp_end": int(r["amplicon_end"]),
             "product": int(r["product_size"])} for r in rows]


def verify(seq, p, primer, orient, max_mm):
    """Full-primer alignment at p; returns (n_mismatch, 3'-exact-run) or None."""
    L = len(primer)
    win = seq[p:p + L]
    if len(win) < L:
        return None
    target = primer if orient == "+" else rc(primer)
    mm = 0
    for a, b in zip(win, target):
        if a != b:
            mm += 1
            if mm > max_mm:
                return None
    # 3' end of the primer: last base for '+', first base of the genome match
    # for '-' (there genome[p] pairs with the primer's 3'-terminal base).
    run = 0
    rng = range(L - 1, -1, -1) if orient == "+" else range(L)
    for k in rng:
        if win[k] == target[k]:
            run += 1
        else:
            break
    return mm, run


def build_seed_maps(primers, seed):
    """Sorted query codes plus code -> [(pair_idx, side, orient), ...].

    Lookups are vectorised with searchsorted over the sorted code array, so the
    per-position Python loop only ever runs over actual hits -- iterating every
    window in a 3.2 Gbp genome in Python is not viable.
    """
    entries = {}
    for idx, pr in enumerate(primers):
        for side in ("fwd", "rev"):
            seed_str = pr[side][-seed:]
            entries.setdefault(kmer_code(seed_str), []).append((idx, side, "+"))
            entries.setdefault(kmer_code(rc(seed_str)), []).append((idx, side, "-"))
    qarr = np.array(sorted(entries), dtype=np.uint32)
    return qarr, entries


def scan_contig(seq, primers, seed, max_mm, hits, contig, qarr, entries):
    """Append every primer site with an exact `seed`-base 3' anchor."""
    if len(seq) < seed:
        return
    arr = CODE[np.frombuffer(seq.encode(), dtype=np.uint8)]
    for cs0 in range(0, len(seq) - seed + 1, CHUNK):
        c = arr[cs0:cs0 + CHUNK + seed - 1]
        n = len(c) - seed + 1
        if n <= 0:
            continue
        code = c[:n].astype(np.uint32)
        for j in range(1, seed):
            code = code * np.uint32(4) + c[j:j + n].astype(np.uint32)
        bad = (c > 3).astype(np.int32)
        cs = np.concatenate((np.zeros(1, dtype=np.int32), np.cumsum(bad)))
        valid = (cs[seed:seed + n] - cs[:n]) == 0
        if not valid.any():
            continue
        vpos = np.nonzero(valid)[0]
        sub = code[vpos]
        o = np.searchsorted(qarr, sub)
        np.clip(o, 0, len(qarr) - 1, out=o)
        hm = qarr[o] == sub
        if not hm.any():
            continue
        for g_local, cd in zip(vpos[hm].tolist(), sub[hm].tolist()):
            g = cs0 + g_local
            for idx, side, orient in entries[cd]:
                q = primers[idx][side]
                p = g - (len(q) - seed) if orient == "+" else g
                if p < 0:
                    continue
                v = verify(seq, p, q, orient, max_mm)
                if v:
                    hits[(idx, side, orient)].append((contig, p, v[0], v[1]))


def pair_products(primers, hits, cfg_name, out_rows):
    """Pair facing F/R sites into candidate amplicons, both strand configs."""
    n_off = 0
    for idx, pr in enumerate(primers):
        Lf, Lr = len(pr["fwd"]), len(pr["rev"])
        # configuration 1: plus-strand product.  F on '+', R on '-'.
        by_contig = {}
        for contig, p, mm, run in hits.get((idx, "fwd", "+"), ()):
            by_contig.setdefault(contig, []).append((p, mm, run))
        rmap = {}
        for contig, q, mm, run in hits.get((idx, "rev", "-"), ()):
            rmap.setdefault(contig, []).append((q, mm, run))
        cands = []
        for contig, fl in by_contig.items():
            rl = sorted(rmap.get(contig, ()))
            if not rl:
                continue
            for p, mmf, runf in fl:
                for q, mmr, runr in rl:
                    length = q + Lr - p
                    if length < MIN_PRODUCT or length > MAX_PRODUCT:
                        continue
                    if q < p + Lf:
                        continue
                    cands.append((contig, p, q + Lr, length, mmf, mmr,
                                  runf, runr, "+"))
        # configuration 2: minus-strand product.  F on '-', R on '+'.
        by_contig2 = {}
        for contig, p, mm, run in hits.get((idx, "fwd", "-"), ()):
            by_contig2.setdefault(contig, []).append((p, mm, run))
        rmap2 = {}
        for contig, q, mm, run in hits.get((idx, "rev", "+"), ()):
            rmap2.setdefault(contig, []).append((q, mm, run))
        for contig, fl in by_contig2.items():
            rl = rmap2.get(contig, ())
            if not rl:
                continue
            for p, mmf, runf in fl:
                for q, mmr, runr in rl:
                    length = p + Lf - q
                    if length < MIN_PRODUCT or length > MAX_PRODUCT:
                        continue
                    if p < q + Lr:
                        continue
                    cands.append((contig, q, p + Lf, length, mmf, mmr,
                                  runf, runr, "-"))

        for contig, s, e, length, mmf, mmr, runf, runr, orient in cands:
            # s/e are 0-based half-open (UCSC/BED style, as in designed.json);
            # the deliverable TSV is 1-based closed, hence amp_start - 1.
            on_target = (mmf == 0 and mmr == 0 and contig == pr["chrom"]
                         and s == pr["amp_start"] - 1 and e == pr["amp_end"])
            if not on_target:
                n_off += 1
            out_rows.append({
                "config": cfg_name, "no": pr["no"], "gene": pr["gene"],
                "kind": "on_target" if on_target else "off_target",
                "contig": contig, "start": s, "end": e,
                "product_len": length, "fwd_mm": mmf, "rev_mm": mmr,
                "fwd_3run": runf, "rev_3run": runr, "strand": orient})
    return n_off


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--configs", default="15:3,12:4,10:6",
                    help="comma-separated seed:mismatches settings")
    ap.add_argument("--genome", default=GENOME_ONE,
                    help="directory of header+single-line FASTA files")
    ap.add_argument("--max-contigs", type=int, default=0,
                    help="debug: only scan the first N contigs")
    ap.add_argument("--out-prefix", default=os.path.join(BASE, "results", "mismatch_pcr"),
                    help="output path prefix ('.tsv' / '.json' appended)")
    args = ap.parse_args()
    configs = []
    for c in args.configs.split(","):
        s, m = c.split(":")
        configs.append((int(s), int(m)))

    genome_dir = args.genome
    primers = load_primers()
    files = sorted(f for f in os.listdir(genome_dir) if f.endswith(".fa"))
    if args.max_contigs:
        files = files[:args.max_contigs]
    scanned_contigs = {f[:-3] for f in files}
    print(f"{len(primers)} pairs, {len(configs)} configs, {len(files)} contigs "
          f"from {genome_dir}", file=sys.stderr)

    out_rows = []
    summary = {"configs": [], "product_range": [MIN_PRODUCT, MAX_PRODUCT],
               "on_target_found": {}, "scanned_files": len(files),
               "genome_dir": genome_dir}
    total_bp = 0
    t0 = time.time()
    for seed, max_mm in configs:
        hits = {}
        for idx in range(len(primers)):
            for side in ("fwd", "rev"):
                for orient in ("+", "-"):
                    hits[(idx, side, orient)] = []
        scanned = 0
        qarr, entries = build_seed_maps(primers, seed)
        for fi, fn in enumerate(files):
            contig = fn[:-3] if fn.endswith(".fa") else fn
            with open(os.path.join(genome_dir, fn)) as f:
                f.readline()
                seq = f.read().strip().upper()
            scanned += len(seq)
            scan_contig(seq, primers, seed, max_mm, hits, contig,
                        qarr, entries)
            if (fi + 1) % 60 == 0:
                n = sum(len(v) for v in hits.values())
                print(f"  seed={seed} mm={max_mm}: {fi+1}/{len(files)} contigs, "
                      f"{scanned/1e9:.2f} Gbp, {n} sites, "
                      f"{time.time()-t0:.0f}s", file=sys.stderr)
        total_bp = scanned

        cfg = f"seed{seed}_mm{max_mm}"
        n_off = pair_products(primers, hits, cfg, out_rows)
        per_primer = [sum(len(hits[(i, s, o)]) for s in ("fwd", "rev")
                          for o in ("+", "-")) for i in range(len(primers))]
        summary["configs"].append({
            "config": cfg, "seed": seed, "max_mismatches": max_mm,
            "sites_total": sum(per_primer),
            "sites_max_per_pair": max(per_primer),
            "sites_median_per_pair": sorted(per_primer)[len(per_primer) // 2],
            "putative_offtarget_products": n_off})
        print(f"  seed={seed} mm={max_mm}: {sum(per_primer)} sites total, "
              f"{n_off} putative off-target products", file=sys.stderr)

    # per-pair on-target recovery, using the strictest config
    strict = f"seed{configs[0][0]}_mm{configs[0][1]}"
    for pr in primers:
        if pr["chrom"] not in scanned_contigs:
            # reference incomplete for this locus: not a failure, just untested
            summary["on_target_found"][pr["no"]] = {
                "gene": pr["gene"], "n": None, "coord_match": None,
                "note": f"{pr['chrom']} not in scanned set"}
            continue
        got = [r for r in out_rows if r["config"] == strict
               and r["no"] == pr["no"] and r["kind"] == "on_target"]
        summary["on_target_found"][pr["no"]] = {
            "gene": pr["gene"], "n": len(got),
            "coord_match": bool(got and got[0]["start"] == pr["amp_start"] - 1
                                and got[0]["end"] == pr["amp_end"]),
            "product_len": got[0]["product_len"] if got else None,
            "expected_len": pr["product"]}

    cols = ["config", "no", "gene", "kind", "contig", "start", "end",
            "product_len", "fwd_mm", "rev_mm", "fwd_3run", "rev_3run", "strand"]
    with open(args.out_prefix + ".tsv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=cols, delimiter="\t")
        w.writeheader()
        for r in out_rows:
            w.writerow(r)

    summary["scanned_bp"] = total_bp
    summary["elapsed_s"] = round(time.time() - t0, 1)
    with open(args.out_prefix + ".json", "w") as f:
        json.dump(summary, f, indent=1, ensure_ascii=False)

    print(f"done in {summary['elapsed_s']}s; "
          f"{sum(1 for r in out_rows if r['kind']=='off_target')} off-target "
          f"candidates in {args.out_prefix}.tsv")
    for c in summary["configs"]:
        print(f"  {c['config']}: {c['sites_total']} sites, "
              f"{c['putative_offtarget_products']} off-target products")
    bad = [k for k, v in summary["on_target_found"].items()
           if v["n"] is not None and (v["n"] != 1 or not v["coord_match"])]
    untested = [k for k, v in summary["on_target_found"].items()
                if v["n"] is None]
    if untested:
        print(f"NOTE: on-target control not testable for {len(untested)} pairs "
              f"(contig not scanned)", file=sys.stderr)
    if bad:
        print(f"WARNING: on-target control failed for pairs {bad}",
              file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    main()
