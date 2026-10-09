#!/usr/bin/env python3
"""Verify the reference under data/ really is UCSC hg38, base for base.

Why this exists: the panel's coordinates and every specificity count are
expressed in hg38, but the assembled working reference comes from two file sets
(GRCh38 no_alt chr1-22 held locally, and chrX fetched from Ensembl because
UCSC's bulk host rate-limited us mid-transfer).  Mixing references is exactly
the failure mode this project cannot afford -- a silent chromosome-level
mismatch would make every off-target count wrong without failing anything.

So the assembly is not trusted because of where it came from; it is checked
against UCSC directly, through the REST API, which is a different service from
the one that blocked us:

  R1  every chromosome's length equals UCSC's published hg38 length
  R2  random windows on every chromosome match UCSC hg38 exactly
  R3  all 24 panel loci (amplicon span) match UCSC hg38 exactly

R2/R3 are sampled, not exhaustive: an exact whole-genome diff is what the
blocked bulk download was for.  The sampling density and what it does and does
not rule out are reported alongside the result.

Writes results/reference_check.json.

Usage: python3 scripts/20_verify_reference.py [--windows 2] [--span 10000]
"""
import argparse
import json
import os
import random
import sys
import time
import urllib.error
import urllib.request

BASE = os.environ.get("PRIMER_DESIGN_BASE") or os.path.dirname(
    os.path.dirname(os.path.abspath(__file__)))
GENOME = os.path.join(BASE, "data", "hg38_oneline")
API = "https://api.genome.ucsc.edu/getData/sequence"
API_LIST = "https://api.genome.ucsc.edu/list/chromosomes?genome=hg38"


def api_get(url, retries=5):
    last = None
    for attempt in range(retries):
        try:
            with urllib.request.urlopen(url, timeout=90) as r:
                return json.loads(r.read().decode())
        except urllib.error.HTTPError as e:
            last = f"HTTP {e.code}"
            # 429 = over the rate limit; hammering harder only extends the block
            time.sleep(30 * (attempt + 1) if e.code == 429 else 2 ** attempt)
        except Exception as e:                          # noqa: BLE001
            last = repr(e)
            time.sleep(2 ** attempt)
    raise RuntimeError(f"{url} failed: {last}")


def api_seq(chrom, start, end):
    d = api_get(f"{API}?genome=hg38;chrom={chrom};start={start};end={end}")
    return d["dna"].upper()


def load(chrom):
    with open(os.path.join(GENOME, f"{chrom}.fa")) as f:
        f.readline()
        return f.read().strip().upper()


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--windows", type=int, default=2,
                    help="random windows checked per chromosome")
    ap.add_argument("--span", type=int, default=10000)
    ap.add_argument("--seed", type=int, default=20261008)
    args = ap.parse_args()

    if not os.path.isdir(GENOME):
        sys.exit(f"no reference at {GENOME}")
    chroms = sorted(f[:-3] for f in os.listdir(GENOME) if f.endswith(".fa"))
    print(f"reference: {GENOME} ({len(chroms)} contigs)")

    print("R1 lengths vs UCSC hg38")
    try:
        sizes = api_get(API_LIST)["chromosomes"]
    except RuntimeError as e:
        sys.exit(f"FATAL: UCSC API unreachable, cannot verify: {e}")
    len_bad = []
    for c in chroms:
        # header + single sequence line: subtract the header and both newlines
        with open(os.path.join(GENOME, f"{c}.fa")) as fh:
            hlen = len(fh.readline().rstrip("\n"))
        local = os.path.getsize(os.path.join(GENOME, f"{c}.fa")) - hlen - 2
        want = sizes.get(c)
        if want is None or local != want:
            len_bad.append(f"{c}: local {local} vs UCSC {want}")
    print(f"  [{'PASS' if not len_bad else 'FAIL'}] {len(chroms)} 条染色体长度与 UCSC hg38 一致"
          + (f" — {'; '.join(len_bad[:3])}" if len_bad else ""))

    print("R2 random windows vs UCSC hg38")
    rng = random.Random(args.seed)
    n_ok = n_bad = 0
    bad_detail = []
    compared_bp = 0
    for c in chroms:
        seq = load(c)
        for _ in range(args.windows):
            if len(seq) <= args.span:
                continue
            s = rng.randrange(0, len(seq) - args.span)
            ref = api_seq(c, s, s + args.span)
            if ref == seq[s:s + args.span]:
                n_ok += 1
                compared_bp += len(ref)
            else:
                n_bad += 1
                diff = sum(1 for a, b in zip(ref, seq[s:s + args.span]) if a != b)
                bad_detail.append(f"{c}:{s}-{s + args.span} ({diff} diffs)")
    print(f"  [{'PASS' if not n_bad else 'FAIL'}] 抽检 {n_ok + n_bad} 个 {args.span} bp 窗口"
          f"（{compared_bp / 1e3:.0f} kb 逐位比对）"
          + (f" — {'; '.join(bad_detail[:3])}" if bad_detail else ""))

    print("R3 panel loci vs UCSC hg38")
    import csv
    with open(os.path.join(BASE, "results", "primer_design_results.tsv")) as f:
        rows = list(csv.DictReader(f, delimiter="\t"))
    loc_bad = []
    cache = {}
    for r in rows:
        c = r["chrom"]
        if c not in cache:
            cache[c] = load(c)
        seq = cache[c]
        a0, a1 = int(r["amplicon_start"]) - 1, int(r["amplicon_end"])
        ref = api_seq(c, a0, a1)
        got = seq[a0:a1]
        if ref != got:
            loc_bad.append(f"#{r['no']} {r['gene']} {c}:{a0}-{a1}")
        elif not got.startswith(r["fwd"]):
            loc_bad.append(f"#{r['no']} {r['gene']} 引物位点不符")
    print(f"  [{'PASS' if not loc_bad else 'FAIL'}] 24 个位点扩增子区间与 UCSC hg38 及交付坐标一致"
          + (f" — {'; '.join(loc_bad[:3])}" if loc_bad else ""))

    out = {"reference": GENOME, "contigs": chroms,
           "R1_lengths_ok": not len_bad, "R1_detail": len_bad,
           "R2_windows_ok": n_ok, "R2_windows_bad": n_bad,
           "R2_compared_bp": compared_bp, "R2_span": args.span,
           "R2_detail": bad_detail,
           "R3_loci_ok": not loc_bad, "R3_detail": loc_bad,
           "sampling_note": (
               f"{n_ok + n_bad} windows x {args.span} bp = {compared_bp / 1e6:.2f} Mb "
               f"({100 * compared_bp / 3.088e9:.3f}% of the 3.088 Gbp assembly) "
               "plus all 24 loci; sampled, not an exhaustive whole-genome diff")}
    with open(os.path.join(BASE, "results", "reference_check.json"), "w") as f:
        json.dump(out, f, indent=1, ensure_ascii=False)

    print()
    if len_bad or n_bad or loc_bad:
        print("REFERENCE CHECK: FAILED")
        sys.exit(1)
    print(f"REFERENCE CHECK: PASSED — {out['sampling_note']}")


if __name__ == "__main__":
    main()
