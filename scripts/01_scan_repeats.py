#!/usr/bin/env python3
"""Scan hg38 simpleRepeat (TRF) annotations genome-wide for GC-rich tandem repeats.

Targets:
  - dinucleotide (GC)n / (CG)n, copy number >= 20
  - trinucleotide pure-GC repeats ((GCC)n family: GCC/CCG/CGC/GGC/CGG/GCG), copy number >= 20
Outputs candidates with bucket assignment (20-25 / 25-30 / >30).
"""
import json
import time
import urllib.request

CHROMS = [f"chr{i}" for i in range(1, 23)] + ["chrX", "chrY"]
OUT = "/home/keith/projects/primer_design/data/candidates.tsv"


def fetch_track(chrom):
    url = f"https://api.genome.ucsc.edu/getData/track?genome=hg38;track=simpleRepeat;chrom={chrom}"
    for attempt in range(3):
        try:
            with urllib.request.urlopen(url, timeout=120) as r:
                return json.load(r)
        except Exception as e:
            print(f"  retry {attempt + 1} for {chrom}: {e}")
            time.sleep(5)
    return None


def bucket(copy_num):
    if 20 <= copy_num < 25:
        return "20-25"
    if 25 <= copy_num < 30:
        return "25-30"
    return ">30"


def is_gc_dinuc(seq):
    return len(seq) == 2 and set(seq) <= {"G", "C"}


def is_gc_trinuc(seq):
    return len(seq) == 3 and set(seq) <= {"G", "C"}


def main():
    rows = []
    for chrom in CHROMS:
        print(f"scanning {chrom} ...", flush=True)
        data = fetch_track(chrom)
        if not data or "simpleRepeat" not in data:
            print(f"  no data for {chrom}")
            continue
        for tr in data["simpleRepeat"]:
            period = tr.get("period", 0)
            seq = tr.get("sequence", "").upper()
            copy_num = tr.get("copyNum", 0)
            per_match = tr.get("perMatch", 0)
            if copy_num < 20 or per_match < 90:
                continue
            motif = None
            if period == 2 and is_gc_dinuc(seq):
                motif = "GC"
            elif period == 3 and is_gc_trinuc(seq):
                motif = "GCC"
            if motif is None:
                continue
            start, end = tr["chromStart"], tr["chromEnd"]
            tract_len = end - start
            # integer copies consistent with annotated copy number
            n_copies = tract_len / len(seq)
            if n_copies < 20:
                continue
            rows.append({
                "chrom": chrom,
                "start": start,
                "end": end,
                "motif": motif,
                "unit": seq,
                "copies": round(n_copies, 1),
                "copyNum_trf": copy_num,
                "perMatch": per_match,
                "bucket": bucket(copy_num),
            })
        time.sleep(1)

    rows.sort(key=lambda r: (r["motif"], r["bucket"], r["chrom"], r["start"]))
    with open(OUT, "w") as f:
        f.write("chrom\tstart\tend\tmotif\tunit\tcopies\tcopyNum_trf\tperMatch\tbucket\n")
        for r in rows:
            f.write(f"{r['chrom']}\t{r['start']}\t{r['end']}\t{r['motif']}\t{r['unit']}\t"
                    f"{r['copies']}\t{r['copyNum_trf']}\t{r['perMatch']}\t{r['bucket']}\n")

    # summary
    from collections import Counter
    c = Counter((r["motif"], r["bucket"]) for r in rows)
    print("\n=== candidate counts (motif, bucket) ===")
    for k in sorted(c):
        print(k, c[k])
    print(f"total candidates: {len(rows)} -> {OUT}")


if __name__ == "__main__":
    main()
