#!/usr/bin/env python3
"""Method D: genome-wide specificity rescan with numpy rolling-hash coding.

Completely independent of the design-time grep -F pass (different tool,
different algorithm). Scans every contig of the reference (same file set as
the grep pass) for:
  - all 3'-end 15-mers of the delivered primers (both strands) -> end15 counts
  - full primer matches, recovered by extending each 15-mer hit -> full counts

Base-4 coding of 15-mers is injective (4^15 < 2^32), windows spanning non-ACGT
bases are excluded via a cumsum mask, and every hit is verified by direct
string comparison.

Writes results/method_D_genomescan.json
"""
import json
import os
import sys
import time

import numpy as np

BASE = "/home/keith/projects/primer_design"
GENOME_ONE = os.path.join(BASE, "data", "hg38_oneline")
K = 15
CHUNK = 20_000_000

RC = str.maketrans("ACGT", "TGCA")


def rc(s):
    return s.translate(RC)[::-1]


CODE = np.full(256, 4, dtype=np.uint8)
for i, ch in enumerate("ACGT"):
    CODE[ord(ch)] = i
    CODE[ord(ch.lower())] = i


def kmer_code(s):
    c = 0
    for ch in s:
        v = int(CODE[ord(ch)])
        if v > 3:
            raise ValueError("kmer contains non-ACGT")
        c = c * 4 + v
    return c


def main():
    # ---- collect primers from the deliverable ----
    with open(os.path.join(BASE, "results", "designed.json")) as f:
        raw = [it for it in json.load(f) if it.get("status") in ("ok", "ok_check")]
    primers = []  # (no, side, seq)
    for it in raw:
        r = it["locus"]
        for side, key in (("fwd", "left_seq"), ("rev", "right_seq")):
            primers.append({"no": len(primers) // 2 + 1, "side": side,
                            "seq": it["pair"][key], "chrom": r["chrom"]})

    # pattern code -> list of (primer_idx, orientation)
    qmap = {}
    for idx, pr in enumerate(primers):
        km = pr["seq"][-K:]
        for orient, s in (("+", km), ("-", rc(km))):
            qmap.setdefault(kmer_code(s), []).append((idx, orient))
    codes_arr = np.array(sorted(qmap.keys()), dtype=np.uint64)
    print(f"{len(primers)} primers, {len(qmap)} distinct 15-mer patterns", file=sys.stderr)

    counts15 = {}   # (idx, chrom) -> hits
    counts_full = {}  # (idx, chrom) -> hits

    files = sorted(os.listdir(GENOME_ONE))
    t0 = time.time()
    scanned_bp = 0
    for fi, fn in enumerate(files):
        chrom = fn.split(".")[0]
        with open(os.path.join(GENOME_ONE, fn)) as f:
            f.readline()
            seq = f.read().strip()
        if len(seq) < K:
            continue
        arr = CODE[np.frombuffer(seq.encode(), dtype=np.uint8)]
        scanned_bp += len(seq)
        for cs0 in range(0, len(arr) - K + 1, CHUNK):
            c = arr[cs0:cs0 + CHUNK + K - 1]
            n = len(c) - K + 1
            if n <= 0:
                continue
            code = c[:n].astype(np.uint32)
            for j in range(1, K):
                code = code * np.uint32(4) + c[j:j + n].astype(np.uint32)
            # windows containing non-ACGT are invalid
            nonacgt = (c > 3).astype(np.int32)
            cs = np.concatenate((np.zeros(1, dtype=np.int32), np.cumsum(nonacgt)))
            valid = (cs[K:K + n] - cs[:n]) == 0
            pos = np.nonzero(valid)[0]
            if pos.size == 0:
                continue
            pc = code[pos]
            order = np.searchsorted(codes_arr, pc)
            order_c = np.clip(order, 0, len(codes_arr) - 1)
            hitmask = codes_arr[order_c] == pc
            if not hitmask.any():
                continue
            hpos = pos[hitmask]
            hcodes = pc[hitmask]
            for g_local, cd in zip(hpos.tolist(), hcodes.tolist()):
                g = cs0 + g_local
                for idx, orient in qmap[cd]:
                    pr = primers[idx]
                    L = len(pr["seq"])
                    counts15[(idx, chrom)] = counts15.get((idx, chrom), 0) + 1
                    # extend to a full-primer match?
                    if orient == "+":
                        p = g - (L - K)
                        if p >= 0 and seq[p:p + L] == pr["seq"]:
                            counts_full[(idx, chrom)] = counts_full.get((idx, chrom), 0) + 1
                    else:
                        p = g
                        if p + L <= len(seq) and seq[p:p + L] == rc(pr["seq"]):
                            counts_full[(idx, chrom)] = counts_full.get((idx, chrom), 0) + 1
        if (fi + 1) % 40 == 0:
            print(f"  {fi + 1}/{len(files)} contigs, {scanned_bp / 1e9:.2f} Gbp, "
                  f"{time.time() - t0:.0f}s", file=sys.stderr)

    # ---- assemble per-primer results ----
    rows = []
    for it_idx, it in enumerate(raw):
        r = it["locus"]
        row = {"no": it_idx + 1, "id": r["id"], "chrom": r["chrom"]}
        for side in ("fwd", "rev"):
            idx = (it_idx + 1 - 1) * 2 + (0 if side == "fwd" else 1)
            total15 = sum(v for (i, c), v in counts15.items() if i == idx)
            own15 = counts15.get((idx, r["chrom"]), 0)
            totalf = sum(v for (i, c), v in counts_full.items() if i == idx)
            ownf = counts_full.get((idx, r["chrom"]), 0)
            row[side] = {"genome_full_hits": totalf,
                         "full_off": totalf - ownf,
                         "end15_hits": total15,
                         "end15_off": total15 - own15}
        rows.append(row)

    out = {"scanned_files": len(files), "scanned_bp": scanned_bp,
           "elapsed_s": round(time.time() - t0, 1), "rows": rows}
    with open(os.path.join(BASE, "results", "method_D_genomescan.json"), "w") as f:
        json.dump(out, f, indent=1)
    print(f"done in {time.time() - t0:.0f}s: {scanned_bp / 1e9:.2f} Gbp scanned", file=sys.stderr)
    n_off = sum(1 for r in rows for s in ("fwd", "rev") if r[s]["full_off"])
    print(f"primers with full-length off-target hits: {n_off}/{len(primers)}", file=sys.stderr)


if __name__ == "__main__":
    main()
