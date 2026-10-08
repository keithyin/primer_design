#!/usr/bin/env python3
"""Pick the most specific primer pair per locus from primer3 candidates.

One genome-wide grep pass over all candidate primers (full seq + 3' 15-mer,
both strands), then per locus choose the pair minimizing off-target risk.

Writes results/designed.json
"""
import json
import os
import subprocess
import sys

BASE = "/home/keith/projects/primer_design"
GENOME_ONE = os.path.join(BASE, "data", "hg38_oneline")
RC = str.maketrans("ACGT", "TGCA")


def rc(s):
    return s.translate(RC)[::-1]


def main():
    with open(os.path.join(BASE, "results", "candidates.json")) as f:
        raw = json.load(f)

    genome_files = sorted(os.listdir(GENOME_ONE))

    def grep_counts(patterns):
        pat_file = os.path.join(BASE, "data", "pat_tmp.txt")
        with open(pat_file, "w") as f:
            f.write("\n".join(patterns) + "\n")
        counts, by_chrom = {}, {}
        cmd = ["grep", "-oFf", pat_file] + [os.path.join(GENOME_ONE, g) for g in genome_files]
        proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.DEVNULL,
                                text=True, bufsize=1 << 20)
        for line in proc.stdout:
            fname, _, match = line.strip().rpartition(":")
            if not match:
                continue
            chrom = os.path.basename(fname).split(".")[0]
            counts[match] = counts.get(match, 0) + 1
            by_chrom[(match, chrom)] = by_chrom.get((match, chrom), 0) + 1
        proc.wait()
        return counts, by_chrom

    full_pats, end_pats = set(), set()
    for item in raw:
        if item.get("status") != "ok":
            continue
        for p in item["candidates"]:
            for s in (p["left_seq"], p["right_seq"]):
                full_pats.add(s)
                full_pats.add(rc(s))
                end_pats.add(s[-15:])
                end_pats.add(rc(s[-15:]))
    full_pats, end_pats = sorted(full_pats), sorted(end_pats)
    print(f"{len(full_pats)} full patterns, {len(end_pats)} 3'-15mer patterns", file=sys.stderr)

    print("grep pass: full sequences ...", file=sys.stderr)
    fc, fbc = grep_counts(full_pats)
    print("grep pass: 3'-15mers ...", file=sys.stderr)
    ec, ebc = grep_counts(end_pats)

    out = []
    for item in raw:
        if item.get("status") != "ok":
            out.append({"locus": item["locus"], "pair": None, "status": "failed"})
            continue
        locus = item["locus"]
        chrom = locus["chrom"]
        best, best_score = None, None
        for p in item["candidates"]:
            side = {}
            for sd, sq in (("left", p["left_seq"]), ("right", p["right_seq"])):
                tf = fc.get(sq, 0) + fc.get(rc(sq), 0)
                of_ = fbc.get((sq, chrom), 0) + fbc.get((rc(sq), chrom), 0)
                te = ec.get(sq[-15:], 0) + ec.get(rc(sq[-15:]), 0)
                oe = ebc.get((sq[-15:], chrom), 0) + ebc.get((rc(sq[-15:]), chrom), 0)
                side[sd] = {"full_off": tf - of_, "end15_off": te - oe}
            fo = side["left"]["full_off"] + side["right"]["full_off"]
            eo = side["left"]["end15_off"] + side["right"]["end15_off"]
            expansion = locus.get("expansion_row", False)
            prod_pen = p["product_size"] / (400.0 if expansion else 1000.0)
            score = min(fo, 5) * 400 + min(eo, 60) * 8 + p["penalty"] + prod_pen * 20
            p2 = dict(p)
            p2["full_off"] = fo
            p2["end15_off"] = eo
            p2["left_full_off"] = side["left"]["full_off"]
            p2["left_end15_off"] = side["left"]["end15_off"]
            p2["right_full_off"] = side["right"]["full_off"]
            p2["right_end15_off"] = side["right"]["end15_off"]
            if best_score is None or score < best_score:
                best_score = score
                best = p2
        if best is None:
            out.append({"locus": locus, "pair": None, "status": "failed"})
            continue
        st = "ok"
        if best["full_off"] > 0 or best["end15_off"] > 10:
            st = "ok_check"
        out.append({"locus": locus, "pair": best, "status": st})
        print(f"{locus['chrom']}:{locus['start']} {locus['motif']}x{locus['n_units']} -> "
              f"prod={best['product_size']}bp full_off={best['full_off']} "
              f"end15_off={best['end15_off']} [{st}]", file=sys.stderr)

    with open(os.path.join(BASE, "results", "designed.json"), "w") as f:
        json.dump(out, f, indent=1)
    n_ok = sum(1 for o in out if o["status"] == "ok")
    print(f"wrote results/designed.json ({n_ok} clean, {len(out) - n_ok} flagged)", file=sys.stderr)


if __name__ == "__main__":
    main()
