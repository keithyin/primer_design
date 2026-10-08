#!/usr/bin/env python3
"""Genome-wide specificity check for designed primers.

Builds a pattern file (full primer sequences + 3'-end 15-mers, both strands),
runs GNU grep -oF over one-line-per-chromosome genome FASTA, and counts hits.

Writes results/specificity.tsv
"""
import json
import os
import subprocess
import sys

BASE = "/home/keith/projects/primer_design"
GENOME_ONE = os.path.join(BASE, "data", "hg38_oneline")
REVCOMP = str.maketrans("ACGT", "TGCA")


def rc(s):
    return s.translate(REVCOMP)[::-1]


def build_genome_oneline():
    """Flatten UCSC chromFa into one-line fasta files for grep."""
    src = os.path.join(BASE, "data", "hg38", "chromFa")
    os.makedirs(GENOME_ONE, exist_ok=True)
    for fn in sorted(os.listdir(src)):
        if not fn.endswith(".fa"):
            continue
        out = os.path.join(GENOME_ONE, fn)
        if os.path.exists(out):
            continue
        with open(os.path.join(src, fn)) as f, open(out, "w") as o:
            name = f.readline().strip()[1:].split()[0]
            seq = f.read().replace("\n", "").upper()
            o.write(f">{name}\n{seq}\n")
        print(f"flattened {name}", file=sys.stderr)


def collect_primers():
    with open(os.path.join(BASE, "results", "designed.json")) as f:
        raw = json.load(f)
    primers = []
    for item in raw:
        if item.get("status") not in ("ok", "ok_check"):
            continue
        c = item["locus"]
        p = item["pair"]
        locus_id = f"{c['chrom']}:{c['start']}_{c['motif']}"
        for side in ("left", "right"):
            seq = p[f"{side}_seq"]
            primers.append({"locus": locus_id, "side": side, "seq": seq})
    return primers


def main():
    build_genome_oneline()
    primers = collect_primers()
    print(f"{len(primers)} primers to check", file=sys.stderr)

    full_pats, end_pats = [], []
    for i, pr in enumerate(primers):
        pr["full_key"] = f"F{i}"
        pr["end_key"] = f"E{i}"
        full_pats.append(pr["seq"])
        full_pats.append(rc(pr["seq"]))
        end_pats.append(pr["seq"][-15:])
        end_pats.append(rc(pr["seq"][-15:]))

    pat_dir = os.path.join(BASE, "data")
    with open(os.path.join(pat_dir, "pat_full.txt"), "w") as f:
        f.write("\n".join(full_pats) + "\n")
    with open(os.path.join(pat_dir, "pat_end.txt"), "w") as f:
        f.write("\n".join(end_pats) + "\n")

    genome_files = sorted(os.listdir(GENOME_ONE))

    def count_matches(pat_file):
        counts, by_chrom = {}, {}
        cmd = ["grep", "-oFf", pat_file] + [os.path.join(GENOME_ONE, g) for g in genome_files]
        proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, text=True, bufsize=1 << 20)
        for line in proc.stdout:
            # with multiple files grep prefixes "path:match"; keep both parts
            fname, _, match = line.strip().rpartition(":")
            if not match:
                continue
            chrom = os.path.basename(fname).split(".")[0]
            counts[match] = counts.get(match, 0) + 1
            by_chrom[(match, chrom)] = by_chrom.get((match, chrom), 0) + 1
        proc.wait()
        return counts, by_chrom

    print("counting full-length matches ...", file=sys.stderr)
    full_counts, full_by_chrom = count_matches(os.path.join(pat_dir, "pat_full.txt"))
    print("counting 3'-15mer matches ...", file=sys.stderr)
    end_counts, end_by_chrom = count_matches(os.path.join(pat_dir, "pat_end.txt"))

    out = os.path.join(BASE, "results", "specificity.tsv")
    with open(out, "w") as f:
        f.write("locus\tside\tprimer\tfull_hits\tfull_hits_other_chrom\tend15_hits\tend15_hits_other_chrom\n")
        for i, pr in enumerate(primers):
            chrom = pr["locus"].split(":")[0]
            seq_f, seq_r = pr["seq"], rc(pr["seq"])
            fh = full_counts.get(seq_f, 0) + full_counts.get(seq_r, 0)
            fh_off = fh - (full_by_chrom.get((seq_f, chrom), 0) + full_by_chrom.get((seq_r, chrom), 0))
            eh = end_counts.get(seq_f[-15:], 0) + end_counts.get(seq_r[-15:], 0)
            eh_off = eh - (end_by_chrom.get((seq_f[-15:], chrom), 0) + end_by_chrom.get((seq_r[-15:], chrom), 0))
            f.write(f"{pr['locus']}\t{pr['side']}\t{pr['seq']}\t{fh}\t{fh_off}\t{eh}\t{eh_off}\n")
    print(f"wrote {out}", file=sys.stderr)


if __name__ == "__main__":
    main()
