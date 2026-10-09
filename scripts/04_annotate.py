#!/usr/bin/env python3
"""Annotate repeat loci with gene context via Ensembl REST (GRCh38).

Writes results/annotation.tsv
"""
import json
import os
import sys
import time
import urllib.request

BASE = os.environ.get("PRIMER_DESIGN_BASE") or os.path.dirname(
    os.path.dirname(os.path.abspath(__file__)))


def ens_request(url, retries=3):
    for i in range(retries):
        try:
            req = urllib.request.Request(url, headers={"Content-Type": "application/json",
                                                       "User-Agent": "primer-design-panel/1.0"})
            with urllib.request.urlopen(req, timeout=30) as r:
                return json.load(r)
        except Exception as e:
            print(f"  retry {i + 1}: {e}", file=sys.stderr)
            time.sleep(3)
    return None


def annotate(chrom, start, end):
    """0-based half-open -> 1-based closed for Ensembl."""
    s1, e1 = start + 1, end
    region = f"{chrom}:{s1}-{e1}"
    genes = ens_request(f"https://rest.ensembl.org/overlap/region/human/{region}?feature=gene")
    if genes is None:
        return {"region": region, "gene": ".", "biotype": ".", "context": "annotation-failed"}
    if genes:
        # prefer protein_coding, smallest overlap
        def gkey(g):
            return (0 if g.get("biotype") == "protein_coding" else 1,
                    (min(end, g["end"]) - max(start + 1, g["start"])))
        g = sorted(genes, key=gkey)[0]
        symbol = g.get("external_name") or g["gene_id"]
        return {"region": region, "gene": symbol,
                "biotype": g.get("biotype", "."), "context": "gene-body", "gene_id": g["gene_id"]}
    # no overlapping gene: find nearest within 100 kb
    lo = max(1, s1 - 100000)
    hi = e1 + 100000
    genes = ens_request(f"https://rest.ensembl.org/overlap/region/human/{chrom}:{lo}-{hi}?feature=gene")
    if not genes:
        return {"region": region, "gene": ".", "biotype": ".", "context": "intergenic(>100kb)"}
    best, bestd = None, None
    for g in genes:
        d = max(0, s1 - g["end"], g["start"] - e1)
        if bestd is None or d < bestd:
            best, bestd = g, d
    symbol = best.get("external_name") or best["gene_id"]
    direction = "upstream" if best["start"] > e1 else "downstream"
    return {"region": region, "gene": symbol, "biotype": best.get("biotype", "."),
            "context": f"intergenic({direction} {bestd}bp of {symbol})", "gene_id": best["gene_id"]}


def main():
    with open(os.path.join(BASE, "results", "designed.json")) as f:
        raw = json.load(f)
    out = os.path.join(BASE, "results", "annotation.tsv")
    with open(out, "w") as f:
        f.write("chrom\tstart\tend\tmotif\tregion\tgene\tbiotype\tcontext\n")
        for item in raw:
            c = item["locus"]
            ann = annotate(c["chrom"], c["start"], c["end"])
            print(f"{c['chrom']}:{c['start']}-{c['end']} -> {ann['gene']} ({ann['context']})", file=sys.stderr)
            f.write(f"{c['chrom']}\t{c['start']}\t{c['end']}\t{c['motif']}\t"
                    f"{ann['region']}\t{ann['gene']}\t{ann['biotype']}\t{ann['context']}\n")
            time.sleep(0.3)
    print(f"wrote {out}", file=sys.stderr)


if __name__ == "__main__":
    main()
