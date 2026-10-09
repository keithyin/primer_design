#!/usr/bin/env python3
"""Design flanking primers for a fixed panel of repeat loci (data/panel_loci.json).

Each panel entry: {chrom, start, end, motif, unit, n_units, exact_units,
interrupt_units, purity, tract_len, bucket_label, note}
"""
import functools
import json
import os
import sys

import primer3

BASE = os.environ.get("PRIMER_DESIGN_BASE") or os.path.dirname(
    os.path.dirname(os.path.abspath(__file__)))
CHROM_FA = os.path.join(BASE, "data", "hg38", "chromFa")

RC = str.maketrans("ACGT", "TGCA")


def revcomp(s):
    return s.translate(RC)[::-1]


@functools.lru_cache(maxsize=6)
def load_chrom(chrom):
    with open(os.path.join(CHROM_FA, f"{chrom}.fa")) as f:
        f.readline()
        return f.read().replace("\n", "").upper()

GLOBAL_ARGS = {
    "PRIMER_OPT_SIZE": 22, "PRIMER_MIN_SIZE": 20, "PRIMER_MAX_SIZE": 28,
    "PRIMER_OPT_TM": 60.0, "PRIMER_MIN_TM": 57.0, "PRIMER_MAX_TM": 64.0,
    "PRIMER_MIN_GC": 35.0, "PRIMER_OPT_GC": 50.0, "PRIMER_MAX_GC": 65.0,
    "PRIMER_SALT_MONOVALENT": 50.0, "PRIMER_SALT_DIVALENT": 1.5,
    "PRIMER_DNA_CONC": 50.0,
    # --- oligo self-structure constraints -------------------------------
    # primer3-py runs with thermodynamic oligo alignment ON, so the plain
    # PRIMER_MAX_SELF_ANY / _SELF_END keys are INERT: only the _TH family is
    # consulted. PRIMER_MAX_HAIRPIN does not exist as a key at all (only
    # PRIMER_MAX_HAIRPIN_TH does) and unknown keys are silently ignored.
    # The values below are primer3's own defaults (libprimer3.c:
    # max_self_any_th / max_self_end_th / max_hairpin_th / pair_compl_*_th
    # are all 47.0), stated explicitly so the config says what it means.
    # Physical equivalents (fitted against calc_homodimer/calc_hairpin):
    #   SELF_ANY_TH 47 ~ homodimer dG -12.5 kcal/mol
    #   HAIRPIN_TH  47 ~ hairpin  dG  -0.8 kcal/mol
    "PRIMER_MAX_SELF_ANY_TH": 47.0, "PRIMER_MAX_SELF_END_TH": 47.0,
    "PRIMER_MAX_HAIRPIN_TH": 47.0,
    # pair-level (hetero-dimer) bounds; previously unset, i.e. these defaults
    "PRIMER_PAIR_MAX_COMPL_ANY_TH": 47.0, "PRIMER_PAIR_MAX_COMPL_END_TH": 47.0,
    "PRIMER_MAX_END_STABILITY": 6.0,
    "PRIMER_MAX_POLY_X": 5,
    "PRIMER_NUM_RETURN": 10,
    "PRIMER_PICK_LEFT_PRIMER": 1, "PRIMER_PICK_RIGHT_PRIMER": 1,
    "PRIMER_PICK_INTERNAL_OLIGO": 0,
}


def product_tiers(locus):
    if locus.get("expansion_row"):
        # short products to leave headroom for patient-length expansions
        return [[[250, 400]], [[200, 550]], [[200, 900]], [[200, 1600]], [[200, 2000]]]
    return [[[250, 900]], [[200, 1000]], [[200, 1600]], [[200, 2000]]]


def gc_fallback(ga):
    ga = dict(ga)
    ga.update({"PRIMER_MAX_GC": 75.0, "PRIMER_MIN_GC": 30.0,
               "PRIMER_MIN_TM": 56.0, "PRIMER_MAX_POLY_X": 6})
    return ga


def design(locus):
    chrom = locus["chrom"]
    seq_all = load_chrom(chrom)
    start, end = locus["start"], locus["end"]

    def attempt(tiers, ga_over=None):
        for flank in (800, 500, 350, 250):
            t_start = max(0, start - flank)
            t_end = min(len(seq_all), end + flank)
            tmpl = seq_all[t_start:t_end]
            if len(tmpl) < (end - start) + 200:
                continue
            off = start - t_start
            tract_offset = max(0, off - 3)
            tract_len = min(len(tmpl), off + (end - start) + 3) - tract_offset
            # widen exclusion over adjacent period-2 simple repeats (CA runs etc.)
            wo, wl = widen_exclusion(tmpl, tract_offset, tract_offset + tract_len)
            tract_len = wl - wo
            tract_offset = wo
            seq_args = {
                "SEQUENCE_ID": locus.get("id", f"{chrom}_{start}"),
                "SEQUENCE_TEMPLATE": tmpl,
                # amplicon MUST contain the repeat tract...
                "SEQUENCE_TARGET": [(tract_offset, tract_len)],
                # ...and primers must not bind inside it (incl. adjacent simple repeats)
                "SEQUENCE_EXCLUDED_REGION": [(tract_offset, tract_len)],
            }
            ga = dict(GLOBAL_ARGS)
            if ga_over:
                ga.update(ga_over)
            all_pairs = []
            for rng in tiers:
                ga["PRIMER_PRODUCT_SIZE_RANGE"] = rng
                res = primer3.bindings.design_primers(seq_args, ga)
                if res.get("PRIMER_PAIR_NUM_RETURNED", 0) > 0:
                    for i in range(res["PRIMER_PAIR_NUM_RETURNED"]):
                        lp, ll = res[f"PRIMER_LEFT_{i}"]
                        rp, rl = res[f"PRIMER_RIGHT_{i}"]
                        # primer3 convention (verified empirically):
                        # LEFT (start,len): start = 0-based leftmost base
                        # RIGHT (start,len): start = 0-based position of the 3'-most base
                        # product = rp + 1 - lp
                        prod = rp + 1 - lp
                        # standard path caps at 1 kb; the GC-rich fallback may use
                        # the user's 2 kb hard cap if 1 kb failed entirely
                        limit = 2000 if ga_over else 1000
                        if prod > limit:
                            continue
                        all_pairs.append({
                            "left_seq": res[f"PRIMER_LEFT_{i}_SEQUENCE"],
                            "right_seq": res[f"PRIMER_RIGHT_{i}_SEQUENCE"],
                            "left_pos": [lp, ll], "right_pos": [rp, rl],
                            "left_tm": res[f"PRIMER_LEFT_{i}_TM"],
                            "right_tm": res[f"PRIMER_RIGHT_{i}_TM"],
                            "left_gc": res[f"PRIMER_LEFT_{i}_GC_PERCENT"],
                            "right_gc": res[f"PRIMER_RIGHT_{i}_GC_PERCENT"],
                            "penalty": res[f"PRIMER_PAIR_{i}_PENALTY"],
                            "product_size": prod,
                            "amplicon_start": t_start + lp,
                            "amplicon_end": t_start + rp + 1,
                            "product_range_used": rng[0],
                        })
            if all_pairs:
                # dedupe by primer pair, keep cheapest penalty per product
                seen = set()
                uniq = []
                for p in sorted(all_pairs, key=lambda x: x["penalty"]):
                    key = (p["left_seq"], p["right_seq"])
                    if key in seen:
                        continue
                    seen.add(key)
                    uniq.append(p)
                return uniq[:10]
        return None

    pairs = attempt(product_tiers(locus))
    if pairs is None:
        pairs = attempt(product_tiers(locus), ga_over={"PRIMER_MAX_GC": 75.0,
                                                       "PRIMER_MIN_GC": 30.0,
                                                       "PRIMER_MIN_TM": 56.0,
                                                       "PRIMER_MAX_POLY_X": 6})
    return pairs


def dinuc_run_penalty(s):
    """Penalize primers containing tandem dinuc units (CA/GC/TG...) in the
    3'-most 12 bases - the main slippage / non-specific-priming risk here."""
    s = s[-12:]
    pen = 0
    for i in range(len(s) - 3):
        if s[i:i + 2] == s[i + 2:i + 4]:
            pen += 1
    return pen


def widen_exclusion(seq, lo, hi):
    """Expand [lo,hi) while the flanks continue period-2 simple repeats
    (CA/TG/GC/AC/GT dinuc runs adjacent to the target tract)."""
    while lo >= 4 and seq[lo - 2:lo] == seq[lo - 4:lo - 2]:
        lo -= 2
    while hi + 4 <= len(seq) and seq[hi:hi + 2] == seq[hi + 2:hi + 4]:
        hi += 2
    return lo, hi


def pick_best(res):
    best, best_score = None, None
    for i in range(res["PRIMER_PAIR_NUM_RETURNED"]):
        lt, rt = res[f"PRIMER_LEFT_{i}_TM"], res[f"PRIMER_RIGHT_{i}_TM"]
        penalty = res[f"PRIMER_PAIR_{i}_PENALTY"]
        ms = (dinuc_run_penalty(res[f"PRIMER_LEFT_{i}_SEQUENCE"])
              + dinuc_run_penalty(res[f"PRIMER_RIGHT_{i}_SEQUENCE"]))
        score = penalty + (abs(lt - 60) + abs(rt - 60)) * 2 + ms * 8
        if best_score is None or score < best_score:
            best_score = score
            best = {
                "left_seq": res[f"PRIMER_LEFT_{i}_SEQUENCE"],
                "right_seq": res[f"PRIMER_RIGHT_{i}_SEQUENCE"],
                "left_pos": res[f"PRIMER_LEFT_{i}"],
                "right_pos": res[f"PRIMER_RIGHT_{i}"],
                "left_tm": lt, "right_tm": rt,
                "left_gc": res[f"PRIMER_LEFT_{i}_GC_PERCENT"],
                "right_gc": res[f"PRIMER_RIGHT_{i}_GC_PERCENT"],
                "penalty": penalty,
                "microsat_penalty": ms,
            }
    return best


def main():
    with open(os.path.join(BASE, "data", "panel_loci.json")) as f:
        loci = json.load(f)
    out = []
    ok = fail = 0
    for i, locus in enumerate(loci):
        try:
            pairs = design(locus)
        except Exception as e:
            pairs = None
            locus["error"] = repr(e)
        if pairs:
            ok += 1
            out.append({"locus": locus, "candidates": pairs, "status": "ok"})
            print(f"[{i + 1}/{len(loci)}] OK  {locus['chrom']}:{locus['start']}-{locus['end']} "
                  f"{locus['motif']} x{locus['n_units']}  {len(pairs)} candidate pairs", file=sys.stderr)
        else:
            fail += 1
            out.append({"locus": locus, "candidates": [], "status": "failed"})
            print(f"[{i + 1}/{len(loci)}] FAIL {locus['chrom']}:{locus['start']} "
                  f"{locus.get('error', 'no pair')}", file=sys.stderr)
    with open(os.path.join(BASE, "results", "candidates.json"), "w") as f:
        json.dump(out, f, indent=1)
    print(f"\ndone: {ok} ok, {fail} failed -> results/candidates.json", file=sys.stderr)


if __name__ == "__main__":
    main()
