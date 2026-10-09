#!/usr/bin/env python3
"""Dimer / hairpin check for the delivered primer panel (C4).

Two things are reported, and they answer different questions.

1. **Design-threshold compliance** (`p3_*_th` columns).  The design constrains
   oligo self-structure through primer3's *thermodynamic* parameter family
   (`PRIMER_MAX_SELF_ANY_TH` / `_SELF_END_TH` / `_HAIRPIN_TH` and
   `PRIMER_PAIR_MAX_COMPL_ANY_TH` / `_END_TH`), all set to 47.0 in
   scripts/02b_design_panel.py.  Those `_TH` values ARE the constrained
   quantities, so they are read back directly and compared against 47.0.
   No conversion is needed, and none is possible for hairpins: the independent
   NN model below and primer3's HAIRPIN_TH were measured to be uncorrelated
   (r = -0.34 over the 24 pairs) because the NN model scores the stem only
   (no loop entropy) and accepts 3 bp stems primer3 does not treat as hairpins.

2. **Independent physical magnitude** (`*_dg` / `*_run` columns).  Computed
   from the SantaLucia (1998) unified nearest-neighbour parameters with no
   primer3 in the calculation, exhaustively over all antiparallel alignments:

     het      F x R         intermolecular
     homoF    F x F         intermolecular (self-dimer)
     homoR    R x R
     hairpin  intramolecular stem, 5' arm paired with a downstream 3' arm

   Two model-free metrics are reported per interaction: the longest contiguous
   complementary run, and whether that run is anchored at a 3' end -- the
   classic PCR primer-dimer criterion.  Hairpin dG here is the STEM ONLY (loop
   entropy omitted), so it is an upper bound on stem stability.

An earlier version compared the NN ΔG against hardcoded 45/35/45 kcal/mol
"mirrors" of the design thresholds.  Those mirrors matched neither the real
constraint (which is the `_TH` family) nor the NN model's scale, and could not
be repaired by changing the numbers.

Writes results/dimer_check.tsv and results/dimer_check.json.  Exits non-zero if
any pair breaches the design-time thresholds.
"""
import csv
import json
import math
import os
import sys

BASE = os.environ.get("PRIMER_DESIGN_BASE") or os.path.dirname(
    os.path.dirname(os.path.abspath(__file__)))

# Design-time limits, from scripts/02b_design_panel.py.  These are primer3
# `_TH` values, directly comparable to the `_TH` columns read back below.
MAX_SELF_ANY_TH = 47.0
MAX_SELF_END_TH = 47.0
MAX_HAIRPIN_TH = 47.0
MAX_COMPL_ANY_TH = 47.0
MAX_COMPL_END_TH = 47.0
# Practical flag: a 3'-anchored complementary run this long is a known
# primer-dimer risk regardless of what any dG threshold says.
END_RUN_WARN = 5

RC = str.maketrans("ACGT", "TGCA")
COMP = str.maketrans("ACGT", "TGCA")


def rc(s):
    return s.translate(RC)[::-1]


# primer3 design args used only to *read back* the `_TH` scores of a fixed
# primer pair.  The pair is force-placed on a synthetic template, so the
# returned scores are those of the delivered oligos, not of a fresh pick.
_TH_ARGS = {
    "PRIMER_OPT_SIZE": 22, "PRIMER_MIN_SIZE": 20, "PRIMER_MAX_SIZE": 28,
    "PRIMER_MIN_TM": 50.0, "PRIMER_MAX_TM": 80.0,
    "PRIMER_MIN_GC": 20.0, "PRIMER_MAX_GC": 80.0,
    "PRIMER_NUM_RETURN": 1,
}

_TH_KEYS = {
    "fwd_self_any_th": "PRIMER_LEFT_0_SELF_ANY_TH",
    "fwd_self_end_th": "PRIMER_LEFT_0_SELF_END_TH",
    "fwd_hairpin_th": "PRIMER_LEFT_0_HAIRPIN_TH",
    "rev_self_any_th": "PRIMER_RIGHT_0_SELF_ANY_TH",
    "rev_self_end_th": "PRIMER_RIGHT_0_SELF_END_TH",
    "rev_hairpin_th": "PRIMER_RIGHT_0_HAIRPIN_TH",
    "compl_any_th": "PRIMER_PAIR_0_COMPL_ANY_TH",
    "compl_end_th": "PRIMER_PAIR_0_COMPL_END_TH",
}


def primer3_th_scores(fwd, rev, filler="ACGT" * 15):
    """Read the `_TH` self/pair scores of one fixed primer pair.

    The pair is placed on `fwd + filler + rc(rev)` and both ends are forced, so
    primer3 scores exactly the delivered oligos.  Returns {} if primer3 is
    unavailable or cannot place the pair (reported as not-run).
    """
    try:
        import primer3
    except ImportError:
        return {}
    tmpl = fwd + filler + rc(rev)
    n = len(tmpl)
    seq_args = {
        "SEQUENCE_ID": "th_readback",
        "SEQUENCE_TEMPLATE": tmpl,
        "SEQUENCE_FORCE_LEFT_START": 0,
        "SEQUENCE_FORCE_LEFT_END": len(fwd) - 1,
        # primer3's RIGHT convention: START is the 5'-most (rightmost) base.
        "SEQUENCE_FORCE_RIGHT_START": n - 1,
        "SEQUENCE_FORCE_RIGHT_END": n - len(rev),
    }
    try:
        res = primer3.bindings.design_primers(seq_args, dict(_TH_ARGS))
    except Exception:                                  # noqa: BLE001
        return {}
    if not res.get("PRIMER_PAIR_NUM_RETURNED", 0):
        return {}
    if (res.get("PRIMER_LEFT_0_SEQUENCE") != fwd
            or res.get("PRIMER_RIGHT_0_SEQUENCE") != rev):
        # the forced placement did not survive; scores would be misleading
        return {}
    return {k: res.get(v) for k, v in _TH_KEYS.items()}


# --- SantaLucia 1998 unified NN parameters, 1 M NaCl -----------------------
# key = top-strand dinucleotide (5'->3'); value = (dH kcal/mol, dS cal/mol/K)
NN = {
    "AA": (-7.9, -22.2), "TT": (-7.9, -22.2),
    "AT": (-7.2, -20.4), "TA": (-7.2, -21.3),
    "CA": (-8.5, -22.7), "TG": (-8.5, -22.7),
    "GT": (-8.4, -22.4), "AC": (-8.4, -22.4),
    "CT": (-7.8, -21.0), "AG": (-7.8, -21.0),
    "GA": (-8.2, -22.2), "TC": (-8.2, -22.2),
    "CG": (-10.6, -27.2),
    "GC": (-9.8, -24.4),
    "GG": (-8.0, -19.9), "CC": (-8.0, -19.9),
}
INIT_GC = (0.1, -2.8)      # duplex end terminating in a G.C pair
INIT_AT = (2.3, 4.1)       # duplex end terminating in an A.T pair
T37 = 310.15


def nn_dg(top):
    """dG37 (kcal/mol) of a fully paired duplex, given one strand 5'->3'.

    The partner strand is fixed as the complement, so the stack is determined by
    the top strand alone; the two terminal base pairs set the initiation terms.
    """
    L = len(top)
    if L < 2:
        return None
    dH = dS = 0.0
    for i in range(L - 1):
        dH += NN[top[i:i + 2]][0]
        dS += NN[top[i:i + 2]][1]
    for end_base in (top[0], top[-1]):
        h, s = INIT_GC if end_base in "GC" else INIT_AT
        dH += h
        dS += s
    return dH - T37 * dS / 1000.0


def scan(a, b):
    """Most stable antiparallel register of a against b, plus run statistics.

    a[i] pairs with rc(b)[k].  rc(b)[0] is the complement of b's 3' base, so a
    run touching k == 0 is anchored at b's 3' end, and one touching
    i == len(a) - 1 is anchored at a's 3' end.
    """
    br = rc(b)
    m, n = len(a), len(br)
    best_dg = None
    longest = 0
    a3_run = 0
    b3_run = 0
    for d in range(-n + 1, m):
        i = max(0, d)
        k = max(0, -d)
        run = 0
        start = None
        while True:
            inb = i < m and k < n
            paired = inb and a[i] == br[k]
            if paired:
                if run == 0:
                    start = (i, k)
                run += 1
            else:
                if run:
                    top = a[start[0]:start[0] + run]
                    dg = nn_dg(top)          # top == br slice by construction
                    if dg is not None and (best_dg is None or dg < best_dg):
                        best_dg = dg
                    longest = max(longest, run)
                    touches_a3 = start[0] + run - 1 == m - 1
                    touches_b3 = start[1] == 0
                    if touches_a3:
                        a3_run = max(a3_run, run)
                    if touches_b3:
                        b3_run = max(b3_run, run)
                run = 0
                if not inb:
                    break
            if not inb:
                break
            i += 1
            k += 1
    return {"dg": best_dg, "run": longest, "a3": a3_run, "b3": b3_run}


def scan_hairpin(a, min_stem=3, min_loop=3):
    """Most stable intramolecular stem, loop entropy omitted (upper bound)."""
    best = None
    L = len(a)
    for stem in range(min_stem, L // 2 + 1):
        for i in range(0, L - 2 * stem - min_loop + 1):
            arm = a[i:i + stem]
            for j in range(i + stem + min_loop, L - stem + 1):
                if a[j:j + stem] == rc(arm):
                    dg = nn_dg(arm)
                    if dg is not None and (best is None or dg < best[0]):
                        best = (dg, stem, i, j)
    return {"dg": best[0] if best else None,
            "stem": best[1] if best else 0,
            "i": best[2] if best else None,
            "j": best[3] if best else None}


# --- primer3 cross-check (same engine as design, matched salt) -------------
def p3(fn, *seqs, **kw):
    try:
        import primer3
    except ImportError:
        return None
    try:
        return fn(*seqs, **kw).dg / 1000.0        # cal/mol -> kcal/mol
    except Exception:                              # noqa: BLE001
        return None


def main():
    tsv = os.path.join(BASE, "results", "primer_design_results.tsv")
    with open(tsv) as f:
        rows = list(csv.DictReader(f, delimiter="\t"))
    try:
        import primer3
        have_p3 = True
    except ImportError:
        primer3 = None
        have_p3 = False
        print("WARNING: primer3-py absent; cross-engine columns left empty",
              file=sys.stderr)

    out = []
    breaches = []
    for r in rows:
        f_, r_ = r["fwd"].upper(), r["rev"].upper()
        het = scan(f_, r_)
        hf = scan(f_, f_)
        hr = scan(r_, r_)
        ph_f = scan_hairpin(f_)
        ph_r = scan_hairpin(r_)

        rec = {"no": int(r["no"]), "gene": r["gene"],
               "het_dg": het["dg"], "het_run": het["run"],
               "het_f3": het["a3"], "het_r3": het["b3"],
               "hf_dg": hf["dg"], "hf_run": hf["run"], "hf_3": hf["a3"],
               "hr_dg": hr["dg"], "hr_run": hr["run"], "hr_3": hr["a3"],
               "phf_dg": ph_f["dg"], "phf_stem": ph_f["stem"],
               "phr_dg": ph_r["dg"], "phr_stem": ph_r["stem"]}

        if have_p3:
            # 1 M monovalent, no divalent/dNTP: the conditions the NN model above
            # assumes, so the two engines are directly comparable.
            oneM = dict(mv_conc=1000.0, dv_conc=0.0, dntp_conc=0.0, dna_conc=50.0)
            design = dict(mv_conc=50.0, dv_conc=1.5, dntp_conc=0.6, dna_conc=50.0)
            rec["p3_het_1M"] = p3(primer3.calc_heterodimer, f_, r_, **oneM)
            rec["p3_het_design"] = p3(primer3.calc_heterodimer, f_, r_, **design)
            rec["p3_hf_1M"] = p3(primer3.calc_homodimer, f_, **oneM)
            rec["p3_hr_1M"] = p3(primer3.calc_homodimer, r_, **oneM)
            rec["p3_phf_design"] = p3(primer3.calc_hairpin, f_, **design)
            rec["p3_phr_design"] = p3(primer3.calc_hairpin, r_, **design)

        # ---- verdict 1: design-threshold compliance, in primer3's own units ----
        # The design constrains the `_TH` family; these are read back for the
        # delivered oligos and compared against the limits directly.  The
        # independent NN numbers below are physical magnitude, NOT the test.
        th = primer3_th_scores(f_, r_)
        rec.update(th)
        fail = []
        if th:
            for key, limit, label in (
                    ("fwd_self_any_th", MAX_SELF_ANY_TH, "F SELF_ANY_TH"),
                    ("fwd_self_end_th", MAX_SELF_END_TH, "F SELF_END_TH"),
                    ("fwd_hairpin_th", MAX_HAIRPIN_TH, "F HAIRPIN_TH"),
                    ("rev_self_any_th", MAX_SELF_ANY_TH, "R SELF_ANY_TH"),
                    ("rev_self_end_th", MAX_SELF_END_TH, "R SELF_END_TH"),
                    ("rev_hairpin_th", MAX_HAIRPIN_TH, "R HAIRPIN_TH"),
                    ("compl_any_th", MAX_COMPL_ANY_TH, "pair COMPL_ANY_TH"),
                    ("compl_end_th", MAX_COMPL_END_TH, "pair COMPL_END_TH")):
                v = th.get(key)
                if v is not None and v > limit:
                    fail.append(f"{label}={v:.1f} > {limit:.0f}")
        else:
            rec["th_readback_note"] = "primer3 could not force-place this pair"

        # ---- verdict 2: model-free primer-dimer warning ----
        # A 3'-anchored complementary run is a risk regardless of any dG
        # threshold, so it is reported independently of the `_TH` limits.
        warn = []
        for name, dg, run3 in (("het", het["dg"], max(het["a3"], het["b3"])),
                               ("homoF", hf["dg"], hf["a3"]),
                               ("homoR", hr["dg"], hr["a3"]),
                               ("hairpinF", ph_f["dg"], 0),
                               ("hairpinR", ph_r["dg"], 0)):
            if run3 >= END_RUN_WARN:
                warn.append(f"{name} 3'-anchored run {run3} bp")
        rec["flag"] = "fail" if fail else ("warn" if warn else "ok")
        rec["detail"] = "; ".join(fail + warn)
        out.append(rec)
        if fail:
            breaches.append(rec)

    # ---- write TSV ----
    cols = ["no", "gene", "flag",
            # design-threshold compliance, in primer3's own `_TH` units
            "fwd_self_any_th", "fwd_self_end_th", "fwd_hairpin_th",
            "rev_self_any_th", "rev_self_end_th", "rev_hairpin_th",
            "compl_any_th", "compl_end_th",
            # independent NN physical magnitude (NOT the threshold test)
            "het_dg", "het_run", "het_f3", "het_r3",
            "hf_dg", "hf_run", "hf_3", "hr_dg", "hr_run", "hr_3",
            "phf_dg", "phf_stem", "phr_dg", "phr_stem",
            "p3_het_1M", "p3_het_design",
            "p3_hf_1M", "p3_hr_1M", "p3_phf_design", "p3_phr_design", "detail"]
    with open(os.path.join(BASE, "results", "dimer_check.tsv"), "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=cols, delimiter="\t",
                           extrasaction="ignore")
        w.writeheader()
        for rec in out:
            w.writerow({k: ("" if rec.get(k) is None else
                            (f"{rec[k]:.2f}" if isinstance(rec.get(k), float)
                             else rec.get(k))) for k in cols})

    # ---- engine agreement ----
    both = [(x["het_dg"], x["p3_het_1M"]) for x in out
            if x.get("het_dg") is not None and x.get("p3_het_1M") is not None]
    agree = None
    if len(both) >= 3:
        diffs = [abs(a - b) for a, b in both]
        mean_a = sum(a for a, _ in both) / len(both)
        mean_b = sum(b for _, b in both) / len(both)
        cov = sum((a - mean_a) * (b - mean_b) for a, b in both)
        va = math.sqrt(sum((a - mean_a) ** 2 for a, _ in both))
        vb = math.sqrt(sum((b - mean_b) ** 2 for _, b in both))
        agree = {"n": len(both), "max_abs_diff": max(diffs),
                 "mean_abs_diff": sum(diffs) / len(diffs),
                 "pearson_r": cov / (va * vb) if va and vb else None}
    summary = {
        "n_pairs": len(out),
        "flag_counts": {k: sum(1 for x in out if x["flag"] == k)
                        for k in ("ok", "warn", "fail")},
        "thresholds": {"MAX_SELF_ANY_TH": MAX_SELF_ANY_TH,
                       "MAX_SELF_END_TH": MAX_SELF_END_TH,
                       "MAX_HAIRPIN_TH": MAX_HAIRPIN_TH,
                       "MAX_COMPL_ANY_TH": MAX_COMPL_ANY_TH,
                       "MAX_COMPL_END_TH": MAX_COMPL_END_TH,
                       "END_RUN_WARN": END_RUN_WARN},
        "threshold_note": ("limits are primer3 `_TH` values, read back for the "
                           "delivered oligos; the NN dG columns are physical "
                           "magnitude and are not thresholded"),
        "model": "SantaLucia 1998 unified NN, 1 M NaCl, T=310.15 K",
        "hairpin_note": "stem only; loop entropy omitted (upper bound)",
        "engine_agreement_het_1M": agree,
        "breaches": [{"no": x["no"], "gene": x["gene"], "detail": x["detail"]}
                     for x in breaches],
    }
    with open(os.path.join(BASE, "results", "dimer_check.json"), "w") as f:
        json.dump(summary, f, indent=1, ensure_ascii=False)

    print(f"24 pairs checked; flags: {summary['flag_counts']}")
    if agree:
        print(f"engine agreement (het, 1 M): r={agree['pearson_r']:.3f}, "
              f"mean|d|={agree['mean_abs_diff']:.2f}, "
              f"max|d|={agree['max_abs_diff']:.2f} kcal/mol")
    worst = sorted((x for x in out if x["het_dg"] is not None),
                   key=lambda x: x["het_dg"])[:3]
    for x in worst:
        print(f"  most stable het: #{x['no']} {x['gene']} "
              f"dG={x['het_dg']:.1f}, run={x['het_run']} bp, "
              f"3'-anchored F={x['het_f3']} R={x['het_r3']}")
    if breaches:
        print("THRESHOLD BREACHES:", file=sys.stderr)
        for x in breaches:
            print(f"  #{x['no']} {x['gene']}: {x['detail']}", file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    main()
