#!/usr/bin/env python3
"""Validation for the mismatch-tolerant pair scan (scripts/18_mismatch_pcr.py).

A "0 off-target products" result is only worth anything if the scanner can be
shown to detect one when it exists.  So this validation is built around
synthetic controls with known answers, and only then applies the scanner to the
real genome.

V1  synthetic positive controls: an engineered off-target pair is found, and the
    mismatch budget behaves as specified (mm2 found everywhere, mm6 only in the
    permissive config)
V2  synthetic negative controls: a pair 5 kb apart, and a wrong-orientation
    pair, are NOT reported
V3  both pairing configurations work (plus-strand and minus-strand products)
V4  on-target recovery against the deliverable: for every pair whose contig is
    present, exactly one perfect product at the delivered coordinates
V5  cross-implementation agreement with method_D: the seed-15 site count can
    never exceed method_D's 3'-15mer hit count, since it is the same seed plus a
    mismatch constraint; and the zero-mismatch site count per primer must be 1,
    independently reproducing full-length genome uniqueness
V6  dimer check reproducibility (scripts/17_dimer_check.py) if its output exists

Genome-dependent checks (V4, V5) are SKIPPED, not failed, when the reference is
absent, and are reported as such.

Usage: python3 scripts/19_validate_mismatch.py
"""
import importlib.util
import json
import os
import random
import sys

BASE = os.environ.get("PRIMER_DESIGN_BASE") or os.path.dirname(
    os.path.dirname(os.path.abspath(__file__)))
RC = str.maketrans("ACGT", "TGCA")

checks = []


def check(name, ok, detail="", skipped=False):
    checks.append((name, ok, detail, skipped))
    tag = "SKIP" if skipped else ("PASS" if ok else "FAIL")
    print(f"  [{tag}] {name}" + (f" — {detail}" if detail else ""))
    return ok


def load_module(fname):
    spec = importlib.util.spec_from_file_location(
        fname[:-3], os.path.join(BASE, "scripts", fname))
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def rc(s):
    return s.translate(RC)[::-1]


def rnd_seq(n, seed=99):
    r = random.Random(seed)
    return "".join(r.choice("ACGT") for _ in range(n))


def find_genome(require_chrX=False):
    """First usable oneline genome: the UCSC one, else the local GRCh38 copy."""
    for cand in ("hg38_oneline", "hg38_local_oneline"):
        d = os.path.join(BASE, "data", cand)
        if not os.path.exists(os.path.join(d, "chr1.fa")):
            continue
        if require_chrX and not os.path.exists(os.path.join(d, "chrX.fa")):
            continue
        return d
    return None


def read_contig(genome_dir, name):
    p = os.path.join(genome_dir, f"{name}.fa")
    if not os.path.exists(p):
        return None
    with open(p) as f:
        f.readline()
        return f.read().strip().upper()


def build_synthetic(primers):
    """Deterministic contigs with known off-target answers, for pair #3."""
    random.seed(1234)
    pr = next(p for p in primers if p["no"] == 3)
    F, R = pr["fwd"], pr["rev"]
    rnd = lambda n: "".join(random.choice("ACGT") for _ in range(n))  # noqa: E731

    def mut(s, idxs):
        l = list(s)
        for i in idxs:
            l[i] = "A" if l[i] != "A" else "C"
        return "".join(l)

    return {
        # detected, 0 mismatches
        "synth_exact": rnd(500) + F + rnd(300) + rc(R) + rnd(500),
        # detected, 2 mismatches in the 5' half (3' 15-mer intact)
        "synth_mm2": rnd(500) + mut(F, [0, 3]) + rnd(300) + rc(R) + rnd(500),
        # only within the permissive budget: 6 mismatches, 3' 15-mer still intact
        "synth_mm6": rnd(500) + mut(F, [0, 1, 2, 3, 4, 5]) + rnd(300) + rc(R) + rnd(500),
        # too far apart -> outside the product window
        "synth_far": rnd(500) + F + rnd(5000) + rc(R) + rnd(500),
        # wrong orientation -> same-direction primers, no product
        "synth_orient": rnd(500) + F + rnd(300) + R + rnd(500),
        # minus-strand product: R left, rc(F) right
        "synth_minus": rnd(500) + R + rnd(300) + rc(F) + rnd(500),
    }, pr


def scan(m18, primers, seqs, seed, mm):
    hits = {(i, s, o): [] for i in range(len(primers))
            for s in ("fwd", "rev") for o in ("+", "-")}
    qarr, entries = m18.build_seed_maps(primers, seed)
    for name, seq in seqs.items():
        m18.scan_contig(seq, primers, seed, mm, hits, name, qarr, entries)
    rows = []
    m18.pair_products(primers, hits, f"seed{seed}_mm{mm}", rows)
    return rows, sum(len(v) for v in hits.values())


def main():
    m18 = load_module("18_mismatch_pcr.py")
    primers = m18.load_primers()
    seqs, _ = build_synthetic(primers)

    print("V1/V2/V3 synthetic controls (pair #3, engineered products)")
    strict, _ = scan(m18, primers, seqs, 15, 3)
    mid, _ = scan(m18, primers, seqs, 12, 4)
    perm, _ = scan(m18, primers, seqs, 10, 6)
    got = {(r["contig"], r["config"]): r for r in strict + mid + perm}

    def found(contig, cfg):
        return got.get((contig, cfg))

    check("V1a synth_exact 在全部三档被检出（0 错配）",
          all(found("synth_exact", c) and
              found("synth_exact", c)["fwd_mm"] == 0 and
              found("synth_exact", c)["rev_mm"] == 0
              for c in ("seed15_mm3", "seed12_mm4", "seed10_mm6")))
    r = found("synth_mm2", "seed15_mm3")
    check("V1b synth_mm2 在严格档被检出且 fwd_mm=2，3' 锚定 18 bp",
          bool(r and r["fwd_mm"] == 2 and r["fwd_3run"] == 18),
          f"fwd_mm={r['fwd_mm']}, 3run={r['fwd_3run']}" if r else "未检出")
    check("V1c synth_mm6 在严格档被拒绝、在宽松档被检出（错配预算生效）",
          found("synth_mm6", "seed15_mm3") is None
          and found("synth_mm6", "seed10_mm6") is not None
          and found("synth_mm6", "seed10_mm6")["fwd_mm"] == 6)
    check("V2a synth_far（相距 5 kb）未被报告为产物",
          found("synth_far", "seed10_mm6") is None)
    check("V2b synth_orient（方向错误）未被报告为产物",
          found("synth_orient", "seed10_mm6") is None)
    rm = found("synth_minus", "seed15_mm3")
    check("V3 负链产物通路（strand='-'）可检出",
          bool(rm and rm["strand"] == "-" and rm["product_len"] == 344),
          f"strand={rm['strand']}, len={rm['product_len']}" if rm else "未检出")

    # ---- genome-dependent checks ----
    gdir = find_genome()
    genome_ok = gdir is not None
    print(f"V4/V5 reference-dependent checks (genome: {gdir})")
    if not genome_ok:
        check("V4 on-target recovery", False,
              "无可用参考基因组，无法执行", skipped=True)
        check("V5 method_D 交叉校验", False, "同上", skipped=True)
    else:
        mm_path = os.path.join(BASE, "results", "mismatch_pcr.json")
        md_path = os.path.join(BASE, "results", "method_D_genomescan.json")
        if not (os.path.exists(mm_path) and os.path.exists(md_path)):
            check("V4 on-target recovery", False,
                  "缺少 mismatch_pcr.json 或 method_D_genomescan.json", skipped=True)
            check("V5 method_D 交叉校验", False, "同上", skipped=True)
        else:
            summary = json.load(open(mm_path))
            md = json.load(open(md_path))
            files = sorted(f for f in os.listdir(gdir) if f.endswith(".fa"))
            scanned = {f[:-3] for f in files}
            rec = summary.get("on_target_found", {})
            testable = {k: v for k, v in rec.items() if v.get("n") is not None}
            bad = {k: v for k, v in testable.items()
                   if v["n"] != 1 or not v["coord_match"]}
            check(f"V4 on-target 回收：{len(testable)} 个可测位点全部唯一且坐标吻合",
                  not bad and len(testable) > 0,
                  f"异常 {list(bad)[:5]}" if bad else
                  f"{len(testable)} 个通过（未测 {'chrX' if 'chrX' not in scanned else '—'}）")

            # V5: seed15 sites <= method_D end15 hits, and exactly one 0-mm site.
            # method_D counts every exact 3'-15mer; this scan counts the same
            # seed intersected with a mismatch budget, so it can only be fewer.
            sites = {c["config"]: c for c in summary.get("configs", [])}
            panel_15 = sum(r["fwd"]["end15_hits"] + r["rev"]["end15_hits"]
                           for r in md["rows"])
            got_panel = sites.get("seed15_mm3", {}).get("sites_total")
            check("V5a seed15 位点总数 ≤ method_D 的 3'15mer 命中总数（同种子+错配约束）",
                  got_panel is not None and got_panel <= panel_15,
                  f"{got_panel} vs {panel_15}")
            on_t = sum(v["n"] for v in testable.values())
            check("V5b 每个可测 pair 的 0 错配位点数 == 1（独立复现全长唯一性）",
                  all(v["n"] == 1 for v in testable.values()),
                  f"{on_t} 个 0 错配位点 / {len(testable)} 个可测 pair")

    # ---- V6 dimer output sanity ----
    print("V6 dimer check")
    dpath = os.path.join(BASE, "results", "dimer_check.json")
    if not os.path.exists(dpath):
        check("V6 dimer_check.json 存在", False, "未生成", skipped=True)
    else:
        d = json.load(open(dpath))
        ag = d.get("engine_agreement_het_1M") or {}
        check("V6a dimer 检查覆盖 24 对且无阈值违例",
              d["n_pairs"] == 24 and d["flag_counts"]["fail"] == 0,
              f"{d['flag_counts']}")
        check("V6b 独立 NN 实现与 primer3 在 1 M 盐下一致（r ≥ 0.9）",
              (ag.get("pearson_r") or 0) >= 0.9,
              f"r={ag.get('pearson_r'):.3f}, mean|Δ|={ag.get('mean_abs_diff'):.2f} kcal/mol"
              if ag.get("pearson_r") else "无")

    # ---- V7 spike-in against a real chromosome background ----
    print("V7 spike-in control on a real chromosome")
    gdir = find_genome()
    if not gdir:
        check("V7 spike-in 对照", False, "无可用参考基因组", skipped=True)
    else:
        chr1 = read_contig(gdir, "chr1")
        if chr1 is None:
            check("V7 spike-in 对照", False, "chr1 不可读", skipped=True)
        else:
            random.seed(99)
            pr5 = next(p for p in primers if p["no"] == 5)
            insert = (pr5["fwd"] + rnd_seq(300) + rc(pr5["rev"]))
            ins_at = 100_000_000
            spiked = chr1[:ins_at] + insert + chr1[ins_at:]
            rows, sites = scan(m18, primers, {"chr1": spiked}, 15, 3)
            off = [r for r in rows if r["kind"] == "off_target"]
            hit = [r for r in off if r["no"] == 5 and r["contig"] == "chr1"
                   and r["start"] == ins_at and r["product_len"] == len(insert)]
            on_t = [r for r in rows if r["kind"] == "on_target"]
            check(f"V7 spike-in：在 chr1 真实背景（{sites} 个位点）中唯一检出人工脱靶",
                  len(hit) == 1 and len(off) == 1,
                  f"检出 {len(off)} 个脱靶，坐标命中 {len(hit)}"
                  f"；同批 on-target {len(on_t)} 个")

    print()
    failed = [c for c in checks if not c[1] and not c[3]]
    skipped = [c for c in checks if c[3]]
    if failed:
        print(f"MISMATCH VALIDATION: FAILED ({len(failed)})")
        for n, _o, d, _s in failed:
            print(f"  - {n}: {d}")
        sys.exit(1)
    print(f"MISMATCH VALIDATION: PASSED "
          f"({len(checks) - len(skipped)} checks, {len(skipped)} skipped)")


if __name__ == "__main__":
    main()
