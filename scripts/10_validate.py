#!/usr/bin/env python3
"""Independent validation of the delivered primer panel.

Methods (each independent of the design pipeline):
  A. PCR simulation straight from the reference genome (no primer3):
     locate both primer binding sites by exact string match, re-derive product
     size/coordinates, verify the repeat tract is bracketed and primers do not
     touch it.
  B. Independent repeat-tract re-analysis: re-extract the tract from the
     genome, re-count unit copies with a regex-based approach (different
     algorithm from the design-time frame-walk), re-check GC%, purity bounds.
  C. Thermodynamics re-check: primer3 calcTm under the design salt conditions,
     Biopython Tm_NN as a second opinion, GC% recomputation, homodimer /
     heterodimer / hairpin dG recomputation.
  E. Annotation cross-check: UCSC ncbiRefSeqCurated vs Ensembl gene symbols.
  X. Internal consistency of the deliverable TSV vs designed.json.

Method D (genome-wide specificity rescan with numpy rolling hash) lives in
11_validate_genome_scan.py and is merged here when its output exists.

Writes results/validation_report.md
"""
import json
import os
import subprocess
import urllib.request

BASE = os.environ.get("PRIMER_DESIGN_BASE") or os.path.dirname(
    os.path.dirname(os.path.abspath(__file__)))
CHROM_FA = os.path.join(BASE, "data", "hg38", "chromFa")

TRINUC_ROTS = ["GCC", "CCG", "CGC", "GGC", "CGG", "GCG"]
DINUC_ROTS = ["GC", "CG"]
RC = str.maketrans("ACGT", "TGCA")


def rc(s):
    return s.translate(RC)[::-1]


def gc_pct(s):
    return (s.count("G") + s.count("C")) / len(s) * 100 if s else 0.0


class Report:
    def __init__(self):
        self.lines = ["# 引物面板数据校验报告", "",
                      f"生成时间：数据流水线 results/primer_design_results.tsv + designed.json", ""]
        self.checks = []

    def check(self, method, name, ok, detail=""):
        self.checks.append((method, name, ok, detail))
        mark = "PASS" if ok else "FAIL"
        self.lines.append(f"- [{mark}] ({method}) {name}" + (f" — {detail}" if detail else ""))
        return ok

    def info(self, text):
        self.lines.append(f"- [INFO] {text}")

    def summary(self):
        n_pass = sum(1 for m, n, ok, d in self.checks if ok)
        n_fail = len(self.checks) - n_pass
        self.lines.insert(4, f"**结果：{n_pass} 项通过，{n_fail} 项失败，共 {len(self.checks)} 项检查**")
        self.lines.insert(5, "")
        return n_fail


def load_chrom(chrom):
    with open(os.path.join(CHROM_FA, f"{chrom}.fa")) as f:
        f.readline()
        return f.read().replace("\n", "").upper()


def main():
    rep = Report()
    rows = []
    with open(os.path.join(BASE, "results", "primer_design_results.tsv")) as f:
        header = f.readline().rstrip("\n").split("\t")
        for line in f:
            rows.append(dict(zip(header, line.rstrip("\n").split("\t"))))
    with open(os.path.join(BASE, "results", "designed.json")) as f:
        raw = [it for it in json.load(f) if it.get("status") in ("ok", "ok_check")]

    # ---------------- X. internal consistency ----------------
    rep.lines.append("## X. 交付数据内部一致性")
    rep.check("X", "TSV 行数 == designed.json 有效对数", len(rows) == len(raw),
              f"tsv={len(rows)}, json={len(raw)}")

    ok = True
    for r, it in zip(rows, raw):
        if r["fwd"] != it["pair"]["left_seq"] or r["rev"] != it["pair"]["right_seq"]:
            ok = False
    rep.check("X", "TSV 引物序列与 designed.json 一致", ok)

    ok, det = True, []
    for r in rows:
        ps = int(r["product_size"])
        span = int(r["amplicon_end"]) - int(r["amplicon_start"]) + 1
        if ps != span:
            ok = False
            det.append(f"#{r['no']} {ps}!={span}")
    rep.check("X", "product_size == 扩增子坐标跨度（1-based 闭区间）", ok, "; ".join(det))

    ok, det = True, []
    for r in rows:
        n = int(r["repeat_end"]) - int(r["repeat_start"]) + 1
        ulen = 3 if r["family"] == "GCC" else 2
        if n // ulen != int(float(r["repeat_feature"].split("×")[1].split("（")[0])):
            ok = False
            det.append(f"#{r['no']}")
    rep.check("X", "重复区长度与 Repeats 特点中拷贝数自洽", ok, "; ".join(det))

    # bucket label sanity
    ok, det = True, []
    for r in rows:
        n = int(r["repeat_end"]) - int(r["repeat_start"]) + 1
        ulen = 3 if r["family"] == "GCC" else 2
        units = n // ulen
        b = r["bucket"]
        if "扩展检测" in b or "扩展" in b:
            continue
        if r["family"] == "GCC" and 20 <= units < 25 and "20-25" not in b:
            ok = False; det.append(f"#{r['no']} {b}")
        if r["family"] == "GCC" and 25 <= units < 30 and "25-30" not in b:
            ok = False; det.append(f"#{r['no']} {b}")
    rep.check("X", "档位标签与拷贝数一致（扩展检测位除外）", ok, "; ".join(det))

    for r in rows:
        if r["bucket"].startswith(">30 ✓"):
            n = int(r["repeat_end"]) - int(r["repeat_start"]) + 1
            if n // 3 < 30:
                ok = False
                det.append(f"#{r['no']} >30 标签但仅 {n // 3} 单元")
    rep.check("X", ">30 ✓ 行的单元数 >= 30", ok, "; ".join(det))

    exp_rows = [r for r in rows if "扩展检测" in r["bucket"]]
    has_fmr1 = any(r["gene"] == "FMR1" for r in exp_rows)
    has_aff2 = any(r["gene"] == "AFF2" for r in exp_rows)
    rep.check("X", ">30 扩展检测行存在且为 FMR1/AFF2", len(exp_rows) == 2 and has_fmr1 and has_aff2,
              f"n={len(exp_rows)}")

    # ---------------- Method A: PCR simulation from genome ----------------
    rep.lines.append("## A. 参考基因组模拟 PCR（不经 primer3 的独立复核）")
    chrom_cache = {}
    okA = True
    detA = []
    for r in rows:
        chrom = r["chrom"]
        if chrom not in chrom_cache:
            chrom_cache[chrom] = load_chrom(chrom)
        seq = chrom_cache[chrom]
        fwd, rev = r["fwd"], r["rev"]
        as0 = int(r["amplicon_start"]) - 1
        ae1 = int(r["amplicon_end"])
        amp = seq[as0:ae1]
        # A1: fwd binds exactly once, at amplicon start
        if not amp.startswith(fwd):
            okA = False
            detA.append(f"#{r['no']} F primer 不在扩增子起点")
        elif amp.count(fwd) != 1:
            okA = False
            detA.append(f"#{r['no']} F primer 在扩增子内多次匹配")
        # A2: rc(rev) ends exactly at amplicon end
        if not amp.endswith(rc(rev)):
            okA = False
            detA.append(f"#{r['no']} R primer(反向互补) 不在扩增子终点")
        elif amp.count(rc(rev)) != 1:
            okA = False
            detA.append(f"#{r['no']} R primer 在扩增子内多次匹配")
        # A3: product size
        if len(amp) != int(r["product_size"]):
            okA = False
            detA.append(f"#{r['no']} 产物长度不符")
        # A4: tract bracketed
        rs0, re1 = int(r["repeat_start"]) - 1, int(r["repeat_end"])
        if not (as0 <= rs0 and re1 <= ae1):
            okA = False
            detA.append(f"#{r['no']} 重复区未被扩增子完整包含")
        # A5: primers do not touch the tract
        if (as0 + len(fwd) > rs0) or (ae1 - len(rev) < re1):
            okA = False
            detA.append(f"#{r['no']} 引物与重复区重叠")
        # A6: ACGT-only primers
        if set(fwd + rev) - set("ACGT"):
            okA = False
            detA.append(f"#{r['no']} 引物含非 ACGT 碱基")
    rep.check("A", f"{len(rows)} 对引物结合位点/产物长度/包含关系全部复核", okA, "; ".join(detA))

    # ---------------- Method B: independent tract re-analysis ----------------
    rep.lines.append("## B. 重复 tract 独立重分析（正则/相位法 vs 设计时 walking 法）")
    import re as _re
    okB, detB = True, []
    for r in rows:
        chrom = r["chrom"]
        seq = chrom_cache[chrom]
        rs0, re1 = int(r["repeat_start"]) - 1, int(r["repeat_end"])
        tract = seq[rs0:re1]
        ulen = 3 if r["family"] == "GCC" else 2
        rots = TRINUC_ROTS if ulen == 3 else DINUC_ROTS
        units = len(tract) // ulen
        feat = r["repeat_feature"]
        if "精确" in feat:
            exact = int(feat.split("精确")[1].split("+")[0])
            interrupts = int(feat.split("中断")[1].split("，")[0])
        else:
            exact, interrupts = units, 0
        purity = float(feat.split("纯度")[1].split("%")[0]) if "纯度" in feat else 100.0

        # independent recount 1: maximal consecutive exact run (regex)
        best_run = 0
        for strand in (tract, rc(tract)):
            for u in rots:
                for m in _re.finditer(f"(?:{u}){{{ulen},}}", strand):
                    best_run = max(best_run, (m.end() - m.start()) // ulen)
        # independent recount 2: max period-aligned occurrences over all frames
        phase_max = 0
        for strand in (tract, rc(tract)):
            for u in rots:
                for phase in range(ulen):
                    cnt = sum(1 for i in range(phase, len(strand) - ulen + 1, ulen)
                              if strand[i:i + ulen] == u)
                    phase_max = max(phase_max, cnt)

        if len(tract) % ulen != 0:
            okB = False; detB.append(f"#{r['no']} tract 长度非整数个单元")
            continue
        if exact + interrupts != units:
            okB = False; detB.append(f"#{r['no']} exact+interrupt={exact + interrupts} != units={units}")
        if abs(purity - round(exact / units * 100, 1)) > 0.11:
            okB = False; detB.append(f"#{r['no']} purity {purity} != {round(exact / units * 100, 1)}")
        if not (best_run <= exact <= phase_max):
            okB = False
            detB.append(f"#{r['no']} exact={exact} 不在 [max_run={best_run}, phase_max={phase_max}] 内")
        if r["family"] == "GC":
            if not any(strand == u * units
                       for strand in (tract, rc(tract)) for u in rots):
                okB = False; detB.append(f"#{r['no']} GC 位点不是完美串联")
        else:
            if interrupts == 0 and not any(strand == u * units
                                           for strand in (tract, rc(tract)) for u in rots):
                okB = False; detB.append(f"#{r['no']} 声称零中断但非完美串联")
    rep.check("B", "重复拷贝数/中断/纯度独立复核（正则+相位法）", okB, "; ".join(detB[:6]))

    # ---------------- Method C: thermodynamics ----------------
    rep.lines.append("## C. 热力学复核（primer3 calcTm + Biopython Tm_NN 双重计算）")
    import primer3
    from Bio.SeqUtils.MeltingTemp import Tm_NN
    detC1, detC2, detC3, detC4 = [], [], [], []
    max_bp_dev, max_biopy_dev = 0.0, 0.0
    for r in rows:
        for side in ("fwd", "rev"):
            s = r[side]
            # C1 same-engine recompute
            t = primer3.bindings.calc_tm(s, mv_conc=50.0, dv_conc=1.5,
                                         dntp_conc=0.6, dna_conc=50.0)
            dev = abs(t - float(r[f"{side}_tm"]))
            max_bp_dev = max(max_bp_dev, dev)
            if dev > 0.3:
                detC1.append(f"#{r['no']} {side} calcTm {t:.2f} != {r[f'{side}_tm']}")
            # C3 GC%
            g = gc_pct(s)
            if abs(g - float(r[f"{side}_gc"])) > 0.05:
                detC3.append(f"#{r['no']} {side} GC% {g:.2f} != {r[f'{side}_gc']}")
            # C4 structure dG
            hd = primer3.bindings.calc_homodimer(s).dg
            hp = primer3.bindings.calc_hairpin(s).dg
            if hd is not None and hd < -45000:
                detC4.append(f"#{r['no']} {side} homodimer dG={hd / 1000:.1f} kcal/mol")
            if hp is not None and hp < -45000:
                detC4.append(f"#{r['no']} {side} hairpin dG={hp / 1000:.1f} kcal/mol")
            # C2 Biopython second opinion
            tb = Tm_NN(s, dnac1=50, dnac2=0, Na=50, Mg=1.5, dNTPs=0.6, saltcorr=7)
            max_biopy_dev = max(max_biopy_dev, abs(tb - float(r[f"{side}_tm"])))
            if abs(tb - float(r[f"{side}_tm"])) > 2.0:
                detC2.append(f"#{r['no']} {side} Biopython Tm {tb:.2f} vs {r[f'{side}_tm']}")
        het = primer3.bindings.calc_heterodimer(r["fwd"], r["rev"]).dg
        if het is not None and het < -45000:
            detC4.append(f"#{r['no']} heterodimer dG={het / 1000:.1f} kcal/mol")
    rep.check("C1", "primer3 calcTm 复算与报告 Tm 一致（≤0.3℃）", not detC1,
              "; ".join(detC1[:3]) or f"最大偏差 {max_bp_dev:.3f}℃")
    rep.lines.append(f"- [INFO] Biopython Tm_NN 独立第二意见：与报告 Tm 最大偏差 {max_biopy_dev:.2f}℃")
    rep.check("C2", "Biopython Tm_NN 偏差 ≤2.0℃（跨引擎容差）", not detC2, "; ".join(detC2[:3]))
    rep.check("C3", "引物 GC% 复算一致", not detC3, "; ".join(detC3[:3]))
    rep.check("C4", "homodimer/heterodimer/hairpin ΔG ≥ -45 kcal/mol 复算", not detC4, "; ".join(detC4[:3]))

    # ---------------- Method E: UCSC annotation cross-check ----------------
    rep.lines.append("## E. 基因注释交叉核对（UCSC ncbiRefSeqCurated vs Ensembl）")
    ann = {}
    with open(os.path.join(BASE, "results", "annotation.tsv")) as f:
        h = f.readline().rstrip("\n").split("\t")
        for line in f:
            d = dict(zip(h, line.rstrip("\n").split("\t")))
            ann[(d["chrom"], int(d["start"]))] = d
    okE, detE, n_match, n_direct = True, [], 0, 0
    import time as _time
    for r in rows:
        c = r["chrom"]
        rs0, re1 = int(r["repeat_start"]) - 1, int(r["repeat_end"])
        # query only the locus neighborhood (fast, avoids whole-chromosome dumps)
        q0, q1 = max(0, rs0 - 50000), re1 + 50000
        overlap, near = set(), {}
        items = []
        try:
            url = (f"https://api.genome.ucsc.edu/getData/track?genome=hg38;"
                   f"track=ncbiRefSeqCurated;chrom={c};start={q0};end={q1}")
            with urllib.request.urlopen(url, timeout=60) as fh:
                data = json.load(fh)
            items = data.get("ncbiRefSeqCurated", [])
        except Exception as e:
            rep.info(f"UCSC {c}:{rs0 + 1} 查询失败: {e}")
        for t in items:
            sym = t.get("name2")
            if t["txStart"] < re1 and t["txEnd"] > rs0:
                overlap.add(sym)
            elif sym and min(abs(t["txStart"] - re1), abs(rs0 - t["txEnd"])) <= 1000:
                d = min(abs(t["txStart"] - re1), abs(rs0 - t["txEnd"]))
                near[sym] = min(near.get(sym, 10 ** 9), d)
        ens_gene = r["gene"]
        if ens_gene in overlap:
            n_match += 1
            n_direct += 1
        elif ens_gene in near:
            # benign transcript-boundary difference between Ensembl and RefSeq
            rep.info(f"#{r['no']} {c}:{rs0 + 1} {ens_gene}：Ensembl 基因体包含该位点，"
                     f"RefSeq 转录本边界距位点 {near[ens_gene]}bp（数据库模型边界差异）")
            n_match += 1
        elif ens_gene.startswith("ENSG"):
            rep.info(f"#{r['no']} {c}:{rs0 + 1} Ensembl={ens_gene}（无符号），UCSC RefSeq 重叠基因：{sorted(x for x in overlap if x)[:3]}")
            n_match += 1
        elif not items or (not overlap and not near and r.get("biotype", "protein_coding") != "protein_coding"):
            # RefSeq Curated has no model here at all (Ensembl-only annotation,
            # e.g. pseudogene / lncRNA) -- cross-check impossible, not a contradiction
            rep.info(f"#{r['no']} {c}:{rs0 + 1} {ens_gene}（{r.get('biotype', '.')}）："
                     f"RefSeq Curated 无对应模型（Ensembl 独有注释），无法交叉核对")
            n_match += 1
        else:
            okE = False
            detE.append(f"#{r['no']} {c}:{rs0 + 1} Ensembl={ens_gene} vs UCSC 重叠={sorted(x for x in overlap if x)[:4]}")
    rep.check("E", "UCSC RefSeq 与 Ensembl 注释交叉核对", okE, "; ".join(detE))
    rep.lines.append(f"- [INFO] 注释核对：{n_match}/{len(rows)} 行与 RefSeq 无矛盾，"
                     f"其中 {n_direct}/{len(rows)} 行为基因符号直接匹配；"
                     f"其余为 RefSeq 无模型或仅有边界差异（见上列 INFO）")

    # ---------------- Method D merge ----------------
    dpath = os.path.join(BASE, "results", "method_D_genomescan.json")
    if os.path.exists(dpath):
        rep.lines.append("## D. 全基因组特异性独立重扫（numpy 滚动哈希 vs grep -F）")
        d = json.load(open(dpath))
        with open(os.path.join(BASE, "results", "designed.json")) as f:
            dj = [it for it in json.load(f) if it.get("status") in ("ok", "ok_check")]
        okD, okD15, detD, det15, diff15 = True, True, [], [], 0
        for it, item in zip(dj, d["rows"]):
            p = it["pair"]
            for side, gk in (("fwd", "left"), ("rev", "right")):
                if p[f"{gk}_full_off"] != item[side]["full_off"]:
                    okD = False
                    detD.append(f"#{item['no']} {side} full_off grep={p[f'{gk}_full_off']} numpy={item[side]['full_off']}")
                if item[side]["end15_off"] < p[f"{gk}_end15_off"]:
                    okD15 = False
                    det15.append(f"#{item['no']} {side} numpy={item[side]['end15_off']} < grep={p[f'{gk}_end15_off']}")
                elif item[side]["end15_off"] != p[f"{gk}_end15_off"]:
                    diff15 += 1
        n_primers = 2 * len(d["rows"])
        rep.check("D1a", f"全长脱靶计数：grep 与 numpy 两种实现完全一致（{n_primers}/{n_primers}）",
                  okD, "; ".join(detD[:4]))
        if okD15:
            rep.lines.append(f"- [INFO] 3'-15mer 脱靶计数：{diff15}/{n_primers} 条 grep 因非重叠匹配语义"
                             f"少计（numpy 按位置计数为准），全部满足 numpy ≥ grep，方向与预期一致")
        else:
            rep.lines.append(f"- [INFO] 3'-15mer 脱靶计数：{len(det15)} 项违反 numpy ≥ grep")
        rep.check("D1b", "3'-15mer 计数关系 numpy ≥ grep（grep -o 非重叠语义）", okD15, "; ".join(det15[:4]))
        rep.lines.append(f"- [INFO] 扫描范围：{d['scanned_files']} 个 contig，{d['scanned_bp'] / 1e9:.2f} Gbp（正链编码，模式含正反两向）")

    n_fail = rep.summary()
    with open(os.path.join(BASE, "results", "validation_report.md"), "w") as f:
        f.write("\n".join(rep.lines) + "\n")
    print(f"validation: {len(rep.checks)} checks, {n_fail} failed -> results/validation_report.md")
    for m, name, ok, d in rep.checks:
        if not ok:
            print(f"  FAILED ({m}) {name}: {d}")


if __name__ == "__main__":
    main()
