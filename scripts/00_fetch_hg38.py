#!/usr/bin/env python3
"""Fetch the UCSC hg38 reference and build the layout the pipeline expects.

The panel's coordinates and all method_D / specificity counts are expressed in
UCSC hg38.  The design session had no full genome locally (the only local copy
was GRCh38 no_alt chr1-chr22, which lacks chrX -- 10 of the 24 loci live on
chrX).  Mixing the two file sets is exactly the reference-mismatch problem this
script removes: every downstream step reads one reference, from one source.

Produces, under {BASE}/data/:
    downloads/hg38.chromFa.tar.gz    tarball (md5 verified against UCSC)
    hg38/chromFa/{name}.fa           455 per-contig FASTA files, as shipped
    hg38_oneline/{name}.fa           header + single upper-case sequence line

Output is written to a fresh directory and never merged with an existing one;
the archive is treated as data, never executed.

Usage:  python3 scripts/00_fetch_hg38.py [--force]
"""
import argparse
import concurrent.futures as cf
import hashlib
import json
import os
import shutil
import sys
import tarfile
import time
import urllib.error
import urllib.request

BASE = os.environ.get("PRIMER_DESIGN_BASE") or os.path.dirname(
    os.path.dirname(os.path.abspath(__file__)))
DATA = os.path.join(BASE, "data")
DEST = os.path.join(DATA, "downloads")
TARBALL = os.path.join(DEST, "hg38.chromFa.tar.gz")

BASE_URL = "https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips"
ARCHIVE = "hg38.chromFa.tar.gz"
CONNECTIONS = 16
# Chunk size is bounded by the per-connection throughput UCSC actually gives
# (single stream measured ~20-30 KB/s): a chunk must finish well inside the
# deadline or every request retries forever.  1 MiB needs ~50 s at that rate,
# comfortably under CHUNK_DEADLINE; 4 MiB did not.
CHUNK = 1 << 20
SOCKET_TIMEOUT = 45       # no bytes at all for this long -> give up on the socket
CHUNK_DEADLINE = 150      # wall clock for one whole chunk, progress or not
RETRIES = 6


def head_len(url):
    req = urllib.request.Request(url, method="HEAD")
    with urllib.request.urlopen(req, timeout=60) as r:
        return int(r.headers["Content-Length"])


def fetch_range(url, start, end, tries=RETRIES):
    """Return bytes [start, end] inclusive, retrying transient failures.

    urlopen's timeout bounds a single socket operation, so a server that
    trickles a few bytes at a time never trips it and the chunk hangs forever
    (observed: the last 164 chunks sat with zero progress).  A wall-clock
    deadline for the whole chunk is what actually bounds a stalled transfer.
    """
    want = end - start + 1
    last = None
    for attempt in range(tries):
        try:
            req = urllib.request.Request(
                url, headers={"Range": f"bytes={start}-{end}"})
            buf = bytearray()
            deadline = time.time() + CHUNK_DEADLINE
            with urllib.request.urlopen(req, timeout=SOCKET_TIMEOUT) as r:
                while len(buf) < want:
                    if time.time() > deadline:
                        raise TimeoutError(
                            f"chunk deadline exceeded at {len(buf)}/{want} bytes")
                    piece = r.read(262144)
                    if not piece:
                        break
                    buf.extend(piece)
            if len(buf) == want:
                return bytes(buf)
            last = f"short read {len(buf)} != {want}"
        except Exception as e:                       # noqa: BLE001
            last = repr(e)
        time.sleep(min(2 ** attempt, 20))
    raise RuntimeError(f"range {start}-{end} failed after {tries} tries: {last}")


def mirror_md5():
    """md5 of the archive, as published by UCSC (md5sum.txt next to bigZips)."""
    try:
        with urllib.request.urlopen(f"{BASE_URL}/md5sum.txt", timeout=60) as r:
            for line in r.read().decode().splitlines():
                parts = line.split()
                if len(parts) == 2 and parts[1].lstrip("./") == ARCHIVE:
                    return parts[0]
    except Exception as e:                           # noqa: BLE001
        print(f"  (md5sum.txt unavailable: {e})", file=sys.stderr)
    return None


def md5_file(path, blk=1 << 20):
    h = hashlib.md5()
    with open(path, "rb") as f:
        for b in iter(lambda: f.read(blk), b""):
            h.update(b)
    return h.hexdigest()


def download(url, dest, connections=CONNECTIONS):
    total = head_len(url)
    state_path = dest + ".state.json"
    state = {}
    if os.path.exists(state_path):
        with open(state_path) as f:
            state = json.load(f)
    if state.get("total") != total or state.get("chunk") != CHUNK:
        # chunk boundaries depend on CHUNK, so a changed size invalidates
        # every recorded index -- restart rather than write into the wrong slots
        state = {"total": total, "chunk": CHUNK, "done": []}
    done = set(state.get("done", []))

    ranges = [(s, min(s + CHUNK - 1, total - 1))
              for s in range(0, total, CHUNK)]
    todo = [(i, s, e) for i, (s, e) in enumerate(ranges) if i not in done]
    print(f"  {ARCHIVE}: {total / 1e6:.1f} MB, {len(ranges)} chunks, "
          f"{len(todo)} to fetch, {connections} connections", file=sys.stderr)

    fd = os.open(dest, os.O_RDWR | os.O_CREAT, 0o644)
    os.ftruncate(fd, total)
    t0 = time.time()
    got = 0
    lock_free = True
    try:
        with cf.ThreadPoolExecutor(max_workers=connections) as ex:
            futs = {ex.submit(fetch_range, url, s, e): (i, s, e)
                    for i, s, e in todo}
            for fut in cf.as_completed(futs):
                i, s, e = futs[fut]
                buf = fut.result()
                os.pwrite(fd, buf, s)
                done.add(i)
                got += len(buf)
                if len(done) % 25 == 0 or len(done) == len(ranges):
                    el = time.time() - t0
                    json.dump({"total": total, "chunk": CHUNK, "done": sorted(done)},
                              open(state_path, "w"))
                    print(f"    {len(done)}/{len(ranges)} chunks, "
                          f"{got / 1e6:.0f} MB this run, "
                          f"{got / max(el, 1) / 1024:.0f} KB/s", file=sys.stderr)
    finally:
        os.close(fd)
        del lock_free
    if len(done) == len(ranges):
        json.dump({"total": total, "chunk": CHUNK, "done": sorted(done)}, open(state_path, "w"))


def build_oneline(chromfa_dir, out_dir):
    """header + one upper-case sequence line, matching the pipeline's reader.

    Upper-casing is required: the specificity scanners compare primer strings
    against the raw text with ==, and hg38.chromFa ships soft-masked (lower
    case) repeats -- exactly the regions this panel targets.
    """
    os.makedirs(out_dir, exist_ok=True)
    n = 0
    for fn in sorted(os.listdir(chromfa_dir)):
        if not fn.endswith(".fa"):
            continue
        name = fn[:-3]
        with open(os.path.join(chromfa_dir, fn)) as f:
            header = f.readline().rstrip("\n")
            seq = "".join(l.strip() for l in f).upper()
        with open(os.path.join(out_dir, fn), "w") as f:
            f.write(f"{header}\n{seq}\n")
        n += 1
    return n


# --- fallback route: the UCSC REST API -------------------------------------
# The bulk host (hgdownload) rate-limits sustained parallel downloads -- it cut
# us off mid-transfer after ~800 MB.  The REST API is a separate service on the
# same assembly, so it backfills whatever the tarball did not deliver, and the
# reference stays single-sourced on UCSC.
UCSC_API = "https://api.genome.ucsc.edu/getData/sequence"
API_SPAN = 5_000_000        # measured: this endpoint serves 5 Mb per request
API_CONNECTIONS = 3     # the REST API 429s quickly above this
PRIMARY = [f"chr{i}" for i in range(1, 23)] + ["chrX", "chrY", "chrM"]


def api_chrom_sizes():
    url = "https://api.genome.ucsc.edu/list/chromosomes?genome=hg38"
    with urllib.request.urlopen(url, timeout=120) as r:
        return json.loads(r.read().decode())["chromosomes"]


def fetch_chrom_via_api(chrom, size, connections=API_CONNECTIONS):
    blocks = [(s, min(s + API_SPAN, size)) for s in range(0, size, API_SPAN)]

    def one(i):
        s, e = blocks[i]
        last = None
        for attempt in range(RETRIES):
            try:
                url = (f"{UCSC_API}?genome=hg38;chrom={chrom}"
                       f";start={s};end={e}")
                with urllib.request.urlopen(url, timeout=180) as r:
                    seq = json.loads(r.read().decode())["dna"]
                if len(seq) != e - s:
                    raise ValueError(f"got {len(seq)} bases, want {e - s}")
                return i, seq.upper()
            except urllib.error.HTTPError as ex:
                last = repr(ex)
                # 429 means we are over the rate limit: back off hard rather
                # than retrying into it, which only extends the block.
                time.sleep(60 * (attempt + 1) if ex.code == 429
                           else min(2 ** attempt, 20))
            except Exception as ex:                   # noqa: BLE001
                last = repr(ex)
                time.sleep(min(2 ** attempt, 20))
        raise RuntimeError(f"api {chrom}:{s}-{e} failed: {last}")

    parts = [None] * len(blocks)
    done = 0
    t0 = time.time()
    with cf.ThreadPoolExecutor(max_workers=connections) as ex:
        for i, seq in ex.map(one, range(len(blocks))):
            parts[i] = seq
            done += 1
            if done % 8 == 0 or done == len(blocks):
                print(f"    api {chrom}: {done}/{len(blocks)} blocks, "
                      f"{time.time() - t0:.0f}s", file=sys.stderr)
    return "".join(parts)


def install_chrom(chrom, seq, chromfa, oneline):
    for d in (chromfa, oneline):
        os.makedirs(d, exist_ok=True)
        with open(os.path.join(d, f"{chrom}.fa"), "w") as f:
            f.write(f">{chrom}\n{seq}\n")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--force", action="store_true",
                    help="re-download and re-extract even if already present")
    ap.add_argument("--connections", type=int, default=CONNECTIONS)
    ap.add_argument("--skip-tarball", action="store_true",
                    help="go straight to the REST API route")
    ap.add_argument("--api-chroms", default="",
                    help="comma list to (re)fetch via the REST API; default is "
                         "whichever primary chromosomes are missing locally")
    args = ap.parse_args()

    chromfa = os.path.join(DATA, "hg38", "chromFa")
    oneline = os.path.join(DATA, "hg38_oneline")
    url = f"{BASE_URL}/{ARCHIVE}"
    os.makedirs(DEST, exist_ok=True)
    os.makedirs(chromfa, exist_ok=True)
    os.makedirs(oneline, exist_ok=True)

    if not args.skip_tarball:
        want = mirror_md5()
        if want:
            print(f"UCSC md5sum.txt: {ARCHIVE} = {want}", file=sys.stderr)
        try:
            if args.force or not os.path.exists(TARBALL) or (
                    want and md5_file(TARBALL) != want):
                download(url, TARBALL, args.connections)
            got = md5_file(TARBALL)
            if want and got != want:
                raise RuntimeError(f"md5 mismatch: got {got}, want {want}")
            print(f"  md5 {'VERIFIED' if want else 'computed'}: {got}",
                  file=sys.stderr)

            tmp = os.path.join(DATA, "hg38.extract.tmp")
            if os.path.exists(tmp):
                shutil.rmtree(tmp)
            os.makedirs(tmp)
            with tarfile.open(TARBALL) as tf:
                members = [m for m in tf.getmembers()
                           if m.isfile() and m.name.endswith(".fa")]
                for m in members:
                    m.name = os.path.basename(m.name)      # flatten chromFa/
                try:
                    tf.extractall(tmp, members=members, filter="data")
                except TypeError:                          # Python < 3.12
                    tf.extractall(tmp, members=members)
            if os.path.isdir(chromfa):
                shutil.rmtree(chromfa)
            shutil.move(tmp, chromfa)
            print(f"  extracted {len(members)} contigs -> {chromfa}",
                  file=sys.stderr)
        except Exception as e:                             # noqa: BLE001
            print(f"  tarball route unusable ({e}); using the REST API route",
                  file=sys.stderr)

    # ---- REST API backfill ----
    sizes = api_chrom_sizes()
    targets = (args.api_chroms.split(",") if args.api_chroms
               else [c for c in PRIMARY
                     if not os.path.exists(os.path.join(chromfa, f"{c}.fa"))])
    for c in targets:
        if not c:
            continue
        if os.path.exists(os.path.join(chromfa, f"{c}.fa")) and not args.force:
            continue
        size = sizes.get(c)
        if not size:
            print(f"  skip {c}: unknown size", file=sys.stderr)
            continue
        print(f"  {c}: fetching {size / 1e6:.1f} Mb via UCSC REST API",
              file=sys.stderr)
        seq = fetch_chrom_via_api(c, size, args.connections)
        if len(seq) != size:
            sys.exit(f"FATAL: {c} assembled length {len(seq)} != {size}")
        install_chrom(c, seq, chromfa, oneline)
        print(f"    {c}: {len(seq)} bp installed", file=sys.stderr)

    n = build_oneline(chromfa, oneline)
    total_bp = sum(os.path.getsize(os.path.join(oneline, f)) - 1
                   for f in os.listdir(oneline))
    print(f"  oneline: {n} contigs, {total_bp / 1e9:.3f} Gbp", file=sys.stderr)
    have = [c for c in PRIMARY if os.path.exists(os.path.join(oneline, f"{c}.fa"))]
    missing = [c for c in PRIMARY if c not in have]
    if missing:
        print(f"  WARNING: primary chromosomes missing: {', '.join(missing)}",
              file=sys.stderr)
    else:
        print("  all 24 primary chromosomes present", file=sys.stderr)
    print("OK", file=sys.stderr)


if __name__ == "__main__":
    main()
