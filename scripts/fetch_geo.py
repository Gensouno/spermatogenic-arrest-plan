#!/usr/bin/env python3
"""Resumable parallel-range downloader for GEO FTP.

Why this exists
---------------
NCBI GEO FTP throttles a single stream to ~16-35 KB/s. Aggregating N range
requests lifts that to ~100-130 KB/s (measured), so parallel ranges are worth it.

Why it is written this way (both lessons came from real failures)
----------------------------------------------------------------
1. The first version accumulated every chunk in a Python dict and only wrote the
   file at the very end. One transient HTTP 503 on the last chunk threw away
   20 minutes of downloading. This version writes every chunk to its own part
   file the moment it arrives, so partial progress survives.
2. The first chunking used `seq` in bash, which silently produced nothing on
   macOS/BSD for some ranges and yielded a 0-byte file. Chunk arithmetic is now
   done in Python only.

It is therefore safe to re-run: completed parts are detected and skipped, and
the final file is assembled only when every byte is present.

Stdlib only, so it runs under /usr/bin/python3 with no venv.

Usage:
    python3 scripts/fetch_geo.py <url> <dest> [--streams 6] [--chunk 3000000]
"""
import argparse
import os
import sys
import time
import urllib.error
import urllib.request
from concurrent.futures import ThreadPoolExecutor, as_completed

CHUNK = 3_000_000
STREAMS = 6
RETRIES = 6


def head_size(url: str) -> int | None:
    req = urllib.request.Request(url, method="HEAD")
    try:
        with urllib.request.urlopen(req, timeout=60) as r:
            n = r.headers.get("Content-Length")
            return int(n) if n else None
    except Exception as e:  # noqa: BLE001
        print(f"  HEAD failed: {e}", file=sys.stderr)
        return None


def fetch_range(url, start, end, tries=RETRIES):
    """Fetch one byte range with exponential backoff. Returns bytes or raises."""
    want = end - start + 1
    last = None
    for k in range(tries):
        try:
            req = urllib.request.Request(url, headers={"Range": f"bytes={start}-{end}"})
            with urllib.request.urlopen(req, timeout=600) as r:
                data = r.read()
            if len(data) == want:
                return data
            last = f"short read {len(data)} != {want}"
        except Exception as e:  # noqa: BLE001
            last = f"{type(e).__name__}: {e}"
        # back off harder for server-side throttling
        time.sleep(min(3 * (k + 1), 30))
    raise RuntimeError(f"range {start}-{end} failed: {last}")


def gzip_ok(path: str, max_bytes: int | None = None) -> bool:
    """True if the file decompresses cleanly (CRCs validated).

    A gzip stream can only be validated by decompressing it, so this does one
    full sequential pass (output is discarded). That costs roughly 1-2 min for
    a 3.7 GB decompressed payload -- trivial next to the ~15 min download, and
    it is the check that would have caught the truncation that corrupted the
    reference matrix.
    """
    import gzip
    try:
        with gzip.open(path, "rb") as fh:
            while True:
                if not fh.read(1 << 22):
                    break
        return True
    except Exception:  # noqa: BLE001
        return False


def download(url, dest, streams=STREAMS, chunk=CHUNK, max_rounds=8):
    size = head_size(url)
    if size is None:
        raise SystemExit(f"cannot determine size of {url}")

    # If the destination is already complete, do not touch it. This is the
    # guard that was missing when a stale .parts directory caused a good file
    # to be overwritten with a truncated reassembly.
    if os.path.exists(dest) and os.path.getsize(dest) == size:
        if not dest.endswith(".gz") or gzip_ok(dest):
            print(f"already complete: {dest} ({size/1e6:.1f} MB) -- skipping")
            return
        print(f"WARNING: {dest} has the right size but is not a valid gzip; "
              f"re-downloading", file=sys.stderr)

    partdir = dest + ".parts"
    os.makedirs(partdir, exist_ok=True)
    ranges = [(o, min(o + chunk - 1, size - 1)) for o in range(0, size, chunk)]
    print(f"size = {size/1e6:.1f} MB | {len(ranges)} chunks of "
          f"{chunk/1e6:.0f} MB | {streams} streams")
    print(f"parts dir: {partdir}  (re-run to resume)")

    # Drop leftover .tmp scratch files: a killed worker can leave one behind and
    # a stale .tmp is never valid input for assembly.
    for j, (a, b) in enumerate(ranges):
        for junk in (os.path.join(partdir, f"part_{j:05d}") + ".tmp",):
            if os.path.exists(junk):
                os.remove(junk)

    t0 = time.time()
    # Multiple passes: a chunk that fails in one pass is retried in the next.
    for round_no in range(1, max_rounds + 1):
        todo = []
        for i, (a, b) in enumerate(ranges):
            p = os.path.join(partdir, f"part_{i:05d}")
            want = b - a + 1
            if os.path.exists(p) and os.path.getsize(p) == want:
                continue
            todo.append((i, a, b, p))
        if not todo:
            break
        got_before = sum(
            os.path.getsize(os.path.join(partdir, f"part_{i:05d}"))
            for i, _ in enumerate(ranges)
            if os.path.exists(os.path.join(partdir, f"part_{i:05d}"))
        )
        print(f"[round {round_no}] {len(todo)} chunk(s) to fetch "
              f"({got_before/1e6:.1f}/{size/1e6:.1f} MB already on disk)")

        def work(item):
            i, a, b, p = item
            try:
                data = fetch_range(url, a, b)
                tmp = p + ".tmp"
                with open(tmp, "wb") as fh:
                    fh.write(data)
                os.replace(tmp, p)
                return (i, None)
            except Exception as e:  # noqa: BLE001
                return (i, str(e))

        with ThreadPoolExecutor(max_workers=streams) as ex:
            futs = [ex.submit(work, it) for it in todo]
            done = 0
            for f in as_completed(futs):
                i, err = f.result()
                done += 1
                if err:
                    print(f"  chunk {i} FAILED: {err}", file=sys.stderr)
                else:
                    got = sum(
                        os.path.getsize(os.path.join(partdir, f"part_{j:05d}"))
                        for j, _ in enumerate(ranges)
                        if os.path.exists(os.path.join(partdir, f"part_{j:05d}"))
                    )
                    el = max(time.time() - t0, 1e-9)
                    print(f"  [{done}/{len(todo)}] {got/1e6:.1f} MB "
                          f"{got/1024/el:.0f} KB/s avg", flush=True)

    # Assemble
    missing = [i for i, (a, b) in enumerate(ranges)
               if not os.path.exists(os.path.join(partdir, f"part_{i:05d}"))
               or os.path.getsize(os.path.join(partdir, f"part_{i:05d}")) != b - a + 1]
    if missing:
        raise SystemExit(f"still missing {len(missing)} chunk(s): {missing[:5]}... "
                         f"re-run this command to resume")

    tmp = dest + ".assembling"
    with open(tmp, "wb") as out:
        for i, (a, b) in enumerate(ranges):
            with open(os.path.join(partdir, f"part_{i:05d}"), "rb") as fh:
                out.write(fh.read())
    if os.path.getsize(tmp) != size:
        raise SystemExit(f"size mismatch: got {os.path.getsize(tmp)}, want {size}")
    if dest.endswith(".gz") and not gzip_ok(tmp):
        raise SystemExit(
            f"{dest} assembled to the right size but is not a valid gzip "
            f"stream; refusing to install it. Delete {partdir} and re-run.")
    os.replace(tmp, dest)
    el = time.time() - t0
    print(f"OK -> {dest}  ({size/1e6:.1f} MB in {el:.0f}s, "
          f"{size/1024/max(el,1e-9):.0f} KB/s)")
    # Remove the parts directory once the destination is known-good. Keeping it
    # is what let a later run rebuild a truncated file over a good one.
    import shutil
    shutil.rmtree(partdir, ignore_errors=True)
    print(f"removed {partdir}")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("url")
    ap.add_argument("dest")
    ap.add_argument("--streams", type=int, default=STREAMS)
    ap.add_argument("--chunk", type=int, default=CHUNK)
    ap.add_argument("--rounds", type=int, default=8)
    a = ap.parse_args()
    download(a.url, a.dest, a.streams, a.chunk, a.rounds)


if __name__ == "__main__":
    main()
