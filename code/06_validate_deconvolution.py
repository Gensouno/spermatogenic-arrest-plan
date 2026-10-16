#!/usr/bin/env python3
"""
06_validate_deconvolution.py -- end-to-end validation of the deconvolution step
on a REAL whole-testis dataset, with the project's positive control.

WHY THIS EXISTS
---------------
The pipeline has one make-or-break assumption: that the GSE112393 single-cell
reference can recover a *biologically meaningful* cell-composition shift from
bulk whole-testis RNA-seq. If it cannot, Act II has no foundation.

So before committing to the full plan we ran the smallest possible real test:
take one whole-testis KO/WT dataset that we already have on disk, deconvolve it,
and check the THE PROJECT'S OWN POSITIVE CONTROL --

    a germ-cell-arrest knockout should show a DECREASE in the germ-cell
    (Spermatogonia / Meiotic / Postmeiotic) fractions and a CONSEQUENT
    INCREASE in the Somatic fraction and/or a disappearance of the post-
    meiotic compartment.

If the estimated composition does not move in that direction, the reference is
not usable and the whole three-act design must be rethought. This was written
*after* the The deconvolution spec (docs/03-null-model.md) flagged that an un-normalised
signature inflates the Somatic fraction by up to 0.26 -- hence the CPM signature
is used here.

WHAT IT DOES
------------
Per sample: linear-space NNLS on log2(CPM+1) values,
    min_p || log2(S_cpm + 1)^T p - log2(CPM_bulk + 1) ||^2,  p >= 0,  sum(p)=1
Reported per sample: proportions + reconstruction R^2.
Then the WT->KO delta in each cell type is printed, which IS the positive control.

Usage
-----
    MPLCONFIGDIR="$PWD/.mplcache" .venv/bin/python code/06_validate_deconvolution.py \
        --sig data/ref/signature_GSE112393_cpm.csv.gz \
        --bulk data/raw/GSE236388/GSE236388_XAP5_counts.txt.gz \
        --wt-pattern WT --ko-pattern cKO
"""

import argparse
import re
import sys

import numpy as np
import pandas as pd
from scipy.optimize import nnls


def load_bulk(path):
    X = pd.read_csv(path, sep="\t", index_col=0)
    dup = X.index.duplicated().sum()
    if dup:
        # Gene symbols legitimately collapse to the same symbol; sum their counts
        # rather than silently keeping the first (which loses expression).
        X = X.groupby(level=0).sum()
    return X


def deconvolve(S, X):
    cpm = X / X.sum(axis=0) * 1e6
    lcpm = np.log2(cpm + 1)
    common = [g for g in S.columns if g in lcpm.index]
    A = np.log2(S[common].to_numpy() + 1).T          # G x K
    Y = lcpm.loc[common].to_numpy()                  # G x N
    props, r2s = [], []
    for j in range(Y.shape[1]):
        y = Y[:, j]
        p, _ = nnls(A, y)
        if p.sum() > 0:
            p = p / p.sum()
        props.append(p)
        r2s.append(1 - ((y - A @ p) ** 2).sum() / max(((y - y.mean()) ** 2).sum(), 1e-12))
    P = pd.DataFrame(np.array(props).T, index=S.index, columns=X.columns)
    return P, np.array(r2s), len(common)


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--sig", required=True)
    ap.add_argument("--bulk", required=True)
    ap.add_argument("--wt-pattern", default="WT")
    ap.add_argument("--ko-pattern", default="KO")
    ap.add_argument("--out", default=None)
    a = ap.parse_args()

    S = pd.read_csv(a.sig, index_col=0)
    S = S.loc[:, ~S.columns.duplicated()]
    print(f"signature {a.sig}: {S.shape[0]} cell types x {S.shape[1]} genes (CPM)")

    X = load_bulk(a.bulk)
    print(f"bulk {a.bulk}: {X.shape[0]} genes x {X.shape[1]} samples")

    P, r2, ncommon = deconvolve(S, X)
    print(f"shared genes used: {ncommon}\n")

    print(f"{'sample':32s} {'R2':>6s}  " + "  ".join(f"{c[:6]:>6s}" for c in S.index))
    for s in P.columns:
        print(f"{s:32s} {r2[list(P.columns).index(s)]:6.3f}  "
              + "  ".join(f"{v:6.3f}" for v in P[s]))

    wt = [c for c in P.columns if re.search(a.wt_pattern, c, re.I)]
    ko = [c for c in P.columns if re.search(a.ko_pattern, c, re.I)]
    wt = [c for c in wt if c not in ko]
    print(f"\nWT samples ({len(wt)}): {wt}")
    print(f"KO samples ({len(ko)}): {ko}")

    if wt and ko:
        print(f"\n=== POSITIVE CONTROL: Delta P = mean(KO) - mean(WT) ===")
        d = P[ko].mean(axis=1) - P[wt].mean(axis=1)
        for ct, v in d.items():
            arrow = "DOWN" if v < -0.02 else ("UP" if v > 0.02 else "flat")
            print(f"  {ct:16s} WT={P.loc[ct, wt].mean():6.3f}  KO={P.loc[ct, ko].mean():6.3f}  "
                  f"delta={v:+7.3f}  {arrow}")
        germ = [c for c in S.index if c.lower() not in ("somatic",)]
        gd = d[germ].sum()
        print(f"\n  sum of germ-cell fractions: delta = {gd:+.3f}")
        print(f"  Somatic (somatic)          : delta = {d.get('Somatic', np.nan):+.3f}")
        ok = gd < 0
        print(f"\n  POSITIVE CONTROL {'PASS' if ok else 'FAIL'}: "
              f"germ-cell fraction {'decreased' if ok else 'did NOT decrease'} under KO")

        if ok and gd < -0.05:
            print("  => the reference recovers a biologically sensible composition shift.")
            print("     Act II (composition null model) has a working foundation.")
        elif ok:
            print("  => direction correct but small; check whether the knockout is expected")
            print("     to produce a large block before tuning anything.")

    print(f"\nreconstruction R^2: min={r2.min():.3f} median={np.median(r2):.3f} max={r2.max():.3f}")
    print(f"  (runbook acceptance threshold: median > 0.6)")

    if a.out:
        P.to_csv(a.out, float_format="%.6g")
        print(f"wrote {a.out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
