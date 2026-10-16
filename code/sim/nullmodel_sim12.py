#!/usr/bin/env python3
"""v12 FINAL: two remaining specifics.
 A. pseudobulk-counts signature S causes depth/RNA-content bias in NNLS --
    quantify it and test the CPM fix.
 B. the compShare ratio is unstable when cov(lfc_A,lfc_B) ~ 0; give the stable
    alternative and its decision rule.
"""
import numpy as np
from scipy.optimize import nnls

GLYPH = "=" * 74
G, K = 4000, 4
BASE = np.array([0.15, 0.10, 0.20, 0.55])
CONC, NW, NK, NOISE = 250.0, 4, 4, 0.18
_r = np.random.default_rng(11)
base_expr = _r.gamma(2.0, 50.0, size=G)
S_clean = np.tile(base_expr[:, None], (1, K))
for _k in range(K):
    _i = _r.choice(G, 200, replace=False)
    _a = 6.0 * (1 + .6 * _k)
    S_clean[_i, _k] *= _r.uniform(_a * .6, _a * 1.6, 200)

# Different cell types have different total RNA / library depth. In pseudobulk
# COUNTS this shows up as a per-cell-type scale factor.
DEPTH = np.array([0.6, 1.0, 1.8, 3.2])     # somatic low, postmeiotic high
S_depth = S_clean * DEPTH[None, :]

SH = {"pre": np.array([+.02, +.10, -.04, -.08]),
      "mei": np.array([+.09, +.02, +.05, -.16]),
      "post": np.array([+.06, +.04, -.06, -.04])}


def make_props(rng, n, shift=None):
    b = BASE + (shift if shift is not None else 0.)
    return rng.dirichlet(np.clip(b, 1e-3, None) * CONC, size=n).T


def deconv(X, Sref):
    P = np.zeros((K, X.shape[1]))
    for j in range(X.shape[1]):
        P[:, j], _ = nnls(Sref, X[:, j])
    t = P.sum(0)
    return np.divide(P, t, out=np.zeros_like(P), where=t > 0)


def cpm(S):
    """columns scaled to a common total (1e6), preserving within-column shape."""
    return S / S.sum(axis=0, keepdims=True) * 1e6


def cov(a, b):
    a = np.asarray(a, float) - np.mean(a); b = np.asarray(b, float) - np.mean(b)
    return float(a @ b / len(a))


print(GLYPH); print("A. pseudobulk-counts S -> depth bias, and the CPM fix"); print(GLYPH)
print(f"{'signature S':>18} {'bulk X':>14} {'MAE(P)':>9} {'somatic bias':>13} {'postmeio bias':>14}")
print("-" * 72)
combos = [
    ("counts", S_depth, "absolute"),
    ("counts", S_depth, "cpm"),
    ("cpm", cpm(S_depth), "absolute"),
    ("cpm", cpm(S_depth), "cpm"),
    ("clean", S_clean, "absolute"),
]
for tag, Sa, xmode in combos:
    maes, bias_s, bias_p = [], [], []
    for rep in range(20):
        r = np.random.default_rng(200 + rep)
        Pt = np.hstack([make_props(r, NW), make_props(r, NK, SH["mei"])])
        X = (S_clean @ Pt) * np.exp(NOISE * r.standard_normal((G, 8)))
        if xmode == "cpm":
            X = cpm(X)
        Pn = deconv(X, Sa)
        maes.append(np.abs(Pn - Pt).mean())
        bias_s.append(np.mean(Pn[0] - Pt[0]))
        bias_p.append(np.mean(Pn[3] - Pt[3]))
    print(f"{tag:>18} {xmode:>14} {np.nanmean(maes):>9.4f} "
          f"{np.nanmean(bias_s):>+13.4f} {np.nanmean(bias_p):>+14.4f}")
print()
print("  Reading: when the atlas is in raw COUNTS and the bulk too, NNLS absorbs")
print("  the per-cell-type depth into the proportion -> systematic bias. CPM on")
print("  BOTH sides is the fix, and it must be applied to the SAME gene set.")

print(); print(GLYPH); print("B. the compShare ratio: instability and the stable alternative"); print(GLYPH)
T = _r.normal(0, 1, G); T -= T.mean()
S_ref = cpm(S_depth * np.exp(0.35 * _r.standard_normal(S_depth.shape)))


def one(node, a, seed):
    r = np.random.default_rng(seed)
    Pw = make_props(r, NW); Pk = make_props(r, NK, SH[node])
    Xw = S_clean @ Pw; Xk = (S_clean @ Pk) * np.exp(0.5 * a * T)[:, None]
    Xw *= np.exp(NOISE * r.standard_normal(Xw.shape))
    Xk *= np.exp(NOISE * r.standard_normal(Xk.shape))
    Xw, Xk = cpm(Xw), cpm(Xk)
    Pn = deconv(np.hstack([Xw, Xk]), S_ref)
    pw, pk = Pn[:, :NW].mean(1), Pn[:, NW:].mean(1)
    lfc = np.log((Xk.mean(1) + 1) / (Xw.mean(1) + 1))
    c = np.log((S_ref @ pk) + 1) - np.log((S_ref @ pw) + 1)
    return lfc, c


print(f"{'a':>6} {'n_models':>9} {'mean compShare':>15} {'sd':>9} "
      f"{'P(cov_obs<=0)':>14}")
print("-" * 58)
for a in (0.0, 0.20, 0.40):
    for nmod in (6, 8):
        vals, zero_frac = [], []
        for rep in range(30):
            ms = [one(list(SH)[i % 3], a, 4000 + rep * 60 + i) for i in range(nmod)]
            for i in range(nmod):
                for j in range(i + 1, nmod):
                    tot = cov(ms[i][0], ms[j][0])
                    zero_frac.append(abs(tot) < 1e-6)
                    vals.append(cov(ms[i][1], ms[j][1]) / tot if abs(tot) > 1e-12 else np.nan)
        print(f"{a:>6.2f} {nmod:>9} {np.nanmean(vals):>15.3f} {np.nanstd(vals):>9.3f} "
              f"{np.mean(zero_frac):>14.3f}")
print()
print("  At a=0 the true cross-model logFC covariance is ~0, so the RATIO blows up")
print("  (mean -3, sd 16 as seen in v11). The ratio is only meaningful when the")
print("  observed cross-model similarity is itself established.")
print()
print("  DECISION RULE: report compShare only for model pairs whose cross-model")
print("  logFC correlation is significantly > 0 by the model-mismatch permutation")
print("  null; for the rest report the covariance DECOMPOSITION as two numbers")
print("  (cov_comp, cov_resid) without dividing.")
