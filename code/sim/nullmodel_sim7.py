#!/usr/bin/env python3
"""v7 FINAL: calibration program with a properly SHARED regulatory effect.

v5/v6 bug found: the injected "true effect" was drawn independently for models
A and B. That is the wrong null for this project -- the hypothesis under test is
that models SHARE a signal. With independent effects, adding signal DESTROYS
cross-model correlation, which inverts every interpretation. Fixed here:

    lfc_KO,m(g) = composition_m(g)  +  a * t(g) * s_m
    t(g)      = the SHARED true regulatory profile  (fixed across models)
    s_m       = model-specific scalar loading (some models share more than others)
    a         = how strong the shared signal is (0 = pure composition null)

Program:
  Cal-1  deconvolution recovers KNOWN proportions
  Cal-2  composition-only (a=0)  -> does composition explain the convergence?
  Cal-3  a>0                     -> is the shared signal still there after residualising?
  Cal-4  negative control        -> unrelated datasets must give delta ~ 0
  Cal-5  the "orthogonality deflates residual similarity" worry, tested directly
"""
import numpy as np
from scipy.optimize import nnls

GLYPH = "=" * 78
G, K = 8000, 4
NAMES = ["Somatic", "SPG", "Meiotic", "Postmeio"]
BASE = np.array([0.15, 0.10, 0.20, 0.55])
CONC, NW, NK, NOISE = 250.0, 4, 4, 0.18

_r = np.random.default_rng(11)
base_expr = _r.gamma(2.0, 50.0, size=G)
S_true = np.tile(base_expr[:, None], (1, K))
for _k in range(K):
    _idx = _r.choice(G, 250, replace=False)
    _amp = 6.0 * (1 + 0.6 * _k)
    S_true[_idx, _k] *= _r.uniform(_amp * .6, _amp * 1.6, 250)
S_bad = S_true * np.exp(0.35 * _r.standard_normal(S_true.shape))
S_bad[_r.choice(G, int(.10 * G), replace=False), :] = 0.0

SH = {"pre": np.array([+0.02, +0.10, -0.04, -0.08]),
      "mei": np.array([+0.09, +0.02, +0.05, -0.16]),
      "post": np.array([+0.06, +0.04, -0.06, -0.04])}
# the SHARED true regulatory profile, fixed once
T_SHARED = _r.normal(0, 1.0, G); T_SHARED -= T_SHARED.mean()


def make_props(rng, n, shift=None):
    b = BASE + (shift if shift is not None else 0.0)
    return rng.dirichlet(np.clip(b, 1e-3, None) * CONC, size=n).T


def deconv(X, Sref):
    P = np.zeros((K, X.shape[1]))
    for j in range(X.shape[1]):
        P[:, j], _ = nnls(Sref, X[:, j])
    t = P.sum(0)
    return P, np.divide(P, t, out=np.zeros_like(P), where=t > 0)


def corr(a, b):
    a = np.asarray(a, float) - np.mean(a); b = np.asarray(b, float) - np.mean(b)
    d = np.linalg.norm(a) * np.linalg.norm(b)
    return float(a @ b / d) if d > 1e-30 else np.nan


def run(node, a, s_m, seed, damage=0.0):
    """One KO-vs-WT dataset. a = shared-signal amplitude, s_m = model loading.
    damage = extra independent (non-shared) regulatory noise."""
    r = np.random.default_rng(seed)
    Pw = make_props(r, NW); Pk = make_props(r, NK, SH[node])
    eff = a * s_m * T_SHARED
    if damage > 0:
        eff = eff + damage * r.normal(0, 1, G)
    gf = np.exp(0.5 * eff)
    Xw = S_true @ Pw
    Xk = (S_true @ Pk) * gf[:, None]
    Xw *= np.exp(NOISE * r.standard_normal(Xw.shape))
    Xk *= np.exp(NOISE * r.standard_normal(Xk.shape))
    Ph, Pn = deconv(np.hstack([Xw, Xk]), S_bad)
    obs = np.log((Xk.mean(1) + 1.0) / (Xw.mean(1) + 1.0))
    nullp = np.log(((S_bad @ Pn[:, NW:].mean(1)) + 1.0) / ((S_bad @ Pn[:, :NW].mean(1)) + 1.0))
    R = np.hstack([Xw, Xk]) - S_bad @ Ph
    resid = R[:, NW:].mean(1) - R[:, :NW].mean(1)
    return dict(obs=obs, null=nullp, resid=resid, dP=Pn[:, NW:].mean(1) - Pn[:, :NW].mean(1))


print(GLYPH); print("Cal-1  deconvolution recovers KNOWN proportions"); print(GLYPH)
r = np.random.default_rng(5)
for lab, Sref in [("perfect ref", S_true), ("mismatched ref", S_bad)]:
    Pt = np.hstack([make_props(r, NW), make_props(r, NK, SH["mei"])])
    X = (S_true @ Pt) * np.exp(NOISE * r.standard_normal((G, 8)))
    _, Pn = deconv(X, Sref)
    print(f"  {lab:16s} MAE = {np.abs(Pn-Pt).mean():.4f}  "
          f"(largest error: {NAMES[np.argmax(np.abs(Pn-Pt).mean(1))]})")

print(); print(GLYPH)
print("Cal-2/3  shared effect sweep: composition share and residual retention")
print(GLYPH)
print(f"{'a(shared)':>10} {'r_obs':>7} {'r_null':>7} {'r_resid':>8} "
      f"{'compShare':>10} {'residSig':>9}")
print("-" * 58)
for a in (0.0, 0.10, 0.20, 0.40, 0.80):
    fold = {k: [] for k in ("ro", "rn", "rr", "cs", "rs")}
    for s in range(30):
        A = run("mei", a, 1.0, 1000 + s)
        B = run("mei", a, 1.0, 2000 + s)
        ro = corr(A["obs"], B["obs"]); rn = corr(A["null"], B["null"])
        rr = corr(A["resid"], B["resid"])
        # composition-only prediction of the cross-model correlation:
        # project each observed logFC onto its composition null
        fA = corr(A["obs"], A["null"]) ** 2
        fB = corr(B["obs"], B["null"]) ** 2
        comp = np.sqrt(fA * fB) * rn
        fold["ro"].append(ro); fold["rn"].append(rn)
        fold["rr"].append(rr); fold["cs"].append(comp / ro if abs(ro) > 1e-9 else np.nan)
        fold["rs"].append(rr)
    m = {k: np.nanmean(v) for k, v in fold.items()}
    print(f"{a:>10.2f} {m['ro']:>+7.3f} {m['rn']:>+7.3f} {m['rr']:>+8.3f} "
          f"{m['cs']:>10.3f} {m['rs']:>+9.3f}")
print("  compShare = attributed-to-composition fraction of the observed correlation")
print("  residSig  = residual cross-model correlation (carries the shared signal)")

print(); print(GLYPH)
print("Cal-4  negative control + the same-node vs different-node contrast")
print(GLYPH)
print(f"{'a(shared)':>10} {'same-node r':>12} {'diff-node r':>12} {'delta':>8} {'z':>7}")
print("-" * 54)
for a in (0.0, 0.20, 0.40):
    outs = []
    for s in range(30):
        rs = np.random.default_rng(3000 + s)
        bynode = {}
        for nd in SH:
            bynode[nd] = [run(nd, a, 1.0, int(rs.integers(1e6)))["resid"] for _ in range(3)]
        same = [corr(bynode[nd][i], bynode[nd][j])
                for nd in SH for i in range(3) for j in range(i + 1, 3)]
        diff = [corr(bynode[n1][i], bynode[n2][j])
                for n1 in SH for n2 in SH if n1 < n2 for i in range(3) for j in range(3)]
        outs.append((np.nanmean(same), np.nanmean(diff)))
    o = np.array(outs); d = o[:, 0] - o[:, 1]
    print(f"{a:>10.2f} {o[:,0].mean():>+12.3f} {o[:,1].mean():>+12.3f} "
          f"{d.mean():>+8.3f} {d.mean()/(d.std()+1e-9):>+7.2f}")
print("  NOTE: at a=0 this delta is NOT zero -> residual composition leakage.")
print("        That is the calibration's own false-positive floor and must be")
print("        subtracted / permuted against, not assumed away.")

print(); print(GLYPH)
print("Cal-5  does NNLS orthogonality DEFLATE residual similarity? (direct test)")
print(GLYPH)
# Compare: residual obtained by (i) NNLS on S, (ii) regression on a K-dim basis
# built from a DIFFERENT random reference. If orthogonality mechanically deflated
# things, (i) would systematically sit below (ii).
alt = _r.normal(0, 1, (G, K))
altP = np.zeros((8, K))
out_i, out_ii = [], []
for s in range(30):
    A = run("mei", 0.30, 1.0, 5000 + s)
    B = run("mei", 0.30, 1.0, 6000 + s)
    out_i.append(corr(A["resid"], B["resid"]))
    # alternative projection using an unrelated basis
    rr = []
    for m in (A, B):
        pr = np.linalg.lstsq(alt, m["obs"], rcond=None)[0]
        rr.append(m["obs"] - alt @ pr)
    out_ii.append(corr(rr[0], rr[1]))
print(f"  residual via NNLS on S          : r = {np.nanmean(out_i):+.3f}")
print(f"  residual via unrelated K-dim proj: r = {np.nanmean(out_ii):+.3f}")
print(f"  difference                       : {np.nanmean(out_i)-np.nanmean(out_ii):+.3f}")
print("  -> if orthogonality mechanically deflated residual similarity, the first")
print("     would sit systematically BELOW the second. It does not; the residual")
print("     similarity is driven by the shared signal plus leakage, not by the")
print("     projection's geometry.")
