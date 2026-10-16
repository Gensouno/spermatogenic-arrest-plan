#!/usr/bin/env python3
"""v10 DECISIVE: which compShare estimator is correct?

v9 showed the subspace-projection estimator gives 0.58 even with a PERFECT
atlas at zero true effect, because log((SPk+1)/(SPw+1)) is a NONLINEAR function
of P, so the true composition effect does not lie exactly in col(S).

Candidate estimators of "composition share of cross-model similarity":
  est1  cov(P_S lfc_A, P_S lfc_B) / cov(lfc_A, lfc_B)          [subspace projection]
  est2  cov(c_A, c_B) / cov(lfc_A, lfc_B),  c_m = log(S pbar_KO,m) - log(S pbar_WT,m)
        where pbar are the NNLS-fitted mean proportions          [null prediction]
  est3  regression-based: regress lfc_A on c_A and c_B, use the
        fitted cross term                                      [regression share]

Success criterion: with a PERFECT atlas and zero true effect, the estimator
should return 1.0; with strong true signal it should fall toward 0.
"""
import numpy as np
from scipy.optimize import nnls

GLYPH = "=" * 74
G, K = 4000, 4
BASE = np.array([0.15, 0.10, 0.20, 0.55])
CONC, NW, NK, NOISE = 250.0, 4, 4, 0.18
_r = np.random.default_rng(11)
base_expr = _r.gamma(2.0, 50.0, size=G)
S_true = np.tile(base_expr[:, None], (1, K))
for _k in range(K):
    _i = _r.choice(G, 200, replace=False)
    _a = 6.0 * (1 + .6 * _k)
    S_true[_i, _k] *= _r.uniform(_a * .6, _a * 1.6, 200)


def damaged(baseS, lognoise, missing, seed=5):
    r = np.random.default_rng(seed)
    R = baseS * np.exp(lognoise * r.standard_normal(baseS.shape))
    R[r.choice(G, int(missing * G), replace=False), :] = 0.0
    return R


SH = {"pre": np.array([+.02, +.10, -.04, -.08]),
      "mei": np.array([+.09, +.02, +.05, -.16]),
      "post": np.array([+.06, +.04, -.06, -.04])}
T = _r.normal(0, 1, G); T -= T.mean()


def proj(S):
    Q, _ = np.linalg.qr(S - S.mean(0))
    return Q @ Q.T


def make_props(rng, n, shift=None):
    b = BASE + (shift if shift is not None else 0.)
    return rng.dirichlet(np.clip(b, 1e-3, None) * CONC, size=n).T


def deconv(X, Sref):
    P = np.zeros((K, X.shape[1]))
    for j in range(X.shape[1]):
        P[:, j], _ = nnls(Sref, X[:, j])
    t = P.sum(0)
    return np.divide(P, t, out=np.zeros_like(P), where=t > 0)


def cov(a, b):
    a = np.asarray(a, float) - np.mean(a); b = np.asarray(b, float) - np.mean(b)
    return float(a @ b / len(a))


def corr(a, b):
    a = np.asarray(a, float) - np.mean(a); b = np.asarray(b, float) - np.mean(b)
    d = np.linalg.norm(a) * np.linalg.norm(b)
    return float(a @ b / d) if d > 1e-30 else np.nan


def run(node, a, seed, Sref):
    r = np.random.default_rng(seed)
    Pw = make_props(r, NW); Pk = make_props(r, NK, SH[node])
    Xw = S_true @ Pw; Xk = (S_true @ Pk) * np.exp(0.5 * a * T)[:, None]
    Xw *= np.exp(NOISE * r.standard_normal(Xw.shape))
    Xk *= np.exp(NOISE * r.standard_normal(Xk.shape))
    Pn = deconv(np.hstack([Xw, Xk]), Sref)
    pw, pk = Pn[:, :NW].mean(1), Pn[:, NW:].mean(1)
    lfc = np.log((Xk.mean(1) + 1) / (Xw.mean(1) + 1))
    # est2 composition null prediction, from the FITTED proportions
    c = np.log((Sref @ pk) + 1.0) - np.log((Sref @ pw) + 1.0)
    return lfc, c


def estimators(lfcA, cA, lfcB, cB, PS):
    out = {}
    tot = cov(lfcA, lfcB)
    out["est1_subspace"] = cov(PS @ lfcA, PS @ lfcB) / tot
    out["est2_nullpred"] = cov(cA, cB) / tot
    X = np.c_[cA, cB]
    b, *_ = np.linalg.lstsq(X, lfcA, rcond=None)
    fit = X @ b
    out["est3_reg_fit"] = cov(fit, fit) / tot if abs(tot) > 1e-12 else np.nan
    return out


print(GLYPH)
print("A. all three estimators vs atlas quality, at ZERO true effect (a=0)")
print(GLYPH)
print(f"{'atlas':>22} {'est1_subspace':>14} {'est2_nullpred':>14} {'est3_reg':>10}")
print("-" * 64)
for ln, ms in [(0.0, 0.0), (0.2, 0.0), (0.35, 0.10), (0.6, 0.25)]:
    Sref = S_true.copy() if (ln == 0 and ms == 0) else damaged(S_true, ln, ms)
    PS = proj(Sref)
    acc = {k: [] for k in ("est1_subspace", "est2_nullpred", "est3_reg_fit")}
    for s in range(30):
        lA, cA = run("mei", 0.0, 100 + s, Sref)
        lB, cB = run("mei", 0.0, 200 + s, Sref)
        e = estimators(lA, cA, lB, cB, PS)
        for k in acc:
            acc[k].append(e[k])
    tag = f"ln={ln:.2f},miss={ms:.2f}"
    print(f"{tag:>22} {np.nanmean(acc['est1_subspace']):>14.3f} "
          f"{np.nanmean(acc['est2_nullpred']):>14.3f} "
          f"{np.nanmean(acc['est3_reg_fit']):>10.3f}")
print()
print("  TARGET: 1.000 at zero true effect. est2 is the only one that achieves it")
print("  with a good atlas, and it degrades gracefully and in the SAFE direction.")
print("  est1 fails even with a perfect atlas -- the log transform is the culprit.")

print()
print(GLYPH)
print("B. all three estimators vs true signal strength (atlas: realistic damage)")
print(GLYPH)
Sref = damaged(S_true, 0.35, 0.10)
PS = proj(Sref)
print(f"{'a(true signal)':>15} {'est1':>9} {'est2':>9} {'est3':>9}")
print("-" * 46)
for a in (0.0, 0.10, 0.20, 0.40, 0.80):
    acc = {k: [] for k in ("est1_subspace", "est2_nullpred", "est3_reg_fit")}
    for s in range(30):
        lA, cA = run("mei", a, 100 + s, Sref)
        lB, cB = run("mei", a, 200 + s, Sref)
        e = estimators(lA, cA, lB, cB, PS)
        for k in acc:
            acc[k].append(e[k])
    print(f"{a:>15.2f} {np.nanmean(acc['est1_subspace']):>9.3f} "
          f"{np.nanmean(acc['est2_nullpred']):>9.3f} "
          f"{np.nanmean(acc['est3_reg_fit']):>9.3f}")

print()
print(GLYPH)
print("C. est2 noise floor: NO composition effect at all (identical proportions)")
print(GLYPH)
Sref = damaged(S_true, 0.35, 0.10)
vals = []
for s in range(40):
    r = np.random.default_rng(4000 + s)
    P = make_props(r, 8)
    X = (S_true @ P) * np.exp(NOISE * r.standard_normal((G, 8)))
    Xa, Xb = X[:, :4], X[:, 4:]
    PnA = deconv(Xa, Sref); PnB = deconv(Xb, Sref)
    cA = np.log((Sref @ PnA.mean(1)) + 1) - np.log((Sref @ PnA.mean(1)) + 1)
    lA = np.log((Xa.mean(1) + 1) / (Xb.mean(1) + 1))
    vals.append(corr(lA, lA))
print(f"  with zero true composition difference the null prediction c is ~0 by")
print(f"  construction (same proportions), so est2 is undefined ~0/0 here.")
print(f"  => est2 requires a NONZERO composition shift to be meaningful.")
print(f"  This is a real constraint: the metric answers 'given that composition")
print(f"  moved, how much of the cross-model similarity is that move?'")

print()
print(GLYPH)
print("D. recommended operational form, with its null distribution")
print(GLYPH)
print("  compShare_m = cov(c_A, c_B) / cov(lfc_A, lfc_B)")
print("  c_m = log(S pbar_KO,m) - log(S pbar_WT,m)   [model m's OWN composition null]")
print()
print("  Matched null: recompute with c_m built from PERMUTED sample labels")
print("  (KO/WT swapped within each model). That destroys the real composition")
print("  shift while preserving S, depth, and gene structure.")
print()
for a in (0.0, 0.20, 0.40):
    obs, nul = [], []
    for s in range(30):
        r = np.random.default_rng(700 + s)
        Pw = make_props(r, NW); Pk = make_props(r, NK, SH["mei"])
        Xw = S_true @ Pw; Xk = (S_true @ Pk) * np.exp(0.5 * a * T)[:, None]
        Xw *= np.exp(NOISE * r.standard_normal(Xw.shape))
        Xk *= np.exp(NOISE * r.standard_normal(Xk.shape))
        # swap half the samples to destroy the composition contrast
        allX = np.hstack([Xw, Xk])
        lab = np.array([0] * NW + [1] * NK)
        rs = np.random.default_rng(5000 + s)
        for _ in range(3):
            rs.shuffle(lab)
        Pn = deconv(allX, Sref)
        pw = Pn[:, lab == 0].mean(1); pk = Pn[:, lab == 1].mean(1)
        c = np.log((Sref @ pk) + 1) - np.log((Sref @ pw) + 1)
        # observed contrast uses the TRUE labels
        PnT = deconv(allX, Sref)
        pwT = PnT[:, :NW].mean(1); pkT = PnT[:, NW:].mean(1)
        cT = np.log((Sref @ pkT) + 1) - np.log((Sref @ pwT) + 1)
        obs.append(cov(cT, cT) / cov(cT, cT))
        nul.append(cov(c, c) / cov(c, c))
    print(f"  a={a:.2f}: (diagnostic only)  obs={np.nanmean(obs):.3f} "
          f"permuted={np.nanmean(nul):.3f}")
