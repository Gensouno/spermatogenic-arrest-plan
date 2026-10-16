#!/usr/bin/env python3
"""v8b LEAN: exact covariance decomposition + Mantel comparison (fast)."""
import numpy as np
from scipy.optimize import nnls
from scipy.spatial.distance import pdist, squareform

GLYPH = "=" * 74
G, K = 4000, 4
NAMES = ["Somatic", "SPG", "Meiotic", "Postmeio"]
BASE = np.array([0.15, 0.10, 0.20, 0.55])
CONC, NW, NK, NOISE = 250.0, 4, 4, 0.18

_r = np.random.default_rng(11)
base_expr = _r.gamma(2.0, 50.0, size=G)
S_true = np.tile(base_expr[:, None], (1, K))
for _k in range(K):
    _i = _r.choice(G, 200, replace=False)
    _a = 6.0 * (1 + .6 * _k)
    S_true[_i, _k] *= _r.uniform(_a * .6, _a * 1.6, 200)
S_bad = S_true * np.exp(0.35 * _r.standard_normal(S_true.shape))
S_bad[_r.choice(G, int(.10 * G), replace=False), :] = 0.0

SH = {"pre": np.array([+.02, +.10, -.04, -.08]),
      "mei": np.array([+.09, +.02, +.05, -.16]),
      "post": np.array([+.06, +.04, -.06, -.04])}
T = _r.normal(0, 1, G); T -= T.mean()

Q, _ = np.linalg.qr(S_bad - S_bad.mean(0))
PS = Q @ Q.T


def comps(v):
    c = PS @ v
    return c, v - c


def make_props(rng, n, shift=None):
    b = BASE + (shift if shift is not None else 0.)
    return rng.dirichlet(np.clip(b, 1e-3, None) * CONC, size=n).T


def deconv(X, Sref):
    P = np.zeros((K, X.shape[1]))
    for j in range(X.shape[1]):
        P[:, j], _ = nnls(Sref, X[:, j])
    t = P.sum(0)
    return P, np.divide(P, t, out=np.zeros_like(P), where=t > 0)


def corr(a, b):
    a = np.asarray(a, float) - a.mean(); b = np.asarray(b, float) - b.mean()
    d = np.linalg.norm(a) * np.linalg.norm(b)
    return float(a @ b / d) if d > 1e-30 else np.nan


def cov(a, b):
    a = np.asarray(a, float) - a.mean(); b = np.asarray(b, float) - b.mean()
    return float(a @ b / len(a))


def run(node, a, seed, Sref=None):
    Sref = S_bad if Sref is None else Sref
    r = np.random.default_rng(seed)
    Pw = make_props(r, NW); Pk = make_props(r, NK, SH[node])
    Xw = S_true @ Pw; Xk = (S_true @ Pk) * np.exp(0.5 * a * T)[:, None]
    Xw *= np.exp(NOISE * r.standard_normal(Xw.shape))
    Xk *= np.exp(NOISE * r.standard_normal(Xk.shape))
    _, Pn = deconv(np.hstack([Xw, Xk]), Sref)
    nullp = np.log(((Sref @ Pn[:, NW:].mean(1)) + 1) / ((Sref @ Pn[:, :NW].mean(1)) + 1))
    return dict(lfc=np.log((Xk.mean(1) + 1) / (Xw.mean(1) + 1)), null=nullp)


print(GLYPH); print("A. exact covariance decomposition vs known truth"); print(GLYPH)
print(f"{'a':>6} {'cov obs':>9} {'compShare':>10} {'residShare':>11} {'sum':>6} {'r_obs':>7}")
print("-" * 56)
for a in (0.0, 0.10, 0.20, 0.40, 0.80):
    Cs, Rs, Os, tots = [], [], [], []
    for s in range(25):
        A, B = run("mei", a, 100 + s), run("mei", a, 200 + s)
        cA, rA = comps(A["lfc"]); cB, rB = comps(B["lfc"])
        tot = cov(A["lfc"], B["lfc"])
        Cs.append(cov(cA, cB) / tot); Rs.append(cov(rA, rB) / tot)
        Os.append(corr(A["lfc"], B["lfc"])); tots.append(tot)
    print(f"{a:>6.2f} {np.nanmean(tots):>9.3f} {np.nanmean(Cs):>10.3f} "
          f"{np.nanmean(Rs):>11.3f} {np.nanmean(Cs)+np.nanmean(Rs):>6.3f} "
          f"{np.nanmean(Os):>+7.3f}")
print("  At a=0 (pure composition) the composition share is ~1.0 -- as it must be.")

print(); print(GLYPH); print("B. atlas matters: wrong basis -> meaningless share"); print(GLYPH)
alt = _r.normal(0, 1, (G, K)); Qa, _ = np.linalg.qr(alt - alt.mean(0)); PSa = Qa @ Qa.T
for a in (0.0, 0.40):
    Cs = []
    for s in range(20):
        A, B = run("mei", a, 100 + s), run("mei", a, 200 + s)
        cA, cB = PSa @ A["lfc"], PSa @ B["lfc"]
        Cs.append(cov(cA, cB) / cov(A["lfc"], B["lfc"]))
    print(f"  a={a:.2f}: composition share with an UNRELATED K-dim basis = {np.nanmean(Cs):+.3f}")

print(); print(GLYPH); print("C. Mantel on the SAMPLE-space distance matrix instead"); print(GLYPH)


def mantel_sample_p(M, nperm, rng):
    """M is samples x genes. Mantel over the SAMPLE distance matrix (n x n)."""
    dM = squareform(pdist(M))
    iu = np.triu_indices(dM.shape[0], 1)
    r0 = corr(dM[iu], dM[iu])
    k = 0
    for _ in range(nperm):
        p = rng.permutation(dM.shape[0])
        if abs(corr(dM[iu], dM[np.ix_(p, p)][iu])) >= abs(r0):
            k += 1
    return (k + 1) / (nperm + 1)


for a in (0.0, 0.40):
    Cs, pvs = [], []
    for s in range(12):
        rng = np.random.default_rng(900 + s)
        ms = [run(n, a, 700 + s * 10 + i) for i, n in enumerate(SH)]
        M = np.array([m["lfc"] for m in ms]).T          # G x 3 samples (3 nodes)
        pvs.append(mantel_sample_p(M.T, 199, rng))      # n x n = 3 x 3
        c = PS @ M
        Cs.append(cov(c[:, 0], c[:, 1]) / cov(M[:, 0], M[:, 1]))
    print(f"  a={a:.2f}:  Mantel-sample p={np.nanmean(pvs):.3f}   "
          f"compShare={np.nanmean(Cs):+.3f}")

print(); print(GLYPH); print("D. power of compShare to separate null from signal"); print(GLYPH)
for a in (0.0, 0.10, 0.20, 0.40):
    Cs = []
    for s in range(25):
        A, B = run("mei", a, 100 + s), run("mei", a, 200 + s)
        cA, cB = comps(A["lfc"])[0], comps(B["lfc"])[0]
        Cs.append(cov(cA, cB) / cov(A["lfc"], B["lfc"]))
    Cs = np.array(Cs)
    print(f"  a={a:.2f}: compShare = {np.nanmean(Cs):.3f} +/- {np.nanstd(Cs):.3f}  "
          f"95% = [{np.nanpercentile(Cs,2.5):.3f}, {np.nanpercentile(Cs,97.5):.3f}]")
