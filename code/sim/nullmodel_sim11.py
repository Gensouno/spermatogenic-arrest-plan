#!/usr/bin/env python3
"""v11: settle the remaining design choices with n = 6-8 datasets (realistic).

1. Power of compShare with only 6-8 models.
2. Which permutation scheme gives a valid null.
3. delta-P: raw vs CLR vs ILR, with K=4.
4. Unknown-cell handling (merge vs drop) -- bias on the somatic estimate.
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
S_ref = S_true * np.exp(0.35 * _r.standard_normal(S_true.shape))
S_ref[_r.choice(G, int(.10 * G), replace=False), :] = 0.0

SH = {"pre": np.array([+.02, +.10, -.04, -.08]),
      "mei": np.array([+.09, +.02, +.05, -.16]),
      "post": np.array([+.06, +.04, -.06, -.04])}
T = _r.normal(0, 1, G); T -= T.mean()


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


def one(node, a, seed):
    r = np.random.default_rng(seed)
    Pw = make_props(r, NW); Pk = make_props(r, NK, SH[node])
    Xw = S_true @ Pw; Xk = (S_true @ Pk) * np.exp(0.5 * a * T)[:, None]
    Xw *= np.exp(NOISE * r.standard_normal(Xw.shape))
    Xk *= np.exp(NOISE * r.standard_normal(Xk.shape))
    Pn = deconv(np.hstack([Xw, Xk]), S_ref)
    pw, pk = Pn[:, :NW].mean(1), Pn[:, NW:].mean(1)
    lfc = np.log((Xk.mean(1) + 1) / (Xw.mean(1) + 1))
    c = np.log((S_ref @ pk) + 1) - np.log((S_ref @ pw) + 1)
    return lfc, c


def comp_share(lfcA, cA, lfcB, cB):
    return cov(cA, cB) / cov(lfcA, lfcB)


print(GLYPH); print("1. POWER with only n models (the real constraint)"); print(GLYPH)
# all pairs among n models -> many pairs, but they share models (not independent)
print(f"{'n models':>9} {'a=0':>16} {'a=0.20':>16} {'a=0.40':>16}")
print("-" * 60)
for nmod in (6, 8, 12):
    row = []
    for a in (0.0, 0.20, 0.40):
        vals = []
        for rep in range(20):
            models = []
            for i in range(nmod):
                node = list(SH)[i % 3]
                models.append(one(node, a, 3000 + rep * 100 + i))
            ps = [comp_share(models[i][0], models[i][1], models[j][0], models[j][1])
                  for i in range(nmod) for j in range(i + 1, nmod)]
            vals.append(np.nanmean(ps))
        row.append(f"{np.nanmean(vals):.2f}+/-{np.nanstd(vals):.2f}")
    print(f"{nmod:>9} {row[0]:>16} {row[1]:>16} {row[2]:>16}")
print()
print("  n=6: sd of the mean compShare is large -> only a coarse claim is supportable.")
print("  Report compShare as a RANGE with a permutation interval, never 2 decimals.")

print(); print(GLYPH); print("2. which permutation scheme?"); print(GLYPH)
a = 0.20
print(f"{'scheme':>34} {'null compShare':>16}")
print("-" * 52)
# scheme A: permute gene labels of one model's null prediction
# scheme B: permute model labels (pair observed with unrelated model's null)
# scheme C: permute KO/WT sample labels then refit
nulA, nulB = [], []
for rep in range(25):
    rng = np.random.default_rng(6000 + rep)
    ms = [one(list(SH)[i % 3], a, 7000 + rep * 50 + i) for i in range(6)]
    # B: pair model i's lfc with model j's null prediction (i != j) -- destroys
    #    the model-matched composition linkage while keeping both marginal dists
    for i in range(6):
        for j in range(6):
            if i == j:
                continue
            nulB.append(comp_share(ms[i][0], ms[j][1], ms[j][0], ms[i][1]))
    # A: shuffle genes of each c before pairing
    for i in range(6):
        for j in range(i + 1, 6):
            cA = rng.permutation(ms[i][1]); cB = rng.permutation(ms[j][1])
            nulA.append(comp_share(ms[i][0], cA, ms[j][0], cB))
print(f"  A gene-shuffle null      : compShare = {np.nanmean(nulA):+.3f} "
      f"+/- {np.nanstd(nulA):.3f}  (n={len(nulA)})")
print(f"  B model-mismatch null    : compShare = {np.nanmean(nulB):+.3f} "
      f"+/- {np.nanstd(nulB):.3f}  (n={len(nulB)})")
print()
print("  A is INVALID: shuffling genes destroys the marker structure that makes")
print("  c a valid composition prediction, so the null is far too wide.")
print("  B is the right null: it asks 'is the composition shift MATCHED to this")
print("  model specifically, or would any model's shift fit equally well?'")

print(); print(GLYPH); print("3. delta-P space: raw vs CLR vs ILR (K=4)"); print(GLYPH)


def clr(p):
    p = np.clip(np.asarray(p, float), 1e-9, None)
    return np.log(p) - np.log(p).mean()


def ilr_basis(K):
    B = np.zeros((K, K - 1))
    for j in range(K - 1):
        B[:j + 1, j] = 1.0 / (j + 1)
        B[j + 1, j] = -1.0
        B[:, j] *= np.sqrt((j + 1) / (j + 2))
    return B


B = ilr_basis(K)
raww, clrw, ilrw = [], [], []
for rep in range(300):
    rng = np.random.default_rng(8000 + rep)
    d1 = rng.dirichlet(np.ones(K)) - 1.0 / K
    d2 = rng.dirichlet(np.ones(K)) - 1.0 / K
    raww.append(np.corrcoef(d1, d2)[0, 1])
    clrw.append(np.corrcoef(clr(np.abs(d1) + 1e-6), clr(np.abs(d2) + 1e-6))[0, 1])
    ilrw.append(np.corrcoef(B.T @ d1, B.T @ d2)[0, 1])
print(f"  INDEPENDENT shifts, correlated across cell types:")
print(f"    raw delta-P (K=4 values) : mean r = {np.nanmean(raww):+.3f} sd {np.nanstd(raww):.3f}")
print(f"    CLR                      : mean r = {np.nanmean(clrw):+.3f} sd {np.nanstd(clrw):.3f}")
print(f"    ILR (K-1=3 coordinates)  : mean r = {np.nanmean(ilrw):+.3f} sd {np.nanstd(ilrw):.3f}")
print(f"  With K-1 = 3 df the sd is ~{np.nanstd(ilrw):.2f}; any single-pair delta-P")
print(f"  correlation is nearly uninformative. USE ILR and report a permutation")
print(f"  interval; never interpret one pairwise delta-P correlation.")

print(); print(GLYPH); print("4. Unknown cells: merge into Postmeiotic vs DROP"); print(GLYPH)
# Unknown = 2205 of 37241 cells = 5.9%. If Unknown is really late germ cells,
# dropping it removes them from BOTH the atlas and (implicitly) the fit.
frac_unknown = 2205 / 37241
print(f"  Unknown fraction of atlas = {frac_unknown:.3f}")
# simulate: true 5th compartment 'Unknown' that is actually late germ cells
S5 = S_true.copy()
S5[:, 3] = S5[:, 3] * np.exp(0.35 * _r.standard_normal(G))   # perturbed postmeio
S_atlas_drop = S_ref.copy()
S_atlas_merge = S_ref.copy(); S_atlas_merge[:, 3] = 0.5 * (S_ref[:, 3] + S5[:, 3])
for tag, Sa in [("drop", S_atlas_drop), ("merge", S_atlas_merge)]:
    mae = []
    for rep in range(20):
        r = np.random.default_rng(9000 + rep)
        P = make_props(r, 8)
        X = (S_true @ P) * np.exp(NOISE * r.standard_normal((G, 8)))
        Pn = deconv(X, Sa)
        mae.append(np.abs(Pn - P).mean())
    print(f"  {tag:6s}: mean |P_hat - P_true| = {np.nanmean(mae):.4f}")
print()
print("  Dropping Unknown is DEFENSIBLE only if Unknown is a technical artefact")
print("  (doublets/ambient). GSE112393 labels 2205 cells Unknown -- if they are")
print("  genuine late germ cells, dropping them biases the Postmeiotic proportion")
print("  DOWN and shifts the deficit onto Somatic. MUST be tested as a sensitivity")
print("  arm, not assumed.")
