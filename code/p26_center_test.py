"""
p26_center_test.py -- T=8 coded, single SNR point.  Is the residual soft gap of eq. (26) a
list-size effect or a centering effect?  Receivers:
  exact          full posterior
  P26-D1 (eta)   paper: lists centered on the greedy decision
  P26-D3 (eta)   same lists, centered on the D3 (coordinate-descent, near-ML) decision
  P26-cells(eta) union of the lists of all T forced-cell candidates (T x larger list)
Usage: python3 p26_center_test.py rho_db ncw [eta]
"""
import numpy as np, time, sys, os
from cubesplit_gap import *

rho_db = float(sys.argv[1]); NCW = int(sys.argv[2]); ETA = int(sys.argv[3]) if len(sys.argv) > 3 else 16
T, B0, N, NBLK = 8, 1, 4, 60
rng = np.random.default_rng(5)
cs = CubeSplit(T, B0); C, bits = cs.codebook(); Phi = codebook_features(C)
n = NBLK * cs.B
ld = LDPC(n, 3, 6, rng); perm = rng.permutation(n); iperm = np.argsort(perm)
fn = f"lists_T{T}_B0{B0}_eta{ETA}.npy"
if os.path.exists(fn):
    lists = np.load(fn)
else:
    t0 = time.time(); lists, nf = precompute_lists_sorted(cs, C, bits, ETA, log=lambda m: print(m, flush=True)); np.save(fn, lists)
    print(f"lists: {time.time()-t0:.0f}s fallbacks {nf}", flush=True)
rho = 10 ** (rho_db / 10); kap = kappa_of(rho, T)


def p26_from_index(k, Y, lists_k):
    """LLRs from the lists of the codewords in lists_k (n,B,2,eta) ; lse over the lists."""
    n_ = len(k); B, em = lists_k.shape[1], lists_k.shape[3]
    cand = C[:, lists_k.reshape(n_, -1)]
    z = np.einsum('ntm,tnc->nmc', Y.conj(), cand)
    a = kap * (z.real ** 2 + z.imag ** 2).sum(1).reshape(n_, B, 2, em)
    return np.clip(lse(a[:, :, 0, :], 2) - lse(a[:, :, 1, :], 2), -30, 30)


RX = ['exact', 'P26-D1', 'P26-D3', 'P26-cells']
err = {r: 0 for r in RX}; ncw = 0; mass = {'D1': [], 'D3': [], 'cells': []}
batch = 50
t0 = time.time()
while ncw < NCW:
    u = rng.integers(2, size=(batch, ld.k)).astype(np.uint8); x = ld.encode(u)
    xb = x[:, perm].reshape(-1, cs.B); cell, dig = cs.from_bits(xb); k = cs.index(cell, dig)
    X = cs.symbols(cell, dig); Y = channel(rng, X, N, rho); G = gram(Y)
    m = ml_metrics_fast(G, Phi).astype(np.float64)
    L = {'exact': exact_llr(m, bits, kap)}
    p = np.exp(kap * (m - m.max(1, keepdims=True))); p /= p.sum(1, keepdims=True)
    c1, d1 = decode_D1(cs, G); k1 = cs.index(c1, d1)
    c2, d2, mus, digs = decode_D2(cs, G, Y); c3, d3, _, _ = refine_D3(cs, G, c2, d2); k3 = cs.index(c3, d3)
    nb = len(k)
    for name, kk in (('D1', k1), ('D3', k3)):
        lk = lists[kk]
        L['P26-' + name] = p26_from_index(kk, Y, lk)
        # posterior mass captured by the candidate set (diagnostic)
        idx = lk.reshape(nb, -1)
        mass[name].append(np.mean([p[i, np.unique(idx[i])].sum() for i in range(nb)]))
    # union over the T forced-cell candidates: concatenate their lists along eta axis
    kc = np.stack([cs.index(np.full(nb, i), digs[:, i]) for i in range(T)], 1)      # (nb,T)
    lk = np.concatenate([lists[kc[:, i]] for i in range(T)], axis=3)                 # (nb,B,2,T*eta)
    L['P26-cells'] = p26_from_index(k, Y, lk)
    idx = lk.reshape(nb, -1)
    mass['cells'].append(np.mean([p[i, np.unique(idx[i])].sum() for i in range(nb)]))
    for r in RX:
        xh = ld.decode(L[r].reshape(batch, -1)[:, iperm], 50)
        err[r] += int(np.any(xh != x, 1).sum())
    ncw += batch
    print(f"  {ncw} cw: " + "  ".join(f"{r} {err[r]/ncw:.4f}" for r in RX) + f"   ({time.time()-t0:.0f}s)", flush=True)
print(f"T={T} B0={B0} N={N} rho={rho_db} dB eta={ETA}  ncw={ncw}")
print("BLER: " + "  ".join(f"{r} {err[r]/ncw:.4f}" for r in RX))
print("mean posterior mass inside candidate set: " + "  ".join(f"{a} {np.mean(v):.3f}" for a, v in mass.items()))
print(f"candidate set sizes: D1/D3 {cs.B*2*ETA}, cells {T*cs.B*2*ETA}")
