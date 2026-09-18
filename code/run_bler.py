"""
run_bler.py -- coded BLER: rate-1/2 (3,6)-regular LDPC (Gallager, sum-product 50 it), random
bit interleaver, exact posterior LLR (full codebook) vs the paper's approximate LLR (26)
(eta nearest per bit, precomputed; only when K <= 2^14) vs D3s (structured on-the-fly list).
Usage: python3 run_bler.py T B0 N nblk_per_codeword [ncw_min] [ncw_max] [seed] [eta]
"""
import numpy as np, time, json, sys
from cubesplit_gap import *

T, B0, N, NBLK = int(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3]), int(sys.argv[4])
NCW_MIN = int(sys.argv[5]) if len(sys.argv) > 5 else 400
NCW_MAX = int(sys.argv[6]) if len(sys.argv) > 6 else 1500
seed = int(sys.argv[7]) if len(sys.argv) > 7 else 5
ETAS = [int(v) for v in sys.argv[8].split(',')] if len(sys.argv) > 8 else [4, 8, 16]
rng = np.random.default_rng(seed)

cs = CubeSplit(T, B0)
C, bits = cs.codebook()
Phi = codebook_features(C)
n = NBLK * cs.B
assert n % 6 == 0, "n must be divisible by d_c = 6"
ld = LDPC(n, 3, 6, rng)
perm = rng.permutation(n)
iperm = np.argsort(perm)
RX = ['exact'] + [f'P26e{e}' for e in ETAS]
t0 = time.time()
lists, nfall = precompute_lists_sorted(cs, C, bits, max(ETAS), log=lambda m: print(m, flush=True))
print(f"P26 lists eta_max={max(ETAS)}: {time.time()-t0:.0f}s, fallback rows {nfall}", flush=True)
print(f"=== CS({T},{B0}) N={N}  B={cs.B} K={cs.K}  LDPC n={n} k={ld.k} rate={ld.k/n:.3f}  blocks/cw={NBLK}  receivers={RX}  seed={seed}", flush=True)


def run_point(rho_db, ncw_min, ncw_max, batch=100, min_err=30):
    rho = 10 ** (rho_db / 10)
    kap = kappa_of(rho, T)
    err = {r: 0 for r in RX}
    hard_err = {'D1': 0, 'D3': 0}     # uncoded symbol errors (diagnostic)
    ncw = 0
    nblocks = 0
    while ncw < ncw_max and (ncw < ncw_min or min(err.values()) < min_err):
        u = rng.integers(2, size=(batch, ld.k)).astype(np.uint8)
        x = ld.encode(u)
        xb = x[:, perm].reshape(-1, cs.B)                       # (batch*NBLK, B)
        cell, dig = cs.from_bits(xb)
        k = cs.index(cell, dig)
        X = cs.symbols(cell, dig)
        Y = channel(rng, X, N, rho)
        G = gram(Y)
        L = {}
        L['exact'] = exact_llr_chunked(G, Phi, bits, kap)
        Lp, c1, d1 = soft_P26_multi(cs, G, Y, C, lists, ETAS, kap)
        for e in ETAS:
            L[f'P26e{e}'] = Lp[e]
        hard_err['D1'] += int(np.sum(cs.index(c1, d1) != k))
        hard_err['D3'] += int(np.sum(ml_metrics_fast(G, Phi).argmax(1) != k))   # ML symbol errors (diagnostic)
        for r in RX:
            Lc = L[r].reshape(batch, -1)[:, iperm]              # de-interleave
            xh = ld.decode(Lc, 50)
            err[r] += int(np.any(xh != x, axis=1).sum())
        ncw += batch
        nblocks += batch * NBLK
    return {r: err[r] / ncw for r in RX}, ncw, {a: v / nblocks for a, v in hard_err.items()}


# ---- probe (or explicit grid lo:hi:step as argv[9])
probe = {}
if len(sys.argv) > 9:
    lo, hi, step = [float(v) for v in sys.argv[9].split(':')]
else:
    t0 = time.time()
    for s in range(-4, 40, 2):
        b, _, _ = run_point(s, 40, 40, batch=40, min_err=0)
        probe[s] = b
        if max(b.values()) < 0.05:
            break
    ps = sorted(probe)
    lo = max([s for s in ps if probe[s]['exact'] >= 0.9] or [ps[0]])
    hi = ps[-1]
    print(f"probe ({time.time()-t0:.0f}s): " + "  ".join(f"{s}:" + "/".join(f"{probe[s][r]:.2f}" for r in RX) for s in ps), flush=True)
    step = 0.5
snrs = [round(lo + step * i, 2) for i in range(int(round((hi - lo) / step)) + 1)]
print(f"fine grid {snrs[0]}..{snrs[-1]} dB step {step}", flush=True)

res = {}
print(f"{'rho dB':>7} " + " ".join(f"{'BLER ' + r:>11}" for r in RX) + f" {'ncw':>5}  {'SER D1':>8} {'SER ML':>8}  time")
for s in snrs:
    t0 = time.time()
    b, ncw, he = run_point(s, NCW_MIN, NCW_MAX)
    res[s] = {'bler': b, 'ncw': ncw, 'ser_hard': he}
    print(f"{s:>7} " + " ".join(f"{b[r]:>11.4f}" for r in RX) + f" {ncw:>5}  {he['D1']:>8.5f} {he['D3']:>8.5f}  {time.time()-t0:5.0f}s", flush=True)
    if max(b.values()) < 3e-3:
        break
snrs = [s for s in snrs if s in res]

print("\nSNR (dB, per channel use) at target BLER, log-linear interpolation; gap = receiver - exact")
for tgt in (1e-1, 1e-2):
    ex = snr_at(tgt, snrs, [res[s]['bler']['exact'] for s in snrs])
    row = []
    for r in RX:
        v = snr_at(tgt, snrs, [res[s]['bler'][r] for s in snrs])
        row.append(f"{r}={v:6.2f} (+{v-ex:5.2f})" if (v is not None and ex is not None) else f"{r}=  n/a")
    print(f"  BLER={tgt:g}: " + "  ".join(row))

json.dump({'T': T, 'B0': B0, 'N': N, 'B': cs.B, 'K': cs.K, 'n': n, 'k': ld.k, 'nblk': NBLK, 'etas': ETAS, 'seed': seed,
           'snrs': snrs, 'res': {str(s): res[s] for s in snrs}, 'probe': {str(s): probe[s] for s in probe}},
          open(f"bler_T{T}_B0{B0}_N{N}_eta{max(ETAS)}.json", "w"), indent=1)
