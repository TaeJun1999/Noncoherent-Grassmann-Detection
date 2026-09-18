"""
run_ser.py -- uncoded symbol-error-rate sweep: exhaustive exact ML vs the Cube-split
efficient demappers on the same codebook.  Output: JSON + tables to stdout.
Usage: python3 run_ser.py T B0 N [nblocks] [seed]
"""
import numpy as np, time, json, sys
from cubesplit_gap import *

T, B0, N = int(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3])
NB = int(sys.argv[4]) if len(sys.argv) > 4 else 4000
seed = int(sys.argv[5]) if len(sys.argv) > 5 else 11
rng = np.random.default_rng(seed)
RX = ['ML', 'D1', 'D1r', 'D2', 'D3']

cs = CubeSplit(T, B0)
C, bits = cs.codebook()
Phi = codebook_features(C)
print(f"=== CS({T},{B0})  T={T} B0={B0} N={N}  B={cs.B} bits/block  K={cs.K}  blocks/point={NB}  seed={seed}", flush=True)


def run_point(rho_db, nblocks, batch=500, want=RX):
    rho = 10 ** (rho_db / 10)
    err = {r: 0 for r in want}
    berr = {r: 0 for r in want}
    cellerr = 0
    done = 0
    while done < nblocks:
        n = min(batch, nblocks - done)
        k = rng.integers(cs.K, size=n)
        cell, dig = cs.from_index(k)
        X = cs.symbols(cell, dig)
        Y = channel(rng, X, N, rho)
        G = gram(Y)
        tb = bits[k]
        dec = {}
        if 'ML' in want:
            dec['ML'] = ml_metrics_fast(G, Phi).argmax(1)
        if 'D1' in want:
            c1, d1 = decode_D1(cs, G)
            dec['D1'] = cs.index(c1, d1)
            cellerr += int(np.sum(c1 != cell))
        if 'D1r' in want:
            c1, d1 = decode_D1r(cs, G)
            dec['D1r'] = cs.index(c1, d1)
        if 'D2' in want or 'D3' in want:
            c2, d2, mus, digs = decode_D2(cs, G, Y)
            dec['D2'] = cs.index(c2, d2)
            if 'D3' in want:
                c3, d3, mf, gm = refine_D3(cs, G, c2, d2)
                dec['D3'] = cs.index(c3, d3)
        for r in want:
            wrong = dec[r] != k
            err[r] += int(wrong.sum())
            berr[r] += int(np.sum(bits[dec[r]] != tb))
        done += n
    out = {r: err[r] / nblocks for r in want}
    ber = {r: berr[r] / (nblocks * cs.B) for r in want}
    return out, ber, cellerr / nblocks


# ---- coarse probe to locate the useful SNR range
t0 = time.time()
probe = {}
for s in range(-4, 44, 4):
    ser, _, _ = run_point(s, 200, want=['ML', 'D1'])
    probe[s] = ser
    if ser['D1'] < 2e-3 and ser['ML'] < 2e-3:
        break
ps = sorted(probe)
lo = max([s for s in ps if probe[s]['ML'] >= 0.4] or [ps[0]])
hi = min([s for s in ps if probe[s]['D1'] <= 2e-3] or [ps[-1]])
print(f"probe ({time.time()-t0:.0f}s): " + "  ".join(f"{s}:{probe[s]['ML']:.3f}/{probe[s]['D1']:.3f}" for s in ps))
snrs = list(range(lo, hi + 1, 2))
print(f"fine grid {snrs[0]}..{snrs[-1]} dB step 2  (per-channel-use SNR rho; block SNR = rho*T = +{10*np.log10(T):.1f} dB)", flush=True)

res = {}
print(f"{'rho dB':>7} " + " ".join(f"{r:>9}" for r in RX) + f" {'cell-err':>9}  {'BER ML':>8} {'BER D1':>8}  time")
for s in snrs:
    t0 = time.time()
    ser, ber, ce = run_point(s, NB)
    res[s] = {'ser': ser, 'ber': ber, 'cellerr_D1': ce}
    print(f"{s:>7} " + " ".join(f"{ser[r]:>9.5f}" for r in RX) + f" {ce:>9.5f}  {ber['ML']:>8.5f} {ber['D1']:>8.5f}  {time.time()-t0:5.0f}s", flush=True)

print("\nSNR (dB, per channel use) at target SER, log-linear interpolation; gap = receiver - ML")
for tgt in (1e-1, 1e-2, 1e-3):
    row = []
    ml = snr_at(tgt, snrs, [res[s]['ser']['ML'] for s in snrs])
    for r in RX:
        v = snr_at(tgt, snrs, [res[s]['ser'][r] for s in snrs])
        row.append(f"{r}={v:6.2f} (+{v-ml:5.2f})" if (v is not None and ml is not None) else f"{r}=  n/a")
    print(f"  SER={tgt:g}: " + "  ".join(row))

# ---- per-block cost (wall clock, this container, 1 core) for reference
rho = 10 ** (snrs[len(snrs) // 2] / 10)
n = 256
k = rng.integers(cs.K, size=n); cell, dig = cs.from_index(k); X = cs.symbols(cell, dig)
Y = channel(rng, X, N, rho); G = gram(Y)
tm = {}
t0 = time.time(); ml_metrics_fast(G, Phi).argmax(1); tm['ML'] = time.time() - t0
t0 = time.time(); decode_D1(cs, G); tm['D1'] = time.time() - t0
t0 = time.time(); decode_D1r(cs, G); tm['D1r'] = time.time() - t0
t0 = time.time(); c2, d2, mus, digs = decode_D2(cs, G, Y); tm['D2'] = time.time() - t0
t0 = time.time(); refine_D3(cs, G, c2, d2); tm['D3'] = tm['D2'] + time.time() - t0
print("\nwall time per block (ms, batched numpy, 1 core; ML excludes G=YY^H):  " +
      "  ".join(f"{r}={tm[r]/n*1e3:.3f}" for r in RX))
c_d1 = T * T * N + 10 * T * T
macs = {'ML': cs.K * T * T, 'D1': c_d1, 'D2': c_d1 + T * (T * N + T), 'D3(1 sweep)': c_d1 + T * (T * N + T) + (T - 1) * (cs.Q ** 2 * 4 + T)}
print("complex-MAC-equivalent per block (analytic; ML = K*T^2 quadratic forms, D1 = Gram + power/eig ~10T^2):  " +
      "  ".join(f"{a}={v:.3g}" for a, v in macs.items()))

json.dump({'T': T, 'B0': B0, 'N': N, 'B': cs.B, 'K': cs.K, 'nblocks': NB, 'seed': seed,
           'snrs': snrs, 'res': {str(s): res[s] for s in snrs}, 'probe': {str(s): probe[s] for s in probe},
           'walltime_ms_per_block': {r: tm[r] / n * 1e3 for r in RX}},
          open(f"ser_T{T}_B0{B0}_N{N}.json", "w"), indent=1)
