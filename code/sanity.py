import numpy as np, time
from cubesplit_gap import *

rng = np.random.default_rng(0)

# ---------- (1) constellation properties + Lemma 2 (paper) for CS(T,1)
for T in (2, 4):
    cs = CubeSplit(T, 1)
    C, bits = cs.codebook(); C = C.T
    cell, dig = cs.from_index(np.arange(cs.K))
    assert np.allclose(np.sum(np.abs(C) ** 2, 1), 1.0)
    assert np.all(np.argmax(np.abs(C), 1) == cell)
    corr = np.abs(C.conj() @ C.T) ** 2
    np.fill_diagonal(corr, 0.0)
    dmin = np.sqrt(1 - corr.max())
    dnn = np.sqrt(1 - corr.max(1))
    m = probit(0.75)
    c = (1 - np.exp(-m ** 2)) / (1 + np.exp(-m ** 2))
    dlem = np.sqrt(1 - abs(1 - (1 + 1j) / (1 / c + T - 1)) ** 2)
    print(f"CS({T},1): K={cs.K} B={cs.B}  dmin computed {dmin:.6f}  Lemma-2 {dlem:.6f}  "
          f"NN-distance spread {dnn.max()-dnn.min():.2e}")

# ---------- (2) bit round trip
cs = CubeSplit(4, 3)
k = rng.integers(cs.K, size=1000)
cell, dig = cs.from_index(k)
b = cs.bits(cell, dig)
c2, d2 = cs.from_bits(b)
assert np.all(c2 == cell) and np.all(d2 == dig)
assert np.all(cs.index(cell, dig) == k)
print("bit round trip ok")

# ---------- (3) noiseless recovery by every decoder
for (T, B0) in [(4, 1), (4, 2), (4, 3), (8, 1)]:
    cs = CubeSplit(T, B0)
    n = 2000
    k = rng.integers(cs.K, size=n)
    cell, dig = cs.from_index(k)
    X = cs.symbols(cell, dig)
    N = 4
    H = (rng.standard_normal((n, N)) + 1j * rng.standard_normal((n, N))) / np.sqrt(2)
    Y = X[:, :, None] * H[:, None, :]
    G = gram(Y)
    ok = {}
    c1, d1 = decode_D1(cs, G); ok['D1'] = np.mean(cs.index(c1, d1) == k)
    c1, d1 = decode_D1r(cs, G); ok['D1r'] = np.mean(cs.index(c1, d1) == k)
    c2_, d2_, mus, digs = decode_D2(cs, G, Y); ok['D2'] = np.mean(cs.index(c2_, d2_) == k)
    c3, d3, mf, gm = refine_D3(cs, G, c2_, d2_); ok['D3'] = np.mean(cs.index(c3, d3) == k)
    assert np.allclose(mf, metric_of(Y, cs.symbols(c3, d3)), rtol=1e-8)
    print(f"noiseless CS({T},{B0}) K={cs.K}: " + "  ".join(f"{a}={v:.4f}" for a, v in ok.items()))

# ---------- (4) exhaustive ML: float32 codebook vs float64 argmax agreement, and timing
cs = CubeSplit(4, 3)
t0 = time.time(); C, bits = cs.codebook(); print(f"codebook CS(4,3) K={cs.K} built in {time.time()-t0:.1f}s")
n = 64
k = rng.integers(cs.K, size=n)
cell, dig = cs.from_index(k)
X = cs.symbols(cell, dig)
Y = channel(rng, X, 4, 10 ** (18 / 10))
t0 = time.time(); m64 = ml_metrics(Y, C); t64 = time.time() - t0
t0 = time.time(); m32 = ml_metrics(Y, C.astype(np.complex64)); t32 = time.time() - t0
print(f"ML metrics: argmax agree {np.mean(m64.argmax(1) == m32.argmax(1)):.4f}, "
      f"max rel diff {np.max(np.abs(m64-m32)/np.abs(m64).max()):.2e}; time/block f64 {t64/n*1e3:.1f} ms, f32 {t32/n*1e3:.1f} ms")
# refined D3 metric must never exceed exhaustive max (sanity) and D3 hit-rate
G = gram(Y)
c2_, d2_, mus, digs = decode_D2(cs, G, Y)
c3, d3, mf, gm = refine_D3(cs, G, c2_, d2_)
print(f"  D3 metric <= ML max: {np.all(mf <= m64.max(1) + 1e-9)}; D3==ML argmax {np.mean(cs.index(c3,d3)==m64.argmax(1)):.3f}")

# ---------- (5) soft outputs: exact vs D3s vs P26 at moderate SNR (T=4,B0=2)
cs = CubeSplit(4, 2)
C, bits = cs.codebook()
n = 400
k = rng.integers(cs.K, size=n)
cell, dig = cs.from_index(k)
X = cs.symbols(cell, dig)
rho = 10 ** (14 / 10)
Y = channel(rng, X, 4, rho)
kap = kappa_of(rho, 4)
Lex = exact_llr(ml_metrics(Y, C), bits, kap)
G = gram(Y)
L3, _, _ = soft_D3(cs, G, Y, kap)
t0 = time.time(); lists = precompute_lists(cs, C, bits, 8); print(f"P26 lists (eta=8) for K={cs.K}: {time.time()-t0:.1f}s")
Lp, _, _ = soft_P26(cs, G, Y, C, lists, kap)
tb = bits[k].astype(float)
sign_ok = lambda L: np.mean((L > 0) == (tb == 0))
print(f"  LLR sign agreement with sent bits: exact {sign_ok(Lex):.4f}  D3s {sign_ok(L3):.4f}  P26 {sign_ok(Lp):.4f}")
print(f"  corr(exact, D3s) {np.corrcoef(Lex.ravel(), L3.ravel())[0,1]:.4f}   corr(exact, P26) {np.corrcoef(Lex.ravel(), Lp.ravel())[0,1]:.4f}")

# ---------- (6) LDPC
ld = LDPC(1008, 3, 6, rng)
print(f"LDPC n={ld.n} m={ld.m} rank={len(ld.piv)} k={ld.k} rate={ld.k/ld.n:.3f}")
u = rng.integers(2, size=(50, ld.k)).astype(np.uint8)
x = ld.encode(u)
assert not ((x.astype(int) @ ld.H.T.astype(int)) % 2).any()
for ebn0 in (1.5, 2.5, 3.5):
    s2 = 1 / (2 * (ld.k / ld.n) * 10 ** (ebn0 / 10))
    y = (1 - 2 * x.astype(float)) + np.sqrt(s2) * rng.standard_normal(x.shape)
    L = 2 * y / s2
    t0 = time.time(); xh = ld.decode(L, 50); dt = time.time() - t0
    print(f"  BPSK/AWGN Eb/N0={ebn0} dB: BLER {np.mean((xh != x).any(1)):.3f}  ({dt:.1f}s for 50 cw)")
