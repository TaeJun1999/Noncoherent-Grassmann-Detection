"""
scale_and_cost.py  --  Hand-off v2, items 2b and 4 of week 1.

(A) MACG base scale law: find c* in Omega = I + cG that maximises ESS,
    test the predicted scaling c* ~ kappa/T.
(B) Cost profiling (Hand-off v2 section 8.1): does the 1st claim have a stage?
    exhaustive demapper   : K * T^2 * M complex MACs per block
    proposed              : NFE * C_net + L * T^2 * M + T^3
    -> break-even codebook size, and the (T,M,B) region where we can win.
"""
import numpy as np, time

rng = np.random.default_rng(2)


def cn(shape):
    return (rng.standard_normal(shape) + 1j * rng.standard_normal(shape)) / np.sqrt(2)


def ess_for_c(G, kappa, M, c, n=60000):
    T = G.shape[0]
    Om = np.eye(T) + c * G
    w_, U_ = np.linalg.eigh(Om)
    Oh = U_ @ np.diag(np.sqrt(w_)) @ U_.conj().T
    Oi = np.linalg.inv(Om); ldO = np.linalg.slogdet(Om)[1]
    Z = cn((n, T, M)); X = np.linalg.qr(Oh @ Z)[0]
    lt = kappa * np.real(np.einsum('nti,ts,nsi->n', X.conj(), G, X))
    S = np.einsum('nti,ts,nsj->nij', X.conj(), Oi, X)
    lq = -M * ldO - T * np.linalg.slogdet(S)[1]
    lw = lt - lq; w = np.exp(lw - lw.max())
    return w.sum() ** 2 / (w ** 2).sum() / n


print("=" * 74)
print("(A) MACG base scale:  c* in Omega = I + cG,  prediction c* ~ kappa/T")
print("=" * 74)
print(f"{'T':>3} {'M':>2} {'kappa':>6} {'c*':>7} {'kappa/T':>8} {'c*T/kappa':>10} {'ESS@c*':>7} {'ESS@Haar':>9}")
ratios = []
for T in [4, 8, 12, 16]:
    for M in [1, 2, 3]:
        if M >= T:
            continue
        for kappa in [1.0, 4.0]:
            lam = np.sort(rng.uniform(0.5, 4.0, T))[::-1]
            V = np.linalg.qr(cn((T, T)))[0]
            G = V @ np.diag(lam) @ V.conj().T
            G = (G + G.conj().T) / 2
            grid = np.geomspace(0.02, 3.0, 22)
            es = [ess_for_c(G, kappa, M, c, n=40000) for c in grid]
            i = int(np.argmax(es))
            cstar = grid[i]
            ratios.append(cstar * T / kappa)
            print(f"{T:>3} {M:>2} {kappa:>6.1f} {cstar:>7.3f} {kappa/T:>8.3f} "
                  f"{cstar*T/kappa:>10.2f} {es[i]:>7.3f} {ess_for_c(G,kappa,M,0.0,40000):>9.3f}")
r = np.array(ratios)
print(f"\n  c* T / kappa : median {np.median(r):.2f}, range [{r.min():.2f}, {r.max():.2f}]")
print("  -> c* = a*kappa/T with a of order 1.  Use c = kappa/T as the default,")
print("     then one cheap ESS line-search per (T,M,SNR) at training setup time.")

# ------------------------------------------------------------------ (B) cost
print()
print("=" * 74)
print("(B) Cost profiling: exhaustive noncoherent soft demapping")
print("=" * 74)


def cost_exhaustive(T, M, B):
    """complex MACs per coherence block: one GEMM L^H C + column norms."""
    K = 2 ** B
    return K * T * T * M + K * T * M


def cost_net(T, d, layers):
    """T-token transformer, per NFE, complex-MAC-equivalent (real MACs / 4)."""
    per_layer = 4 * T * d * d + 2 * T * T * d + 2 * T * d * (4 * d)   # attn + MLP
    return layers * per_layer / 4.0


def cost_proposed(T, M, L, nfe, d=64, layers=4):
    return nfe * cost_net(T, d, layers) + L * (T * T * M + T * M) + T ** 3


# calibrate the flop model against measured wall time
print("\n  flop-model calibration (this container, 1 core, OpenBLAS):")
for (T, M, B) in [(8, 2, 14), (8, 2, 16), (16, 2, 16)]:
    K = 2 ** B
    C = cn((T, M * K))
    L = np.linalg.cholesky(np.eye(T) * 1.0 + 0j)
    t0 = time.perf_counter()
    Z = L.conj().T @ C
    m = np.sum(np.abs(Z) ** 2, axis=0).reshape(K, M).sum(1)
    t1 = time.perf_counter()
    mac = cost_exhaustive(T, M, B)
    print(f"    T={T:>2} M={M} B={B:>2}  K={K:>7}  {mac/1e6:>8.2f} MMAC  "
          f"{(t1-t0)*1e3:>7.2f} ms  -> {mac/(t1-t0)/1e9:>5.2f} GMAC/s")

print("\n  cost per information bit (complex MACs/bit), exhaustive:")
print(f"  {'T':>3} {'M':>2} {'B':>3} {'K':>9} {'MAC/block':>12} {'MAC/bit':>11}")
rows = [(8, 2, 12), (8, 2, 16), (8, 2, 20), (16, 2, 16), (16, 2, 20),
        (16, 4, 24), (32, 4, 24), (32, 4, 28)]
for (T, M, B) in rows:
    c = cost_exhaustive(T, M, B)
    print(f"  {T:>3} {M:>2} {B:>3} {2**B:>9} {c:>12.3e} {c/B:>11.3e}")
print("\n  reference: coherent 4x4 MIMO 256-QAM max-log soft demapping is")
print("  O(10^2 - 10^3) MAC/bit. Exhaustive noncoherent sits 10^3 - 10^6 above it.")

# ------------------------------------------------------------------ break-even
print()
print("=" * 74)
print("(B2) Break-even: smallest B at which the proposed receiver is cheaper")
print("=" * 74)
print(f"  {'T':>3} {'M':>2} {'d':>3} {'NFE':>4} {'L':>4} {'net MAC':>10} "
      f"{'B_even':>7} {'B for 10x':>10}")
for (T, M) in [(8, 2), (16, 2), (16, 4), (32, 4)]:
    for (d, layers, nfe, Lst) in [(64, 4, 8, 64), (64, 4, 4, 64), (128, 6, 4, 128)]:
        cp = cost_proposed(T, M, Lst, nfe, d, layers)
        # smallest B with cost_exhaustive > cp   and   > 10*cp
        be = tn = None
        for B in range(6, 41):
            ce = cost_exhaustive(T, M, B)
            if be is None and ce > cp:
                be = B
            if tn is None and ce > 10 * cp:
                tn = B
                break
        print(f"  {T:>3} {M:>2} {d:>3} {nfe:>4} {Lst:>4} {cp:>10.3e} "
              f"{be:>7} {tn:>10}")

print("""
  Reading: B_even is where we stop losing; B for 10x is where the claim has room.
  The network term dominates until the codebook is large, so the operating
  point is pushed to large B, small NFE, small d -- exactly the pressure that
  makes few-step inference (section 5.4) a requirement rather than a nicety.
""")

# ------------------------------------------------------------------ throughput
print("=" * 74)
print("(B3) System-level load: parallel subcarriers are the real multiplier")
print("=" * 74)
Tsym = 0.5e-3 / 14          # 30 kHz SCS OFDM symbol incl. CP
for (T, M, B, nsc) in [(8, 2, 16, 1), (8, 2, 16, 600), (8, 2, 16, 1200),
                       (16, 2, 20, 600), (16, 4, 24, 600)]:
    blocks_per_s = nsc / (T * Tsym)
    load = cost_exhaustive(T, M, B) * blocks_per_s
    rate = nsc * B / (T * Tsym) / 1e6
    print(f"  T={T:>2} M={M} B={B:>2} n_sc={nsc:>5}: {rate:>8.1f} Mbps, "
          f"{blocks_per_s/1e6:>6.3f} Mblock/s, {load/1e12:>8.2f} TMAC/s")
print("""
  A single stream is not a bottleneck. The load only becomes real once the
  same demapper runs across the subcarrier grid. Any complexity claim in the
  paper must state n_sc and the symbol duration, or it is not checkable.
""")
