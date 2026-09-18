"""
gap_experiment.py  --  is there an accuracy gap for a cheap receiver to lose?

The complexity claim only has room if cheap detectors are actually worse.
Cheapest meaningful shortcut = throw away the spectrum and match subspaces
(exactly what Hand-off v2 section 1.2 says is lossy). Measure the cost.

  ML       : argmax_i  tr(X_i^H G X_i)                 [exact, uses spectrum]
  subspace : argmax_i  tr(X_i^H U U^H X_i)             [chordal to span(Y)]
  MMSE-ish : argmax_i  tr(X_i^H U diag(g) U^H X_i)     [span(Y) + flat weights]

Reports symbol error rate vs SNR and the SNR gap at a fixed SER.
"""
import numpy as np

rng = np.random.default_rng(7)


def cn(shape, var=1.0):
    return np.sqrt(var / 2) * (rng.standard_normal(shape) + 1j * rng.standard_normal(shape))


def haar(T, M, n=1):
    Q, _ = np.linalg.qr(cn((n, T, M)))
    return Q


T, M, N, B = 8, 2, 4, 10
K = 2 ** B
n_blocks = 4000
codebook = haar(T, M, K)                       # unstructured random codebook
Cflat = codebook.transpose(1, 0, 2).reshape(T, K * M)   # T x (K*M)

print(f"T={T} M={M} N={N} B={B} (K={K}), {n_blocks} blocks per SNR point")
print(f"{'SNR dB':>7} {'SER ML':>10} {'SER subsp':>11} {'SER spec-only':>14}")

snr_db_list = [0, 4, 8, 12, 16, 20]
res = {}
for snr_db in snr_db_list:
    rho, s2 = 10 ** (snr_db / 10), 1.0
    err_ml = err_sub = err_spec = 0
    for _ in range(n_blocks):
        i = rng.integers(K)
        X = codebook[i]
        Y = X @ cn((M, N), rho) + cn((T, N), s2)
        G = Y @ Y.conj().T

        # exact ML metric
        Lc = np.linalg.cholesky(G + 1e-12 * np.eye(T))
        Z = Lc.conj().T @ Cflat
        m_ml = np.sum(np.abs(Z) ** 2, axis=0).reshape(K, M).sum(1)

        # subspace-only: top-min(N,M) left singular vectors of Y, equal weight
        U, s, _ = np.linalg.svd(Y, full_matrices=False)
        r = min(N, M)
        Us = U[:, :r]
        Zs = Us.conj().T @ Cflat
        m_sub = np.sum(np.abs(Zs) ** 2, axis=0).reshape(K, M).sum(1)

        # full left basis, singular values replaced by their mean (spectrum flattened)
        gflat = np.full(len(s), s.mean())
        Zf = (U * gflat) .conj().T @ Cflat
        m_spec = np.sum(np.abs(Zf) ** 2, axis=0).reshape(K, M).sum(1)

        err_ml += int(np.argmax(m_ml) != i)
        err_sub += int(np.argmax(m_sub) != i)
        err_spec += int(np.argmax(m_spec) != i)

    ser = (err_ml / n_blocks, err_sub / n_blocks, err_spec / n_blocks)
    res[snr_db] = ser
    print(f"{snr_db:>7} {ser[0]:>10.4f} {ser[1]:>11.4f} {ser[2]:>14.4f}")

# SNR gap at a target SER, by log-linear interpolation
def snr_at(target, idx):
    xs = np.array(snr_db_list, float)
    ys = np.array([max(res[s][idx], 1e-6) for s in snr_db_list])
    ly = np.log10(ys)
    for k in range(len(xs) - 1):
        if ly[k] >= np.log10(target) >= ly[k + 1]:
            t = (np.log10(target) - ly[k]) / (ly[k + 1] - ly[k])
            return xs[k] + t * (xs[k + 1] - xs[k])
    return None

print()
for tgt in [1e-1, 1e-2]:
    a, b, c = snr_at(tgt, 0), snr_at(tgt, 1), snr_at(tgt, 2)
    if a is not None and b is not None:
        print(f"  SER={tgt:g}:  ML {a:.2f} dB,  subspace-only {b:.2f} dB  "
              f"-> gap {b-a:.2f} dB")
    if a is not None and c is not None:
        print(f"  SER={tgt:g}:  ML {a:.2f} dB,  flat-spectrum {c:.2f} dB  "
              f"-> gap {c-a:.2f} dB")
