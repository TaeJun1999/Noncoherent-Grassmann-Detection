"""
normalizer_and_sampler.py  --  Hand-off v2, items 2 and 3 of week 1.

(1) Normalizer  Z(kappa G) = int_{G(M,T)} exp(kappa tr(X^H G X)) [dX] = 1F1~(M;T;kappa G)
    - exact closed form for M=1 (divided difference of exp / Hermite-Genocchi)
    - MACG importance-sampling estimator for general M, validated against it
(2) Gold-standard posterior sampler
    - independence Metropolis-Hastings with MACG(I + cG) proposal
    - validated for T=2, M=1 against exact 1-D quadrature (CP^1)
"""
import numpy as np
from scipy.integrate import quad

rng = np.random.default_rng(1)


def cn(shape):
    return (rng.standard_normal(shape) + 1j * rng.standard_normal(shape)) / np.sqrt(2)


# ------------------------------------------------------------------ M = 1 exact
def F11_M1_exact(a):
    """1F1~(1;T;diag(a)) = E_{u~unif simplex}[exp(<a,u>)]  (distinct a_t)."""
    a = np.asarray(a, float)
    T = len(a)
    from math import factorial
    tot = 0.0
    mx = a.max()
    for t in range(T):
        den = np.prod([a[t] - a[s] for s in range(T) if s != t])
        tot += np.exp(a[t] - mx) / den
    return factorial(T - 1) * tot, mx      # (mantissa, log-offset)


# ------------------------------------------------------------------ MACG tools
def macg_sample(Omega_half, n, T, M):
    Z = cn((n, T, M))
    Q, _ = np.linalg.qr(Omega_half @ Z)
    return Q


def log_macg_density(X, Oinv, T, M, logdetO):
    """log p_MACG(X) = -M logdet(Omega) - T logdet(X^H Omega^-1 X)."""
    S = np.einsum('nti,ts,nsj->nij', X.conj(), Oinv, X)
    _, ld = np.linalg.slogdet(S)
    return -M * logdetO - T * ld


def log_target(X, G, kappa):
    return kappa * np.real(np.einsum('nti,ts,nsi->n', X.conj(), G, X))


def normalizer_IS(G, kappa, M, n=200000, c=None):
    """Self-normalised-free IS estimate of log Z, plus ESS."""
    T = G.shape[0]
    lam = np.linalg.eigvalsh(G)[::-1]
    if c is None:                       # heuristic: match curvature at the mode
        c = 2 * kappa
    Omega = np.eye(T) + c * G
    w_, U_ = np.linalg.eigh(Omega)
    Ohalf = U_ @ np.diag(np.sqrt(w_)) @ U_.conj().T
    Oinv = np.linalg.inv(Omega)
    logdetO = np.linalg.slogdet(Omega)[1]

    X = macg_sample(Ohalf, n, T, M)
    lt = log_target(X, G, kappa)
    lq = log_macg_density(X, Oinv, T, M, logdetO)
    lw = lt - lq
    m = lw.max()
    w = np.exp(lw - m)
    logZ = m + np.log(w.mean())
    ess = w.sum() ** 2 / (w ** 2).sum()
    return logZ, ess / n, c


# ------------------------------------------------------------------ validation 1
print("=" * 72)
print("(1) Normalizer: MACG importance sampling vs exact M=1 closed form")
print("=" * 72)
print(f"{'T':>3} {'kappa':>7} {'log Z exact':>13} {'log Z  IS':>13} {'rel err':>10} {'ESS/n':>7}")
for T in [2, 4, 8, 12]:
    for kappa in [0.5, 2.0, 8.0]:
        lam = np.sort(rng.uniform(0.3, 5.0, T))[::-1]
        V = np.linalg.qr(cn((T, T)))[0]
        G = V @ np.diag(lam) @ V.conj().T
        G = (G + G.conj().T) / 2
        mant, off = F11_M1_exact(kappa * lam)
        logZ_exact = np.log(mant) + off
        logZ_is, ess, c = normalizer_IS(G, kappa, M=1, n=120000)
        rel = abs(np.expm1(logZ_is - logZ_exact))
        print(f"{T:>3} {kappa:>7.2f} {logZ_exact:>13.6f} {logZ_is:>13.6f} "
              f"{rel:>10.2e} {ess:>7.3f}")

# ------------------------------------------------------------------ general M
print()
print("=" * 72)
print("(1b) General M: IS estimate of log Z, ESS as a function of the MACG scale c")
print("=" * 72)
T, M, kappa = 8, 3, 2.0
lam = np.sort(rng.uniform(0.3, 5.0, T))[::-1]
V = np.linalg.qr(cn((T, T)))[0]
G = V @ np.diag(lam) @ V.conj().T
G = (G + G.conj().T) / 2
print(f"T={T} M={M} kappa={kappa}  lam={np.round(lam,3)}")
best = None
for c in [0.0, 0.5, 1.0, 2.0, 4.0, 8.0]:
    logZ, ess, _ = normalizer_IS(G, kappa, M, n=120000, c=c)
    tag = "  <- Haar baseline" if c == 0.0 else ""
    print(f"  c={c:>4.1f}   log Z = {logZ:>10.5f}   ESS/n = {ess:.4f}{tag}")
    if best is None or ess > best[1]:
        best = (c, ess)
print(f"  best c = {best[0]} (ESS/n = {best[1]:.4f})")

# ------------------------------------------------------------------ sampler
print()
print("=" * 72)
print("(2) Gold-standard sampler: independence MH with MACG proposal")
print("=" * 72)


def mh_sampler(G, kappa, M, n_iter=60000, burn=5000, c=2.0, thin=1):
    T = G.shape[0]
    Omega = np.eye(T) + c * G
    w_, U_ = np.linalg.eigh(Omega)
    Ohalf = U_ @ np.diag(np.sqrt(w_)) @ U_.conj().T
    Oinv = np.linalg.inv(Omega)
    logdetO = np.linalg.slogdet(Omega)[1]

    prop = macg_sample(Ohalf, n_iter, T, M)
    lt = log_target(prop, G, kappa)
    lq = log_macg_density(prop, Oinv, T, M, logdetO)
    lw = lt - lq

    idx = np.zeros(n_iter, dtype=int)
    cur = 0
    acc = 0
    u = np.log(rng.uniform(size=n_iter))
    for i in range(1, n_iter):
        if u[i] < lw[i] - lw[cur]:
            cur = i
            acc += 1
        idx[i] = cur
    keep = idx[burn::thin]
    return prop[keep], acc / n_iter


# --- exact reference on CP^1  (T=2, M=1): |x_1|^2 ~ Uniform(0,1) under Haar
T, M, kappa = 2, 1, 3.0
lam = np.array([4.0, 0.7])
G = np.diag(lam).astype(complex)
Zex, _ = quad(lambda u: np.exp(kappa * (lam[0] * u + lam[1] * (1 - u))), 0, 1)
Eu_ex = quad(lambda u: u * np.exp(kappa * (lam[0] * u + lam[1] * (1 - u))), 0, 1)[0] / Zex
Eu2_ex = quad(lambda u: u**2 * np.exp(kappa * (lam[0] * u + lam[1] * (1 - u))), 0, 1)[0] / Zex

S, acc = mh_sampler(G, kappa, M, n_iter=200000, burn=20000, c=2.0)
u_mc = np.abs(S[:, 0, 0]) ** 2
print(f"  CP^1 check (T=2,M=1,kappa={kappa},lam={lam})   acceptance = {acc:.3f}")
print(f"    E[|x1|^2] : exact {Eu_ex:.6f}   MH {u_mc.mean():.6f}   "
      f"err {abs(u_mc.mean()-Eu_ex):.2e}")
print(f"    E[|x1|^4] : exact {Eu2_ex:.6f}   MH {(u_mc**2).mean():.6f}   "
      f"err {abs((u_mc**2).mean()-Eu2_ex):.2e}")

# --- cross-check MH against IS on a larger case
T, M, kappa = 8, 3, 2.0
lam = np.sort(rng.uniform(0.3, 5.0, T))[::-1]
V = np.linalg.qr(cn((T, T)))[0]
G = V @ np.diag(lam) @ V.conj().T
G = (G + G.conj().T) / 2
S, acc = mh_sampler(G, kappa, M, n_iter=150000, burn=20000, c=2.0)
r_mh = np.sort(np.mean(np.sum(np.abs(np.einsum('ts,ntm->nsm', V.conj(), S)) ** 2, axis=2), axis=0))[::-1]

Omega = np.eye(T) + 2.0 * G
w_, U_ = np.linalg.eigh(Omega)
Ohalf = U_ @ np.diag(np.sqrt(w_)) @ U_.conj().T
Xis = macg_sample(Ohalf, 300000, T, M)
lw = log_target(Xis, G, kappa) - log_macg_density(Xis, np.linalg.inv(Omega), T, M,
                                                  np.linalg.slogdet(Omega)[1])
wv = np.exp(lw - lw.max())
r_is = (wv[:, None] * np.sum(np.abs(np.einsum('ts,ntm->nsm', V.conj(), Xis)) ** 2, axis=2)).sum(0) / wv.sum()
print(f"\n  T={T} M={M} kappa={kappa}: acceptance = {acc:.3f}")
print(f"    E[r] MH  = {np.round(r_mh,4)}")
print(f"    E[r] SNIS= {np.round(r_is,4)}")
print(f"    max|diff| = {np.max(np.abs(r_mh - r_is)):.2e}   (sum must be M={M}: "
      f"{r_mh.sum():.4f} / {r_is.sum():.4f})")
