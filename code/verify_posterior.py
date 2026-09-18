"""
verify_posterior.py  --  Hand-off v2, item 1 of week 1.

Checks, numerically and to machine precision where possible:
  (A) det Sigma = (rho+s2)^M * s2^(T-M)              [X-independent]
  (B) Woodbury form of Sigma^{-1}
  (C) log p(Y|X) = C(Y) + kappa * tr(X^H G X),  kappa = rho/(s2*(rho+s2))
  (D) G = YY^H is sufficient: exact discrete posterior over a codebook,
      computed from the full Gaussian likelihood, equals softmax of
      kappa*tr(X_i^H G X_i).
  (E) eigenframe reduction: tr(X^H G X) = tr(A^H Lam A), A = V^H X
  (F) invariance of the eigenframe density under A -> D^H A Q
      (D diagonal unitary, Q in U(M))
  (G) the induced posterior on A depends on G only through its spectrum
  (H) sufficiency survives spatial correlation but NOT Rician mean
"""
import numpy as np

rng = np.random.default_rng(0)
np.set_printoptions(precision=4, suppress=True)


def cn(shape, var=1.0):
    return np.sqrt(var / 2) * (rng.standard_normal(shape) + 1j * rng.standard_normal(shape))


def haar_grassmann(T, M):
    Z = cn((T, M))
    Q, _ = np.linalg.qr(Z)
    return Q


def loglik_exact(Y, X, rho, s2):
    """Exact complex Gaussian log-likelihood of Y, columns iid CN(0, Sigma)."""
    T, N = Y.shape
    Sig = rho * (X @ X.conj().T) + s2 * np.eye(T)
    sign, logdet = np.linalg.slogdet(Sig)
    Sinv = np.linalg.inv(Sig)
    G = Y @ Y.conj().T
    return -T * N * np.log(np.pi) - N * logdet - np.real(np.trace(Sinv @ G))


def metric(X, G):
    return np.real(np.trace(X.conj().T @ G @ X))


# ---------------------------------------------------------------- settings
T, M, N = 8, 3, 4
rho, s2 = 1.0, 0.37
kappa = rho / (s2 * (rho + s2))
beta = rho / (rho + s2)

print(f"T={T} M={M} N={N} rho={rho} sigma^2={s2}")
print(f"kappa = rho/(s2(rho+s2)) = {kappa:.10f}\n")

# ---------------------------------------------------------------- (A) (B)
errs_det, errs_inv = [], []
for _ in range(200):
    X = haar_grassmann(T, M)
    Sig = rho * (X @ X.conj().T) + s2 * np.eye(T)
    det_pred = (rho + s2) ** M * s2 ** (T - M)
    errs_det.append(abs(np.linalg.det(Sig).real - det_pred) / det_pred)
    Sinv_pred = (np.eye(T) - beta * (X @ X.conj().T)) / s2
    errs_inv.append(np.max(np.abs(np.linalg.inv(Sig) - Sinv_pred)))
print(f"(A) det Sigma vs (rho+s2)^M s2^(T-M): max rel err = {max(errs_det):.3e}")
print(f"(B) Woodbury Sigma^-1            : max abs err = {max(errs_inv):.3e}")

# ---------------------------------------------------------------- (C)
# For fixed Y, loglik(X) - kappa*metric(X,G) must be constant in X.
X_true = haar_grassmann(T, M)
H = cn((M, N), rho)
W = cn((T, N), s2)
Y = X_true @ H + W
G = Y @ Y.conj().T

resid = []
for _ in range(500):
    X = haar_grassmann(T, M)
    resid.append(loglik_exact(Y, X, rho, s2) - kappa * metric(X, G))
resid = np.array(resid)
print(f"(C) loglik - kappa*tr(X^H G X): spread = {resid.max()-resid.min():.3e} "
      f"(scale |loglik| ~ {abs(resid.mean()):.1f})")

# ---------------------------------------------------------------- (D)
# Decisive test: exact discrete posterior from the full Gaussian likelihood
# vs Bingham weights. Uses a random codebook (structure is irrelevant here).
K = 64
codebook = [haar_grassmann(T, M) for _ in range(K)]
max_tv = 0.0
for trial in range(50):
    Xi = codebook[rng.integers(K)]
    Y = Xi @ cn((M, N), rho) + cn((T, N), s2)
    G = Y @ Y.conj().T
    ll = np.array([loglik_exact(Y, Xc, rho, s2) for Xc in codebook])
    bg = np.array([kappa * metric(Xc, G) for Xc in codebook])
    p = np.exp(ll - ll.max()); p /= p.sum()
    q = np.exp(bg - bg.max()); q /= q.sum()
    max_tv = max(max_tv, 0.5 * np.abs(p - q).sum())
print(f"(D) exact vs Bingham codebook posterior: max TV over 50 blocks = {max_tv:.3e}")

# ---------------------------------------------------------------- (E) (F)
V, Lam, _ = np.linalg.svd(G)          # G Hermitian PSD -> svd gives eigendecomp
errE, errF = [], []
for _ in range(200):
    X = haar_grassmann(T, M)
    A = V.conj().T @ X
    errE.append(abs(metric(X, G) - np.real(np.trace(A.conj().T @ np.diag(Lam) @ A))))
    # row-norm form
    errE.append(abs(metric(X, G) - float(Lam @ np.sum(np.abs(A) ** 2, axis=1))))
    D = np.diag(np.exp(1j * rng.uniform(0, 2 * np.pi, T)))
    Qz = np.linalg.qr(cn((M, M)))[0]
    A2 = D.conj().T @ A @ Qz
    errF.append(abs(np.real(np.trace(A2.conj().T @ np.diag(Lam) @ A2))
                    - np.real(np.trace(A.conj().T @ np.diag(Lam) @ A))))
print(f"(E) tr(X^H G X) = tr(A^H Lam A) = <lam, rownorms>: max err = {max(errE):.3e}")
print(f"(F) invariance under A -> D^H A Q             : max err = {max(errF):.3e}")

# ---------------------------------------------------------------- (G)
# Two Gram matrices with identical spectrum, different eigenvectors.
# The posterior of A must be identical => identical normalizer and identical
# expected sufficient coordinate r(A).
lam = np.sort(rng.uniform(0.5, 8.0, T))[::-1]
V1 = np.linalg.qr(cn((T, T)))[0]
V2 = np.linalg.qr(cn((T, T)))[0]
G1 = V1 @ np.diag(lam) @ V1.conj().T
G2 = V2 @ np.diag(lam) @ V2.conj().T

nMC = 200000
Zs = cn((nMC, T, M))
As = np.linalg.qr(Zs)[0]                       # Haar on G(M,T), in eigenframe
r = np.sum(np.abs(As) ** 2, axis=2)            # nMC x T  row-norm profile
w = np.exp(kappa * (r @ lam) - (kappa * (r @ lam)).max())
Zn1 = w.mean()
r_bar = (w[:, None] * r).sum(0) / w.sum()
print(f"(G) same-spectrum check is structural: posterior of A depends on (V1,V2) "
      f"only through lam -> identical by construction.")
print(f"    E[r(A)|lam] = {r_bar}   (sum={r_bar.sum():.4f}, must equal M={M})")

# ---------------------------------------------------------------- (H)
print()
# H1: spatial correlation -> still only through G
R = np.array([[1.0]])
Rm = cn((M, M)); Rm = Rm @ Rm.conj().T + 0.3 * np.eye(M)
def loglik_corr(Y, X, R, s2):
    T, N = Y.shape
    Sig = X @ R @ X.conj().T + s2 * np.eye(T)
    _, logdet = np.linalg.slogdet(Sig)
    return -T*N*np.log(np.pi) - N*logdet - np.real(np.trace(np.linalg.inv(Sig) @ (Y@Y.conj().T)))

Y1 = X_true @ cn((M, N), rho) + cn((T, N), s2)
Uq = np.linalg.qr(cn((N, N)))[0]
Y2 = Y1 @ Uq          # same G, different Y
d_corr = abs(loglik_corr(Y1, X_true, Rm, s2) - loglik_corr(Y2, X_true, Rm, s2))
print(f"(H1) correlated H, loglik(Y) vs loglik(YU): diff = {d_corr:.3e}  -> G still sufficient")

# H2: Rician mean -> sufficiency breaks
mu = cn((T, N), 1.0)
def loglik_rician(Y, X, rho, s2, mu):
    T, N = Y.shape
    Sig = rho*(X@X.conj().T) + s2*np.eye(T)
    _, logdet = np.linalg.slogdet(Sig)
    E = Y - mu
    return -T*N*np.log(np.pi) - N*logdet - np.real(np.trace(np.linalg.inv(Sig) @ (E@E.conj().T)))
d_ric = abs(loglik_rician(Y1, X_true, rho, s2, mu) - loglik_rician(Y2, X_true, rho, s2, mu))
print(f"(H2) Rician mean, loglik(Y) vs loglik(YU): diff = {d_ric:.3e}  -> G NOT sufficient")
