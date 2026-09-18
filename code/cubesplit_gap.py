"""
cubesplit_gap.py  --  gate experiment (00_STATE.md): structured-codebook efficient
demapper vs exact ML, on the SAME codebook.  M = 1 (Cube-split, arXiv 1905.08745).

Conventions follow the paper:  Y = sqrt(rho T) x h^T + Z,  ||x|| = 1,  h ~ CN(0,I_N),
Z ~ CN(0,1).  rho = SNR per channel use per receive antenna.  Block SNR = rho*T.
Exact posterior weight (project R1 / paper eq. 24):  exp(kappa ||Y^H x||^2),
kappa = rho T / (1 + rho T).

Receivers
  ML   : exhaustive argmax_k ||Y^H c_k||^2                         [ground truth]
  D1   : paper's greedy decoder (Sec. IV-A): top left singular vector u,
         cell = argmax|u_j|, inverse map, per-real-dimension quantization
  D1r  : same but ratios taken from one row of G = YY^H (no SVD)    [even cheaper]
  D2   : greedy candidate in EVERY cell (T candidates) + exact metric rescoring
  D3   : D2 + coordinate-descent on the per-coordinate Q^2 grid with exact metric
Soft outputs
  exact: full posterior bit marginals over the codebook
  P26  : paper's eq. (26): greedy x_hat + precomputed per-bit eta-nearest lists
  D3s  : cell bits from the T cell candidates, coordinate bits from the Q^2
         per-coordinate grid at the D3 point (log-sum, exact metric)
"""
import numpy as np
from scipy.special import erfinv, erf


# ----------------------------------------------------------------------------- maps
def probit(p):
    return np.sqrt(2.0) * erfinv(2.0 * p - 1.0)


def ncdf(x):
    return 0.5 * (1.0 + erf(x / np.sqrt(2.0)))


def cs_map(w):
    """eq. (12): Gaussian-domain w -> unit disc t."""
    e = np.exp(-0.5 * np.abs(w) ** 2)
    return np.sqrt((1.0 - e) / (1.0 + e)) * w / np.maximum(np.abs(w), 1e-300)


def cs_unmap(t):
    """eq. (14): unit disc t -> w."""
    s = np.clip(np.abs(t) ** 2, 0.0, 1.0 - 1e-12)
    r = np.sqrt(2.0 * np.log((1.0 + s) / (1.0 - s)))
    return r * t / np.maximum(np.abs(t), 1e-300)


def _bits(v, nb):
    v = np.asarray(v, np.int64)
    return ((v[:, None] >> (nb - 1 - np.arange(nb))[None, :]) & 1).astype(np.uint8)


def _from_bits(b):
    nb = b.shape[1]
    return (b.astype(np.int64) << (nb - 1 - np.arange(nb))[None, :]).sum(1)


# ----------------------------------------------------------------------------- constellation
class CubeSplit:
    def __init__(self, T, B0):
        self.T, self.B0 = T, B0
        self.nc = int(round(np.log2(T)))
        assert 2 ** self.nc == T, "T must be a power of 2 for the bit labeling"
        self.Q = 2 ** B0
        self.grid = (2 * np.arange(self.Q) + 1) / (2.0 * self.Q)      # eq. (9)
        self.gw = probit(self.grid)
        self.nd = 2 * (T - 1)                                          # real coordinates
        self.B = self.nc + self.nd * B0
        self.K = 2 ** self.B
        self.Kloc = 2 ** (self.nd * B0)
        q = np.arange(self.Q)
        self.gray = q ^ (q >> 1)
        self.igray = np.argsort(self.gray)
        ia = np.repeat(q, self.Q)
        ib = np.tile(q, self.Q)
        self.tgrid = cs_map(self.gw[ia] + 1j * self.gw[ib])            # (Q^2,) complex grid of one coordinate
        self.gdig = np.stack([ia, ib], 1)                              # (Q^2, 2) digit pairs
        self.gbits = np.concatenate([_bits(self.gray[ia], B0), _bits(self.gray[ib], B0)], 1)  # (Q^2, 2B0)

    # --- indexing / labeling
    def positions(self, cell):
        m = np.arange(self.T - 1)[None, :]
        return m + (m >= cell[:, None])

    def index(self, cell, digits):
        loc = np.zeros(len(cell), np.int64)
        for m in range(self.nd):
            loc += digits[:, m].astype(np.int64) * (self.Q ** m)
        return cell.astype(np.int64) * self.Kloc + loc

    def from_index(self, k):
        k = np.asarray(k, np.int64)
        cell = k // self.Kloc
        loc = k % self.Kloc
        digits = np.stack([(loc // (self.Q ** m)) % self.Q for m in range(self.nd)], 1)
        return cell, digits

    def bits(self, cell, digits):
        cb = _bits(cell, self.nc)
        db = _bits(self.gray[digits].reshape(-1), self.B0).reshape(len(cell), self.nd * self.B0)
        return np.concatenate([cb, db], 1)

    def from_bits(self, bits):
        cell = _from_bits(bits[:, :self.nc])
        g = _from_bits(bits[:, self.nc:].reshape(-1, self.B0)).reshape(bits.shape[0], self.nd)
        return cell, self.igray[g]

    def symbols(self, cell, digits):
        w = self.gw[digits[:, 0::2]] + 1j * self.gw[digits[:, 1::2]]
        t = cs_map(w)
        x = np.ones((len(cell), self.T), complex)
        np.put_along_axis(x, self.positions(cell), t, axis=1)
        return x / np.sqrt(1.0 + np.sum(np.abs(t) ** 2, 1))[:, None]

    def codebook(self):
        cell, digits = self.from_index(np.arange(self.K))
        return np.ascontiguousarray(self.symbols(cell, digits).T), self.bits(cell, digits)   # C (T,K), bits (K,B)

    # --- greedy step 2 of the paper, from any unit-norm estimate u  (n,T)
    def quantize(self, u, cell=None):
        if cell is None:
            cell = np.argmax(np.abs(u), 1)
        pos = self.positions(cell)
        t = np.take_along_axis(u, pos, 1) / np.take_along_axis(u, cell[:, None], 1)
        w = cs_unmap(t)
        digits = np.empty((len(u), self.nd), np.int64)
        digits[:, 0::2] = np.clip(np.floor(ncdf(w.real) * self.Q), 0, self.Q - 1)
        digits[:, 1::2] = np.clip(np.floor(ncdf(w.imag) * self.Q), 0, self.Q - 1)
        return cell, digits


# ----------------------------------------------------------------------------- channel
def channel(rng, X, N, rho):
    """X (n,T) unit-norm rows -> Y (n,T,N)."""
    n, T = X.shape
    H = (rng.standard_normal((n, N)) + 1j * rng.standard_normal((n, N))) / np.sqrt(2)
    Z = (rng.standard_normal((n, T, N)) + 1j * rng.standard_normal((n, T, N))) / np.sqrt(2)
    return np.sqrt(rho * T) * X[:, :, None] * H[:, None, :] + Z


def kappa_of(rho, T):
    return rho * T / (1.0 + rho * T)


def gram(Y):
    return Y @ Y.conj().transpose(0, 2, 1)


def metric_of(Y, X):
    """||Y^H x||^2 for matched rows: Y (n,T,N), X (n,T) -> (n,)"""
    z = np.einsum('ntm,nt->nm', Y.conj(), X)
    return (z.real ** 2 + z.imag ** 2).sum(1)


# ----------------------------------------------------------------------------- exact ML
def ml_metrics(Y, C, chunk_elems=1.2e7):
    """All-codebook metrics (n,K) = ||Y^H c_k||^2.  C (T,K).  Chunked over blocks."""
    n, T, N = Y.shape
    K = C.shape[1]
    out = np.empty((n, K), np.float64)
    nb = max(1, int(chunk_elems // (N * K)))
    for s in range(0, n, nb):
        Yh = Y[s:s + nb].conj().transpose(0, 2, 1).reshape(-1, T).astype(C.dtype)
        Z = Yh @ C
        out[s:s + nb] = (Z.real.astype(np.float64) ** 2 + Z.imag.astype(np.float64) ** 2).reshape(-1, N, K).sum(1)
    return out


def exact_llr(metrics, bits, kappa, clip=30.0):
    """Exact posterior bit LLRs, L = log P(b=0) - log P(b=1).  bits (K,B) uint8."""
    a = kappa * metrics
    a -= a.max(1, keepdims=True)
    p = np.exp(a)
    p /= p.sum(1, keepdims=True)
    p1 = p @ bits.astype(np.float64)
    p0 = 1.0 - p1
    L = np.log(np.maximum(p0, 1e-300)) - np.log(np.maximum(p1, 1e-300))
    return np.clip(L, -clip, clip)


# ----------------------------------------------------------------------------- greedy family
def top_vec(G):
    """principal eigenvector of G = Y Y^H (n,T,T) == top left singular vector of Y."""
    _, V = np.linalg.eigh(G)
    return V[:, :, -1]


def decode_D1(cs, G):
    return cs.quantize(top_vec(G))


def decode_D1r(cs, G):
    """row-ratio variant: cell = argmax G_jj, t_j = G_{j,cell} / G_{cell,cell}."""
    d = np.real(np.einsum('ntt->nt', G))
    cell = np.argmax(d, 1)
    u = np.take_along_axis(G, cell[:, None, None].repeat(G.shape[1], 1), 2)[:, :, 0]   # column `cell` of G
    u = u / np.sqrt(np.sum(np.abs(u) ** 2, 1))[:, None]
    return cs.quantize(u, cell=cell)


def decode_D2(cs, G, Y):
    """greedy candidate forced in every cell, exact rescoring.  Returns cell, digits, X, per-cell metrics (n,T), per-cell digits."""
    u = top_vec(G)
    n, T = u.shape
    mus = np.empty((n, T))
    digs = np.empty((n, T, cs.nd), np.int64)
    for i in range(T):
        ci = np.full(n, i)
        _, di = cs.quantize(u, cell=ci)
        Xi = cs.symbols(ci, di)
        mus[:, i] = metric_of(Y, Xi)
        digs[:, i] = di
    cell = np.argmax(mus, 1)
    digits = np.take_along_axis(digs, cell[:, None, None], 1)[:, 0]
    return cell, digits, mus, digs


def refine_D3(cs, G, cell, digits, max_sweeps=3):
    """coordinate descent over the Q^2 grid of each coordinate, exact metric, starting at (cell,digits).
    Returns cell, digits, final metric (n,), per-coordinate grid metrics at the final point (n,T-1,Q^2)."""
    n, T = G.shape[:2]
    ar = np.arange(n)
    pos = cs.positions(cell)                                    # (n,T-1)
    t = cs.tgrid[digits[:, 0::2] * cs.Q + digits[:, 1::2]]     # (n,T-1) current coordinate values
    u = np.ones((n, T), complex)
    np.put_along_axis(u, pos, t, axis=1)
    Gu = np.einsum('nts,ns->nt', G, u)
    quad = np.real(np.einsum('nt,nt->n', u.conj(), Gu))
    s = 1.0 + np.sum(np.abs(t) ** 2, 1)
    g = cs.tgrid[None, :]
    g2 = np.abs(cs.tgrid) ** 2
    digits = digits.copy()
    for sweep in range(max_sweeps):
        changed = False
        for m in range(T - 1):
            j = pos[:, m]
            told = t[:, m]
            Gjj = np.real(G[ar, j, j])
            Guj = Gu[ar, j]
            quad_v = quad - 2 * np.real(np.conj(told) * Guj) + np.abs(told) ** 2 * Gjj
            Gvj = Guj - told * Gjj
            s_v = s - np.abs(told) ** 2
            mg = (quad_v[:, None] + 2 * np.real(np.conj(g) * Gvj[:, None]) + g2[None, :] * Gjj[:, None]) / (s_v[:, None] + g2[None, :])
            best = np.argmax(mg, 1)
            gnew = cs.tgrid[best]
            upd = best != (digits[:, 2 * m] * cs.Q + digits[:, 2 * m + 1])
            if upd.any():
                changed = True
            digits[:, 2 * m] = cs.gdig[best, 0]
            digits[:, 2 * m + 1] = cs.gdig[best, 1]
            Gu += (gnew - told)[:, None] * G[ar, :, j]
            quad = quad_v + 2 * np.real(np.conj(gnew) * Gvj) + np.abs(gnew) ** 2 * Gjj
            s = s_v + np.abs(gnew) ** 2
            t[:, m] = gnew
        if not changed:
            break
    # final pass: grid metrics at the converged point (for soft output)
    grid_m = np.empty((n, T - 1, cs.Q * cs.Q))
    for m in range(T - 1):
        j = pos[:, m]
        told = t[:, m]
        Gjj = np.real(G[ar, j, j])
        Guj = Gu[ar, j]
        quad_v = quad - 2 * np.real(np.conj(told) * Guj) + np.abs(told) ** 2 * Gjj
        Gvj = Guj - told * Gjj
        s_v = s - np.abs(told) ** 2
        grid_m[:, m] = (quad_v[:, None] + 2 * np.real(np.conj(g) * Gvj[:, None]) + g2[None, :] * Gjj[:, None]) / (s_v[:, None] + g2[None, :])
    return cell, digits, quad / s, grid_m


def lse(a, axis):
    m = a.max(axis, keepdims=True)
    return (m + np.log(np.sum(np.exp(a - m), axis, keepdims=True))).squeeze(axis)


def soft_D3(cs, G, Y, kappa, clip=30.0):
    """D3s soft output.  Returns LLR (n,B), hard cell, hard digits."""
    cell, digits, mus, digs = decode_D2(cs, G, Y)
    cell, digits, mfin, grid_m = refine_D3(cs, G, cell, digits)
    n = len(cell)
    mus = mus.copy()
    mus[np.arange(n), cell] = mfin                       # refined metric for the winning cell
    a = kappa * mus                                      # (n,T)
    cb = _bits(np.arange(cs.T), cs.nc).astype(bool)      # (T,nc)
    L = np.empty((n, cs.B))
    for b in range(cs.nc):
        m0 = np.where(cb[None, :, b], -np.inf, a)
        m1 = np.where(cb[None, :, b], a, -np.inf)
        L[:, b] = lse(m0, 1) - lse(m1, 1)
    gb = cs.gbits.astype(bool)                            # (Q^2, 2B0)
    ag = kappa * grid_m                                   # (n,T-1,Q^2)
    for q in range(2 * cs.B0):
        m0 = np.where(gb[None, None, :, q], -np.inf, ag)
        m1 = np.where(gb[None, None, :, q], ag, -np.inf)
        Lq = lse(m0, 2) - lse(m1, 2)                      # (n,T-1)
        L[:, cs.nc + q::2 * cs.B0] = Lq
    return np.clip(L, -clip, clip), cell, digits


# ----------------------------------------------------------------------------- paper eq. (26)
def precompute_lists(cs, C, bits, eta, chunk=512, log=None):
    """For each codeword c, bit j, value b: indices of the eta nearest (chordal) codewords with bit j = b.
    Returns int32 array (K, B, 2, eta)."""
    K, B = bits.shape
    lists = np.empty((K, B, 2, eta), np.int32)
    Ch = C.conj().T                                     # (K,T)
    mask1 = bits.astype(bool)
    for s in range(0, K, chunk):
        corr = np.abs(Ch[s:s + chunk] @ C) ** 2          # (chunk,K) |c^H c'|^2, larger = closer
        for j in range(B):
            for b in (0, 1):
                v = np.where(mask1[:, j] == bool(b), corr, -1.0)
                idx = np.argpartition(-v, eta - 1, axis=1)[:, :eta]
                lists[s:s + chunk, j, b] = idx
        if log is not None and (s // chunk) % 8 == 0:
            log(f"    lists {s + chunk}/{K}")
    return lists


def soft_P26(cs, G, Y, C, lists, kappa, clip=30.0):
    """paper's approximate LLR (26) with the paper's greedy x_hat."""
    cell, digits = decode_D1(cs, G)
    k = cs.index(cell, digits)
    n = len(k)
    L = lists[k]                                         # (n,B,2,eta)
    B, eta = L.shape[1], L.shape[3]
    cand = C[:, L.reshape(n, -1)]                         # (T, n, B*2*eta)
    z = np.einsum('ntm,tnc->nmc', Y.conj(), cand)         # (n,N,C)
    met = (z.real ** 2 + z.imag ** 2).sum(1).reshape(n, B, 2, eta)
    a = kappa * met
    out = lse(a[:, :, 0, :], 2) - lse(a[:, :, 1, :], 2)
    return np.clip(out, -clip, clip), cell, digits


# ----------------------------------------------------------------------------- LDPC
def make_ldpc(n, dv, dc, rng):
    """Gallager (dv,dc)-regular parity-check matrix, n divisible by dc."""
    assert n % dc == 0
    rb = n // dc
    m = rb * dv
    H = np.zeros((m, n), np.uint8)
    base = np.zeros((rb, n), np.uint8)
    for r in range(rb):
        base[r, r * dc:(r + 1) * dc] = 1
    for b in range(dv):
        perm = np.arange(n) if b == 0 else rng.permutation(n)
        H[b * rb:(b + 1) * rb] = base[:, perm]
    return H


def gf2_rref(H):
    A = H.copy()
    m, n = A.shape
    piv = []
    r = 0
    for c in range(n):
        if r >= m:
            break
        rows = np.nonzero(A[r:, c])[0]
        if len(rows) == 0:
            continue
        p = r + rows[0]
        if p != r:
            A[[r, p]] = A[[p, r]]
        others = np.nonzero(A[:, c])[0]
        others = others[others != r]
        A[others] ^= A[r]
        piv.append(c)
        r += 1
    return A[:r], np.array(piv)


class LDPC:
    def __init__(self, n, dv, dc, rng):
        self.H = make_ldpc(n, dv, dc, rng)
        self.n, self.dv, self.dc = n, dv, dc
        R, piv = gf2_rref(self.H)
        self.R, self.piv = R, piv
        self.free = np.setdiff1d(np.arange(n), piv)
        self.k = len(self.free)
        self.m = self.H.shape[0]
        # edge structure (regular)
        self.E_var = np.stack([np.nonzero(self.H[c])[0] for c in range(self.m)])     # (m,dc)
        flat = np.arange(self.m * self.dc).reshape(self.m, self.dc)
        V = [[] for _ in range(n)]
        for c in range(self.m):
            for s_ in range(self.dc):
                V[self.E_var[c, s_]].append(flat[c, s_])
        self.V_edge = np.array(V)                                                      # (n,dv)

    def encode(self, u):
        """u (nc,k) -> x (nc,n)"""
        x = np.zeros((u.shape[0], self.n), np.uint8)
        x[:, self.free] = u
        x[:, self.piv] = (u.astype(np.int64) @ self.R[:, self.free].T.astype(np.int64)) % 2
        return x

    def decode(self, L, iters=50):
        """sum-product, L (nc,n) with positive -> bit 0.  Returns hard decisions (nc,n)."""
        nc = L.shape[0]
        L = np.clip(L, -30, 30)
        Mvc = L[:, self.E_var]
        xhat = (L < 0).astype(np.uint8)
        for it in range(iters):
            Tn = np.tanh(0.5 * np.clip(Mvc, -30, 30))
            left = np.cumprod(Tn, 2)
            right = np.cumprod(Tn[:, :, ::-1], 2)[:, :, ::-1]
            loo = np.ones_like(Tn)
            loo[:, :, 1:] *= left[:, :, :-1]
            loo[:, :, :-1] *= right[:, :, 1:]
            Mcv = 2.0 * np.arctanh(np.clip(loo, -1 + 1e-12, 1 - 1e-12))
            Ltot = L + Mcv.reshape(nc, -1)[:, self.V_edge].sum(2)
            Mvc = Ltot[:, self.E_var] - Mcv
            xhat = (Ltot < 0).astype(np.uint8)
            synd = xhat[:, self.E_var].sum(2) % 2
            if not synd.any():
                break
        return xhat


# ----------------------------------------------------------------------------- helpers
def snr_at(target, snrs, sers):
    """log-linear interpolation of the SNR at which the error rate crosses `target`."""
    xs = np.asarray(snrs, float)
    ly = np.log10(np.maximum(np.asarray(sers, float), 1e-9))
    lt = np.log10(target)
    for k in range(len(xs) - 1):
        if ly[k] >= lt >= ly[k + 1] and ly[k] > ly[k + 1]:
            return xs[k] + (lt - ly[k]) / (ly[k + 1] - ly[k]) * (xs[k + 1] - xs[k])
    return None


# ----------------------------------------------------------------------------- fast exact ML via real features
def codebook_features(C):
    """Phi (K, T^2) float32 s.t.  c_k^H G c_k = Phi[k] . gvec(G)  for Hermitian G."""
    T, K = C.shape
    iu = np.triu_indices(T, 1)
    P = np.einsum('tk,sk->tsk', C.conj(), C)            # (T,T,K)
    diag = np.real(np.einsum('ttk->tk', P)).T              # (K,T)
    cross = P[iu[0], iu[1], :].T                           # (K, T(T-1)/2)
    return np.concatenate([diag, cross.real, cross.imag], 1).astype(np.float32)


def gvec(G):
    """G (n,T,T) Hermitian -> (n,T^2) float32 matching codebook_features."""
    n, T, _ = G.shape
    iu = np.triu_indices(T, 1)
    diag = np.real(np.einsum('ntt->nt', G))
    cross = G[:, iu[0], iu[1]]
    return np.concatenate([diag, 2 * cross.real, -2 * cross.imag], 1).astype(np.float32)


def ml_metrics_fast(G, Phi, nb=64):
    """(n,K) float32 metrics c_k^H G c_k for all codewords, batched sgemm."""
    n = G.shape[0]
    out = np.empty((n, Phi.shape[0]), np.float32)
    for s in range(0, n, nb):
        out[s:s + nb] = gvec(G[s:s + nb]) @ Phi.T
    return out


def exact_llr_chunked(G, Phi, bits, kappa, chunk=256, clip=30.0):
    """memory-safe exact posterior LLRs over many blocks."""
    out = np.empty((G.shape[0], bits.shape[1]))
    for s in range(0, G.shape[0], chunk):
        m = ml_metrics_fast(G[s:s + chunk], Phi, nb=64).astype(np.float64)
        out[s:s + chunk] = exact_llr(m, bits, kappa, clip)
    return out


def precompute_lists_sorted(cs, C, bits, eta_max, chunk=256, topm=2048, log=None):
    """Paper's B_j(c,b): for each codeword, bit j, value b, the eta_max nearest codewords (chordal,
    descending closeness) with bit j = b.  Two-stage: top-`topm` overall per row, then per-(j,b)
    selection; rows with an unfilled (j,b) list fall back to the full row.  Returns (K,B,2,eta_max) int32
    and the number of fallback rows."""
    K, B = bits.shape
    T = C.shape[0]
    lists = np.empty((K, B, 2, eta_max), np.int32)
    Ch = np.ascontiguousarray(C.conj().T).astype(np.complex64)
    C32 = C.astype(np.complex64)
    mask1 = bits.astype(bool)
    nfall = 0
    for s in range(0, K, chunk):
        corr = np.abs(Ch[s:s + chunk] @ C32) ** 2                       # (r,K)
        r = corr.shape[0]
        if topm < K:
            top = np.argpartition(-corr, topm - 1, axis=1)[:, :topm]      # (r,topm)
            ctop = np.take_along_axis(corr, top, 1)
        else:
            top = np.tile(np.arange(K), (r, 1)); ctop = corr
        btop = mask1[top]                                                 # (r,topm,B)
        for j in range(B):
            for b in (0, 1):
                v = np.where(btop[:, :, j] == bool(b), ctop, -1.0)
                idx = np.argsort(-v, axis=1, kind='stable')[:, :eta_max]
                lists[s:s + r, j, b] = np.take_along_axis(top, idx, 1)
                bad = np.take_along_axis(v, idx, 1)[:, -1] < 0            # unfilled list
                if bad.any():
                    for rr in np.nonzero(bad)[0]:
                        nfall += 1
                        vf = np.where(mask1[:, j] == bool(b), corr[rr], -1.0)
                        lists[s + rr, j, b] = np.argsort(-vf, kind='stable')[:eta_max]
        if log is not None and (s // chunk) % 32 == 0:
            log(f"    lists {min(s + chunk, K)}/{K}  fallbacks so far {nfall}")
    return lists, nfall


def soft_P26_multi(cs, G, Y, C, lists, etas, kappa, clip=30.0, chunk=400):
    """paper eq. (26) for several eta at once (lists sorted by closeness). Returns dict eta -> LLR (n,B)."""
    cell, digits = decode_D1(cs, G)
    k = cs.index(cell, digits)
    n = len(k)
    B, em = lists.shape[1], lists.shape[3]
    out = {e: np.empty((n, B)) for e in etas}
    for s0 in range(0, n, chunk):
        kk = k[s0:s0 + chunk]
        m = len(kk)
        L = lists[kk]                                      # (m,B,2,em)
        cand = C[:, L.reshape(m, -1)]                      # (T,m,B*2*em)
        z = np.einsum('ntm,tnc->nmc', Y[s0:s0 + chunk].conj(), cand)
        a = kappa * (z.real ** 2 + z.imag ** 2).sum(1).reshape(m, B, 2, em)
        for e in etas:
            out[e][s0:s0 + chunk] = np.clip(lse(a[:, :, 0, :e], 2) - lse(a[:, :, 1, :e], 2), -clip, clip)
    return out, cell, digits
