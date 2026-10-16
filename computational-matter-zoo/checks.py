"""Numerical checks backing the claims in REPORT.md.

Run: python3 checks.py   (needs numpy only)
Each check prints PASS/FAIL plus the numbers it relied on.
"""
import numpy as np

rng = np.random.default_rng(0)


def check(name, ok, detail=""):
    print(f"[{'PASS' if ok else 'FAIL'}] {name}" + (f"  -- {detail}" if detail else ""))


# ---------------------------------------------------------------------------
# 1. Three-site flux triangle: routing table, isospectrality, gauge invariance
# ---------------------------------------------------------------------------
def H_tri(phi):
    return np.array([[0, 1, np.exp(-1j * phi)],
                     [1, 0, 1],
                     [np.exp(1j * phi), 1, 0]])


def S_matrix(H, E=0.0, gamma=2.0):
    n = H.shape[0]
    G = np.linalg.inv(E * np.eye(n) - H + 1j * gamma / 2 * np.eye(n))
    return np.eye(n) - 1j * gamma * G


rows = {}
for phi in (0.0, np.pi / 2, -np.pi / 2):
    S = S_matrix(H_tri(phi))
    rows[phi] = np.abs(S[:, 0]) ** 2
    print(f"    flux {phi/np.pi:+.2f} pi  -> |S_p1|^2 = {np.round(rows[phi], 3)}"
          f"  spectrum {np.round(np.linalg.eigvalsh(H_tri(phi)), 3)}")
check("flux 0 gives 20/40/40", np.allclose(rows[0.0], [0.2, 0.4, 0.4]))
check("flux +pi/2 routes 1->2", np.allclose(rows[np.pi / 2], [0, 1, 0], atol=1e-12))
check("flux -pi/2 routes 1->3", np.allclose(rows[-np.pi / 2], [0, 0, 1], atol=1e-12))
check("+pi/2 and -pi/2 are isospectral",
      np.allclose(np.linalg.eigvalsh(H_tri(np.pi / 2)), np.linalg.eigvalsh(H_tri(-np.pi / 2))))

# Isospectral by construction: H(-phi) = conj(H(phi)) (time reversal).
check("H(-phi) == conj(H(phi)) for random phi",
      all(np.allclose(H_tri(-p), H_tri(p).conj()) for p in rng.uniform(-np.pi, np.pi, 5)))

# The pair is also a mirror pair: swapping ports 2<->3 maps one to the other (up to gauge).
P = np.eye(3)[[0, 2, 1]]
loop = lambda M: M[0, 1] * M[1, 2] * M[2, 0]
check("swap(2,3) of H(+pi/2) has the loop flux of H(-pi/2)",
      np.isclose(loop(P @ H_tri(np.pi / 2) @ P.T), loop(H_tri(-np.pi / 2))))

# Local rephasing (gauge) changes the matrix, not the routing probabilities.
D = np.diag(np.exp(1j * rng.uniform(0, 2 * np.pi, 3)))
H0 = H_tri(0.7)
S0, S1 = S_matrix(H0), S_matrix(D @ H0 @ D.conj().T)
check("gauge transform leaves |S|^2 unchanged",
      np.allclose(np.abs(S0) ** 2, np.abs(S1) ** 2),
      f"raw matrix change {np.linalg.norm(D @ H0 @ D.conj().T - H0):.2f}")

# ---------------------------------------------------------------------------
# 2. Repertoire dimension = rank of d|S|^2 / d(controls)
# ---------------------------------------------------------------------------
def H_general(t):
    """3-site loop: three real hoppings + one loop flux = 4 controls."""
    t12, t23, t31, phi = t
    return np.array([[0, t12, t31 * np.exp(-1j * phi)],
                     [t12, 0, t23],
                     [t31 * np.exp(1j * phi), t23, 0]])


def repertoire_rank(Hfun, x, eps=1e-6):
    f = lambda y: (np.abs(S_matrix(Hfun(y))) ** 2).ravel()
    J = np.stack([(f(x + eps * e) - f(x - eps * e)) / (2 * eps) for e in np.eye(len(x))], 1)
    s = np.linalg.svd(J, compute_uv=False)
    return int((s > 1e-6 * s[0]).sum()), s


r1, _ = repertoire_rank(lambda y: H_tri(y[0]), np.array([0.3]))
r4, s4 = repertoire_rank(H_general, np.array([1.0, 0.8, 1.2, 0.3]))
check("flux-only triangle: repertoire rank 1", r1 == 1)
check("4-control triangle: repertoire rank 4 <= Birkhoff bound (N-1)^2 = 4",
      r4 <= 4, f"rank {r4}, singular values {np.round(s4, 3)}")

# ---------------------------------------------------------------------------
# 3. Attention logits as a 1-form on the query-key bipartite graph
# ---------------------------------------------------------------------------
m, n = 5, 7
S = rng.normal(size=(m, n))
softmax = lambda X: np.exp(X - X.max(1, keepdims=True)) / np.exp(X - X.max(1, keepdims=True)).sum(1, keepdims=True)
A = softmax(S)
A_shift = softmax(S + rng.normal(size=(m, 1)))           # row shift = softmax gauge
A_key = softmax(S + rng.normal(size=(1, n)))             # key bias = physical
check("softmax invariant to row shifts (gauge)", np.allclose(A, A_shift))
check("softmax NOT invariant to key biases (potential)", not np.allclose(A, A_key))

# Cycle space of K_{m,n} has dim (m-1)(n-1) = dim of double-centered matrices.
proj = np.kron(np.eye(m) - 1 / m, np.eye(n) - 1 / n)
check("dim(double-centered) == dim(cycle space of K_mn)",
      np.linalg.matrix_rank(proj) == (m - 1) * (n - 1), f"{(m-1)*(n-1)}")

# Plaquette flux = log odds ratio of attention, independent of the softmax normalizers.
i, k, j, l = 0, 1, 2, 3
flux_logits = S[i, j] - S[i, l] - S[k, j] + S[k, l]
log_or = np.log(A[i, j] * A[k, l] / (A[i, l] * A[k, j]))
check("4-cycle flux of logits == log odds ratio of attention", np.isclose(flux_logits, log_or))

# QRR-style node-phase transport has zero loop flux (pure gauge).
gam = np.exp(1j * rng.uniform(0, 2 * np.pi, 3))
U = lambda a, b: gam[a] * gam[b].conj()
check("U_ij = g_i conj(g_j) has trivial triangle holonomy",
      np.isclose(U(0, 1) * U(1, 2) * U(2, 0), 1))

# ---------------------------------------------------------------------------
# 4. Agent-diversity metric: Laplacian entropy is inverted; kernel (Vendi) is not
# ---------------------------------------------------------------------------
def spectral_entropy_eff(evals):
    ev = np.clip(evals, 0, None)
    p = ev / ev.sum()
    p = p[p > 0]
    return float(np.exp(-(p * np.log(p)).sum()))


N = 6
W_same = np.ones((N, N)) - np.eye(N)                     # all agents identical
L_same = np.diag(W_same.sum(1)) - W_same
lap_same = spectral_entropy_eff(np.linalg.eigvalsh(L_same))

K_div = np.eye(N)                                         # mutually orthogonal agent embeddings
vendi = lambda K: spectral_entropy_eff(np.linalg.eigvalsh(K / K.shape[0]))
check("Laplacian 'diversity' of 6 IDENTICAL agents is ~5 (near max) -> metric inverted",
      lap_same > 4.9, f"{lap_same:.2f}")
check("Vendi score of 6 identical agents == 1", np.isclose(vendi(np.ones((N, N))), 1.0))
check("Vendi score of 6 orthogonal agents == 6", np.isclose(vendi(K_div), 6.0))
