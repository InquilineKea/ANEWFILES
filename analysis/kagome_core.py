"""Independent numpy implementation of the two models behind the kagome pages.

twisted-kagome-moire-bands.html   -> moire_* functions (commensurate supercell tight binding)
twisted-kagome-weyl-hyperlattice.html (Weyl block) -> weyl_* functions (stacked kagome, 3x3 Bloch)

Written separately from the JavaScript so the two can check each other: the JS uses a real
2n x 2n embedding with Householder + QL; this uses numpy.linalg.eigh on the complex matrix.
Only numpy is required.
"""
from __future__ import annotations

import math
from dataclasses import dataclass, field

import numpy as np

SQ3 = math.sqrt(3.0)
TAU = 2.0 * math.pi


# --------------------------------------------------------------------------- moire supercell
@dataclass
class MoireParams:
    m: int = 3            # commensurate index, theta from cos(theta) = (3m^2+3m+1/2)/(3m^2+3m+1)
    t: float = 1.0        # intralayer hopping (-t convention, flat band at +2t)
    breath: float = 0.0   # t_up = t(1-delta), t_down = t(1+delta)
    flux: float = 0.0     # degrees per triangle
    tp: float = 0.35      # interlayer amplitude
    xi: float = 0.25      # Gaussian decay length of interlayer hopping
    rc: float = 0.8       # cutoff (capped at 4 xi as in the JS)
    V: float = 0.0        # interlayer bias, layer 0 at -V/2
    sx: float = 0.0       # slide of layer 1 (lattice units, before rotation)
    sy: float = 0.0


@dataclass
class Supercell:
    p: MoireParams
    N: int
    theta: float
    S1: np.ndarray
    S2: np.ndarray
    xy: np.ndarray        # (n, 2)
    layer: np.ndarray     # (n,)
    sub: np.ndarray       # (n,)
    hop_i: np.ndarray = field(repr=False)
    hop_j: np.ndarray = field(repr=False)
    hop_amp: np.ndarray = field(repr=False)   # complex amplitude of H[j, i]
    hop_d: np.ndarray = field(repr=False)     # displacement r_j - r_i (image)

    @property
    def n(self) -> int:
        return len(self.layer)

    def reciprocal(self):
        S1, S2 = self.S1, self.S2
        det = S1[0] * S2[1] - S1[1] * S2[0]
        b1 = TAU * np.array([S2[1], -S2[0]]) / det
        b2 = TAU * np.array([-S1[1], S1[0]]) / det
        return b1, b2

    def high_symmetry(self):
        b1, b2 = self.reciprocal()
        return {"G": np.zeros(2), "K": (2 * b1 + b2) / 3, "M": (b1 + b2) / 2}


def _rot(v, a):
    c, s = math.cos(a), math.sin(a)
    return np.array([c * v[0] - s * v[1], s * v[0] + c * v[1]])


def moire_theta(m: int) -> float:
    N = 3 * m * m + 3 * m + 1
    return math.acos((3 * m * m + 3 * m + 0.5) / N)


def moire_build(p: MoireParams) -> Supercell:
    m = p.m
    N = 3 * m * m + 3 * m + 1
    th = moire_theta(m)
    S1 = _rot([m + (m + 1) / 2, (m + 1) * SQ3 / 2], -th / 2)
    S2 = _rot(S1, math.pi / 3)
    Sinv = np.linalg.inv(np.column_stack([S1, S2]))
    r = 0.5 * (1 + p.breath)
    subs = np.array([[0, 0], [r, 0], [r / 2, r * SQ3 / 2]])
    R = int(math.ceil(np.hypot(*S1))) * 2 + 3
    ii, jj = np.meshgrid(np.arange(-R, R + 1), np.arange(-R, R + 1), indexing="ij")
    base = np.stack([ii + jj / 2, jj * SQ3 / 2], -1).reshape(-1, 2)
    xy, lay, sub = [], [], []
    for l, ang in enumerate((-th / 2, th / 2)):
        shift = np.array([p.sx, p.sy]) if l else np.zeros(2)
        c, s = math.cos(ang), math.sin(ang)
        Rm = np.array([[c, -s], [s, c]])
        o = (base + shift) @ Rm.T
        f = o @ Sinv.T
        keep = np.all((f >= -1e-9) & (f < 1 - 1e-9), axis=1)
        cells = o[keep]
        if len(cells) != N:
            raise RuntimeError(f"layer {l}: {len(cells)} cells, expected {N}")
        for q in range(3):
            xy.append(cells + subs[q] @ Rm.T)
            lay.append(np.full(N, l))
            sub.append(np.full(N, q))
    xy = np.concatenate(xy)
    lay = np.concatenate(lay)
    sub = np.concatenate(sub)
    n = len(lay)

    shifts = np.array([u * S1 + v * S2 for u in (-1, 0, 1) for v in (-1, 0, 1)])
    d = xy[None, :, None, :] + shifts[None, None, :, :] - xy[:, None, None, :]  # (i, j, s, 2)
    dist = np.linalg.norm(d, axis=-1)
    same = lay[:, None] == lay[None, :]
    cyc = ((sub[:, None] == 0) & (sub[None, :] == 1)) | ((sub[:, None] == 1) & (sub[None, :] == 2)) | (
        (sub[:, None] == 2) & (sub[None, :] == 0))
    tu, td = p.t * (1 - p.breath), p.t * (1 + p.breath)
    ph = math.radians(p.flux) / 3
    I, J, A, D = [], [], [], []
    # intralayer, oriented A->B->C->A, both triangle types
    lengths = ((r, tu),) if abs(r - 0.5) < 1e-9 else ((r, tu), (1 - r, td))  # equal lengths: count once
    for length, amp in lengths:
        mask = same[:, :, None] & cyc[:, :, None] & (np.abs(dist - length) < 1e-6)
        i, j, s = np.nonzero(mask)
        I.append(i); J.append(j); D.append(d[i, j, s]); A.append(np.full(len(i), -amp * np.exp(1j * ph)))
    rc = min(p.rc, 4 * p.xi + 0.05)
    if p.tp != 0:
        mask = (lay[:, None, None] == 0) & (lay[None, :, None] == 1) & (dist < rc)
        i, j, s = np.nonzero(mask)
        I.append(i); J.append(j); D.append(d[i, j, s])
        A.append((p.tp * np.exp(-dist[i, j, s] ** 2 / (2 * p.xi ** 2))).astype(complex))
    return Supercell(p, N, th, S1, S2, xy, lay, sub,
                     np.concatenate(I), np.concatenate(J), np.concatenate(A), np.concatenate(D))


def moire_hk(sc: Supercell, k) -> np.ndarray:
    n = sc.n
    H = np.zeros((n, n), complex)
    val = sc.hop_amp * np.exp(1j * (sc.hop_d @ np.asarray(k, float)))
    np.add.at(H, (sc.hop_j, sc.hop_i), val)
    np.add.at(H, (sc.hop_i, sc.hop_j), np.conj(val))
    H[np.diag_indices(n)] += np.where(sc.layer == 1, sc.p.V / 2, -sc.p.V / 2)
    return H


def moire_bands(sc: Supercell, ks) -> np.ndarray:
    return np.array([np.linalg.eigvalsh(moire_hk(sc, k)) for k in ks])


def moire_path(sc: Supercell, nseg: int = 14):
    hs = sc.high_symmetry()
    pts = [hs["G"], hs["K"], hs["M"], hs["G"]]
    ks = []
    for a in range(3):
        for s in range(nseg + (a == 2)):
            u = s / nseg
            ks.append(pts[a] + u * (pts[a + 1] - pts[a]))
    return np.array(ks)


def moire_grid(sc: Supercell, nk: int):
    b1, b2 = sc.reciprocal()
    return np.array([((i + .5) / nk) * b1 + ((j + .5) / nk) * b2 for i in range(nk) for j in range(nk)])


def registry(sc: Supercell) -> np.ndarray:
    """AA-ness of each site: moire envelope of the two layer densities, in [-1.5, 3]."""
    th = sc.theta
    b = [np.array([TAU, -TAU / SQ3]), np.array([0, 2 * TAU / SQ3])]
    b.append(-b[0] - b[1])
    G0 = [_rot(v, -th / 2) for v in b]
    G1 = [_rot(v, th / 2) for v in b]
    sh = _rot([sc.p.sx, sc.p.sy], th / 2)
    f = np.zeros(sc.n)
    for g0, g1 in zip(G0, G1):
        f += np.cos(sc.xy @ (g0 - g1) + g1 @ sh)
    return f


def manifold_stats(sc: Supercell, E: np.ndarray, tol: float = 1e-3) -> dict:
    """E: (nk, n) sorted energies. Top-band manifold = top 2N bands."""
    lo, hi = E.min(0), E.max(0)
    n, N = sc.n, sc.N
    m0 = n - 2 * N
    widths = hi[m0:] - lo[m0:]
    gaps = lo[m0 + 1:] - hi[m0:-1]
    return {
        "manifold_width": float(hi[-1] - lo[m0]),
        "narrowest_width": float(widths.min()),
        "narrowest_band": int(m0 + widths.argmin()),
        "pieces": int(1 + np.sum(gaps > tol * sc.p.t)),
        "gap_below": float(lo[m0] - hi[m0 - 1]),
    }


def isolated_band_scores(E: np.ndarray) -> list[dict]:
    """For every band: width, gap to the nearest neighbour, flatness ratio gap/width."""
    lo, hi = E.min(0), E.max(0)
    out = []
    for b in range(E.shape[1]):
        below = lo[b] - hi[b - 1] if b > 0 else np.inf
        above = lo[b + 1] - hi[b] if b < E.shape[1] - 1 else np.inf
        gap = min(below, above)
        w = hi[b] - lo[b]
        out.append({"band": b, "width": float(w), "gap": float(gap), "ratio": float(gap / max(w, 1e-12)), "E": float((lo[b] + hi[b]) / 2)})
    return out


def ldos_window(sc: Supercell, k, E0: float, win: float):
    w, v = np.linalg.eigh(moire_hk(sc, k))
    sel = np.abs(w - E0) <= win
    wt = (np.abs(v[:, sel]) ** 2).sum(1)
    return wt, int(sel.sum())


def localisation_metrics(sc: Supercell, wt: np.ndarray, reg: np.ndarray | None = None) -> dict:
    if reg is None:
        reg = registry(sc)
    s = wt.sum()
    if s <= 0:
        return {"participation": float("nan"), "aa_enhancement": float("nan"), "layer1_share": float("nan")}
    aa = reg > 1.5
    return {
        "participation": float(s * s / (sc.n * (wt ** 2).sum())),
        "aa_enhancement": float((wt[aa].sum() / s) / max(aa.mean(), 1e-12)),
        "layer1_share": float(wt[sc.layer == 1].sum() / s),
    }


# --------------------------------------------------------------------------- Weyl block
@dataclass
class WeylParams:
    t: float = 1.0
    breath: float = 0.0
    flux: float = 90.0   # degrees per triangle
    mA: float = 0.0
    tz: float = 0.0
    dz: float = 0.8


def weyl_h(p: WeylParams, k1, k2, kz) -> np.ndarray:
    """Vectorised 3x3 Bloch Hamiltonian; k1 = k.a1, k2 = k.a2. Returns (..., 3, 3)."""
    k1, k2, kz = np.broadcast_arrays(np.asarray(k1, float), np.asarray(k2, float), np.asarray(kz, float))
    tu, td, f = p.t * (1 - p.breath), p.t * (1 + p.breath), math.radians(p.flux) / 3
    H = np.zeros(k1.shape + (3, 3), complex)
    H[..., 1, 0] = -tu * np.exp(1j * f) - td * np.exp(1j * (f + k1))
    H[..., 2, 1] = -tu * np.exp(1j * f) - td * np.exp(1j * (f - (k1 - k2)))
    H[..., 0, 2] = -tu * np.exp(1j * f) - td * np.exp(1j * (f - k2))
    H[..., 0, 1] = np.conj(H[..., 1, 0])
    H[..., 1, 2] = np.conj(H[..., 2, 1])
    H[..., 2, 0] = np.conj(H[..., 0, 2])
    c = np.cos(kz)
    H[..., 0, 0] = p.mA + 2 * (p.tz + p.dz) * c
    H[..., 1, 1] = 2 * (p.tz - p.dz / 2) * c
    H[..., 2, 2] = H[..., 1, 1]
    return H


def weyl_chern(p: WeylParams, kz: float, band: int, nk: int = 24):
    """Fukui-Hatsugai-Suzuki Chern number of one band on an nk x nk grid; also min gap of that band."""
    g = TAU * np.arange(nk) / nk
    K1, K2 = np.meshgrid(g, g, indexing="ij")
    w, v = np.linalg.eigh(weyl_h(p, K1, K2, kz))
    u = v[..., :, band]
    U1 = np.einsum("ijc,ijc->ij", u.conj(), np.roll(u, -1, 0))
    U2 = np.einsum("ijc,ijc->ij", u.conj(), np.roll(u, -1, 1))
    F = np.angle(U1 * np.roll(U2, -1, 0) * np.conj(np.roll(U1, -1, 1)) * np.conj(U2))
    gaps = []
    if band > 0:
        gaps.append((w[..., band] - w[..., band - 1]).min())
    if band < 2:
        gaps.append((w[..., band + 1] - w[..., band]).min())
    return F.sum() / TAU, float(min(gaps))


def weyl_gap_scan(p: WeylParams, lower_band: int, n: int = 48, refine: bool = True):
    """Brute-force 3D search for touchings between bands lower_band and lower_band+1.

    Returns a list of local minima (k1, k2, kz, gap), refined by coordinate pattern search.
    This is the ground truth the Chern-jump detector is compared against.
    """
    g = TAU * np.arange(n) / n - math.pi
    K1, K2, KZ = np.meshgrid(g, g, g, indexing="ij")
    w = np.linalg.eigvalsh(weyl_h(p, K1, K2, KZ))
    gap = w[..., lower_band + 1] - w[..., lower_band]
    is_min = np.ones_like(gap, bool)
    for ax in range(3):
        for sh in (1, -1):
            is_min &= gap <= np.roll(gap, sh, ax)
    idx = np.argwhere(is_min & (gap < 0.25))
    found = []
    for a, b, c in idx:
        x = np.array([g[a], g[b], g[c]])
        fx = gap[a, b, c]
        if refine:
            st = TAU / n
            while st > 1e-7:
                moved = False
                for dim in range(3):
                    for sgn in (1, -1):
                        y = x.copy(); y[dim] += sgn * st
                        wy = np.linalg.eigvalsh(weyl_h(p, *y))
                        fy = wy[lower_band + 1] - wy[lower_band]
                        if fy < fx:
                            x, fx, moved = y, fy, True
                if not moved:
                    st /= 2
        found.append((*[(v + math.pi) % TAU - math.pi for v in x], float(fx)))
    # merge duplicates (periodic distance)
    uniq = []
    for f in sorted(found, key=lambda z: z[3]):
        if all(np.linalg.norm([(f[i] - u[i] + math.pi) % TAU - math.pi for i in range(3)]) > 0.05 for u in uniq):
            uniq.append(f)
    return uniq


def weyl_detect_like_js(p: WeylParams, gap_sel: str = "12", nk: int = 18, nz: int = 48, gap_tol: float = 0.01):
    """Re-implementation of the page's detector: C(kz) on nz+1 planes, transitions skipped where gapped < tol."""
    band = 0 if gap_sel == "12" else 2
    sg = 1 if band == 0 else -1
    kz = -math.pi + TAU * np.arange(nz + 1) / nz
    res = [weyl_chern(p, z, band, nk) for z in kz]
    C = [round(sg * c) for c, _ in res]
    G = [g for _, g in res]
    transitions = []
    for q in range(nz):
        if C[q] != C[q + 1] and G[q] >= gap_tol and G[q + 1] >= gap_tol:
            transitions.append({"kz_lo": kz[q], "kz_hi": kz[q + 1], "dC": C[q + 1] - C[q]})
    n_nodes = sum(abs(t["dC"]) for t in transitions)
    return {"kz": kz, "C": C, "gap": G, "transitions": transitions, "n_nodes": n_nodes,
            "ungapped_planes": int(sum(g < gap_tol for g in G))}


def weyl_berry_charge(p: WeylParams, k1: float, k2: float, kz: float, band: int, r: float = 0.03,
                      nt: int = 10, nph: int = 20) -> int:
    """Berry flux of one band out of a small sphere around a point, in units of 2*pi."""
    th = math.pi * np.arange(nt + 1) / nt
    ph = TAU * np.arange(nph) / nph
    T, F = np.meshgrid(th, ph, indexing="ij")
    _, v = np.linalg.eigh(weyl_h(p, k1 + r * np.sin(T) * np.cos(F), k2 + r * np.sin(T) * np.sin(F), kz + r * np.cos(T)))
    u = v[..., :, band]
    U1 = np.einsum("ijc,ijc->ij", u[:-1].conj(), u[1:])                       # theta links
    U2 = np.einsum("ijc,ijc->ij", u.conj(), np.roll(u, -1, 1))                 # phi links
    flux = np.angle(U1 * U2[1:] * np.conj(np.roll(U1, -1, 1)) * np.conj(U2[:-1]))
    return int(round(flux.sum() / TAU))


def weyl_detect_v2(p: WeylParams, gap_sel: str = "12", nk: int = 18, nz: int = 48, gap_tol: float = 0.01,
                   node_tol: float = 1e-3):
    """Mirror of the page's fixed detector (after the AI-Scientist-v3 red team).

    Search windows: Chern jumps between gapped planes, and runs of nearly gapless planes (which the
    old detector skipped). In each window, gap minima on a 36x36 grid are refined in 3D; a point
    counts only if its gap falls below node_tol*t, and it is a Weyl node only if its Berry charge
    is nonzero. Jumps with no converged touching are reported as Chern-grid flicker.
    """
    band = 0 if gap_sel == "12" else 2
    gl = 0 if band == 0 else 1
    old = weyl_detect_like_js(p, gap_sel, nk, nz, gap_tol)
    kz, C, G = old["kz"], old["C"], old["gap"]
    sg = 1 if band == 0 else -1

    def gap(x):
        w = np.linalg.eigvalsh(weyl_h(p, *x))
        return w[gl + 1] - w[gl]

    windows = []
    for q in range(nz):
        if C[q] != C[q + 1] and G[q] >= gap_tol and G[q + 1] >= gap_tol:
            lo, hi = kz[q], kz[q + 1]
            for _ in range(7):
                mid = (lo + hi) / 2
                if round(sg * weyl_chern(p, mid, band, nk)[0]) == C[q]:
                    lo = mid
                else:
                    hi = mid
            windows.append((kz[q], kz[q + 1], [(lo + hi) / 2], True))
    q = 0
    while q <= nz:
        if G[q] < gap_tol:
            r = q
            while r < nz and G[r + 1] < gap_tol:
                r += 1
            lo, hi = kz[max(q - 1, 0)], kz[min(r + 1, nz)]
            windows.append((lo, hi, list(np.arange(lo, hi + 1e-9, TAU / nz / 2)), False))
            q = r
        q += 1
    M = 36
    g = TAU * np.arange(M) / M
    K1, K2 = np.meshgrid(g, g, indexing="ij")
    nodes, neutral, found, flicker, extended = [], [], [], 0, False
    for lo, hi, slices, jump in windows:
        cand = []
        for kc in slices:
            w = np.linalg.eigvalsh(weyl_h(p, K1, K2, kc))
            gg = w[..., gl + 1] - w[..., gl]
            is_min = np.ones_like(gg, bool)
            for a in (-1, 0, 1):
                for b in (-1, 0, 1):
                    if a or b:
                        is_min &= gg <= np.roll(np.roll(gg, a, 0), b, 1)
            for i, j in np.argwhere(is_min & (gg <= 0.08)):
                cand.append((gg[i, j], g[i], g[j], kc))
        cand = sorted(cand)[: 6 if jump else 24]
        accepted = 0
        for _, a1, a2, kc in cand:
            x, fx, st = np.array([a1, a2, kc]), None, 0.15
            fx = gap(x)
            for _ in range(400):
                if st <= 1e-6:
                    break
                moved = False
                for d in range(3):
                    for sgn in (1, -1):
                        y = x.copy(); y[d] += sgn * st
                        if d == 2:
                            y[2] = min(max(y[2], lo - .05), hi + .05)
                        fy = gap(y)
                        if fy < fx:
                            x, fx, moved = y, fy, True
                if not moved:
                    st /= 2
            if fx > node_tol * p.t:
                continue
            x = (x + math.pi) % TAU - math.pi
            if any(np.linalg.norm((f - x + math.pi) % TAU - math.pi) < .05 for f in found):
                continue
            found.append(x); accepted += 1
            chi = weyl_berry_charge(p, *x, gl)
            (nodes if chi else neutral).append({"k": x.tolist(), "gap": float(fx), "chi": chi})
        if jump and not accepted:
            flicker += 1
        if not jump and accepted > 6:
            extended = True
    return {"n_nodes": len(nodes), "nodes": nodes, "neutral": neutral, "flicker": flicker, "extended": extended}
