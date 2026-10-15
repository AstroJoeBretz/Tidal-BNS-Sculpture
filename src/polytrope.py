"""n=1 (Gamma=2) Newtonian polytrope: f/p-mode overlaps, Love numbers, shape numbers h_l.

The stellar model behind the sculpture stars (docs/math.md sec 5): `static_love` and `modes`
give the dimensionless $k_l$, $h_l$, $I_l$ and $\\omega_{a0}/\\omega_0$ that
`tools/make_star_modes.py` scales to each star. Units G=M=R=1.

For N=0 (neutrally stratified) and Gamma = polytropic index the displacement is a
potential flow, xi = grad(chi)/w^2 with chi = 2K rho' + Phi'. Two coupled radial ODEs
for (chi, Phi') are solved by Chebyshev collocation; modes are normalized as in
Yu et al. 2022 Eq. 3:  2 w^2 int rho |xi|^2 d^3x = E0 = M^2/R.

Static (w -> 0) tide is also solved exactly: psi = Phi' + Phi_tid obeys a Helmholtz
equation, psi = A j_l(pi r/R) Y_lm, matched to the exterior -> k_l^Love, h_l = 1 + k_l^Love.
Paper convention: k_l = k_l^Love / 2  and  k_l = 4 pi I^2 / (2l+1)  (Yu+2024 Eq. 43).

Usage: python src/polytrope.py [l ...]      (default: 2 3) prints the mode table and sum rules.
"""
import sys

import numpy as np
from numpy.polynomial import chebyshev as Ch
from scipy.linalg import eig
from scipy.special import spherical_jn


def static_love(l):
    """Exact n=1 static tidal response. Returns (k_Love, k_paper, h)."""
    j = spherical_jn(l, np.pi)
    jp = spherical_jn(l, np.pi, derivative=True)
    y = np.pi * jp / j
    k_love = (l - y) / (l + 1 + y)
    return k_love, k_love / 2, 1 + k_love


def modes(l, nmodes=8, N=61):
    s = (-1) ** l                                     # parity of chi(r) ~ r^l
    x = np.cos(np.pi * np.arange(N + 1) / N)
    c = np.ones(N + 1); c[0] = c[-1] = 2; c *= (-1) ** np.arange(N + 1)
    X = np.tile(x, (N + 1, 1)).T
    D = np.outer(c, 1 / c) / (X - X.T + np.eye(N + 1))
    D -= np.diag(D.sum(axis=1))
    pos = np.arange(0, (N + 1) // 2); neg = N - pos   # N odd: no node at r=0
    D1 = D[np.ix_(pos, pos)] + s * D[np.ix_(pos, neg)]
    DD = D @ D
    D2 = DD[np.ix_(pos, pos)] + s * DD[np.ix_(pos, neg)]
    r = x[pos]; n = len(r)

    rho_c = np.pi / 4
    rho = rho_c * np.sin(np.pi * r) / (np.pi * r)
    drho = rho_c * (np.pi * r * np.pi * np.cos(np.pi * r) - np.pi * np.sin(np.pi * r)) / (np.pi * r) ** 2
    twoK = 4.0 / np.pi
    ll = l * (l + 1)
    In = np.eye(n)
    Lap = D2 + 2 * np.diag(1 / r) @ D1 - ll * np.diag(1 / r ** 2)

    A = np.zeros((2 * n, 2 * n)); B = np.zeros((2 * n, 2 * n))
    A[:n, :n] = -(np.diag(rho) @ Lap + np.diag(drho) @ D1)
    B[:n, :n] = In / twoK; B[:n, n:] = -In / twoK
    A[n:, :n] = -4 * np.pi * In / twoK
    A[n:, n:] = Lap + 4 * np.pi * In / twoK
    A[n, :] = 0; B[n, :] = 0                           # surface BC (r=R): Phi'_r + (l+1) Phi'/R = 0 (exterior)
    A[n, n:] = D1[0, :] + (l + 1) * In[0, :]

    w, V = eig(A, B)
    ok = np.isfinite(w) & (np.abs(w.imag) < 1e-8) & (w.real > 0.1)
    idx = np.where(ok)[0]; idx = idx[np.argsort(w[idx].real)]

    def interp(f, rr, par):
        full = np.empty(N + 1); full[pos] = f; full[neg] = par * f
        return Ch.chebval(rr, Ch.chebfit(x, full, N))

    rr = np.linspace(1e-6, 1, 40001)
    rh = rho_c * np.sin(np.pi * rr) / (np.pi * rr)
    out = []
    for k in idx[:nmodes]:
        w2 = w[k].real; chi = V[:n, k].real
        xr = interp(D1 @ chi / w2, rr, -s); xh = interp(chi / (w2 * r), rr, -s)   # xi_r, xi_h have parity -s
        sc = np.sqrt(1.0 / (2 * w2 * np.trapezoid(rh * (xr ** 2 + ll * xh ** 2) * rr ** 2, rr)))
        xr *= sc; xh *= sc
        I = np.trapezoid(rh * (xr * l * rr ** (l - 1) + xh * ll * rr ** (l - 1)) * rr ** 2, rr)
        out.append((np.sqrt(w2), I, xr[-1], xh[-1]))
    return out


def report(l):
    k_love, k_pap, h = static_love(l)
    print(f"\n=== l = {l} ===")
    print(f"exact static:  k_Love = {k_love:.5f}   k_paper = k_Love/2 = {k_pap:.5f}   h_l = 1 + k_Love = {h:.5f}")
    print(" n   w/w0       I        xi_r(R)/R   2 xi_r I   cum     4pi/(2l+1) I^2   cum")
    cum_h = cum_k = 0.0
    for nn, (w, I, xr, xh) in enumerate(modes(l)):
        cum_h += 2 * xr * I
        cum_k += 4 * np.pi / (2 * l + 1) * I ** 2
        print(f"{nn:2d} {w:8.4f} {I:10.5f} {xr:10.4f} {2 * xr * I:10.4f} {cum_h:8.4f} {4 * np.pi / (2 * l + 1) * I ** 2:12.5f} {cum_k:8.5f}")
    f = modes(l, 1)[0]
    print(f"f-mode: I = {abs(f[1]):.4f}, xi_r(R)/R = {abs(f[2]):.4f}, calibration h/(2I) = {h / (2 * abs(f[1])):.4f}")
    print(f"sum rules: sum 4pi/(2l+1) I^2 = {cum_k:.5f} (k_paper {k_pap:.5f});  sum 2 xi_r I = {cum_h:.4f} (h_l {h:.4f})")


if __name__ == "__main__":
    for l in ([int(a) for a in sys.argv[1:]] or [2, 3]):
        report(l)
