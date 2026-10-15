"""Retarded per-star Newtonian potential and the carved height field.

docs/design.md sec 3. Each star is treated as a point mass on its exact circular
orbit; the field at (x, y, Z_EVAL_MM) is the difference between the
self-consistently retarded potential and the same-instant (unretarded)
potential. Near each star this difference vanishes (light-crossing time
across a near-zone distance is negligible against the orbital period, so
retarded ~ instantaneous); far away it captures the propagation-delay
("radiative") content, including a trailing spiral, without ever collapsing
the two-body source to a point at the origin the way a point-quadrupole
model would.
"""
import numpy as np

import params as P


def _star_geometry(b):
    """(distance from origin, orbital angle at t0=0) for each star.

    Matches binary.solve()'s placement: heavy star on the -x axis (angle
    pi), light star on the +x axis (angle 0).
    """
    d_heavy, theta_heavy0 = abs(b.x1_mm), np.pi
    d_light, theta_light0 = b.x2_mm, 0.0
    m_heavy_kg = P.M1_MSUN * P.M_SUN
    m_light_kg = P.M2_MSUN * P.M_SUN
    return (d_heavy, theta_heavy0, m_heavy_kg), (d_light, theta_light0, m_light_kg)


def mm_to_m(r_mm, S_mm_per_km):
    """Convert a scaled sculpture-mm distance back to physical metres."""
    return r_mm * 1000.0 / S_mm_per_km


def solve_retarded_r(x, y, z_eval, d_star, theta0, k_half, n_iter=P.RETARDED_ITERS):
    """Self-consistent retarded distance from (x, y, z_eval) to a star on a
    circular orbit of radius d_star, at angle theta0 at t0 = 0.

    Fixed-point iteration on phi = k_half * r (the retardation phase, with
    k_half = Omega/c at scale): the star's retarded angular position is
    theta0 - phi (it was further back along +phi_hat at the earlier time
    t0 - phi/Omega). Contracts geometrically since d(k_half*r)/dphi ~
    k_half*d_star << 1 here; a fixed iteration count is used rather than a
    convergence check since the contraction factor is the same everywhere
    on the grid (docs/design.md sec 3).
    """
    xi = d_star * np.cos(theta0)
    yi = d_star * np.sin(theta0)
    r = np.sqrt((x - xi) ** 2 + (y - yi) ** 2 + z_eval ** 2)
    phi = k_half * r
    for _ in range(n_iter):
        xi = d_star * np.cos(theta0 - phi)
        yi = d_star * np.sin(theta0 - phi)
        r = np.sqrt((x - xi) ** 2 + (y - yi) ** 2 + z_eval ** 2)
        phi = k_half * r
    return r


def instantaneous_r(x, y, z_eval, d_star, theta0):
    """Distance to a star at its t0 = 0 position -- no retardation."""
    xi = d_star * np.cos(theta0)
    yi = d_star * np.sin(theta0)
    return np.sqrt((x - xi) ** 2 + (y - yi) ** 2 + z_eval ** 2)


def delta_phi(x, y, z_eval, b, k_half):
    """(retarded - instantaneous) Newtonian potential, summed over both
    stars, in SI (m^2/s^2). This isolates the propagation-delay content:
    it vanishes where retardation is negligible (near each star) and
    approaches the retarded potential itself where the instantaneous
    quadrupole-like piece has decayed away (far zone)."""
    (d_h, th_h, m_h), (d_l, th_l, m_l) = _star_geometry(b)
    S = b.S_mm_per_km

    total_ret = 0.0
    total_inst = 0.0
    for d_star, theta0, m_kg in [(d_h, th_h, m_h), (d_l, th_l, m_l)]:
        r_ret = solve_retarded_r(x, y, z_eval, d_star, theta0, k_half)
        r_inst = instantaneous_r(x, y, z_eval, d_star, theta0)
        total_ret += -P.G * m_kg / mm_to_m(r_ret, S)
        total_inst += -P.G * m_kg / mm_to_m(r_inst, S)
    return total_ret - total_inst


def height_field(x, y, b, c_eff=P.C_EFF):
    """z(x, y) in mm, z in [-A_pp, 0]. docs/design.md sec 4.3.

    No edge taper: the field is carved all the way to the cut edge of the
    panel.

    c_eff: the wave speed used for retardation (sec 3.5). Defaults to
    P.C_EFF (= C_EFF_FACTOR * C; pass c_eff=params.C to recover the real
    speed of light). Changes only the retardation phase
    (binary.retardation_k_half), not the plate/star geometry, which stays
    fixed at whatever scale the real C produced.

    The normalization N is the maximum |u_raw| over the points passed in, so z is only the
    plate's relief when (x, y) covers the whole plate (plate.sample_grid). To read the relief
    at a point, interpolate the plate grid (plate.relief_at) instead.
    """
    import binary

    k_half = binary.retardation_k_half(b, c_eff)
    u_raw = delta_phi(x, y, P.Z_EVAL_MM, b, k_half) / P.C ** 2

    N = np.max(np.abs(u_raw))
    u = np.clip(u_raw / N, -1.0, 1.0)
    return -P.A_PP_MM / 2 * (1.0 - u)


if __name__ == "__main__":
    import binary
    import plate

    b = binary.solve()
    X, Y, _, _ = plate.sample_grid()
    edge = np.zeros(X.shape, bool)
    edge[0, :] = edge[-1, :] = edge[:, 0] = edge[:, -1] = True
    for label, c_eff in [("real c", P.C), (f"{P.C_EFF_FACTOR:g} c", P.C_EFF)]:
        k_half = binary.retardation_k_half(b, c_eff)
        u_raw = delta_phi(X, Y, P.Z_EVAL_MM, b, k_half) / P.C ** 2
        N = np.max(np.abs(u_raw))
        print(f"{label}: k_half = {k_half:.5f} /mm  peak|u_raw| = {N:.4e}  "
              f"edge mean|u|/peak = {np.mean(np.abs(u_raw[edge])) / N:.3f}")
