"""Checks of the equations in docs/math.md against the production data and code.

verification/math_walkthrough.py re-implements the chain orbit -> amplitudes -> surface in float64 with
numpy/scipy only; these tests tie it to data/default/star_modes.toml (written by tools/make_star_modes.py) and to
src/tidal_shape.py.
"""
import tomllib
from pathlib import Path

import numpy as np
import pytest
from scipy.special import sph_harm_y

import math_walkthrough as MW
import polytrope as poly
import tidal_shape as TS

ROOT = Path(__file__).resolve().parent.parent
RAW = tomllib.load(open(ROOT / "data" / "default" / "star_modes.toml", "rb"))


@pytest.fixture(scope="module", params=["heavy", "light"])
def star(request):
    return next(s for s in RAW["star"] if s["name"] == request.param)


@pytest.fixture(scope="module")
def chains():
    return {s["name"]: MW.chain(s) for s in RAW["star"]}


def unit_vectors(n=300, seed=3):
    v = np.random.default_rng(seed).normal(size=(n, 3))
    return v / np.linalg.norm(v, axis=1, keepdims=True)


def test_independent_chain_reproduces_the_generated_amplitudes(chains):
    for name, c in chains.items():
        for r in c["modes"]:
            rel = abs(r["C"] - r["C_toml"]) / abs(r["C_toml"])
            # non-resonant modes agree to float64 round-off; the two prograde modes differ only by the
            # code's trapezoid-rule accumulation of the dynamical-tide integral F(u)
            limit = 1e-3 if r["m"] == r["l"] else 1e-12
            assert rel < limit, (name, r["l"], r["m"], rel)


def test_dynamical_tide_is_small_but_nonzero_at_the_snapshot(chains):
    for name, c in chains.items():
        for r in c["modes"]:
            if r["m"] == r["l"]:
                ratio = abs(r["b_dyn"]) / abs(r["C"])
                assert 1e-3 < ratio < 2e-2, (name, r["l"], ratio)       # 0.3 - 1.2 % of |C|
            else:
                assert r["b_dyn"] == 0


def test_resummation_reproduces_eq18_in_the_point_mass_limit():
    # b0 = w V / Delta, b1 = i db0/dt / Delta, b_eq = b0 / (1 - b1/b0)  =>  b_eq = w Delta V / (Delta^2 - i kappa)
    l, m, w_a0, Om_s = 2, 2, 1.0, 0.0
    omega, omega_dot = 0.30, 1e-3
    V = lambda w: 0.01 * w ** (2.0 * (l + 1) / 3.0)                  # V ~ r^-(l+1), r ~ omega^(-2/3)
    Delta = lambda w: w_a0 + m * Om_s - m * w
    b0 = lambda w: w_a0 * V(w) / Delta(w)
    h = 1e-6
    db0_dt = omega_dot * (b0(omega + h) - b0(omega - h)) / (2 * h)
    b1 = 1j * db0_dt / Delta(omega)
    b_resummed = b0(omega) / (1 - b1 / b0(omega))
    dlogr = -(2.0 / 3.0) * omega_dot / omega
    kappa = m * omega_dot - Delta(omega) * (l + 1) * dlogr
    b_eq18 = Delta(omega) * w_a0 * V(omega) / (Delta(omega) ** 2 - 1j * kappa)
    assert b_resummed == pytest.approx(b_eq18, rel=1e-8)
    assert np.angle(b_eq18) == pytest.approx(np.arctan(kappa / Delta(omega) ** 2), rel=1e-10)


def test_surface_equals_equilibrium_tide_times_transfer_function(chains):
    """dR/R = h_l (M'/M)(R/r)^(l+1) Re sum_a T_a W_lm Y_lm(n), T_a = C_a / V_a: the overlap I cancels."""
    n = unit_vectors()
    theta, phi = np.arccos(n[:, 2]), np.arctan2(n[:, 1], n[:, 0])
    for star_raw in RAW["star"]:
        c = chains[star_raw["name"]]
        ratio = star_raw["Mc_msun"] / star_raw["M_msun"]
        modes = [TS.ModeAmp(r["l"], r["m"], r["C"], r["I"], r["h"]) for r in c["modes"]]   # chain amplitudes
        lhs = TS.delta_r_over_R(n, modes)
        rhs = np.zeros(len(n))
        for r in c["modes"]:
            T = r["C"] / r["V"]
            rhs += r["h"] * ratio * c["r_over_R"] ** (-(r["l"] + 1)) * np.real(
                T * MW.W(r["l"], r["m"]) * sph_harm_y(r["l"], r["m"], theta, phi))
        np.testing.assert_allclose(lhs, rhs, rtol=0, atol=1e-13)


def test_l2_closed_form_in_transfer_functions(chains):
    """e_xx = h2 eps_t [3/8 T_c + 1/4 Re T_0], e_yy = h2 eps_t [-3/8 T_c + 1/4 Re T_0], e_zz = -1/2 h2 eps_t Re T_0,
    e_xy = -3/8 h2 eps_t T_s, with T_c + i T_s = T_{+2} + conj(T_{-2}) and eps_t = (M'/M)(R/r)^3."""
    for star_raw in RAW["star"]:
        c = chains[star_raw["name"]]
        r = {(x["l"], x["m"]): x for x in c["modes"]}
        eps_t = star_raw["Mc_msun"] / star_raw["M_msun"] * c["r_over_R"] ** -3
        h2 = r[(2, 2)]["h"]
        T = {k: r[k]["C"] / r[k]["V"] for k in [(2, 2), (2, 0), (2, -2)]}
        Tcs = T[(2, 2)] + np.conj(T[(2, -2)])
        q = TS.quadrupole_form(r[(2, 2)]["C"], r[(2, -2)]["C"], r[(2, 0)]["C"], r[(2, 2)]["I"], h2)
        assert q["e_xx"] == pytest.approx(h2 * eps_t * (3 / 8 * Tcs.real + 0.25 * T[(2, 0)].real), rel=1e-12)
        assert q["e_yy"] == pytest.approx(h2 * eps_t * (-3 / 8 * Tcs.real + 0.25 * T[(2, 0)].real), rel=1e-12)
        assert q["e_zz"] == pytest.approx(h2 * eps_t * (-0.5 * T[(2, 0)].real), rel=1e-12)
        assert q["e_xy"] == pytest.approx(-3 / 8 * h2 * eps_t * Tcs.imag, rel=1e-12)


def test_W_N_products():
    # W_22 N_2 = 3/8 and W_20 N_0 = -1/4 make the adiabatic l=2 surface exactly P_2(n_x)
    assert MW.W(2, 2) * TS.N2 == pytest.approx(3 / 8, rel=1e-14)
    assert MW.W(2, 0) * TS.N0 == pytest.approx(-1 / 4, rel=1e-14)


def test_resonance_root_and_closed_form(star):
    M1, M2 = star["M_msun"] * MW.MS, star["Mc_msun"] * MW.MS
    X1, X2 = M1 / (M1 + M2), M2 / (M1 + M2)
    tMt = MW.G * (M1 + M2) / MW.C_LIGHT ** 3
    Om_s = 2 * np.pi * star["spin_hz"]
    for md in star["mode"]:
        l, m = md["l"], md["m"]
        if m != l:
            continue
        A = star["A_nl_l2"] if l == 2 else 0.0
        w_res = MW.resonance(l, m, md["omega_a0_rad_s"], Om_s, A, X1, X2, tMt)
        D, _ = MW.delta_a(w_res, l, m, md["omega_a0_rad_s"], Om_s, A, X1, X2, tMt)
        assert abs(D) / md["omega_a0_rad_s"] < 1e-12                       # exact root of Delta_a = 0
        if l == 2:                                                         # get_ba_l2's 1PN closed form
            closed = MW.resonance(l, m, md["omega_a0_rad_s"], Om_s, A, X1, X2, tMt, closed_form_l2=True)
            assert abs(closed - w_res) / w_res < 0.03


def test_v11_margin_uses_newtonian_delta_but_conclusion_survives_gr(star):
    M1, M2 = star["M_msun"] * MW.MS, star["Mc_msun"] * MW.MS
    X1, X2 = M1 / (M1 + M2), M2 / (M1 + M2)
    orb = MW.orbit(M1, M2, np.pi * star["f_gw_hz"])
    Om_s = 2 * np.pi * star["spin_hz"]
    for md in star["mode"]:
        l, m, w0 = md["l"], md["m"], md["omega_a0_rad_s"]
        if m <= 0:
            continue
        wcor = w0 - (m / l) * Om_s
        newt = (wcor + m * Om_s - m * np.pi * star["f_gw_hz"]) / wcor          # V11 definition (Newtonian)
        D_gr, _ = MW.delta_a(np.pi * star["f_gw_hz"], l, m, w0, Om_s, star["A_nl_l2"] if l == 2 else 0.0, X1, X2, orb["tMt"])
        assert np.sign(D_gr) == np.sign(newt)            # never flips the sign in the final configuration
        if m == l:                                       # the resonant (prograde) modes: margin shifts by < 0.06
            assert abs(D_gr / wcor - newt) < 0.06


def test_polytrope_shape_numbers():
    k_love2, k2, h2 = poly.static_love(2)
    assert h2 == pytest.approx(1.51982, abs=1e-5) and k2 == pytest.approx((15 - np.pi ** 2) / (2 * np.pi ** 2), rel=1e-12)
    assert poly.static_love(3)[2] == pytest.approx(1.21291, abs=1e-5)
    for l in (2, 3):                                     # mode sum rule: (4 pi/(2l+1)) sum I^2 = k_l
        k_exact = poly.static_love(l)[1]
        ks = sum(4 * np.pi / (2 * l + 1) * I ** 2 for (_, I, _, _) in poly.modes(l, 8))
        assert ks == pytest.approx(k_exact, abs=2e-5)
