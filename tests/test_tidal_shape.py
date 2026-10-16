"""Tests V1-V12 of docs/neutron_star_shape.md (V3 and V9 are checked by src/polytrope.py and
verification/verify_amplitudes.py)."""
import dataclasses
import subprocess
import sys
import warnings
from pathlib import Path

import numpy as np
import pytest
from scipy.special import eval_legendre, sph_harm_y

import binary
import params as P
import star as STAR
import tidal_shape as TS

ROOT = Path(__file__).resolve().parent.parent
MODES_TOML = ROOT / "data" / "default" / "star_modes.toml"
RNG = np.random.default_rng(12345)


def unit_vectors(n=200):
    v = RNG.normal(size=(n, 3))
    return v / np.linalg.norm(v, axis=1, keepdims=True)


def W(l, m):
    return 4.0 * np.pi / (2 * l + 1) * np.real(sph_harm_y(l, m, np.pi / 2, 0.0))


def adiabatic_modes(l, mass_ratio, r_over_R, I_abs, h, flip_I=False):
    """C_a = V_a = (M'/M) W_lm I_a (R/r)^(l+1) with I_{l,-m} = (-1)^m I_{l,m} (V1)."""
    out = []
    for m in range(-l, l + 1, 2):
        I = I_abs * ((-1) ** abs(m) if m < 0 else 1.0)
        if flip_I:
            I = -I
        V = mass_ratio * r_over_R ** (-(l + 1)) * W(l, m) * I
        out.append(TS.ModeAmp(l, m, complex(V), I, h))
    return out


def random_modes(l, scale=0.02):
    modes = []
    I_abs = RNG.uniform(0.1, 0.4)
    for m in range(-l, l + 1, 2):
        I = I_abs * ((-1) ** abs(m) if m < 0 else 1.0)
        modes.append(TS.ModeAmp(l, m, complex(*RNG.normal(scale=scale, size=2)), I, 1.3))
    return modes


@pytest.fixture(scope="module")
def models():
    assert MODES_TOML.exists(), "run  python tools/make_star_modes.py  first"
    return TS.load_star_modes(MODES_TOML)


# --------------------------------------------------------------------- V1 adiabatic limit

@pytest.mark.parametrize("l,flip", [(2, False), (2, True), (3, False), (3, True)])
def test_v1_adiabatic_limit(l, flip):
    ratio, r_over_R, I_abs, h = 1.37, 2.9, 0.27, 1.43
    n = unit_vectors()
    modes = adiabatic_modes(l, ratio, r_over_R, I_abs, h, flip_I=flip)
    got = TS.delta_r_over_R(n, modes)
    expected = h * ratio * r_over_R ** (-(l + 1)) * eval_legendre(l, n[:, 0])   # companion on +x
    np.testing.assert_allclose(got, expected, atol=1e-12, rtol=0)


def test_v1_l2_plus_l3_superpose():
    n = unit_vectors()
    m2 = adiabatic_modes(2, 0.8, 2.7, 0.3, 1.5)
    m3 = adiabatic_modes(3, 0.8, 2.7, 0.24, 1.2)
    both = TS.delta_r_over_R(n, m2 + m3)
    np.testing.assert_allclose(both, TS.delta_r_over_R(n, m2) + TS.delta_r_over_R(n, m3), atol=1e-14)


# --------------------------------------------------------------------- V2 no l=0 component

@pytest.mark.parametrize("X_lag", [1.0, 3.0, 0.0])
def test_v2_zero_mean(X_lag):
    modes = random_modes(2) + random_modes(3)
    n, w = TS.gauss_sphere(96, 192)
    mean = TS.sphere_mean(TS.delta_r_over_R(n, modes, X_lag), w, 192)
    assert abs(mean) < 1e-12


# --------------------------------------------------------------------- V4 symmetry

def test_v4_real_amplitudes_mirror_symmetry():
    modes = []
    for md in random_modes(2) + random_modes(3):
        modes.append(dataclasses.replace(md, C=complex(md.C.real, 0.0)))
    n = unit_vectors()
    f = TS.delta_r_over_R
    np.testing.assert_allclose(f(n, modes), f(n * [1, -1, 1], modes), atol=1e-13)   # y -> -y
    np.testing.assert_allclose(f(n, modes), f(n * [1, 1, -1], modes), atol=1e-13)   # z -> -z


def test_v4_l2_invariant_under_pi_rotation_about_z():
    n = unit_vectors()
    modes = random_modes(2)                                  # complex amplitudes: still true for any l=2 pattern
    np.testing.assert_allclose(TS.delta_r_over_R(n, modes), TS.delta_r_over_R(n * [-1, -1, 1], modes), atol=1e-13)


# --------------------------------------------------------------------- V5 closed form vs harmonics

def _q2(Ca, Cb, Cc, I):
    return [TS.ModeAmp(2, 2, Ca, I, 1.0), TS.ModeAmp(2, 0, Cc, I, 1.0), TS.ModeAmp(2, -2, Cb, I, 1.0)]


@pytest.mark.parametrize("seed", range(5))
def test_v5_quadratic_form_matches_harmonics(seed):
    rng = np.random.default_rng(seed)
    Ca, Cb, Cc = (complex(*rng.normal(scale=0.03, size=2)) for _ in range(3))
    I, h2 = rng.choice([-1, 1]) * rng.uniform(0.15, 0.4), rng.uniform(1.1, 1.8)
    modes = [dataclasses.replace(md, h=h2) for md in _q2(Ca, Cb, Cc, I)]
    q = TS.quadrupole_form(Ca, Cb, Cc, I, h2)
    n = unit_vectors(400)
    form = (q["e_xx"] * n[:, 0] ** 2 + q["e_yy"] * n[:, 1] ** 2 + q["e_zz"] * n[:, 2] ** 2
            + 2 * q["e_xy"] * n[:, 0] * n[:, 1])
    np.testing.assert_allclose(form, TS.delta_r_over_R(n, modes), rtol=0, atol=1e-13)
    assert q["e_xx"] + q["e_yy"] + q["e_zz"] == pytest.approx(0.0, abs=1e-14)       # traceless
    assert q["psi_deg"] == pytest.approx(TS.lag_angles(modes)[(2, 2)], abs=1e-10)


# --------------------------------------------------------------------- V6 principal axes

def test_v6_principal_axes_follow_lag_angle():
    rng = np.random.default_rng(7)
    Ca, Cb, Cc = (complex(*rng.normal(scale=0.03, size=2)) for _ in range(3))
    Ca, Cb = complex(0.05, 0.004), complex(0.012, 0.0006)             # a typical lagged l=2 pattern
    I, h2 = 0.32, 1.52
    q = TS.quadrupole_form(Ca, Cb, Cc, I, h2)
    A = np.array([[q["e_xx"], q["e_xy"], 0.0], [q["e_xy"], q["e_yy"], 0.0], [0.0, 0.0, q["e_zz"]]])
    w, v = np.linalg.eigh(A)
    # z is an exact eigenvector
    assert np.isclose(abs(A @ [0, 0, 1] @ [0, 0, 1]), abs(q["e_zz"]))
    assert np.allclose(A @ [0, 0, 1], q["e_zz"] * np.array([0, 0, 1]))
    # in-plane eigenvector of the larger in-plane eigenvalue sits at phi = -psi (mod pi)
    wi, vi = np.linalg.eigh(A[:2, :2])
    long_axis = np.degrees(np.arctan2(vi[1, -1], vi[0, -1])) % 180.0
    assert long_axis == pytest.approx((-q["psi_deg"]) % 180.0, abs=1e-8)
    # and the equatorial profile of the actual surface peaks there too
    phi = np.radians(np.linspace(0, 360, 360001)[:-1])
    n = np.stack([np.cos(phi), np.sin(phi), np.zeros_like(phi)], axis=1)
    modes = [dataclasses.replace(md, h=h2) for md in _q2(Ca, Cb, Cc, I)]
    prof = TS.delta_r_over_R(n, modes)
    phi_max = np.degrees(phi[np.argmax(prof)])
    assert (phi_max % 180.0) == pytest.approx((-q["psi_deg"]) % 180.0, abs=2e-3)


# --------------------------------------------------------------------- V7 meshes and gap

@pytest.fixture(scope="module")
def stars(models):
    b = binary.solve()
    out = {}
    for name, R_mm in (("heavy", b.R1_mm), ("light", b.R2_mm)):
        fn = models[name].radius_fn(R_mm)
        out[name] = (fn, STAR.build_star(fn), R_mm)
    return b, out


def test_v7_meshes_watertight_positive_radius(stars):
    b, s = stars
    for name, (fn, mesh, R_mm) in s.items():
        assert mesh.is_watertight and mesh.is_winding_consistent, name
        assert np.linalg.norm(mesh.vertices, axis=1).min() > 0.0
        assert fn(unit_vectors(2000)).min() > 0.0


def test_v7_gap_between_stars_positive(stars):
    from scipy.spatial import cKDTree
    b, s = stars
    vh = s["heavy"][1].vertices + [b.x1_mm, 0, 0]
    vl = s["light"][1].vertices + [b.x2_mm, 0, 0]
    d, _ = cKDTree(vh).query(vl)
    assert d.min() > 0.5                                       # mm on the sculpture (a km-scale gap of ~3 km = 2.4 mm)
    tip_h = s["heavy"][0](np.array([[1.0, 0, 0]]))[0]
    tip_l = s["light"][0](np.array([[-1.0, 0, 0]]))[0]
    line_gap = b.a_mm - tip_h - tip_l
    assert line_gap > 0.0
    # the 3-D gap can only be smaller than the line-of-centres gap; the vertex estimate exceeds the
    # true distance by at most the mesh spacing effect (~h^2/R), and here the closest approach IS the line
    assert d.min() <= line_gap + 0.01
    assert line_gap - d.min() < 0.05


def test_v7_mesh_volume_matches_volume_drift(models):
    b = binary.solve()
    for name, R_mm in (("heavy", b.R1_mm), ("light", b.R2_mm)):
        m = models[name]
        bare = STAR.build_star(m.radius_fn(R_mm), socket=False)
        drift_mesh = bare.volume / (4.0 / 3.0 * np.pi * R_mm ** 3) - 1.0
        drift = TS.volume_drift(m.rel_radius_star_frame)
        assert drift_mesh == pytest.approx(drift, abs=2e-3), name


def test_v7_socket_axis_is_the_pole(stars):
    b, s = stars
    for name, (fn, mesh, R_mm) in s.items():
        r_pole = STAR.polar_radius(fn)
        # the bore removes the pole point; the lowest material is its rim, within a fraction of a mm of
        # -r_pole (a sphere would put it rho^2/(2 r_pole) = 0.18 mm up; the pole region is flatter
        # than a sphere, and for the light star the rim is ~0.015 mm *below* the pole point)
        assert -r_pole - 0.05 < mesh.bounds[0][2] < -r_pole + 0.25, name
        low = mesh.vertices[mesh.vertices[:, 2] < mesh.bounds[0][2] + 0.01]
        assert np.hypot(low[:, 0], low[:, 1]).max() < P.SOCKET_DIAMETER_MM, name   # lowest material hugs the bore
        # the blind end of the bore is SOCKET_DEPTH_MM above the pole point and below the equator
        rho = np.hypot(mesh.vertices[:, 0], mesh.vertices[:, 1])
        bore = mesh.vertices[(rho < P.SOCKET_DIAMETER_MM / 2 + 1e-3) & (mesh.vertices[:, 2] < 3.0)]
        assert bore[:, 2].max() == pytest.approx(-r_pole + P.SOCKET_DEPTH_MM, abs=1e-3), name
        assert bore[:, 2].max() < 0.0, name


# --------------------------------------------------------------------- V8 lag sign (Paper 2 Eq. 18)

def _eq18_amplitudes(omega_dot, omega=0.3, wa0=1.0, I=0.3, ratio=0.01, l=2):
    """Paper 2 Eq. 18 for m=+-2 and m=0 along a point-mass inspiral (dimensionless units)."""
    dlogr = -(2.0 / 3.0) * omega_dot / omega
    out = {}
    for m in (2, 0, -2):
        V = ratio * W(2, m) * I
        D = wa0 - m * omega
        kappa = m * omega_dot - D * (l + 1) * dlogr
        out[m] = D * wa0 * V / (D ** 2 - 1j * kappa)
    return out, I


def test_v8_lag_sign_and_adiabatic_limit():
    psis = []
    for s in (0.0, 1e-4, 1e-3, 4e-3):
        C, I = _eq18_amplitudes(s)
        q = TS.quadrupole_form(C[2], C[-2], C[0], I, 1.5)
        psis.append(q["psi_deg"])
    assert psis[0] == pytest.approx(0.0, abs=1e-12)           # omega_dot -> 0: no lag
    assert all(p > 0 for p in psis[1:])                       # inspiral: bulge trails the companion
    assert psis[1] < psis[2] < psis[3]                        # and grows with the decay rate


# --------------------------------------------------------------------- V10 X_lag

def test_v10_xlag_one_is_identity_and_amplitude_preserved():
    modes = random_modes(2, 0.03) + random_modes(3, 0.01)    # arbitrary complex amplitudes (|D| check only)
    n = unit_vectors()
    np.testing.assert_array_equal(TS.delta_r_over_R(n, modes), TS.delta_r_over_R(n, modes, X_lag=1.0))
    D1, _ = TS.combined_amplitudes(modes, 1.0)
    for X_lag in (0.0, 2.0, 5.0):
        DX, _ = TS.combined_amplitudes(modes, X_lag)
        for k in D1:
            assert abs(DX[k]) == pytest.approx(abs(D1[k]), rel=1e-12)


def _lagged_adiabatic(a=0.05, b=0.02):
    """Physically signed adiabatic amplitudes (V_{-m} has the sign of V_m) with a small lag: C -> V (1 + i a)
    for the prograde modes and V (1 + i b) for the retrograde ones."""
    modes = []
    for md in adiabatic_modes(2, 0.8, 2.7, 0.3, 1.5) + adiabatic_modes(3, 0.8, 2.7, 0.24, 1.2):
        f = 1.0 + 1j * (a if md.m > 0 else b if md.m < 0 else 0.0)
        modes.append(dataclasses.replace(md, C=md.C * f))
    return modes


def test_v10_adiabatic_limit_has_zero_lag_for_every_harmonic():
    # includes (3,3), whose Y_33 is negative on the line of centres: the reference phase must be removed
    for l, psi in TS.lag_angles(adiabatic_modes(2, 0.8, 2.7, 0.3, 1.5) + adiabatic_modes(3, 0.8, 2.7, 0.24, 1.2)).items():
        assert psi == pytest.approx(0.0, abs=1e-10), l


def test_v10_xlag_scales_the_lag_angle():
    modes = _lagged_adiabatic()
    base = TS.lag_angles(modes)
    assert set(base) == {(2, 2), (3, 1), (3, 3)} and all(0.0 < p < 10.0 for p in base.values())
    for X_lag in (0.5, 2.0, 3.0):
        scaled = TS.lag_angles(modes, X_lag)
        for (l, m), psi in base.items():
            assert abs(X_lag * psi) < 180.0 / m
            assert scaled[(l, m)] == pytest.approx(X_lag * psi, rel=1e-9)


def test_v10_xlag_rotates_the_pattern_without_changing_amplitude():
    modes = _lagged_adiabatic()
    D1, _ = TS.combined_amplitudes(modes, 1.0)
    D3, _ = TS.combined_amplitudes(modes, 3.0)
    n, w = TS.gauss_sphere(64, 128)
    r1 = TS.delta_r_over_R(n, modes, 1.0)
    r3 = TS.delta_r_over_R(n, modes, 3.0)
    # amplitude of every harmonic is untouched, so the rms of the surface is identical
    assert TS.sphere_mean(r1 ** 2, w, 128) == pytest.approx(TS.sphere_mean(r3 ** 2, w, 128), rel=1e-10)
    for k in D1:
        assert abs(D3[k]) == pytest.approx(abs(D1[k]), rel=1e-12)


def test_v10_real_amplitudes_unaffected_by_xlag():
    modes = [dataclasses.replace(md, C=complex(md.C.real, 0.0)) for md in random_modes(2) + random_modes(3)]
    n = unit_vectors()
    np.testing.assert_allclose(TS.delta_r_over_R(n, modes, 7.0), TS.delta_r_over_R(n, modes, 1.0), atol=1e-14)


# --------------------------------------------------------------------- V11 regime gate

def test_v11_final_configuration(models):
    with pytest.warns(UserWarning, match=r"mode \(l=3, m=\+3\) is close to resonance"):
        rep_h = TS.check_regime(models["heavy"])               # one known warning (Delta_3/omega = 0.13)
    assert len(rep_h["warnings"]) == 1 and not rep_h["errors"]
    with warnings.catch_warnings():
        warnings.simplefilter("error")
        rep_l = TS.check_regime(models["light"])               # clean
    assert not rep_l["warnings"] and not rep_l["errors"]


def test_v11_past_resonance_is_an_error(models):
    bad = dataclasses.replace(models["light"], spin_hz=-500.0)   # the rejected configuration
    with pytest.raises(TS.RegimeError, match="past resonance"):
        TS.check_regime(bad)
    with pytest.warns(UserWarning):
        TS.check_regime(bad, allow_post_resonance=True)


def test_v11_amplitude_and_volume_warnings(models):
    big = dataclasses.replace(models["light"], X=2.0)
    rep = TS.regime_report(big)
    assert any("max|dR/R|" in w for w in rep["warnings"])
    assert any("|dV/V|" in w for w in rep["warnings"])


# --------------------------------------------------------------------- V12 spin oblateness

def test_v12_off_is_identical_to_tide_only(models):
    m = models["heavy"]
    n = unit_vectors()
    a = TS.surface_radius(n, 1.0, m.modes, m.X, m.X_lag, None)
    b = TS.surface_radius(n, 1.0, m.modes, m.X, m.X_lag, 0.0)
    np.testing.assert_array_equal(a, b)
    np.testing.assert_allclose(a, 1.0 + m.X * TS.delta_r_over_R(n, m.modes, m.X_lag), atol=1e-15)


def test_v12_axisymmetric_about_z_and_independent_of_spin_sign():
    n = unit_vectors()
    eps = 0.05
    rot = np.array([[np.cos(0.7), -np.sin(0.7), 0], [np.sin(0.7), np.cos(0.7), 0], [0, 0, 1]])
    np.testing.assert_allclose(TS.spin_oblateness(n, eps), TS.spin_oblateness(n @ rot.T, eps), atol=1e-15)
    assert TS.eps_omega(12.0, 1.6, +500.0, 1.52) == TS.eps_omega(12.0, 1.6, -500.0, 1.52)
    assert TS.spin_oblateness(np.array([[0, 0, 1.0]]), eps)[0] == pytest.approx(-eps)       # pole
    assert TS.spin_oblateness(np.array([[1.0, 0, 0]]), eps)[0] == pytest.approx(+eps / 2)   # equator


def test_v12_maclaurin_flattening():
    # uniform-density star: h2 = 5/2  =>  R_eq - R_pol = (5/4) q R, the classical small-rotation result
    R_km, M, spin = 12.0, 1.6, 50.0
    eps = TS.eps_omega(R_km, M, spin, 2.5)
    q = (2 * np.pi * spin) ** 2 * (R_km * 1e3) ** 3 / (P.G * M * P.M_SUN)
    diff = (TS.spin_oblateness(np.array([[1.0, 0, 0]]), eps) - TS.spin_oblateness(np.array([[0, 0, 1.0]]), eps))[0]
    assert diff == pytest.approx(1.25 * q, rel=1e-12)


def test_v12_final_configuration_sizes(models):
    assert models["heavy"].eps_omega == pytest.approx(0.0407, abs=5e-4)
    assert models["light"].eps_omega == pytest.approx(0.0098, abs=2e-4)


# --------------------------------------------------------------------- loader, frames, modules

def test_loader_final_configuration(models):
    assert set(models) == {"heavy", "light"}
    assert models["heavy"].spin_hz == -500.0 and models["light"].spin_hz == -200.0
    for m in models.values():
        assert m.spin_oblateness and m.X == 1.0 and m.X_lag == 1.0
        assert m.f_gw_hz == pytest.approx(binary.solve().f_gw, abs=1e-6)
        assert [(md.l, md.m) for md in m.modes] == [(2, 2), (2, 0), (2, -2), (3, 3), (3, 1), (3, -1), (3, -3)]
        I3 = {md.m: md.I for md in m.modes if md.l == 3}
        assert I3[-3] == -I3[3] and I3[-1] == -I3[1]            # I_{l,-m} = (-1)^m I_{l,m}
        assert all(md.omega_a0 > 0 for md in m.modes)
    assert models["heavy"].companion_sign == +1 and models["light"].companion_sign == -1


def test_mode_validation():
    with pytest.raises(ValueError, match="l\\+m odd"):
        TS.ModeAmp(3, 2, 0.0, 0.2, 1.2)
    with pytest.raises(ValueError, match="nonzero overlap"):
        TS.ModeAmp(2, 2, 0.0, 0.0, 1.5)
    with pytest.raises(ValueError, match="duplicate"):
        TS.combined_amplitudes([TS.ModeAmp(2, 2, 0.1, 0.3, 1.5), TS.ModeAmp(2, 2, 0.2, 0.3, 1.5)])
    with pytest.raises(ValueError, match="inconsistent h"):
        TS.combined_amplitudes([TS.ModeAmp(2, 2, 0.1, 0.3, 1.5), TS.ModeAmp(2, 0, 0.2, 0.3, 1.4)])


def test_star_frame_mapping_is_an_involution_and_flips_odd_l(models):
    n = unit_vectors()
    for sign in (+1, -1):
        np.testing.assert_array_equal(TS.star_frame_to_world(TS.star_frame_to_world(n, sign), sign), n)
    light = models["light"]
    fn = light.radius_fn(1.0)
    tip_world = fn(np.array([[-1.0, 0, 0]]))[0]                 # light star: companion on -x
    tip_star = light.rel_radius_star_frame(np.array([[1.0, 0, 0]]))[0]
    assert tip_world == pytest.approx(tip_star, rel=1e-14)


def test_generator_refuses_the_rejected_configuration(tmp_path):
    out = tmp_path / "bad.toml"
    r = subprocess.run([sys.executable, str(ROOT / "tools" / "make_star_modes.py"),
                        "--spin-light", "-500", "--out", str(out)], capture_output=True, text=True)
    assert r.returncode == 1 and not out.exists()
    assert "past resonance" in r.stdout
