"""Independent float64 re-implementation of the chain  orbit -> mode amplitudes -> surface.

Companion to docs/math.md: every equation of that document, written out step by
step with numpy/scipy only (no tidal_waveform import), so the numbers in the worked example can be
reproduced and the production code (src/tidal_waveform.py + src/tidal_shape.py) can be checked against it.

Input is data/default/star_modes.toml (the default configuration's reference amplitudes; masses, radii, spin, snapshot, I, h, omega_a0, A_nl) -- *not* its amplitudes C,
which are what this script recomputes and compares with.

Run:  python verification/math_walkthrough.py [--star heavy|light|both]
"""
import argparse
import sys
import tomllib
from pathlib import Path

import numpy as np
from scipy.integrate import quad
from scipy.special import sph_harm_y

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "src"))

# cgs constants -- the same values src/tidal_waveform.py uses
G, C_LIGHT, MS = 6.67259e-08, 29979245800.0, 1.989e+33
W_REF = np.pi * 1000.0          # reference GW frequency of the nonlinear shift (tidal_waveform default)


def W(l, m):
    """W_lm = 4 pi Y_lm(pi/2, 0)/(2l+1)."""
    return 4.0 * np.pi / (2 * l + 1) * np.real(sph_harm_y(l, m, np.pi / 2, 0.0))


def gr_factors(x, X1, X2, gr=1.0):
    """Redshift factor zz and frame-dragging factor fd = omega_FD/omega (2PN expansions of S+21 eqs 4.7-4.8)."""
    zz = 1.0 + (x / 2 * X2 * (X1 - 3.0) + x ** 2 / 24 * X2 * (5 * X1 * X1 * X2 - 9 * X1 - 6 * X1 * X2 - 27)) * gr
    fd = (x / 2 * X2 * (X1 + 3.0) - x ** 2 / 24 * X2 * (X1 * X1 * X2 - 45 * X1 + 30 * X1 * X2 - 27)) * gr
    return zz, fd


def orbit(M1, M2, omega):
    """Quantities of the point-mass quasi-circular inspiral at orbital angular frequency omega [rad/s]."""
    Mt = M1 + M2
    eta = M1 * M2 / Mt ** 2
    tMt = G * Mt / C_LIGHT ** 3
    v = (tMt * omega) ** (1 / 3)
    omega_dot = 96.0 / 5 * eta * v ** 11 / tMt ** 2
    return dict(Mt=Mt, eta=eta, tMt=tMt, v=v, x=v * v, r=(G * Mt / omega ** 2) ** (1 / 3),
                omega_dot=omega_dot, t_gw=omega / omega_dot, dlogr=-(2.0 / 3.0) * omega_dot / omega)


def delta_a(omega, l, m, w_a0, Om_s, A_nl, X1, X2, tMt):
    """Delta_a(omega) of the code, with GR and the nonlinear shift."""
    zz, fd = gr_factors((tMt * omega) ** (2 / 3), X1, X2)
    w_a = w_a0 - (m / l) * Om_s
    return zz * (w_a + m * Om_s) - m * omega + m * fd * omega + A_nl * (omega / W_REF) ** 2 * w_a0, zz


def equilibrium(l, m, w_a0, V, omega, orb, Om_s, A_nl, X1, X2):
    """Paper 2 Eq. 18 with GR / nonlinear corrections: returns b_eq and its ingredients."""
    D, zz = delta_a(omega, l, m, w_a0, Om_s, A_nl, X1, X2, orb["tMt"])
    kappa = m * orb["omega_dot"] - D * (l + 1) * orb["dlogr"]
    b = zz * D * w_a0 * V / (D ** 2 - 1j * kappa)
    return dict(b=b, Delta=D, kappa=kappa, zz=zz)


def resonance(l, m, w_a0, Om_s, A_nl, X1, X2, tMt, closed_form_l2=False, n_iter=60):
    """Orbital frequency at which Delta_a = 0 for the prograde mode (m = l).

    closed_form_l2 reproduces get_ba_l2's 1PN closed form; otherwise a fixed-point iteration (get_ba_l3)."""
    w_a = w_a0 - (m / l) * Om_s
    if closed_form_l2:
        wN = w_a / 2 + Om_s
        x = (tMt * wN) ** (2 / 3)
        eta = X1 * X2
        return (1 + eta * x + A_nl * (wN / W_REF) ** 2) * wN
    w = (w_a + m * Om_s) / m
    for _ in range(n_iter):
        zz, fd = gr_factors((tMt * w) ** (2 / 3), X1, X2)
        w = (zz * (w_a + m * Om_s) + A_nl * (w / W_REF) ** 2 * w_a0) / (m * (1 - fd))
    return w


def dynamical(l, m, w_a0, I, omega, omega_res, M1, M2, R, Om_s, X1, X2, closed_form_l2):
    """Paper 2 Eqs. 23, 24, 30, 31: resonantly excited component at orbital frequency omega (exact F(u))."""
    orb = orbit(M1, M2, omega)
    res = orbit(M1, M2, omega_res)
    tMt, eta = orb["tMt"], orb["eta"]
    tau = 5.0 / 256 / eta * tMt * (1 / res["v"] ** 8 - 1 / orb["v"] ** 8)               # Eq. 26
    phi = 1.0 / 32 / eta * (1 / res["v"] ** 5 - 1 / orb["v"] ** 5)                       # Eq. 31
    scale = np.sqrt(m / 2.0 * res["omega_dot"])
    u = scale * tau                                                                      # Eq. 25
    eps = np.sqrt((2.0 / m) * res["omega_dot"] / omega_res ** 2)
    a1, a2 = (2.0 * l - 9) / 6 * eps, -2.0 / 3 * (l + 1.0) * eps
    f = lambda t: 2j * (1 + a1 * t) * np.exp(-1j * t * t) / (1 + a2 * t + 2j * t * t) ** 2   # Eq. 24 integrand
    lo = -80.0                                                                           # tail beyond |u|=80 is < 1e-8
    re = quad(lambda t: f(t).real, lo, u, limit=4000)[0]
    im = quad(lambda t: f(t).imag, lo, u, limit=4000)[0]
    F = re + 1j * im
    zz_r, _ = gr_factors(res["x"], X1, X2)
    V_r = (M2 / M1) * (R / res["r"]) ** (l + 1) * W(l, m) * I
    b = zz_r * w_a0 * V_r / scale * F * np.exp(1j * m * (phi - omega_res * tau))         # Eqs. 23, 30
    return dict(b=b, u=u, F=F, tau=tau, phi=phi, a1=a1, a2=a2)


def chain(star):
    """All the numbers for one [[star]] block of star_modes.toml."""
    M1, M2, R = star["M_msun"] * MS, star["Mc_msun"] * MS, star["R_km"] * 1e5
    X1, X2 = M1 / (M1 + M2), M2 / (M1 + M2)
    Om_s = 2 * np.pi * star["spin_hz"]
    omega = np.pi * star["f_gw_hz"]
    orb = orbit(M1, M2, omega)
    zz, fd = gr_factors(orb["x"], X1, X2)
    out = dict(name=star["name"], omega=omega, r_over_R=orb["r"] / R, omega_dot=orb["omega_dot"],
               t_gw=orb["t_gw"], dlogr=orb["dlogr"], zz=zz, fd=fd, x=orb["x"], modes=[])
    for md in star["mode"]:
        l, m, I, w_a0 = md["l"], md["m"], md["I"], md["omega_a0_rad_s"]
        A_nl = star["A_nl_l2"] if l == 2 else 0.0
        V = (M2 / M1) * (R / orb["r"]) ** (l + 1) * W(l, m) * I
        eq = equilibrium(l, m, w_a0, V, omega, orb, Om_s, A_nl, X1, X2)
        row = dict(l=l, m=m, I=I, h=md["h"], V=V, w_a0=w_a0, w_corot=w_a0 - (m / l) * Om_s, **eq,
                   b_dyn=0.0 + 0.0j, C_toml=complex(md["C_re"], md["C_im"]))
        if m == l:                                                                       # prograde: can resonate
            w_res = resonance(l, m, w_a0, Om_s, A_nl, X1, X2, orb["tMt"], closed_form_l2=(l == 2))
            dyn = dynamical(l, m, w_a0, I, omega, w_res, M1, M2, R, Om_s, X1, X2, closed_form_l2=(l == 2))
            row.update(b_dyn=dyn["b"], f_res_gw=w_res / np.pi, u=dyn["u"], F=dyn["F"])
        row["C"] = row["b"] + row["b_dyn"]
        out["modes"].append(row)
    return out


def surface_from_amplitudes(modes, X_lag=1.0):
    """Surface pieces from the recomputed amplitudes: C^ = C/I, D_m, lag psi, closed-form l=2 tensor."""
    import tidal_shape as TS
    ma = [TS.ModeAmp(r["l"], r["m"], r["C"], r["I"], r["h"]) for r in modes]
    c2 = {r["m"]: r["C"] for r in modes if r["l"] == 2}
    I2 = next(r["I"] for r in modes if r["l"] == 2)
    h2 = next(r["h"] for r in modes if r["l"] == 2)
    q = TS.quadrupole_form(c2[2], c2[-2], c2[0], I2, h2)
    return dict(modes=ma, D=TS.combined_amplitudes(ma, X_lag)[0], lag=TS.lag_angles(ma, X_lag), quad=q)


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--star", default="both", choices=["heavy", "light", "both"])
    args = ap.parse_args(argv)
    raw = tomllib.load(open(ROOT / "data" / "default" / "star_modes.toml", "rb"))
    worst = 0.0
    for star in raw["star"]:
        if args.star not in ("both", star["name"]):
            continue
        c = chain(star)
        print(f"\n===== {c['name']} =====  omega = {c['omega']:.2f} rad/s, r/R = {c['r_over_R']:.4f}, "
              f"omega_dot = {c['omega_dot']:.4e} rad/s^2, t_gw = {c['t_gw']*1e3:.3f} ms, dlogr/dt = {c['dlogr']:.3f} /s")
        print(f"   x = (G Mt omega/c^3)^(2/3) = {c['x']:.5f},  zz = {c['zz']:.5f},  omega_FD/omega = {c['fd']:.5f}")
        print("   l  m      V          Delta [1/s]   kappa [1/s^2]   b_eq                      |b_dyn|     C (here)   "
              "rel diff vs toml")
        for r in c["modes"]:
            d = abs(r["C"] - r["C_toml"]) / abs(r["C_toml"])
            worst = max(worst, d)
            print(f"   {r['l']} {r['m']:+d}  {r['V']:+.4e}  {r['Delta']:+.4e}  {r['kappa']:+.4e}  "
                  f"{r['b']:+.5e}{r['b'].imag:+.2e}j  {abs(r['b_dyn']):.2e}  {abs(r['C']):.5e}  {d:.1e}")
        s = surface_from_amplitudes(c["modes"])
        q = s["quad"]
        print(f"   l=2 closed form: e_xx={q['e_xx']:+.5f} e_yy={q['e_yy']:+.5f} e_zz={q['e_zz']:+.5f} "
              f"e_xy={q['e_xy']:+.6f} psi={q['psi_deg']:.3f} deg")
        print("   combined D_m:", {k: f"{v:+.4e}" for k, v in s["D"].items()},
              " lag [deg]:", {k: round(v, 3) for k, v in s["lag"].items()})
    print(f"\nworst relative difference of C vs data/default/star_modes.toml: {worst:.2e}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
