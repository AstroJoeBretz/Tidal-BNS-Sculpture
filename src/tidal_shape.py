"""Tidally deformed neutron-star surface from f-mode amplitudes.

Implements docs/neutron_star_shape.md: complex mode amplitudes C_a of Yu et al.
(arXiv:2211.07002, 2404.00147) -> surface r(theta, phi) via Love-number
calibrated linear displacement (eq. star):

    dR/R = (h_l / I_l) Re sum_a C_a Y_{l m_a}(theta, phi_s)

Conventions (sec 2): star frame with the companion on +x at the snapshot, +z the
orbital angular momentum, Y_lm orthonormal with Condon-Shortley phase
(scipy.special.sph_harm_y), phi_s the azimuth measured in the orbit-corotating
frame. numpy + scipy only: the amplitudes are *data* (star_modes.toml,
written by tools/make_star_modes.py).
"""
from __future__ import annotations

import tomllib
import warnings
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np
from scipy.special import sph_harm_y

import params as P

N2 = 0.25 * np.sqrt(15.0 / (2.0 * np.pi))   # Y_22 = N2 sin^2(theta) e^{2 i phi}
N0 = 0.25 * np.sqrt(5.0 / np.pi)            # Y_20 = N0 (3 cos^2(theta) - 1)

# V11 thresholds (docs/neutron_star_shape.md sec 6)
MARGIN_WARN = 0.15          # Delta_a / omega_a^corot
MAX_DR_WARN = 0.30          # max |dR/R| of the tidal surface
MAX_DV_WARN = 0.05          # |dV/V|


class RegimeError(RuntimeError):
    """V11 error: the snapshot lies outside the linear regime (past resonance)."""


@dataclass(frozen=True)
class ModeAmp:
    """One f-mode amplitude at the snapshot.

    C        complex total amplitude (equilibrium + dynamical tide), paper normalization
    I        *signed* tidal overlap used when C was computed (I_{l,-m} = (-1)^m I_{l,m})
    h        equilibrium shape number h_l = 1 + 2 k_l of this multipole
    omega_a0 non-rotating mode frequency [rad/s]; only needed for the V11 regime gate
    """
    l: int
    m: int
    C: complex
    I: float
    h: float
    omega_a0: float | None = None

    def __post_init__(self):
        if self.l < 2 or abs(self.m) > self.l:
            raise ValueError(f"invalid mode (l={self.l}, m={self.m})")
        if (self.l + self.m) % 2:
            raise ValueError(f"mode (l={self.l}, m={self.m}) has l+m odd: it does not couple "
                             "to a companion in the orbital plane")
        if self.I == 0 or not np.isfinite(self.I):
            raise ValueError(f"mode (l={self.l}, m={self.m}) needs a nonzero overlap I")
        if not self.h > 0:
            raise ValueError(f"mode (l={self.l}, m={self.m}) needs a positive h_l")


@dataclass(frozen=True)
class StarModel:
    """Per-star snapshot configuration (one [[star]] block of star_modes.toml)."""
    name: str
    M_msun: float
    R_km: float
    Mc_msun: float
    companion_sign: int          # +1: companion on +x in the world frame, -1: on -x
    f_gw_hz: float
    spin_hz: float               # negative = anti-aligned with the orbit
    spin_oblateness: bool
    X: float
    X_lag: float
    modes: tuple[ModeAmp, ...]
    info: dict = field(default_factory=dict, compare=False)

    def h_of_l(self, l: int) -> float:
        hs = {md.h for md in self.modes if md.l == l}
        if not hs:
            raise KeyError(f"star {self.name!r} has no l={l} modes")
        return hs.pop()

    @property
    def eps_omega(self) -> float:
        """Centrifugal oblateness amplitude (sec 3.7) for this star's spin."""
        return eps_omega(self.R_km, self.M_msun, self.spin_hz, self.h_of_l(2))

    def rel_radius_star_frame(self, n_star):
        """r/R in the star frame (companion on +x), including X, X_lag and, if enabled, oblateness."""
        eps = self.eps_omega if self.spin_oblateness else None
        return surface_radius(n_star, 1.0, self.modes, self.X, self.X_lag, eps)

    def radius_fn(self, R_mm: float):
        """Callable n_world (..., 3) -> surface radius [mm] in the *world* frame."""
        def fn(n_world):
            return R_mm * self.rel_radius_star_frame(star_frame_to_world(n_world, self.companion_sign))
        return fn


# ---------------------------------------------------------------- geometry helpers

def _unit(n):
    n = np.asarray(n, dtype=float)
    return n / np.linalg.norm(n, axis=-1, keepdims=True)


def _angles(n):
    n = _unit(n)
    theta = np.arccos(np.clip(n[..., 2], -1.0, 1.0))
    phi = np.arctan2(n[..., 1], n[..., 0])
    return theta, phi


def star_frame_to_world(v, companion_sign: int):
    """Map between the star frame (companion on +x) and the world frame.

    +1: identity (companion on +x, the heavy star). -1: rotation by pi about z,
    (x, y, z) -> (-x, -y, z), which preserves the lag direction relative to the
    companion (sec 5). The map is an involution, so it also maps world -> star.
    """
    v = np.asarray(v, dtype=float)
    if companion_sign not in (+1, -1):
        raise ValueError("companion_sign must be +1 or -1")
    if companion_sign > 0:
        return v
    return v * np.array([-1.0, -1.0, 1.0])


def gauss_sphere(n_mu: int = 128, n_phi: int = 256):
    """Unit vectors (n_mu, n_phi, 3) with Gauss-Legendre weights in cos(theta)."""
    mu, w = np.polynomial.legendre.leggauss(n_mu)
    phi = np.linspace(0.0, 2.0 * np.pi, n_phi, endpoint=False)
    s = np.sqrt(1.0 - mu ** 2)[:, None]
    n = np.stack([s * np.cos(phi)[None, :], s * np.sin(phi)[None, :],
                  np.broadcast_to(mu[:, None], (n_mu, n_phi))], axis=-1)
    return n, w


# ---------------------------------------------------------------- mode combination (sec 3.3.1)

def _phase_ref(l: int, mm: int) -> float:
    """sign of Y_lm(pi/2, 0): the zero-lag reference of harmonic (l, m).

    Y_33 is negative on the line of centres (Condon-Shortley), so an adiabatic l=3, m=3
    amplitude D is *negative* real; the lag is measured from this orientation (sec 3.3.1).
    """
    return float(np.sign(np.real(sph_harm_y(l, mm, np.pi / 2, 0.0))))


def combined_amplitudes(modes, X_lag: float = 1.0):
    """Combined complex amplitudes D_m per (l, |m|) pair, with the lag exaggerated by X_lag.

    D_m = C^_{+m} + (-1)^m conj(C^_{-m}),  C^ = C / I_a   (m > 0),
    D_0 = Re C^_0, so that Re[C^_{+m} Y_lm + C^_{-m} Y_l,-m] = Re[D_m Y_lm].
    With D~ = D sign(Y_lm(pi/2, 0)) the equatorial pattern is |D| |Y| cos(m phi_s + arg D~), so
    the lag is psi_m = arg(D~)/m. X_lag rotates the pattern only: D~ -> |D| exp(i X_lag arg D~).
    Returns ({(l, m>=0): D}, {l: h_l}).
    """
    by = {}
    h = {}
    for md in modes:
        key = (md.l, md.m)
        if key in by:
            raise ValueError(f"duplicate mode (l={md.l}, m={md.m})")
        by[key] = md.C / md.I
        if h.setdefault(md.l, md.h) != md.h:
            raise ValueError(f"inconsistent h_l for l={md.l}")
    D = {}
    for (l, mm) in sorted({(l, abs(m)) for (l, m) in by}):
        if mm == 0:
            D[(l, 0)] = complex(by[(l, 0)].real)
            continue
        cp = by.get((l, +mm), 0.0)
        cn = by.get((l, -mm), 0.0)
        d = cp + (-1) ** mm * np.conj(cn)
        if X_lag != 1.0 and d != 0:
            s = _phase_ref(l, mm)
            dt = d * s
            d = s * abs(d) * np.exp(1j * X_lag * np.arctan2(dt.imag, dt.real))
        D[(l, mm)] = complex(d)
    return D, h


def lag_angles(modes, X_lag: float = 1.0):
    """Lag psi_m = arg(D~_m)/m in degrees per (l, m>0); positive = behind the companion."""
    D, _ = combined_amplitudes(modes, X_lag)
    out = {}
    for (l, mm), d in D.items():
        if mm > 0:
            dt = d * _phase_ref(l, mm)
            out[(l, mm)] = float(np.degrees(np.arctan2(dt.imag, dt.real) / mm))
    return out


# ---------------------------------------------------------------- the surface (sec 3, 5)

def delta_r_over_R(n_star, modes, X_lag: float = 1.0):
    """Tidal dR/R at unit vectors n_star (..., 3) in the star frame: eq. (star) with D_m of sec 3.3.1."""
    theta, phi = _angles(n_star)
    D, h = combined_amplitudes(modes, X_lag)
    out = np.zeros(theta.shape)
    for (l, mm), d in D.items():
        out += h[l] * np.real(d * sph_harm_y(l, mm, theta, phi))
    return out


def eps_omega(R_km: float, M_msun: float, spin_hz: float, h2: float) -> float:
    """eps_Omega = (h2/3) Omega^2 R^3 / (G M)  (sec 3.7)."""
    Om = 2.0 * np.pi * spin_hz
    return (h2 / 3.0) * Om ** 2 * (R_km * 1e3) ** 3 / (P.G * M_msun * P.M_SUN)


def spin_oblateness(n_star, eps: float):
    """dR_spin/R = -eps P2(cos theta), theta from the spin axis (+-z). Sign of the spin is irrelevant."""
    nz = _unit(n_star)[..., 2]
    return -eps * 0.5 * (3.0 * nz ** 2 - 1.0)


def surface_radius(n_star, R, modes, X: float = 1.0, X_lag: float = 1.0, eps_omega_value=None):
    """r = R [1 + X dR_tide/R + s_obl dR_spin/R] in the star frame (sec 5). Raises if r <= 0 anywhere."""
    rel = 1.0 + X * delta_r_over_R(n_star, modes, X_lag)
    if eps_omega_value is not None:
        rel = rel + spin_oblateness(n_star, eps_omega_value)
    if np.min(rel) <= 0.0:
        raise ValueError(f"surface radius is not positive (min r/R = {np.min(rel):.3f})")
    return R * rel


def quadrupole_form(Ca, Cb, Cc, I, h2):
    """Closed form of the l=2 tidal surface (sec 3.3).

    dR/R = e_xx nx^2 + e_yy ny^2 + e_zz nz^2 + 2 e_xy nx ny, with S_c = Re(Ca+Cb), S_s = Im(Ca-Cb),
    K = h2/I, psi = 0.5 atan2(S_s, S_c): the long in-plane axis sits at phi_s = -psi (behind the
    companion for psi > 0). Ca, Cb, Cc are the (2,+2), (2,-2), (2,0) amplitudes.
    """
    K = h2 / I
    Sc = (Ca + Cb).real
    Ss = (Ca - Cb).imag
    Rc = Cc.real
    # psi is taken from C/I (sec 3.3.1: D_2 = S_c/I + i S_s/I), so it does not depend on the sign
    # convention of the eigenfunction normalization; for I > 0 this is atan2(S_s, S_c) of sec 3.3.
    s = np.sign(I)
    return {
        "e_xx": K * (N2 * Sc - N0 * Rc),
        "e_yy": K * (-N2 * Sc - N0 * Rc),
        "e_zz": K * 2.0 * N0 * Rc,
        "e_xy": -K * N2 * Ss,
        "psi_deg": float(np.degrees(0.5 * np.arctan2(s * Ss, s * Sc))),
        "rho": float(np.hypot(Sc, Ss)),
    }


def volume_drift(rel_radius_fn, n_mu: int = 128, n_phi: int = 256) -> float:
    """dV/V0 of the star-shaped surface r/R = rel_radius_fn(n): <(r/R)^3> - 1 over the sphere."""
    n, w = gauss_sphere(n_mu, n_phi)
    rel = rel_radius_fn(n)
    return float(np.sum(w[:, None] * rel ** 3) / (2.0 * n_phi) - 1.0)


def sphere_mean(f_values, w, n_phi: int) -> float:
    """Mean of f over the sphere from Gauss sphere samples."""
    return float(np.sum(w[:, None] * f_values) / (2.0 * n_phi))


# ---------------------------------------------------------------- V11: linear-regime gate

def regime_report(star: StarModel, n_mu: int = 96, n_phi: int = 192) -> dict:
    """Linear-regime diagnostics of one star's snapshot (sec 6, V11). Does not raise."""
    Om = 2.0 * np.pi * star.spin_hz
    w = np.pi * star.f_gw_hz
    modes = []
    for md in star.modes:
        if md.omega_a0 is None:
            continue
        w_cor = md.omega_a0 - (md.m / md.l) * Om
        delta = w_cor + md.m * Om - md.m * w
        modes.append({"l": md.l, "m": md.m, "delta_over_omega": float(delta / w_cor)})

    n, wq = gauss_sphere(n_mu, n_phi)
    tide = star.X * delta_r_over_R(n, star.modes, star.X_lag)
    rel = star.rel_radius_star_frame(n)
    dV_tide = float(np.sum(wq[:, None] * (1.0 + tide) ** 3) / (2.0 * n_phi) - 1.0)
    dV_total = float(np.sum(wq[:, None] * rel ** 3) / (2.0 * n_phi) - 1.0)

    errors, warns = [], []
    if not modes:
        warns.append("no mode has omega_a0: cannot evaluate the distance to resonance")
    for r in modes:
        tag = f"(l={r['l']}, m={r['m']:+d})"
        if r["m"] > 0 and r["delta_over_omega"] < 0:
            errors.append(f"{star.name}: mode {tag} is past resonance "
                          f"(Delta/omega = {r['delta_over_omega']:.3f} < 0)")
        elif r["delta_over_omega"] < MARGIN_WARN:
            warns.append(f"{star.name}: mode {tag} is close to resonance "
                         f"(Delta/omega = {r['delta_over_omega']:.3f} < {MARGIN_WARN})")
    max_dr = float(np.max(np.abs(tide)))
    if max_dr > MAX_DR_WARN:
        warns.append(f"{star.name}: max|dR/R| = {max_dr:.3f} > {MAX_DR_WARN} (linear surface is strained)")
    if abs(dV_tide) > MAX_DV_WARN:
        warns.append(f"{star.name}: |dV/V| = {abs(dV_tide):.3f} > {MAX_DV_WARN}")
    return {
        "name": star.name, "modes": modes,
        "min_margin": min((r["delta_over_omega"] for r in modes), default=float("nan")),
        "tide_range": (float(tide.min()), float(tide.max())), "max_abs_dr": max_dr,
        "total_range": (float(rel.min() - 1.0), float(rel.max() - 1.0)),
        "min_over_max_r": float(rel.min() / rel.max()),
        "dV_tide": dV_tide, "dV_total": dV_total,
        "errors": errors, "warnings": warns,
    }


def check_regime(star: StarModel, allow_post_resonance: bool = False, quiet: bool = False) -> dict:
    """Run V11: raise RegimeError for past-resonance modes (unless allowed), warn for the soft limits."""
    rep = regime_report(star)
    if rep["errors"] and not allow_post_resonance:
        raise RegimeError("; ".join(rep["errors"]))
    if not quiet:
        for msg in rep["warnings"] + (rep["errors"] if allow_post_resonance else []):
            warnings.warn(msg, stacklevel=2)
    return rep


# ---------------------------------------------------------------- data file

def load_star_modes(path) -> dict[str, StarModel]:
    """Read a star_modes.toml -> {star name: StarModel}."""
    with open(path, "rb") as fh:
        raw = tomllib.load(fh)
    stars = {}
    for s in raw["star"]:
        modes = tuple(
            ModeAmp(l=int(m["l"]), m=int(m["m"]), C=complex(m["C_re"], m["C_im"]), I=float(m["I"]),
                    h=float(m["h"]), omega_a0=float(m["omega_a0_rad_s"]) if "omega_a0_rad_s" in m else None)
            for m in s["mode"])
        info = {k: v for k, v in s.items() if k not in ("mode",)}
        stars[s["name"]] = StarModel(
            name=s["name"], M_msun=float(s["M_msun"]), R_km=float(s["R_km"]), Mc_msun=float(s["Mc_msun"]),
            companion_sign=int(s["companion_sign"]), f_gw_hz=float(s["f_gw_hz"]),
            spin_hz=float(s["spin_hz"]), spin_oblateness=bool(s["spin_oblateness"]),
            X=float(s["X"]), X_lag=float(s["X_lag"]), modes=modes, info=info)
    return stars
