"""Orbital dynamics and the physical -> model length scale. docs/design.md sec 1-2."""
from dataclasses import dataclass

import numpy as np

import params as P


@dataclass(frozen=True)
class Binary:
    M1: float           # solar masses
    M2: float
    a_km: float
    Omega: float         # rad/s
    f_orb: float          # Hz
    f_gw: float            # Hz
    lambda_gw_km: float
    S_mm_per_km: float     # master scale, sec 1.2
    a_mm: float
    R1_mm: float
    R2_mm: float
    x1_mm: float           # heavy-star centre, world x
    x2_mm: float           # light-star centre, world x


def solve() -> Binary:
    """Keplerian circular orbit + the scale factor that makes one lambda_GW
    SCALE_LAMBDA_MM long (sec 1.2)."""
    M1 = P.M1_MSUN * P.M_SUN
    M2 = P.M2_MSUN * P.M_SUN
    M = M1 + M2
    a = P.A_KM * 1e3

    Omega = np.sqrt(P.G * M / a ** 3)
    f_orb = Omega / (2 * np.pi)
    f_gw = Omega / np.pi
    lambda_gw_km = P.C / f_gw / 1e3

    S = P.SCALE_LAMBDA_MM / lambda_gw_km

    a_mm = P.A_KM * S
    R1_mm = P.R1_KM * S
    R2_mm = P.R2_KM * S

    # barycentre: the heavier star sits closer to the origin, on the -x side
    x1_mm = -(P.M2_MSUN / (P.M1_MSUN + P.M2_MSUN)) * a_mm
    x2_mm = +(P.M1_MSUN / (P.M1_MSUN + P.M2_MSUN)) * a_mm

    return Binary(
        M1=P.M1_MSUN, M2=P.M2_MSUN, a_km=P.A_KM,
        Omega=Omega, f_orb=f_orb, f_gw=f_gw, lambda_gw_km=lambda_gw_km,
        S_mm_per_km=S, a_mm=a_mm, R1_mm=R1_mm, R2_mm=R2_mm,
        x1_mm=x1_mm, x2_mm=x2_mm,
    )


def wavenumber_per_mm(b: Binary) -> float:
    """k = 2*pi / lambda_GW at scale. Equals 2*pi / SCALE_LAMBDA_MM by
    construction (sec 1.2) -- kept as a function of b for traceability."""
    del b  # invariant under the physical parameters by construction
    return 2 * np.pi / P.SCALE_LAMBDA_MM


def retardation_k_half(b: Binary, c_eff: float = P.C) -> float:
    """Omega/c_eff expressed in scaled mm^-1 -- the wavenumber governing each
    star's own retardation phase (docs/design.md sec 3.2).

    Equals wavenumber_per_mm(b)/2 when c_eff == P.C: S is defined (sec 1.2)
    so that exactly one omega_GW = 2*Omega wavelength is
    SCALE_LAMBDA_MM long, and that identity is only true for the real speed
    of light used to derive S. Passing a different c_eff (e.g. to explore a
    hypothetical wave speed, sec 3.5) changes the retardation phase
    without touching the plate/star geometry, which stays fixed at whatever
    S the real C produced.
    """
    return b.Omega * 1000.0 / (b.S_mm_per_km * c_eff)


if __name__ == "__main__":
    b = solve()
    print(f"Omega        = {b.Omega:.1f} rad/s")
    print(f"f_orb        = {b.f_orb:.1f} Hz")
    print(f"f_GW         = {b.f_gw:.1f} Hz")
    print(f"lambda_GW    = {b.lambda_gw_km:.2f} km")
    print(f"S            = {b.S_mm_per_km:.5f} mm/km")
    print(f"a            = {b.a_mm:.3f} mm")
    print(f"R1, R2       = {b.R1_mm:.3f}, {b.R2_mm:.3f} mm")
    print(f"x1, x2       = {b.x1_mm:.3f}, {b.x2_mm:.3f} mm")
    print(f"surface gap  = {b.a_mm - b.R1_mm - b.R2_mm:.3f} mm (undeformed spheres)")
