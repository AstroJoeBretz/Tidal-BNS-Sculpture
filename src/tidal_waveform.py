"""f-mode amplitudes of a tidally driven neutron star along a quasi-circular inspiral.

NumPy/SciPy implementation of the amplitude model of

  Paper 1: Yu, Weinberg, Arras, Kwon & Venumadhav, arXiv:2211.07002
  Paper 2: Yu, Arras & Weinberg, arXiv:2404.00147

for the $l=2$ and $l=3$ f-modes of star 1 driven by a point-mass companion (star 2). Each
function returns the *total* complex amplitude $b_a$ in the star's co-rotating frame: the
equilibrium tide (Paper 2 Eq. 18, including its orbital-decay imaginary part, GR corrections
and the nonlinear frequency shift) plus the dynamical tide excited at the prograde resonance
(Paper 2 Eqs. 23-24, 30-31). The surface built from these amplitudes is
`tidal_shape.py`; the derivation of every step is in `docs/math.md`.

All inputs are cgs (cm, g, s, rad/s) except the GW frequency grid `ff`, which is in Hz.
Negative `OmegaS1` is a spin anti-aligned with the orbit (Paper 2 sign convention).
"""
import warnings

import numpy as np
from scipy.integrate import cumulative_trapezoid
from scipy.special import sph_harm_y

# cgs constants (the values the amplitude model and the worked example in docs/math.md use)
G, c, Ms = 6.67259e-08, 29979245800.0, 1.989e+33

# Mode azimuthal numbers m, in the row order of the returned amplitude arrays. Only modes with
# l+m even couple to a companion in the orbital plane.
MM_L2 = (2, 0, -2)
MM_L3 = (3, 1, -1, -3)

# Default reference GW angular frequency of the nonlinear frequency shift, $\pi\cdot$1 kHz.
W_NL_REF = np.pi * 1000


def _gr_factors(x, X1, X2, gr_flag):
    """Redshift factor and frame-dragging shift for a 1PN-accurate mode frequency.

    Eqs. 4.7 and 4.8 of S+21 (as quoted in Paper 2). `x = (G Mt omega / c^3)^(2/3)`.
    Returns `(zz, fd)` with the frame-dragging frequency shift equal to `fd * omega`.
    """
    zz = 1. + (x / 2. * X2 * (X1 - 3.)
               + x**2 / 24. * X2 * (5 * X1 * X1 * X2 - 9 * X1 - 6 * X1 * X2 - 27)) * gr_flag
    fd = (x / 2. * X2 * (X1 + 3.)
          - x**2 / 24. * X2 * (X1 * X1 * X2 - 45 * X1 + 30 * X1 * X2 - 27)) * gr_flag
    return zz, fd


def _inspiral(ff, M1, M2):
    """Quasi-circular inspiral quantities on the GW frequency grid (Paper 2 Eqs. 66-67)."""
    Mt = M1 + M2
    eta = M1 * M2 / Mt**2
    ww = np.pi * ff
    rr = np.cbrt(G * Mt / ww**2.)
    tMt = G * Mt / c**3
    vv = np.cbrt(tMt * ww)
    dwdt = 96. / 5 * eta * vv**11. / tMt**2.
    dlogrdt = -2. / 3. / (ww / dwdt)
    return dict(Mt=Mt, eta=eta, X1=M1 / Mt, X2=M2 / Mt, ww=ww, rr=rr, tMt=tMt, vv=vv,
                xx=vv * vv, dwdt=dwdt, dlogrdt=dlogrdt)


def _dynamical_tide(orb, wres, m_r, l, zz_r, wa01, Wlm_r, I_r, M1, M2, R1):
    """Dynamical-tide amplitude of the prograde resonant mode on the grid (Paper 2 Eqs. 23-26, 30-31).

    `wres` is the resonance GW angular frequency, `Wlm_r * I_r` the drive overlap of the resonant
    mode. Returns `(ba_dyn, Fu)`; `Fu` is the Fresnel-like integral F(u) of Eq. 24, accumulated
    from `ff[0]` (a warning is issued if the grid starts inside the resonance window).
    """
    eta, tMt, vv = orb["eta"], orb["tMt"], orb["vv"]
    vv_res = np.cbrt(tMt * wres)

    # eqs. 25-26 of 2404.00147
    dwdt_r = 96. / 5 * eta * vv_res**11. / tMt**2.
    tau_r = 5. / 256. / eta * tMt * (1. / vv_res**8. - 1. / vv**8.)
    phi_r = 1. / 32. / eta * (1. / vv_res**5. - 1. / vv**5.)

    uu = np.sqrt(m_r / 2. * dwdt_r) * tau_r
    if uu[0] > -5.:
        warnings.warn(f"ff[0] is within the ({l},{m_r}) resonance window (u={uu[0]:.2f} > -5); "
                      "F(u) is accumulated from ff[0], so start the grid at a lower frequency.")

    # a1, a2 below eq. 25, with s_a = +1 for the prograde mode
    eps = np.sqrt((2. / m_r) * dwdt_r / wres**2.)
    a1 = (2. * l - 9) / 6 * eps
    a2 = -2. / 3 * (l + 1.) * eps

    # eq. 24 of 2404.00147
    dFu = 2.j * (1. + a1 * uu) * np.exp(-1j * uu**2.) / (1 + a2 * uu + 2.j * uu**2.)**2.
    Fu = cumulative_trapezoid(dFu, x=uu, initial=0)

    # eqs. 23, 30 and 31 of 2404.00147; redshift factor at resonance
    rr_r = np.cbrt(G * orb["Mt"] / wres**2.)
    va_r = (M2 / M1) * (R1 / rr_r)**(l + 1.) * Wlm_r * I_r
    ba_dyn = zz_r * wa01 * va_r / np.sqrt(m_r / 2. * dwdt_r) * Fu \
        * np.exp(1j * (m_r * (phi_r - wres * tau_r)))
    return ba_dyn, Fu


def get_ba_l2(ff, M1, M2, R1, lam_M1_5, wa01, A_nl_1=-3e-2, w_nl_ref=W_NL_REF, OmegaS1=0, gr_flag=1):
    """$l=2$ f-mode amplitudes $b_a$ of star 1 (mass `M1`, radius `R1`) with companion `M2`.

    ff:        GW frequency grid [Hz], strictly increasing. The dynamical-tide integral is
               accumulated from `ff[0]`, so start well below the (2,2) resonance.
    lam_M1_5:  dimensionless tidal deformability $\\Lambda_2$ of star 1; the Love number is
               $k_2 = 1.5\\,\\Lambda_2 (G M_1/(R_1 c^2))^5$.
    wa01:      non-rotating $l=2$ f-mode angular frequency [rad/s].
    A_nl_1:    nonlinear frequency-shift coefficient (default: equal-mass 1.3+1.3 Msun polytrope).
    w_nl_ref:  reference GW angular frequency of the nonlinear shift [rad/s].
    OmegaS1:   star-1 spin angular frequency [rad/s]; negative = anti-aligned.
    gr_flag:   1 includes the 1PN redshift / frame-dragging corrections, 0 is Newtonian.

    Returns `(ba1, ba_eq1, va1, Fu1)`. `ba1`, `ba_eq1`, `va1` have shape (3, len(ff)), rows in
    the order `MM_L2 = (2, 0, -2)`; `ba1 = ba_eq1` plus the dynamical tide on the (2,2) mode.
    `va1` is the drive $v_a$ (Paper 2 Eq. 9, phase omitted) and `Fu1` is F(u) of Eq. 24.

    The closed-form resonance frequency used here is the 1PN one, which is 0.5-2% off the
    exact root of $\\Delta_a=0$ once GR and the nonlinear shift are on (`get_ba_l3` solves for
    the root instead).
    """
    ff = np.asarray(ff, dtype=float)
    orb = _inspiral(ff, M1, M2)
    ww, rr, eta, tMt, dwdt, dlogrdt = (orb[k] for k in ("ww", "rr", "eta", "tMt", "dwdt", "dlogrdt"))
    zz1, fd1 = _gr_factors(orb["xx"], orb["X1"], orb["X2"], gr_flag)
    omega_FD_SL1 = fd1 * ww

    ll1 = np.array([2, 2, 2])
    mm1 = np.array(MM_L2)
    wa1 = -mm1 / ll1 * OmegaS1 + wa01           # co-rotating frame frequency

    # eq. 9 of 2404.00147, phase omitted
    Wlm1 = np.array([np.sqrt(3. * np.pi / 10), -np.sqrt(np.pi / 5), np.sqrt(3. * np.pi / 10)])
    rM1 = G * M1 / c**2.
    k2 = 1.5 * lam_M1_5 * (rM1 / R1)**5
    Qa = np.sqrt(k2 * 5 / 4. / np.pi)
    drive = (M2 / M1) * (R1 / rr)**(ll1[0] + 1.) * Qa
    va1 = np.stack([drive * W for W in Wlm1])

    # eq. 18 of 2404.00147 including GR & nonlinear tide. (The m=0 expression is not strictly
    # right for m=0 but that mode is far off resonance.)
    ba_eq1 = np.empty((3, len(ff)), dtype=complex)
    for i in range(3):
        Del = zz1 * (wa1[i] + mm1[i] * OmegaS1) - mm1[i] * ww + mm1[i] * omega_FD_SL1 \
            + A_nl_1 * (ww / w_nl_ref)**2. * wa01
        ba_eq1[i] = zz1 * Del * wa01 * va1[i] \
            / (Del**2. - 1j * (mm1[i] * dwdt - Del * (ll1[i] + 1.) * dlogrdt))

    # closed-form 1PN resonance: Newtonian, linear value, then GR (raises) and NL (lowers)
    wres_N1 = wa1[0] / 2. + OmegaS1
    xx_res = np.cbrt(tMt * wres_N1)**2
    wres1 = (1 + eta * xx_res * gr_flag + A_nl_1 * (wres_N1 / w_nl_ref)**2.) * wres_N1

    ba_dyn, Fu1 = _dynamical_tide(orb, wres1, 2, 2, zz1, wa01, Wlm1[0], Qa, M1, M2, R1)
    ba1 = ba_eq1.copy()
    ba1[0] += ba_dyn
    return ba1, ba_eq1, va1, Fu1


def get_overlaps_l3(k3):
    """Signed $l=3$ f-mode overlaps $I_{3,m}$, in the row order `MM_L3 = (3, 1, -1, -3)`.

    $I_{3,-m} = (-1)^m I_{3,m}$ (Paper 2, below Eq. 11), so $I_{3,\\pm3}$ and $I_{3,\\pm1}$ flip
    sign between $+m$ and $-m$. The surface shape needs $C_a/I_a$ mode by mode, so keep these
    signs with the amplitudes from `get_ba_l3`.
    """
    Q3 = np.sqrt(k3 * 7. / 4. / np.pi)          # |I_3|, from k_l0 = 4 pi I^2 / (2l+1)
    return Q3 * np.array([1., 1., -1., -1.])


def _w_lm_l3():
    # W_lm = 4 pi Y_lm(pi/2, 0) / (2l+1), Condon-Shortley Y_lm (Paper 2 below Eq. 9)
    return np.array([4. * np.pi / 7. * np.real(sph_harm_y(3, m, np.pi / 2., 0.)) for m in MM_L3])


def get_ba_l3(ff, M1, M2, R1, k3, wa01, A_nl_1=0.0, w_nl_ref=W_NL_REF, OmegaS1=0, gr_flag=1):
    """$l=3$ f-mode amplitudes $b_a$ (complex), same conventions as `get_ba_l2`.

    ff:     GW frequency grid [Hz], strictly increasing. Start it well below the (3,3)
            resonance (a warning is issued if u(ff[0]) > -5), because F(u) is accumulated from
            `ff[0]` rather than from u = -infinity.
    k3:     dimensionless $l=3$ Love number, paper convention $k_{l0} = 4\\pi I_l^2/(2l+1)$
            (Paper 2 uses 0.02 for SLy).
    wa01:   non-rotating $l=3$ f-mode angular frequency [rad/s].
    A_nl_1: $l=3$ nonlinear frequency-shift coefficient. Defaults to 0: the -3e-2 used for
            $l=2$ is a calibration for the $l=2$ modes and has not been established for $l=3$.

    Returns `(ba1, ba_eq1, va1, Fu1)`. `ba1`, `ba_eq1`, `va1` have shape (4, len(ff)), rows in
    the order `MM_L3`. `ba1 = ba_eq1` plus the dynamical tide on the (3,3) mode only (the other
    $l=3$ modes never resonate before merger). The resonance frequency is the root of
    $\\Delta_a(\\omega) = 0$ by fixed-point iteration, including GR and the nonlinear shift.
    """
    ff = np.asarray(ff, dtype=float)
    orb = _inspiral(ff, M1, M2)
    ww, rr, tMt, dwdt, dlogrdt = (orb[k] for k in ("ww", "rr", "tMt", "dwdt", "dlogrdt"))
    X1, X2 = orb["X1"], orb["X2"]
    zz1, fd1 = _gr_factors(orb["xx"], X1, X2, gr_flag)
    omega_FD_SL1 = fd1 * ww

    l = 3
    ll1 = np.array([l] * len(MM_L3))
    mm1 = np.array(MM_L3)
    wa1 = -mm1 / ll1 * OmegaS1 + wa01           # co-rotating frame frequency

    # eq. 9 of 2404.00147, phase omitted; I_{l,-m} = (-1)^m I_{l,m}
    Wlm1 = _w_lm_l3()
    Ia1 = get_overlaps_l3(k3)
    va1 = np.stack([(M2 / M1) * (R1 / rr)**(l + 1.) * Wlm1[i] * Ia1[i] for i in range(len(MM_L3))])

    # eq. 18 of 2404.00147 including GR & nonlinear tide
    ba_eq1 = np.empty((len(MM_L3), len(ff)), dtype=complex)
    for i in range(len(MM_L3)):
        Del = zz1 * (wa1[i] + mm1[i] * OmegaS1) - mm1[i] * ww + mm1[i] * omega_FD_SL1 \
            + A_nl_1 * (ww / w_nl_ref)**2. * wa01
        ba_eq1[i] = zz1 * Del * wa01 * va1[i] \
            / (Del**2. - 1j * (mm1[i] * dwdt - Del * (ll1[i] + 1.) * dlogrdt))

    # resonance of the prograde (3,3) mode: root of Delta_a(omega) = 0
    m_r = MM_L3[0]
    wa_r = float(wa1[0])
    wres1 = (wa_r + m_r * OmegaS1) / m_r
    for _ in range(40):
        zz_r, fd_r = _gr_factors((tMt * wres1)**(2. / 3.), X1, X2, gr_flag)
        wres1 = (zz_r * (wa_r + m_r * OmegaS1) + A_nl_1 * (wres1 / w_nl_ref)**2. * wa01) \
            / (m_r * (1. - fd_r))
    zz_r, _ = _gr_factors((tMt * wres1)**(2. / 3.), X1, X2, gr_flag)

    ba_dyn, Fu1 = _dynamical_tide(orb, wres1, m_r, l, zz_r, wa01, Wlm1[0], Ia1[0], M1, M2, R1)
    ba1 = ba_eq1.copy()
    ba1[0] += ba_dyn
    return ba1, ba_eq1, va1, Fu1
