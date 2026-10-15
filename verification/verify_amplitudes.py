"""V9 regression (docs/neutron_star_shape.md sec 6): get_ba_l2 / get_ba_l3 against a time-domain
integration of Paper 2 Eq. 15,  db/dt + i Delta b = i omega_a0 V,  along a point-mass inspiral.

Runs for the sculpture's own stars (tools/make_star_modes.py parameters) with gr_flag = 0 and
A_nl = 0, because the ODE contains neither (the generator's GR / nonlinear terms are
exercised by the closed-form pieces only). The ODE is started from b_eq at the first grid
point. Run:  python verification/verify_amplitudes.py [--spin-heavy -500 --spin-light -200]
"""
import argparse
import sys
from pathlib import Path

import numpy as np
from scipy.integrate import solve_ivp
from scipy.special import sph_harm_y

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT / "tools"))
sys.path.insert(0, str(Path(__file__).resolve().parent))

import binary                          # noqa: E402
import tidal_waveform as jw          # noqa: E402
import make_star_modes as MSM          # noqa: E402


def reference(ff, l, m, OmS, wa0, V_of_w, M1, M2, b0):
    Mt = M1 + M2
    eta = M1 * M2 / Mt ** 2
    tMt = jw.G * Mt / jw.c ** 3
    wa = wa0 - (m / l) * OmS
    w = np.pi * ff
    wdot = lambda ww: 96.0 / 5 * eta * (tMt * ww) ** (11.0 / 3) / tMt ** 2

    def rhs(ww, y):
        b = y[0] + 1j * y[1]
        D = wa + m * OmS - m * ww
        db = (-1j * D * b + 1j * wa0 * V_of_w(ww)) / wdot(ww)
        return [db.real, db.imag]

    sol = solve_ivp(rhs, (w[0], w[-1]), [b0.real, b0.imag], t_eval=w, method="DOP853", rtol=1e-9, atol=1e-15)
    return sol.y[0] + 1j * sol.y[1]


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--spin-heavy", type=float, default=MSM.STARS["heavy"][4])
    ap.add_argument("--spin-light", type=float, default=MSM.STARS["light"][4])
    args = ap.parse_args(argv)
    f_snap = binary.solve().f_gw
    nums = MSM.polytrope_numbers()
    worst = 0.0
    for name, spin in (("heavy", args.spin_heavy), ("light", args.spin_light)):
        M, R_km, Mc, _, _ = MSM.STARS[name]
        M1, M2, R = M * jw.Ms, Mc * jw.Ms, R_km * 1e5
        Mt = M1 + M2
        w0 = np.sqrt(jw.G * M1 / R ** 3)
        C = jw.G * M1 / (R * jw.c ** 2)
        k2, k3 = nums[2]["k"], nums[3]["k"]
        wa2, wa3 = nums[2]["w"] * w0, nums[3]["w"] * w0
        OmS = 2 * np.pi * spin
        ff = np.unique(np.concatenate([np.linspace(60.0, 1.05 * f_snap, 20001), [f_snap]]))
        i = int(np.argmin(np.abs(ff - f_snap)))
        ba2, _, _, _ = jw.get_ba_l2(ff, M1, M2, R, k2 / (1.5 * C ** 5), wa2, A_nl_1=0.0, OmegaS1=OmS, gr_flag=0)
        ba3, _, _, _ = jw.get_ba_l3(ff, M1, M2, R, k3, wa3, A_nl_1=0.0, OmegaS1=OmS, gr_flag=0)
        ba2, ba3 = np.asarray(ba2), np.asarray(ba3)
        I2 = np.sqrt(5 * k2 / (4 * np.pi))
        I3 = np.asarray(jw.get_overlaps_l3(k3), dtype=float)
        print(f"\n{name}: spin {spin:.0f} Hz, snapshot index {i} (f = {ff[i]:.2f} Hz)")
        for l, ms, ba, Is, wa0 in ((2, (2, 0, -2), ba2, [I2] * 3, wa2), (3, jw.MM_L3, ba3, I3, wa3)):
            for k, m in enumerate(ms):
                W = 4 * np.pi / (2 * l + 1) * np.real(sph_harm_y(l, m, np.pi / 2, 0.0))
                V = lambda ww, W=W, I=Is[k], l=l: (M2 / M1) * (R / (jw.G * Mt / ww ** 2) ** (1 / 3)) ** (l + 1) * W * I
                ref = reference(ff, l, m, OmS, wa0, V, M1, M2, ba[k, 0])
                err = abs(ba[k, i] - ref[i]) / abs(ref[i])
                worst = max(worst, err)
                print(f"   l={l} m={m:+d}: |b| = {abs(ba[k, i]):.4e}  ODE {abs(ref[i]):.4e}  rel err at snapshot = {err:.2e}")
    print(f"\nworst relative error at the snapshot: {worst:.2e}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
