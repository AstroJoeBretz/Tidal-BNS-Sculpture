# Tidally Deformed Neutron Star Surface — Shape Specification

Scope: turn the f-mode amplitudes of Yu, Weinberg, Arras, Kwon & Venumadhav
(arXiv:2211.07002, "Paper 1") and Yu, Arras & Weinberg (arXiv:2404.00147, "Paper 2")
into a deformed stellar surface $r(\theta,\phi)$. The derivation of the whole amplitude-to-surface
chain, with a worked example reproduced by `verification/math_walkthrough.py`, is in
[`math.md`](math.md); this document is the specification and its checks.

Out of scope here: STL/mesh generation, sockets, posts and the wood panel
([`design.md`](design.md) §6.3–§7, `src/star.py`); only their *input* comes from this document.

Implemented in `src/tidal_shape.py`, `src/tidal_waveform.py` and `tools/make_star_modes.py`;
results and checks in §9.

---

## 0. Design choices

| Topic | Decision |
|---|---|
| Amplitudes | **Complex** rotating-frame $C_a$ for the $l=2$ f-modes $a=(2,+2)$, $b=(2,-2)$, $c=(2,0)$, at one chosen separation $r/R$ (and $l=3$ modes, §4). They are the **total amplitude** `ba1` $=b_a^{\rm(eq)}+b_a^{\rm(dyn)}$ from `src/tidal_waveform.py`: the equilibrium component of Paper 2 Eq. 18 (its imaginary part comes from orbital decay, §3.3) **plus** the resonantly excited dynamical component (Paper 2 Eqs. 23–24, 30). The dynamical part is nonzero only past a mode's resonance and makes the shape a snapshot (§3.5). |
| Surface eigenfunction | **Calibrated from the Love number** (§3.4). |
| Accuracy | **Linear displacement only**: $r=R+\xi_r(R,\theta,\phi)$. |
| $l=3$ | Not in Paper 1: the $l=3$ f-mode amplitudes come from `get_ba_l3`, the Paper 2 formalism extended to $l=3$; $h_3$, $I_3$ are computed for the $\Gamma=2$ polytrope (§4). |
| D1 exaggeration | **$X=1$** default, kept as a parameter. |
| D2 volume | **Do not rescale**; report the volume drift. |
| D3 stellar model | Use the paper's $\Gamma=2$ polytrope numbers ($k_2\approx0.26$, $I_2\approx0.32$; $h$ values in §4). |
| D4 two stars | A separate amplitude set **per star** (own $M'/M$, $r/R$). |
| D5 complex | Complex $C_a$ accepted; reality of the *surface* is automatic (§3.3). |
| D6 lag exaggeration | **Include $X_{\rm lag}$, default 1.** Multiplies the lag phase only (§3.3.1). |
| Amplitude source | `get_ba_l2` / `get_ba_l3` in `src/tidal_waveform.py`, **total** `ba1` (equilibrium + dynamical tide). |
| Spin (D9) | Anti-aligned, **heavy $-500$ Hz, light $-200$ Hz** (per-star `spin_hz`). Supported through the amplitudes at no cost to the shape pipeline (§3.6). Heavy passes V11 with one warning ($l=3$ margin $0.13<0.15$); light passes cleanly. $-500$ Hz on the light star was rejected (outside the linear regime, §8). |
| D7 spin oblateness | **Optional switch** `spin_oblateness`; library default **off**, formula and checks in §3.7. |
| D10 | **`spin_oblateness = true` for the sculpture build** (any build with nonzero spin), labelled as an addition to the papers' physics. |
| D8 snapshot | **Sculpture's own $f_{\rm GW}=1023.6$ Hz** for both stars ($r/R=2.75$ heavy, $2.64$ light), `design.md` §1.1. |

---

## 1. Inputs (per star)

| Symbol | Meaning | Source |
|---|---|---|
| $R$ | unperturbed radius (mm on the sculpture, via `binary.py`) | `params.py` |
| $I_2$ | $l=2$ f-mode tidal overlap (Paper 1 Eq. 5 / Paper 2 Eq. 10), **the value used when computing $C_a$** | polytrope, 0.3215 (§4) |
| $h_2$ | equilibrium shape number, $1+2k_2$ (§3.4) | polytrope: **1.51982** |
| $C_a,C_b,C_c$ | **complex** total `ba1` (equilibrium + dynamical tide) for $(2,+2)$, $(2,-2)$, $(2,0)$ | `tidal_waveform.get_ba_l2` |
| $(l,m,C,I_l,h_l)$ | extra multipoles ($l=3\ldots$), same meaning; $l+m$ even only | `get_ba_l3`; $I_3,h_3$ in §4 |
| $M$, $M'$ | star and companion mass; $r/R$ follows from $a$ and $f_{\rm GW}$ | `design.md` §1.1 |
| $f_{\rm GW}$ | snapshot frequency, **1023.6 Hz** (D8) | `design.md` §1.1 |
| `spin_hz` | stellar spin, negative = anti-aligned (heavy −500, light −200; D9) | design choice |
| `spin_oblateness` | add the centrifugal $l=2$ term of §3.7 (library default off; sculpture build on, D10) | design choice |
| $X$ | multiplier on the tidal deformation (default **1**) | design choice |
| $X_{\rm lag}$ | multiplier on the lag phase (default **1**, §3.3.1) | design choice |

### 1a. Amplitude generator (`src/tidal_waveform.py`)

| | `get_ba_l2` | `get_ba_l3` |
|---|---|---|
| Inputs (cgs) | `ff, M1, M2, R1, lam_M1_5` ($=\Lambda_2$)`, wa01, A_nl_1, w_nl_ref, OmegaS1, gr_flag` | `ff, M1, M2, R1, k3, wa01, A_nl_1=0, w_nl_ref, OmegaS1, gr_flag` |
| Mode order of the returned arrays | $[m=+2,0,-2]$ | $[m=+3,+1,-1,-3]$ |
| Returns | `ba1` (total), `ba_eq1`, `va1`, `Fu1` | same |
| $I$ | $\sqrt{5k_2/4\pi}$, same sign for all $m$ | `get_overlaps_l3(k3)`, **signed**: $I_{3,-m}=(-1)^m I_{3,m}$ |
| Dynamical tide | on $m=+2$ | on $m=+3$ only |
| Resonance frequency | closed-form 1PN expression | root of $\Delta_a=0$ by fixed-point iteration, GR and nonlinear shift included |

- **Mapping to this spec** ($l=2$): $C_a=$ `ba1[0]`, $C_c=$ `ba1[1]`, $C_b=$ `ba1[2]`.
- **The shape uses $\hat C_a=C_a/I_a$ mode by mode** (§3.3, §4). For $l=3$ the sign of $I_a$ differs
  between $+m$ and $-m$, so `get_overlaps_l3` must be used with `get_ba_l3`.
- `ff` must be strictly increasing and start well below resonance (the dynamical-tide integral
  $F(u)$ starts at `ff[0]`). Both functions warn when $u(\texttt{ff[0]})>-5$.
- Arrays are `complex128`.
- `A_nl_1` for $l=3$ defaults to 0; the $-3\times10^{-2}$ of $l=2$ is an $l=2$ calibration. For the sculpture's $l=2$ modes
  `A_nl_1` is derived per star (§8), not left at the $-3\times10^{-2}$ default.
- `k3` is the paper-convention Love number ($k_{l0}=4\pi I^2/(2l+1)$). From the dimensionless
  deformability $\Lambda_3=\tfrac{2}{15}\,k_3\,(Rc^2/GM)^7$, $k_3=7.5\,\Lambda_3\,(GM/Rc^2)^7$ (the analogue of
  $k_2=1.5\,\Lambda_2\,C^5$ used by `get_ba_l2`).

---

## 2. Conventions

- **Star frame**, origin at the star's centre, right-handed. $+z$ = orbital angular
  momentum. Companion on the **$+x$ axis** at the snapshot ($\theta=\pi/2$, $\phi_s=0$). The
  orbit turns in the $+\phi$ sense.
- **Corotating azimuth** $\phi_s=\phi-\phi_{\rm orb}$. Mode amplitudes obey
  $q_a=b_a\,e^{-im_a\phi_{\rm orb}}$ (Paper 2 Eq. 14 with zero stellar spin; Paper 1
  $C_a=c_a e^{im_a\phi}$), so $q_a\,Y_{lm}(\theta,\phi)=b_a\,Y_{lm}(\theta,\phi_s)$. We write $C_a\equiv b_a$.
- $Y_{lm}$: orthonormal, Condon–Shortley phase (`scipy.special.sph_harm_y`). With it
  $W_{lm}=4\pi\,Y_{lm}(\pi/2,0)/(2l+1)$ gives both papers' $W_{22}=\sqrt{3\pi/10}$, $W_{20}=-\sqrt{\pi/5}$.
- **Nonspinning star** ($\Omega_{\rm spin}=0$): Paper 2's $\Delta_a=\omega_a-m_a\sigma$ with $\sigma=\omega-\Omega_{\rm spin}$
  becomes $\Delta_a=\omega_a-m_a\omega$ ($\omega=\dot\phi$ orbital), the same denominator as Paper 1 Eq. 10.
- Both papers: $G=c=1$; normalization $2\omega^2\!\int\!\rho\,\lvert\boldsymbol\xi\rvert^2\,d^3x=E_0=M^2/R$ (Paper 1 Eq. 3,
  Paper 2 Eq. 8); same $I_a$, same $W_{lm}$, same $U_a$ (Paper 2 Eq. 9 $\equiv$ Paper 1 Eq. 4).
  So $(\star)$ below applies to both. All module quantities are dimensionless.
- **Paper's Love number** $k_l$ is the Hinderer/Damour–Nagar one: $k_{l0}=4\pi I_l^2/(2l+1)$
  (Paper 2 Eq. 43). It is **half** the classical Love number: $k_l^{\rm Love}=2k_l$.

---

## 3. Theory: from $C_a$ to the surface

### 3.1 Real displacement field

Both papers expand $\boldsymbol\xi=\sum_a q_a\boldsymbol\xi_a$ over modes of both frequency signs, with
$(-\omega,-m)$ the complex conjugate of $(\omega,m)$ ($q_{-a}=q_a^*$, $\boldsymbol\xi_{-a}=\boldsymbol\xi_a^*$). With the
normalization above the mode energy is $E_0\lvert q_a\rvert^2$ precisely when the real field is
$\boldsymbol\xi=q_a\boldsymbol\xi_a+\mathrm{c.c.}$, so over positive-frequency modes

$$
\boldsymbol\xi(\mathbf x,t)=\sum_{\omega_a>0}2\,\mathrm{Re}\!\left[q_a(t)\,\boldsymbol\xi_a(\mathbf x)\right],\qquad
\boldsymbol\xi_a=\xi_{r,a}(r)\,Y_{lm}\,\hat{\mathbf r}+\xi_{h,a}(r)\,r\nabla Y_{lm}
$$

The factor 2 is not optional (§6 V3).

### 3.2 Surface radius to linear order

The surface is a material surface; to first order its radial position is
$R+\xi_r(R,\theta,\phi)$. The tangential part $\xi_h$ slides fluid along a sphere and does not
change the shape at this order, so $\xi_h(R)$ is never needed.

$$
\frac{\Delta R(\theta,\phi_s)}{R}=\frac{2}{R}\sum_{\omega_a>0}\xi_{r,a}(R)\,\mathrm{Re}\!\left[C_a\,Y_{lm_a}(\theta,\phi_s)\right]
\tag{3.2}
$$

### 3.3 The $l=2$ shape with complex amplitudes

**Where the imaginary part comes from.** Paper 2 Eq. 18:

$$
b_a^{\rm(eq)}=\frac{\Delta_a\,\omega_{a0}\,V_a}{\Delta_a^{2}-i\,\kappa_a},\qquad
\kappa_a=m_a\dot\omega-\Delta_a\,(l_a+1)\,\frac{\dot r}{r}
$$

Without orbital decay ($\dot\omega=\dot r=0$) this reduces to Paper 1's real
$\omega_{a0}V_a/\Delta_a$. With decay the bracket $\kappa_a$ makes $b_a$ complex: the mode's response
trails the forcing, which is a **tidal lag**. For a point-mass inspiral
($\dot r/r=-\tfrac23\,\dot\omega/\omega$), $\kappa_a=m_a\dot\omega+2\Delta_a\dot\omega/\omega$ is positive for both
$m=+2$ and $m=-2$, and $\mathrm{Im}\,b_a>0$ for both.

**Shape.** Write $N_2=\tfrac14\sqrt{15/2\pi}$, $N_0=\tfrac14\sqrt{5/\pi}$, and $\mathbf n=(n_x,n_y,n_z)$ the unit vector.
$Y_{22}=N_2\sin^2\theta\,e^{2i\phi_s}$, $Y_{2,-2}=Y_{22}^*$, $Y_{20}=N_0(3\cos^2\theta-1)$ (real). Then

$$
\mathrm{Re}\!\left[C_aY_{22}+C_bY_{22}^*\right]=N_2\sin^2\theta\left[S_c\cos2\phi_s-S_s\sin2\phi_s\right],
\qquad S_c=\mathrm{Re}(C_a+C_b),\quad S_s=\mathrm{Im}(C_a-C_b)
$$

$$
\mathrm{Re}\!\left[C_cY_{20}\right]=N_0\,\mathrm{Re}(C_c)\,(3\cos^2\theta-1)\qquad(\mathrm{Im}\,C_c\ \text{drops out})
$$

so, with $K=h_2/I_2$ and using $\sin^2\theta\cos2\phi=n_x^2-n_y^2$, $\sin^2\theta\sin2\phi=2n_xn_y$:

$$
\frac{\Delta R}{R}=K\Big\{N_2S_c\,(n_x^2-n_y^2)-2N_2S_s\,n_xn_y+N_0\,\mathrm{Re}(C_c)\,(2n_z^2-n_x^2-n_y^2)\Big\}
$$

This is a traceless quadratic form: a **triaxial** shape with the $z$ axis principal and
the $x$–$y$ axes rotated by the lag angle. Equatorial profile:

$$
\frac{\Delta R(\phi_s)}{R}\Big|_{\theta=\pi/2}=K\Big\{N_2\,\rho\cos(2\phi_s+\delta)-N_0\,\mathrm{Re}(C_c)\Big\},\qquad
\rho\cos\delta=S_c,\ \ \rho\sin\delta=S_s
$$

$$
\psi\equiv\frac{\delta}{2}=\tfrac12\,\mathrm{atan2}(S_s,S_c),\qquad\text{bulge maximum at }\ \phi_s=-\psi
$$

$\psi>0$ means the long axis trails **behind** the companion (companion at $\phi_s=0$, orbit
turning in $+\phi$). Checks:

- $\psi=0$ and the form reduces to the real-amplitude result (Paper 1) when $\mathrm{Im}\,C_a=\mathrm{Im}\,C_b=0$.
- Real amplitudes give $S_s=0$: mirror symmetry about the $x$–$z$ plane.
- The $m=0$ amplitude's imaginary part encodes fluid velocity ($\dot{\boldsymbol\xi}\propto\omega\,\mathrm{Im}\,q$), not shape,
  because $Y_{20}$ is real; it is correctly ignored here.
- Principal axes are $z$ and the in-plane directions at $\phi_s=-\psi$ and $-\psi+\pi/2$, which is
  also how `star.py` orients the mesh.

**Axisymmetry about the line of centres** requires both
$S_s=0$ and $(C_a+C_b)/\mathrm{Re}\,C_c=-\sqrt6$, the adiabatic ratio $2W_{22}/W_{20}$. At finite
orbital frequency the $m=\pm2$ response is enhanced relative to $m=0$, so neither holds in
general, and the surface is triaxial.

**Size of the lag, from Eq. 18 with a point-mass inspiral** ($1.3+1.3\,M_\odot$, $R=12$ km,
compactness 0.16, $\omega_a=1.2\,\omega_0$ as quoted in Paper 1 (the solved polytrope value used in §8 is $1.227\,\omega_0$),
$\Omega^2=M_t/r^3$): a computed illustration, not a spec value.

| $r/R$ | $\psi$ (bulge behind companion) |
|---|---|
| 4.0 | 0.06° |
| 3.0 | 0.40° |
| 2.75 | 0.75° |
| 2.5 | 1.60° |

So the physical lag is small. It is carried correctly, but it will be nearly invisible on a
Ø17 mm star unless exaggerated (§3.3.1).

### 3.3.1 Lag exaggeration $X_{\rm lag}$ (D6, default 1)

For any $(l,\lvert m\rvert)$ pair define the **combined complex amplitude** (using $\hat C_a=C_a/I_a$)

$$
D_m=\hat C_{+m}+(-1)^m\,\hat C_{-m}^{*},\qquad
\mathrm{Re}\!\left[\hat C_{+m}Y_{lm}+\hat C_{-m}Y_{l,-m}\right]=\mathrm{Re}\!\left[D_m\,Y_{lm}\right]
$$

($Y_{l,-m}=(-1)^mY_{lm}^*$). For $l=2$, $m=2$: $D_2=(C_a+C_b^*)/I_2$, i.e. $S_c+iS_s$
of §3.3 divided by $I_2$ (and $\psi$ of §3.3 is computed from $\hat C=C/I$, i.e. with $S_c,S_s$ multiplied by
$\mathrm{sign}(I)$, so that it does not depend on the sign convention of the eigenfunctions).

**Zero-lag reference.** The lag must be measured from the *adiabatic* orientation, whose amplitude
$D_m$ is real but not always positive: $Y_{33}(\pi/2,0)<0$ (Condon–Shortley), so an adiabatic $(3,3)$
amplitude is negative real. Define $s_{lm}=\mathrm{sign}\,Y_{lm}(\pi/2,0)$ and $\tilde D_m=s_{lm}D_m$. The equatorial
pattern of that harmonic is $\lvert D_m\rvert\,\lvert Y_{lm}\rvert\cos(m\phi_s+\arg\tilde D_m)$, so its
lag is

$$
\psi_m=\frac{\arg\tilde D_m}{m}\qquad(\text{behind the companion for }\psi_m>0),\qquad m=0:\ \text{no lag.}
$$

($s_{22}=s_{31}=+1$, $s_{33}=-1$.) Measuring $\arg D_m$ instead would report a spurious $60^\circ$ lag for an
adiabatic $(3,3)$ mode and rotate it when $X_{\rm lag}\neq1$ (test V10).

The knob rotates the pattern and leaves its amplitude alone:

$$
\tilde D_m\ \longrightarrow\ \lvert D_m\rvert\,\exp\!\big(i\,X_{\rm lag}\,\arg\tilde D_m\big)\qquad(X_{\rm lag}=1\ \text{reproduces the physics})
$$

- It is applied per $(l,\lvert m\rvert)$ pair, so every harmonic's lag is scaled by the same factor.
- It changes the shape only through the $n_xn_y$ cross term ($l=2$) and the equivalent terms
  for $l=3$. $\lvert D_m\rvert$ (hence the semi-axis lengths of the $x$–$y$ pattern) is unchanged.
- $X_{\rm lag}$ and $X$ (§5) are independent. $X_{\rm lag}$ must be reported in any output label.
- For large $X_{\rm lag}\,\psi$ the pattern simply points elsewhere; $\arg\tilde D_m$ is taken on
  $(-\pi,\pi]$ and not unwrapped, so lag that is physically $>180^\circ/m$ aliases (relevant only near
  resonance, §3.6).

### 3.4 Surface eigenfunction by Love-number calibration

The **equilibrium shape number** $h_l$: in the adiabatic limit the surface is an
equipotential of (tidal + induced) potential, so

$$
\frac{\Delta R}{R}\Big|_{\rm adiabatic}=h_l\,\frac{M'}{M}\Big(\frac{R}{r}\Big)^{l+1}P_l(\cos\gamma),\qquad
h_l=1+k_l^{\rm Love}=1+2k_l
$$

($\gamma$ = angle from the companion; $k_l$ = the paper's Love number.)
In the same limit $C_a\to V_a=\frac{M'}{M}\,W_{lm}\,I\,(R/r)^{l+1}$ and
$\sum_m W_{lm}Y_{lm}(\theta,\phi_s)=P_l(\cos\gamma)$, so (3.2) reduces to
$\Delta R/R=\frac{2\xi_r(R)}{R}\,I\,\frac{M'}{M}(R/r)^{l+1}P_l$. Equating:

$$
\frac{\xi_{r,a}(R)}{R}=\frac{h_l}{2I_l}\qquad\Longrightarrow\qquad
\frac{\Delta R}{R}=\frac{h_l}{I_l}\,\mathrm{Re}\sum_aC_a\,Y_{lm_a}(\theta,\phi_s)
\tag{\star}
$$

- $C_a\propto I_l$, and $(\star)$ uses $C_a/I_l$, so the sign convention of the eigenfunction
  normalization cancels (Paper 1 footnote 1).
- Exact adiabatic limit by construction (V1). Dynamics and lag live entirely in $C_a$.
- Approximation: the f-mode alone gives an adiabatic surface sum a few % above the
  full value ($l=2$: 1.551 vs 1.520; $l=3$: 1.256 vs 1.213); the p-modes subtract the
  difference. $(\star)$ folds them into the f-mode, so its *effective* $\xi_r(R)$ is 2% ($l=2$)
  or 3.5% ($l=3$) below the f-mode's true value. P-modes sit far above the driving
  frequency and respond adiabatically, so this is the right place for them.

### 3.5 Dynamical tide: the shape is a snapshot

The amplitudes include the resonantly excited component, `ba1 = ba_eq1 + ba_dyn`. The
equilibrium part is static in the orbit-corotating frame; **`ba_dyn` is not**: it oscillates at the
mode frequency, so after resonance the $m=\pm2$ pattern keeps changing relative to the companion.
Consequences:

- The surface is defined **at one orbital frequency** ($f_{\rm GW}$ or $r/R$) on the generator's grid; the spec's
  input is the pair $(f_{\rm GW},\ \texttt{ba1[:, i]})$. The snapshot is fixed at the sculpture's $f_{\rm GW}=1023.6$ Hz (D8).
- Before resonance `ba_dyn` is small but not zero: at the final snapshot it is 0.3–1.2% of $\lvert C\rvert$ for the
  prograde modes (§9, and `math.md` §3.5); the equilibrium part carries the rest.
- After resonance $\psi$ is not a smooth function of $f_{\rm GW}$: for the SLy-like star of §3.6 at spin $-850$ Hz,
  $\psi_{\rm total}=11^\circ,33^\circ,73^\circ,-28^\circ,54^\circ$ at $f_{\rm GW}=800,1000,1200,1500,1800$ Hz. A sculpture frozen there is
  a legitimate snapshot, but a small change of $f_{\rm GW}$ changes the lag substantially.
- The $F(u)$ integral starts at `ff[0]`, so the grid must begin well below resonance.

### 3.6 Anti-aligned neutron-star spin

Anti-aligned spin needs no change to the shape pipeline: the generator already takes
`OmegaS1 < 0` for anti-aligned spin (Paper 2 sign).

Why the shape pipeline is unchanged: with spin, Paper 2 writes $q_a=b_a\,e^{-im(\phi-\Omega_st)}$
in the star's corotating frame, whose azimuth is $\phi_c=\phi_{\rm inertial}-\Omega_st$. Then
$q_a\,Y_{lm}(\theta,\phi_c)=b_a\,Y_{lm}(\theta,\phi_{\rm inertial}-\phi_{\rm orb})$, which is exactly the orbit-frame
form of §2. Spin enters only through $b_a$ (via $\Delta_a=(\omega_a+m\Omega_s)-m\omega$ and $\omega_a=\omega_{a0}-(m/l)\Omega_s$).

**What anti-aligned spin does to the lag** (`get_ba_l2`, `gr_flag=1`, $\omega_{a0}=2\pi\cdot1.8$ kHz,
$k_2=0.08$, $1.4+1.4\,M_\odot$, $R=11.7$ km; snapshot at $f_{\rm GW}=1023.6$ Hz, the sculpture's value, $r/R=2.8$).
Resonance of the $(2,2)$ mode moves down to $f_{\rm res}=(\omega_{a0}+\Omega_s)/2\pi$. $\Delta_a/\omega_a$ below uses the V11 normalization
($\Delta_a=\omega_a^{\rm corot}+2\Omega_s-2\omega$, $\omega_a^{\rm corot}=\omega_{a0}-\Omega_s$; Newtonian, no GR):

| spin [Hz] | $f_{\rm res}$ [Hz] | $\Delta_a/\omega_a$ at snapshot | $\psi$ (equilibrium) | $\psi$ (total) |
|---|---|---|---|---|
| 0 | 1800 | 0.43 | 1.5° | 1.5° |
| −200 | 1600 | 0.29 | 2.6° | 2.4° |
| −300 | 1500 | 0.23 | 3.6° | 3.9° |
| −400 | 1400 | 0.17 | 5.5° | 6.1° |
| −500 | 1300 | 0.12 | 9.6° | 9.3° |
| −600 | 1200 | 0.07 | 20° | 15° |

Lag grows as the snapshot approaches resonance: $6\times$ at $-500$ Hz, $13\times$ at $-600$ Hz. At fixed spin $-850$ Hz, $\psi_{\rm eq}$ at 500 / 800 Hz is
$0.4^\circ/11.8^\circ$ versus $0.06^\circ/0.4^\circ$ without spin
(about $28\times$ at 800 Hz), and the $(2,2)$ amplitude rises $6\times$ ($\lvert b_{\rm eq}\rvert/\lvert V\rvert$: 11.3 versus 1.8).

**Caveats, in order of importance:**

1. **Do not snapshot past resonance for a clean pattern.** Beyond $f_{\rm res}$ the equilibrium response
   changes sign ($\Delta_a<0$): the bulge toward the companion turns into a *depression*
   (at $-850$ Hz, the $x$-deviation goes $+0.142\to-0.038$ between 800 and 1000 Hz) and $\psi_{\rm eq}\to-90^\circ$.
   Practical window: $\Delta_a/\omega_a\gtrsim0.15$, i.e. spin $\gtrsim-450$ Hz for this SLy-like star at this separation. The
   sculpture's stars have different mode frequencies and tidal strengths; their limits are in §8 (the light
   star is limited to about $-250$ Hz).
2. **Centrifugal oblateness is not in the paper's model.** The amplitudes are first-order in
   $\Omega_s/\omega_1$; the star's rotational flattening is ignored. It is comparable to the tidal bulge at
   large spin and is provided as an optional term (§3.7, D7).
3. **Near resonance Eq. 18 is a first-order truncation**: errors versus the exact ODE reach 5–8% for
   $\Delta_a/\omega_a\approx0.06$–$0.17$ (§6 V9). Fine for a sculpture; do not quote the numbers.
4. **Eigenfunctions are not corrected for rotation**; $(\star)$ holds only to $O(\Omega_s/\omega_1)$.

### 3.7 Optional spin oblateness (D7)

Not part of the papers' model; an additive, clearly labelled extension, enabled by
`spin_oblateness = true` (default **false**).

The centrifugal potential per unit mass is
$\Phi_c=-\tfrac12\Omega_s^2r^2\sin^2\theta=-\tfrac13\Omega_s^2r^2\left[1-P_2(\cos\theta)\right]$; its $l=2$ part is
$+\tfrac13\Omega_s^2r^2P_2(\cos\theta)$. The star's equilibrium response to any static $l=2$ potential is the same shape number
$h_2$ as for the tide (§3.4), so

$$
\frac{\Delta R_{\rm spin}}{R}=-\,\varepsilon_\Omega\,P_2(\cos\theta),\qquad
\varepsilon_\Omega=\frac{h_2}{3}\,\frac{\Omega_s^2R^3}{GM}\qquad(\theta\ \text{from the spin axis, }\pm z)
$$

Equivalently, in the quadratic-form language of §3.3: $e_{xx}\mathrel{+}=\varepsilon_\Omega/2$,
$e_{yy}\mathrel{+}=\varepsilon_\Omega/2$, $e_{zz}\mathrel{-}=\varepsilon_\Omega$.
Properties:

- **Axisymmetric about $z$**, so it is independent of the sign of the spin and of the lag.
  Poles $-\varepsilon_\Omega$, equator $+\varepsilon_\Omega/2$.
- **Check against Maclaurin.** A uniform-density star has $h_2=5/2$ and the formula gives
  $R_{\rm eq}-R_{\rm pol}=\tfrac32\,\varepsilon_\Omega R=\tfrac54\,qR$, $q=\Omega^2R^3/GM$; the classical small-rotation result is
  flattening $5q/4$. For the $n=1$ polytrope ($h_2=1.51982$) it gives $0.76\,q$.
- **$l=0$ part ignored** (a uniform $O(\Omega_s^2)$ radius shift); it appears only through the reported $\Delta V/V$.
- Added to the tidal terms; neither $X$ nor $X_{\rm lag}$ scales it ($X$ exaggerates the tidal deformation only).
- Uses the same $\Omega_s$ as the amplitude generator ($2\pi\times$`spin_hz`).

Size for the §8 stars ($\varepsilon_\Omega\propto\Omega_s^2$): at $-500$ Hz, $0.041$ (heavy) and $0.061$ (light); at $-200$ Hz, $0.0098$ (light).
**In the final configuration (heavy $-500$ Hz, light $-200$ Hz): heavy $0.041$, light $0.0098$.** At $-850$ Hz,
$\approx0.10$ for the $1.4\,M_\odot$ SLy-like star of §3.6 ($h_2=1.16$). Compare the $l=2$ tidal bulge along $x$ at spin 0
(heavy $+0.07$, light $+0.16$): for the heavy star the oblateness is $\sim4\%$ of $R$ at the poles, not negligible.

---

## 4. Higher multipoles ($l=3$)

Same derivation for any $l$:

$$
\frac{\Delta R_l}{R}=\frac{h_l}{I_l}\,\mathrm{Re}\sum_{m_a}C_{l,m_a}\,Y_{lm_a}(\theta,\phi_s),\qquad m_a\equiv l\pmod 2
$$

Only $l+m$ even modes couple to a companion in the orbital plane
($Y_{lm}(\pi/2,0)=0$ otherwise): for $l=3$, $m=\pm1,\pm3$. The amplitudes are complex with the
same lag mechanism (Eq. 18, $l_a=3$); the code uses spherical harmonics directly for
$l\ge3$, no closed Cartesian form is needed.

**Computed for the paper's $\Gamma=2$ ($n=1$) Newtonian polytrope**
(`src/polytrope.py`), paper normalization, $G=M=R=1$:

| | $l=2$ | **$l=3$** | $l=4$ |
|---|---|---|---|
| f-mode $\omega/\omega_0$ | 1.227 | 1.698 | 2.037 |
| f-mode $I_l$ | 0.3215 | **0.2434** | 0.2076 |
| $k_l=\frac{4\pi}{2l+1}\sum I^2$ (mode sum) | 0.25991 | **0.10645** | 0.06024 |
| $k_l$ exact static | 0.25991 | **0.10645** | 0.06024 |
| $h_l$ exact static $=1+2k_l$ | 1.51982 | **1.21291** | 1.12048 |
| $h_l$ mode sum $\sum 2\xi_r(R)\,I/R$ | 1.5198 | **1.2129** | 1.1205 |
| f-mode $\xi_r(R)/R$ (true) | 2.411 | 2.581 | 2.813 |
| $h_l/(2I_l)$ (calibrated, $(\star)$) | 2.363 | **2.492** | 2.699 |

**$h_3=1.21291$, $I_3=0.2434$ ($k_3=0.10645$).** Two independent derivations agree to
4–5 digits: the exact static tide (Helmholtz solution $\psi\propto j_l(\pi r/R)$, matched to the exterior
potential) and the mode sum. This also confirms that $h_l=1+2k_l$ holds for $l=2,3,4$, not just $l=2$.

Sign check: $I_3$ is negative in the script's eigenfunction sign convention and
$I_2$ is negative likewise; only $C/I$ enters $(\star)$, so the sign is irrelevant provided the
$I$ is the one the $C_a$ were computed with. Paper 1 Table 1 and Paper 2 quote
magnitudes.

**Dynamical tide.** For a nonspinning NS Paper 2 notes the $l=3$ f-mode can still be resonantly
excited near the ISCO ($\Omega_{\rm crit}=+2\pi\times5$ Hz). `get_ba_l3` includes the dynamical component on the
$(3,3)$ mode, so the $l=3$ part of the shape is a snapshot in the same sense as §3.5. For the
SLy-like parameters the $(3,3)$ resonance sits at $f_{\rm GW}=1667$ Hz (nonspinning, Newtonian), 1772 Hz
with GR, and 533 Hz at spin $-850$ Hz.

Code requirements: reject $(l,m)$ with $l+m$ odd; require $I_l$ and $h_l$ for every
supplied $l$; default $h_2,h_3$ and the $I_l$ from the table above.

---

## 5. From $\Delta R/R$ to a printable surface

$$
r(\theta,\phi_s)=R\left[\,1+X\,\frac{\Delta R_{\rm tide}(\theta,\phi_s)}{R}+s_{\rm obl}\,\frac{\Delta R_{\rm spin}(\theta)}{R}\right]
$$

- library defaults (D1, D2, D6, D7): $X=1$, $X_{\rm lag}=1$, $s_{\rm obl}=0$, $R_0=R$;
- sculpture build (D10): $X=1$, $X_{\rm lag}=1$, $s_{\rm obl}=1$.

Here $\Delta R_{\rm tide}/R$ is $(\star)$ evaluated with the lag-exaggerated $D_m$ of §3.3.1, and $\Delta R_{\rm spin}/R$ is §3.7
($s_{\rm obl}=1$ when `spin_oblateness = true`).

- $X$ multiplies the deformation (all of $C_a$); $X_{\rm lag}$ multiplies the lag phase only.
- Star-shaped surface; assert $r>0$ everywhere and report $\min r/\max r$.
- At first order no $l=0$ component ($\oint\Delta R\,d\Omega=0$, V2). At second order the enclosed
  volume changes; the code **reports $\Delta V/V$** and does not correct it (D2).
- **Companion on the $-x$ side** (the light star): compute with the companion on $+x$, then map
  $(x,y,z)\to(-x,-y,z)$ (rotation by $\pi$ about $z$). This preserves the lag direction relative
  to the companion. $l=2$ terms are unchanged by the map; odd-$l$ terms flip sign.

---

## 6. Validation

| ID | Check | Expected |
|---|---|---|
| V1 | Adiabatic, real: $C_a=C_b=V_2$, $C_c=V_0$ | $\Delta R/R=h_2\,\frac{M'}{M}(R/r)^3P_2(\cos\gamma)$ to $10^{-12}$ |
| V2 | Mean of $\Delta R$ over sphere | 0, to $10^{-12}$ |
| V3 | $n=1$ polytrope sum rules (`src/polytrope.py`, $l=2,3,4$) | $\frac{4\pi}{2l+1}\sum I^2=k_l$; $\sum2\xi_r(R)I/R=h_l=1+2k_l$; exact-static $h_l$ agrees |
| V4 | Symmetry | real $C_a$: invariant under $y\to-y$, $z\to-z$; $l=2$ invariant under $(x,y)\to(-x,-y)$ |
| V5 | Closed form vs spherical harmonics, **complex** $C_a$ | quadratic form (§3.3) equals $\mathrm{Re}\sum CY$ to $10^{-15}$ (already verified for Eq. 18 amplitudes) |
| V6 | Principal axes | $z$ axis principal; in-plane axes at $-\psi$, $-\psi+\pi/2$; code reports $\psi$ and $(e_{xx},e_{yy},e_{zz})$ |
| V7 | Mesh | $r>0$, tip-to-tip gap $>0$ for the two-star assembly, $\Delta V/V$ reported |
| V8 | **Lag sign** | for Eq. 18 with $\dot\omega>0$ (inspiral): $\psi>0$ and $\psi\to0$ as $\dot\omega\to0$ |
| V9 | **Amplitude generator vs Paper 2 Eq. 15** (below) | `get_ba_l2/l3` vs time-domain ODE |
| V10 | $X_{\rm lag}$ | $X_{\rm lag}=1$ identical to no knob; $\lvert D_m\rvert$ unchanged for any $X_{\rm lag}$; $\psi_m\to X_{\rm lag}\psi_m$ for $\lvert X_{\rm lag}\psi_m\rvert<180^\circ/m$; adiabatic limit has $\psi_m=0$ for **every** $(l,m)$ incl. $(3,3)$ |
| V11 | **Linear-regime gate**, per star, on the snapshot | **error** if any resonant mode has $\Delta_a/\omega_a^{\rm corot}<0$ (past resonance) unless `allow_post_resonance`; **warn** if $\Delta_a/\omega_a^{\rm corot}<0.15$ for any $l=2,3$ mode, if $\max\lvert\Delta R/R\rvert>0.3$, or if $\lvert\Delta V/V\rvert>0.05$; with $\omega_a^{\rm corot}=\omega_{a0}-(m/l)\Omega_s$, $\Delta_a=\omega_a^{\rm corot}+m\Omega_s-m\omega$ (**Newtonian**: no GR / nonlinear terms; they shift the margin by $\le0.05$ for the resonant modes here, see `math.md` §7) |
| V12 | Spin oblateness (§3.7) | off $\Rightarrow$ output identical to the tide-only surface; on $\Rightarrow$ axisymmetric about $z$; Maclaurin check $R_{\rm eq}-R_{\rm pol}=\frac{h_2}{2}\,qR$ with $h_2=5/2$ gives $5q/4$; independent of the sign of the spin |

**V3 evidence** ($l=2$): $I=0.3215$ (paper 0.32), p$_1$ $I=-5.5\times10^{-3}$
(paper $-5.5\times10^{-3}$), $\sum\frac{4\pi}{5}I^2=0.25991=k_2$ exact; $\sum2\xi_r(R)\,I/R=1.5198=1+2k_2$. The
last fixes $h_2=1+2k_2$ and the factor 2 in §3.1.

**V9 evidence.** Paper 2 Eq. 15 is integrated in time along a point-mass
inspiral ($\dot b+i\Delta\,b=i\omega_{a0}V$, DOP853, rtol $10^{-9}$, started from $b_{\rm eq}$ at 60 Hz or 30 Hz,
`gr_flag=0`, `A_nl=0`) and compared with the generator (SLy-like: $1.4+1.4\,M_\odot$, $R=11.7$ km;
$l=2$: $\omega_{a0}=2\pi\cdot1.8$ kHz, $k_2=0.08$; $l=3$: $\omega_{a0}=2\pi\cdot2.5$ kHz, $k_3=0.02$). Max relative error:

| mode | far from resonance ($\le800$ Hz) | worst case, and where |
|---|---|---|
| $l=2$, $m=0,-2$, spin 0 | $\le10^{-4}$ | $8\times10^{-3}$ ($m=0$, 1800 Hz), $10^{-3}$ ($m=-2$) |
| $l=2$, $m=+2$, spin 0 (resonance 1800 Hz) | $4.5\times10^{-4}$ at 800 Hz | 8% at 1500 Hz, just before resonance |
| $l=2$, $m=+2$, spin $-850$ Hz (resonance 950 Hz) | — | 5% at 800 Hz (pre-res.); **0.6–1.9% after resonance** |
| $l=3$, $m=-3,-1,+1$, spin 0 | $\le3\times10^{-4}$ | $2\times10^{-2}$ ($m=+1$, 1800 Hz) |
| $l=3$, $m=+3$, spin 0 (resonance 1667 Hz) | $3\times10^{-4}$ at 800 Hz | 4–6% at 1200–1500 Hz |
| $l=3$, $m=+3$, spin $-850$ Hz (resonance 533 Hz) | — | **$\le0.5\%$ after resonance, $4\times10^{-4}$ by 1800 Hz** |

The residual far from resonance grows as $(\omega_at_{\rm gw})^{-2}$ (e.g. $m=0$, 1800 Hz: $2/(\omega_at_{\rm gw})=0.077$,
squared $6\times10^{-3}$, as seen), which is the truncation of Eq. 18, not a bug. The post-resonance rows are the
test of the dynamical-tide term (amplitude and phase), where the dynamical component is comparable to the
equilibrium one. Further findings:

- **$l=2$ resonance frequency.** `get_ba_l2`'s closed form `wres1` is exact for `gr_flag=0` but differs
  from the exact root of $\Delta_a=0$ by $-2.1\%$ (1882 vs 1923 Hz, nonspinning, GR on), $-0.8\%$
  (spin $-850$), $+2.3\%$ (spin $-850$ with $A_{\rm nl}=-3\times10^{-2}$). Cause: 1PN truncation, and $A_{\rm nl}\,\omega_{a0}$ is
  not divided by $m$ when spin is present. This displaces where the dynamical tide switches
  on by 8–40 Hz, $\sim0.1$–$0.4$ of its width. `get_ba_l3` uses a fixed-point root of $\Delta_a=0$,
  agreeing with `brentq` to all printed digits; `get_ba_l2` keeps the closed form.

**Illustrative real-amplitude shape** (Paper 1 linear response $C_a=V_a\,\omega_a/(\omega_a-m_a\Omega)$,
$\omega_a=1.2\,\omega_0$ as quoted in Paper 1, $\Omega^2=M_t/r^3$, $X=1$, no lag; not the §8 configuration).
Semi-axis deviations $\Delta R/R$ at $x/y/z$:

| $r/R$ | adiabatic $x$ | dynamic $x$ | dynamic $y$ | dynamic $z$ |
|---|---|---|---|---|
| 3.0 | +0.056 | +0.067 | −0.039 | −0.028 |
| 2.75 | +0.073 | +0.093 | −0.057 | −0.037 |
| 2.5 | +0.097 | +0.138 | −0.089 | −0.049 |
| 2.3 | +0.125 | +0.204 | −0.141 | −0.063 |

$y\neq z$ is the triaxiality. At $r/R=2.3$ Paper 1 stops integration at tip contact
($\Delta R\approx0.15R$); linear dynamic $C_a$ give $0.20R$. The paper's nonlinear corrections and
exact contact definition may account for it. Informational, not pass/fail.

---

## 7. Implementation

```
src/tidal_shape.py     pure functions, numpy + scipy only
    ModeAmp(l, m, C: complex, I: signed, h)     dataclass
    StarModel                                   per-star modes + spin, X, X_lag; radius_fn(R_mm)
    delta_r_over_R(n_xyz, modes, X_lag)         general (l,m) via sph_harm_y, complex C, D_m of §3.3.1
    quadrupole_form(C_a, C_b, C_c, I, h2)       closed form (§3.3) -> (e_xx,e_yy,e_zz,e_xy), psi
    spin_oblateness(n_xyz, eps)                 §3.7, returns dR_spin/R
    surface_radius(n_xyz, R, modes, X, X_lag)   §5
    regime_report(...) / check_regime(...)      V11: Delta_a/omega_a, max|dR/R|, dV/V
    volume_drift(...)                           numerical, reported
    star_frame_to_world(vertices, companion_sign)
    load_star_modes(path)                       reads data/star_modes.toml
src/tidal_waveform.py  get_ba_l2 / get_ba_l3: total complex amplitudes along the inspiral (§1a)
src/polytrope.py       Gamma=2 polytrope: k_l, h_l, I_l, omega_a0/omega_0 (§4)
tools/make_star_modes.py    calls get_ba_l2 / get_ba_l3 at the chosen f_gw snapshot (§3.5) and writes the
                            .toml below, mapping [+2,0,-2] -> (C_a,C_c,C_b) and [+3,+1,-1,-3] with
                            signed I (§1a). Uses the §8 stellar model (polytrope k, I, h, omega_a0 scaled
                            per star) and the derived per-star A_nl; refuses a configuration that fails V11.
data/star_modes.toml        per-star blocks: M, R, M', f_gw snapshot (1023.6 Hz), spin_hz (negative =
                            anti-aligned), spin_oblateness, X, X_lag,
                            modes [(l, m, C re, C im, I, h, omega_a0)]
tests/test_tidal_shape.py   V1, V2, V4, V5, V6, V7, V8, V10, V11, V12
verification/verify_amplitudes.py   V9 regression at the sculpture stars
src/star.py                 build_star(radius_fn): mesh bulge axis follows psi
```

Amplitudes are data (the `.toml`), not constants in `params.py`.

---

## 8. Sculpture configuration and its checks

**Final configuration (D8–D10).** Stars (`design.md` §1.1): heavy $M=1.6\,M_\odot$, $R=12.0$ km, companion $1.2\,M_\odot$;
light $M=1.2\,M_\odot$, $R=12.5$ km, companion $1.6\,M_\odot$; $a=33$ km; snapshot $f_{\rm GW}=1023.6$ Hz $\Rightarrow$ $r/R=2.75$
(heavy), $2.64$ (light). Spin, anti-aligned, $\Omega_s=2\pi\times$`spin_hz`: **heavy $-500$ Hz, light $-200$ Hz**.
`gr_flag = 1`, $X=1$, $X_{\rm lag}=1$, `spin_oblateness = true`.

**Stellar model (D3, $\Gamma=2$ polytrope; `src/polytrope.py`).** Scaled to each star through
$\omega_0=\sqrt{GM/R^3}$: $l=2$: $k_2=0.25991$, $I_2=0.32153$, $h_2=1.51982$, $\omega_{a0}=1.227\,\omega_0$; $l=3$: $k_3=0.10645$,
$I_3=0.24341$, $h_3=1.21291$, $\omega_{a0}=1.698\,\omega_0$. ($1.227\,\omega_0$ reproduces Paper 1's $\omega_a\simeq1.95$ kHz for
$1.3\,M_\odot$, 12 km; the paper's "$1.2\,\omega_0$" is a rounding.) $\Lambda_2=k_2/(1.5\,C^5)$ with $C=GM/(Rc^2)$. Resulting
$\omega_{a0}/2\pi$: heavy 2.165 kHz ($l=2$) / 2.996 kHz ($l=3$); light 1.763 kHz / 2.440 kHz.

**Nonlinear frequency shift, derived rather than the default.** The $A_{\rm nl}=-3\times10^{-2}$ default of
`get_ba_l2` is the equal-mass $1.3+1.3\,M_\odot$ polytrope: Paper 1 Eq. 16 gives

$$
\frac{\Delta\omega_a}{\omega_a}=-0.623\,\frac{M'}{M}\,\frac{\Omega^2R^3}{GM_t},\qquad
\sqrt{\pi/5}\,(J_2+4\kappa_2I)=-0.623\ \ (\text{Table 1}),
$$

i.e.

$$
A_{\rm nl}=-0.623\,\frac{M'}{M}\,\frac{\omega_{\rm ref}^2R^3}{GM_t},\qquad\omega_{\rm ref}=\pi\cdot1000\ \mathrm{rad/s};
$$

this evaluates to $-3.04\times10^{-2}$ for that system.
For the sculpture stars: **heavy $-0.0214$, light $-0.0431$**. $A_{\rm nl}$ for $l=3$ is 0.

**Results at the snapshot** ($l=2+l=3$ tidal surface from the generator; $\Delta_a/\omega^{\rm corot}$ for the $(2,2)$ and $(3,3)$ modes;
$\Delta R/R$ range excludes the oblateness term; line-of-centres gap uses the tips toward the companion, a proxy for the
3-D gap of V7):

| star | spin | $\Delta_2/\omega$ | $\Delta_3/\omega$ | $\psi$ ($l=2$) | $l=2$ $x$ tip | $\Delta R/R$ range | $\Delta V/V$ | $\varepsilon_\Omega$ |
|---|---|---|---|---|---|---|---|---|
| heavy | 0 | 0.53 | 0.49 | 0.75° | +0.068 | [−0.042, +0.087] | +0.3% | 0 |
| heavy | −500 Hz | 0.24 | **0.13** | 2.2° | +0.097 | [−0.080, +0.135] | +0.7% | 0.041 |
| light | 0 | 0.42 | 0.37 | 1.7° | +0.161 | [−0.110, +0.207] | +1.7% | 0 |
| light | −500 Hz | **0.11** | **−0.03** | 12.8° | +0.314 | **[−0.450, +0.579]** | **+11.9%** | 0.061 |

The light star at $-500$ Hz triggers V11: its $(3,3)$ resonance is at $f_{\rm GW}\approx961$ Hz (Newtonian; 997 Hz with GR + the nonlinear shift), 6% (2.6%) below the snapshot, and the
dynamical tide dominates its $l=3$ surface ($\pm0.24$). The linear approximation is not valid there.

**Final configuration, computed** (the table above and the scan below show how it was chosen). Surface $=l=2+l=3$ tide
($X=X_{\rm lag}=1$) plus spin oblateness; $\Delta V/V$ and ranges from a $361\times720$ angular grid:

| | heavy ($-500$ Hz) | light ($-200$ Hz) |
|---|---|---|
| $\Delta_2/\omega$, $\Delta_3/\omega$ (V11 threshold 0.15) | 0.24, **0.13 (warn)** | 0.28, 0.19 |
| $(2,2)$ / $(3,3)$ resonance, $f_{\rm GW}$, **Newtonian** (used by V11) | 1665 Hz / 1331 Hz | 1563 Hz / 1361 Hz |
| same, exact root of $\Delta_a=0$ with GR + nonlinear shift | 1616 Hz / 1397 Hz | 1451 Hz / 1431 Hz |
| $\Delta_2/\omega$, $\Delta_3/\omega$ with GR + nonlinear shift | 0.218, 0.137 | 0.223, 0.191 |
| lag $\psi$ ($l=2$, total amplitude) | 2.15° | 3.37° |
| $\varepsilon_\Omega$ (oblateness, polar) | 0.0407 | 0.0098 |
| $\Delta R/R$ range, tide only | [−0.080, +0.135] | [−0.159, +0.267] |
| $\Delta R/R$ range, tide + oblateness | [−0.070, +0.155] | [−0.154, +0.272] |
| $\min r/\max r$ | 0.805 | 0.665 |
| $\Delta V/V$, tide only $\to$ with oblateness | +0.67% → +0.92% | +2.9% → +3.0% |
| tip toward companion, $\Delta R/R$ (tide + oblateness) | +0.154 | +0.271 |

Line-of-centres gap, $a=33$ km: **3.57 km (2.58 mm) tide only, 3.26 km (2.36 mm) with oblateness**, against 4.9 km
for the nonspinning case. V11: heavy passes with one warning ($\Delta_3/\omega$), light passes. The only 3-D check still
missing is V7 on the real meshes.

The light star remains the more deformed one (tip $+0.27$, $\min r/\max r=0.665$), a 27% bulge at spin $-200$ Hz, even
though its spin is the smaller: its mass ratio $M'/M=1.33$ and larger $R/r$ give it $\sim2\times$ the heavy star's tidal
response (spin 0: $+0.16$ vs $+0.07$). At these amplitudes the linear surface is an approximation (second-order
terms are $\sim(0.27)^2$), which is the accuracy chosen in §0; V11 warns beyond $\max\lvert\Delta R/R\rvert>0.3$.

**Light-star spin scan** (heavy fixed at $-500$ Hz, line-of-centres gap with both stars' tidal tips, $a=33$ km):

| light spin | $\Delta_2/\omega$ | $\Delta_3/\omega$ | $\psi$ | $\Delta R/R$ range | $\Delta V/V$ | gap | V11 |
|---|---|---|---|---|---|---|---|
| 0 | 0.42 | 0.37 | 1.7° | [−0.110, +0.207] | +1.7% | 4.3 km | pass |
| −100 | 0.34 | 0.28 | 2.2° | [−0.130, +0.231] | +2.2% | 4.0 km | pass |
| **−200** | 0.28 | 0.19 | 3.4° | [−0.159, +0.267] | +2.9% | 3.6 km | pass |
| −300 | 0.21 | 0.11 | 5.2° | [−0.200, +0.319] | +4.0% | 3.0 km | warn ($l=3$ margin; range $>0.3$) |
| −400 | 0.16 | 0.04 | 7.8° | [−0.268, +0.390] | +6.0% | 2.3 km | warn |
| −500 | 0.11 | −0.03 | 12.8° | [−0.450, +0.579] | +11.9% | 1.3 km | **error** |

If both stars share one spin: $-300$ Hz gives lags $1.3^\circ/5.2^\circ$ (heavy / light), gap 3.3 km; $-400$ Hz $1.6^\circ/7.8^\circ$,
gap 2.5 km. The light star's $\sim2\times$ larger tidal response and lower mode frequencies make it the limiting star.

On the sculpture scale ($S=0.72427$ mm/km) a $1^\circ$ lag moves the ends of a Ø$\sim$18 mm star's long axis by
$R\sin1^\circ\approx0.15$ mm; $\psi=3.4^\circ$ is $\sim0.5$ mm, $12.8^\circ$ is $\sim2$ mm. $X_{\rm lag}$ (§3.3.1) is the knob if this is too small.

---

## 9. Results

Usage: `python tools/make_star_modes.py` (writes `data/star_modes.toml`; refuses configurations that fail V11),
`python -m pytest tests`, `python build_all.py` (STLs into `out/`).

- **Tests:** 48 pass, 37 for the V-list and 11 in `tests/test_math_walkthrough.py` (V1 for $l=2$ and $l=3$ with both
  overlap signs to $10^{-12}$; V2 to $10^{-12}$; V5 to $10^{-13}$; the generator exits 1 and writes nothing for the light
  star at $-500$ Hz).
- **Generated amplitudes** reproduce the §8 table to the printed precision: heavy $\Delta_2/\omega=0.241$,
  $\Delta_3/\omega=0.132$ (the one warning), $\mathrm{d}R/R\in[-0.080,+0.135]$, $\Delta V/V=+0.67\%$;
  light $0.275$, $0.191$, $[-0.159,+0.267]$, $+2.9\%$. Lags per harmonic: heavy $(2,2)$ $2.15^\circ$, $(3,1)$ $0.51^\circ$,
  $(3,3)$ $2.93^\circ$; light $3.37^\circ$, $0.63^\circ$, $2.52^\circ$.
- **Overlaps used.** The generator uses $I_2=\sqrt{5k_2/4\pi}=0.32158$ and $I_3=\sqrt{7k_3/4\pi}=0.24352$ from the exact
  static $k_l$, since that is what `get_ba_l2` / `get_ba_l3` compute internally; the f-mode-only values in §4 / §8
  ($0.32153$, $0.24341$) differ by $0.02\%$ and $0.05\%$. The shape divides by the same $I$, so only $C/I$ matters.
- **Meshes** (subdivision 6): heavy 81 478 faces, light 81 526 faces, both watertight and winding-consistent.
  Bare-mesh volume agrees with `volume_drift` to $2\times10^{-3}$ (mesh discretization).
- **Gap (V7):** line of centres $2.363$ mm; the minimum 3-D vertex distance is $2.363$ mm, so the closest approach
  *is* the line of centres (the $2.58$ mm tide-only figure of §8 is also reproduced). The 3-D gap is
  $\le$ the line-of-centres gap by construction; here they agree to $4\times10^{-7}$ mm.
- **V9 at the sculpture stars** (point-mass limit, `gr_flag=0`, $A_{\rm nl}=0$, relative error at the snapshot against the
  Eq. 15 integration, `verification/verify_amplitudes.py`): $l=2$: $m=0$ $2\text{–}3\times10^{-4}$, $m=-2$ $5\text{–}7\times10^{-5}$,
  $m=+2$ $1.0\%$; $l=3$: $m=\pm1$ $\le5\times10^{-4}$, $m=-3$ $3\text{–}5\times10^{-5}$, $m=+3$ $3.7\%$ (light) / $4.2\%$ (heavy).
  The prograde modes nearest resonance carry the Eq. 18 truncation error ($\Delta/\omega=0.13\text{–}0.28$); at $\approx4\%$
  of an $l=3$ pattern of amplitude $\lesssim0.05$ this is $\lesssim0.002\,R$ ($\approx0.02$ mm).
- **Pipeline:** the full build passes `validate.py` with one warning (the heavy star's $(3,3)$ margin).

**Limits.** The socket is $\varnothing3.375\times6$ mm, measured from the pole point $r(0,0,-1)$ and entirely in the
lower half of the star (`design.md` §6.4); the printed star is the modelled surface plus the bore. The $l=2$ part of
the $(2,2)$ and $l=3$ part of the $(3,3)$ amplitudes of the heavy star are the least accurate inputs ($1\%$ and $4\%$,
above); `get_ba_l2`'s closed-form resonance frequency is the 1PN one (§6).
