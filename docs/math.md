# From Tidal Mode Amplitudes to the Deformed Neutron-Star Surface — Math Specification

Companion to [`neutron_star_shape.md`](neutron_star_shape.md). That document specifies *what* is computed and how it is
validated; this one derives *why the formulas are what they are* and ties every equation to the code.
Every number in the worked example (§6) is reproduced by `verification/math_walkthrough.py` (an independent float64
re-implementation of §1–§4, numpy + scipy only) and checked by `tests/test_math_walkthrough.py`; the polytrope numbers
(§5) come from `src/polytrope.py`, and the other checks are listed in §8.

Sources: Yu, Weinberg, Arras, Kwon & Venumadhav, arXiv:2211.07002 ("Paper 1"); Yu, Arras & Weinberg,
arXiv:2404.00147 ("Paper 2"). Equation numbers refer to those papers. Where a formula is adopted from them
rather than derived here, the text says so.

---

## 0. The pipeline in one page

```
 orbit (M, M', f_GW)  ──►  drive V_a  ──►  amplitude ODE  ──►  b_a = C_a  ──►  C_a / I_a  ──►  ΔR(θ,φ)  ──►  r(n)  ──►  mesh
   §1.1                     §1.2, §2.2       §3.1–3.5            §3.6           §4.3–4.6         §4.5–4.6      §4.8
   (ω, ω̇, ṙ/r)             (W_lm I_a)       (Δ_a, κ_a, F(u))    (complex)      (calibration)    (D_m, ψ)      (+ spin §4.7)
```

The whole construction compresses into two statements.

1. **Amplitudes.** Each mode amplitude is the *adiabatic drive* $V_a$ multiplied by a complex *transfer function*
   $T_a=b_a/V_a$ that carries the finite-frequency enhancement, the tidal lag, spin, and the dynamical tide (§3).
2. **Surface.** The surface is the *equilibrium tide* (shape number $h_l$) with each $(l,m)$ component multiplied by
   its $T_a$ (§4.4):
$$
\boxed{\ \frac{\Delta R}{R}(\theta,\phi_s)=\sum_{l}h_l\,\frac{M'}{M}\Big(\frac{R}{r}\Big)^{l+1}
\mathrm{Re}\sum_{m}T_{lm}\,W_{lm}\,Y_{lm}(\theta,\phi_s)\ }
$$
   The mode overlap $I_a$ cancels. The code carries it anyway because it works with $\hat C_a=C_a/I_a$.

### 0.1 Notation

| Symbol | Meaning |
|---|---|
| $M,R$; $M'$ | deformed star's mass and radius; companion mass |
| $M_t=M+M'$, $\eta=MM'/M_t^2$, $X_1=M/M_t$, $X_2=M'/M_t$ | total mass, symmetric mass ratio, mass fractions |
| $\omega=\pi f_{\rm GW}$, $\phi_{\rm orb}$ | orbital angular frequency and phase |
| $\Omega_s$ | stellar spin rate, **negative = anti-aligned** with the orbit |
| $a=(l_a,m_a,n_a,\pm)$ | a mode; $\bar a$ its conjugate $(−\omega_a,−m_a)$ |
| $\omega_{a0}$, $\omega_a=\omega_{a0}-\frac{m_a}{l_a}\Omega_s$ | non-rotating and corotating-frame mode frequency |
| $I_a$ | tidal overlap (§2.2); $\hat C_a=C_a/I_a$ |
| $V_a$, $b_a\equiv C_a$ | adiabatic drive and slowly varying mode amplitude (§3) |
| $\Delta_a$, $\kappa_a$ | detuning and decay term of the amplitude equation (§3.2) |
| $T_a=C_a/V_a$ | complex transfer function (§3.3) |
| $h_l=1+2k_l$ | equilibrium shape number of the $l$-th multipole (§2.3) |
| $\phi_s=\phi-\phi_{\rm orb}$ | azimuth in the orbit-corotating frame (companion at $\phi_s=0$) |
| $Y_{lm}$, $W_{lm}$, $N_2,N_0$ | orthonormal harmonics (Condon–Shortley), $W_{lm}=\frac{4\pi}{2l+1}Y_{lm}(\tfrac\pi2,0)$, $N_2=\frac14\sqrt{15/2\pi}$, $N_0=\frac14\sqrt{5/\pi}$ |

Units: the generator works in cgs with $G$, $c$ explicit; the papers use $G=c=1$. All quantities that enter the surface
($C_a$, $I_a$, $T_a$, $h_l$, $\Delta R/R$) are dimensionless.

### 0.2 Equation ↔ code map

| Eq. | Content | Code |
|---|---|---|
| (1.1)–(1.4) | orbit, $\dot\omega$, $\dot r/r$ | `tidal_waveform._inspiral`; `math_walkthrough.orbit` |
| (1.6) | $W_{lm}$ | `get_ba_l2` (`Wlm1`, hard-coded $l=2$); `tidal_waveform._w_lm_l3`; `math_walkthrough.W` |
| (2.5)–(2.6) | $V_a=\frac{M'}{M}W_{lm}I_a(R/r)^{l+1}$ | `get_ba_l2`, `get_ba_l3` (`va1`) |
| (2.7)–(2.9) | $k_l=\frac{4\pi}{2l+1}I^2$, $h_l=1+2k_l$ | `get_ba_l2` (`k2`, $k_2\to I_2$); `get_overlaps_l3`; `polytrope.static_love` |
| (3.9), (3.10) | $\Delta_a$ with GR/NL, $z$, $\omega_{\rm FD}$ | `tidal_waveform._gr_factors`; `Del` in `get_ba_l2` / `get_ba_l3`; `math_walkthrough.delta_a`, `gr_factors` |
| (3.6), (3.9) | Eq. 18, $b_{\rm eq}$ | `ba_eq1` in `get_ba_l2` / `get_ba_l3`; `math_walkthrough.equilibrium` |
| (3.11) | $A_{\rm nl}$ | `tools/make_star_modes.py` `PAPER1_NL` (L31), `A_nl2` (L58) |
| (3.12)–(3.14) | dynamical tide | `tidal_waveform._dynamical_tide`; resonance in `get_ba_l2` / `get_ba_l3`; `math_walkthrough.dynamical`, `resonance` |
| (4.9)–(4.11) | $\hat C$, $D_m$, $s_{lm}$, $\psi_m$, $X_{\rm lag}$ | `tidal_shape.combined_amplitudes` L157–189, `lag_angles` L191 |
| (4.10) | $\Delta R/R$ from $D_m$ | `tidal_shape.delta_r_over_R` L204 |
| (4.6)–(4.8) | closed-form $l=2$ tensor, $\psi$ | `tidal_shape.quadrupole_form` L236 |
| (4.13) | spin oblateness | `tidal_shape.eps_omega` L214, `spin_oblateness` L220 |
| (4.14) | $r(\mathbf n)$, frames, $\Delta V/V$ | `tidal_shape.surface_radius` L226, `star_frame_to_world` L121, `StarModel.radius_fn` L100, `volume_drift` L260 |
| (5.1)–(5.4) | polytrope $I,h,\omega/\omega_0$ | `src/polytrope.py`; used by `tools/make_star_modes.py` `polytrope_numbers` (L40), `star_block` (L49) |

---

## 1. The orbit and the tidal drive

### 1.1 Quasi-circular inspiral

The binary is Keplerian at each instant with point-mass quadrupole radiation reaction. With
$t_M\equiv GM_t/c^3$, $v\equiv(t_M\omega)^{1/3}$ and $x\equiv v^2$:

$$
r=\Big(\frac{GM_t}{\omega^2}\Big)^{1/3},\qquad
\dot\omega=\frac{96}{5}\,\eta\,\frac{v^{11}}{t_M^{2}},\qquad
t_{\rm gw}\equiv\frac{\omega}{\dot\omega},\qquad
\frac{\dot r}{r}=-\frac23\,\frac{\dot\omega}{\omega}=-\frac{2}{3\,t_{\rm gw}}
\tag{1.1–1.4}
$$

The last relation follows from $\omega^2r^3=GM_t$ (Paper 2 Eqs. 66–67). Time and phase from a reference frequency
$\omega_r$ follow by integrating $\dot\omega$ (used in §3.5):

$$
\tau(\omega)=\frac{5}{256\,\eta}\,t_M\left(v_r^{-8}-v^{-8}\right),\qquad
\phi(\omega)-\phi_r=\frac{1}{32\,\eta}\left(v_r^{-5}-v^{-5}\right)
\tag{1.5}
$$

(Paper 2 Eqs. 26, 31).
*Heavy star at the snapshot* ($f_{\rm GW}=1023.58$ Hz): $\omega=3215.66\ \mathrm{s^{-1}}$, $r=33.00$ km ($r/R=2.750$),
$x=0.12529$, $\dot\omega=2.7018\times10^5\ \mathrm{s^{-2}}$, $t_{\rm gw}=11.90$ ms, $\dot r/r=-56.01\ \mathrm{s^{-1}}$.

### 1.2 Tidal potential of the companion

The companion sits at radius $r$, colatitude $\pi/2$ and azimuth $\phi_{\rm orb}$. Its potential at a point $\mathbf x$ in the
star expands in Legendre polynomials; the $l=0,1$ terms do not deform the star. With the addition theorem
$P_l(\cos\gamma)=\frac{4\pi}{2l+1}\sum_mY_{lm}(\hat{\mathbf x})Y^*_{lm}(\hat{\mathbf n}_c)$ and
$Y^*_{lm}(\tfrac\pi2,\phi_{\rm orb})=Y_{lm}(\tfrac\pi2,0)\,e^{-im\phi_{\rm orb}}$:

$$
U(\mathbf x)=-GM'\sum_{l\ge2}\frac{|\mathbf x|^{l}}{r^{l+1}}P_l(\cos\gamma)
=-GM'\sum_{l\ge2}\sum_{m}W_{lm}\,\frac{|\mathbf x|^{l}}{r^{l+1}}\,Y_{lm}(\hat{\mathbf x})\,e^{-im\phi_{\rm orb}}
\tag{1.6}
$$

$Y_{lm}(\pi/2,0)=0$ when $l+m$ is odd, so only **$l+m$ even** modes are driven. Values:

| $(l,m)$ | $W_{lm}$ | | $(l,m)$ | $W_{lm}$ |
|---|---|---|---|---|
| $(2,\pm2)$ | $+\sqrt{3\pi/10}=+0.97081$ | | $(3,+3)$, $(3,-3)$ | $-0.74900$, $+0.74900$ |
| $(2,0)$ | $-\sqrt{\pi/5}=-0.79267$ | | $(3,+1)$, $(3,-1)$ | $+0.58017$, $-0.58017$ |

For odd $l$, $W_{l,-m}=-W_{lm}$ ($m$ odd). Two products are used throughout:
$W_{22}N_2=\tfrac38$ and $W_{20}N_0=-\tfrac14$.

---

## 2. Normal modes of the star

### 2.1 Mode decomposition, normalization and reality

The perturbation is expanded in the star's eigenmodes (phase-space decomposition, Paper 1 Eq. 1):

$$
\begin{pmatrix}\boldsymbol\xi\\ \dot{\boldsymbol\xi}\end{pmatrix}=\sum_a c_a(t)\begin{pmatrix}\boldsymbol\xi_a\\ -i\omega_a\boldsymbol\xi_a\end{pmatrix},\qquad
\boldsymbol\xi_a=\xi_{r,a}(r)\,Y_{lm}\,\hat{\mathbf r}+\xi_{h,a}(r)\,r\nabla Y_{lm}
\tag{2.1}
$$

The sum runs over both signs of $\omega_a$. The modes are normalized by (Paper 1 Eq. 3, Paper 2 Eq. 8)

$$
2\omega_a^{2}\int\rho\,\boldsymbol\xi_a^{*}\!\cdot\boldsymbol\xi_b\,d^3x=E_0\,\delta_{ab},\qquad E_0=\frac{GM^2}{R}
\tag{2.2}
$$

and the real field requires $c_{\bar a}=c_a^*$, $\boldsymbol\xi_{\bar a}=\boldsymbol\xi_a^*$, so over positive-frequency modes

$$
\boldsymbol\xi=\sum_{\omega_a>0}2\,\mathrm{Re}\left[c_a\,\boldsymbol\xi_a\right]
\tag{2.3}
$$

The factor 2 is exactly what makes (2.2) the mode energy: with $\boldsymbol\xi=c\,\boldsymbol\xi_a+\mathrm{c.c.}$ the time-averaged kinetic
energy is $\omega^2|c|^2\int\rho|\boldsymbol\xi_a|^2$, and equipartition gives $E=2\omega^2|c|^2\!\int\rho|\boldsymbol\xi_a|^2=E_0|c|^2$.
(The $\sum2\,\xi_rI/R$ sum rule of §2.3 confirms the factor numerically.)

### 2.2 Overlap integral and the drive

The linear coupling of the tide to a mode is its overlap (Paper 1 Eq. 5, Paper 2 Eq. 10):

$$
I_a=\frac{1}{MR^{l}}\int\rho\,\boldsymbol\xi_a^{*}\!\cdot\nabla\!\left(r^{l}Y_{lm}\right)d^3x
\tag{2.4}
$$

Projecting $-\nabla U$ of (1.6) on mode $a$: because $\rho$ is spherical, only the harmonic $(l_a,m_a)$ of (1.6) survives, so

$$
\langle\boldsymbol\xi_a,-\nabla U\rangle=\int\rho\,\boldsymbol\xi_a^*\!\cdot(-\nabla U)\,d^3x
=GM'\,W_{lm}\frac{e^{-im\phi_{\rm orb}}}{r^{l+1}}\;MR^{l}I_a
$$

and dividing by $E_0$ gives the dimensionless drive (Paper 2 Eq. 9)

$$
U_a=\frac{\langle\boldsymbol\xi_a,-\nabla U\rangle}{E_0}=V_a\,e^{-im\phi_{\rm orb}},\qquad
V_a=\frac{M'}{M}\,W_{lm}\,I_a\Big(\frac{R}{r}\Big)^{l+1}
\tag{2.5–2.6}
$$

$V_a$ is real. In the code $V$ is called `va1`; for $l=2$ all three $I_a$ equal $\sqrt{5k_2/4\pi}$, for $l=3$ they are
signed (§3.6).

### 2.3 Quasi-static response, Love numbers and the shape number

The mass multipole induced by the displacement is $Q_{lm}=\int\rho\,\boldsymbol\xi\cdot\nabla(r^lY^*_{lm})\,d^3x=MR^l\sum_aI_ac_a$
(Paper 1 Eq. C1). In the adiabatic limit $c_a\to V_a$, and Paper 2 Eqs. 42–43 give the paper-convention Love number

$$
k_l=\frac{4\pi}{2l+1}\,I_l^{2}
\tag{2.7}
$$

(the factor 2 of $\mathrm{Re}$ in (2.3) is absorbed here). Its relation to the classical Love number follows from the
definitions. With $\Phi_{\rm tid}=\frac1{l!}E_Lx^L$ and $Q_L=-\frac{2}{(2l-1)!!}k_lR^{2l+1}E_L$ (Paper 2 Eq. 41), the induced potential
$\Phi_{\rm ind}=-\frac{(2l-1)!!}{l!}Q_L\,n^L/r^{l+1}$ evaluated at the surface is

$$
\Phi_{\rm ind}(R)=2k_l\,\Phi_{\rm tid}(R)\quad\Longrightarrow\quad k_l^{\rm Love}=2k_l
\tag{2.8}
$$

A fluid star's surface is an equipotential of the total potential. With $g=GM/R^2$ and the static l-th component of (1.6),
$g\,\xi_r+\Phi_{\rm tid}(R)+\Phi_{\rm ind}(R)=0$, so the **equilibrium tide** is

$$
\frac{\xi_r}{R}\Big|_{\rm eq}=h_l\,\frac{M'}{M}\Big(\frac{R}{r}\Big)^{l+1}P_l(\cos\gamma),\qquad
h_l=1+k_l^{\rm Love}=1+2k_l
\tag{2.9}
$$

For the polytrope $h_2=1.51982$, $h_3=1.21291$ (§5.3). The classical limit check: a uniform-density star has $k^{\rm Love}=3/2$ and
$h_2=5/2$, the textbook values.

---

## 3. Dynamics of the mode amplitudes

### 3.1 Amplitude equation in the orbit frame

Projecting the equation of motion $\ddot{\boldsymbol\xi}+\mathcal C\boldsymbol\xi=-\nabla U$ on a mode gives (Paper 2 Eq. 11)

$$
\dot q_a+i\omega_aq_a=i\omega_{a0}\,U_a
\tag{3.1}
$$

For a spinning star the mode frequency in the star's corotating frame is shifted, $\omega_a=\omega_{a0}-\frac{m_a}{l_a}\Omega_s$
(Paper 2 Eq. 7), and the azimuth that multiplies $Y_{lm}$ is $\phi_c=\phi_{\rm inertial}-\Omega_st$.
Define the slowly varying amplitude $b_a=q_a\,e^{im(\phi_{\rm orb}-\Omega_st)}$. Since $U_a=V_ae^{-im(\phi_{\rm orb}-\Omega_st)}$ and
$\dot\phi_{\rm orb}=\omega$,

$$
\dot b_a+i\Delta_a\,b_a=i\omega_{a0}V_a,\qquad
\Delta_a=\omega_a+m_a\Omega_s-m_a\omega
\tag{3.2}
$$

(Paper 2 Eq. 15). $\Delta_a$ is the detuning between the mode and the tidal forcing seen in the star's frame: it vanishes at
**resonance**, $m\omega=\omega_a+m\Omega_s$. Two consequences of the phase factor:

$$
q_a\,Y_{lm}(\theta,\phi_c)=b_a\,Y_{lm}(\theta,\phi_{\rm inertial}-\phi_{\rm orb})=b_a\,Y_{lm}(\theta,\phi_s)
\tag{3.3}
$$

The surface pattern in the orbit-corotating frame depends on the spin only through $b_a$, so **spin needs no change to the
shape code**; and with $b_a\equiv C_a$ the amplitude is static there whenever $b_a$ is constant.

### 3.2 The equilibrium solution (Paper 2 Eq. 18)

Away from resonance $b_a$ varies slowly, so (3.2) is solved by iteration. Setting $\dot b=0$:

$$
b^{(0)}=\frac{\omega_{a0}}{\Delta_a}V_a
\tag{3.4}
$$

The next order is $b^{(1)}=i\dot b^{(0)}/\Delta_a$. Along the inspiral $V_a\propto r^{-(l+1)}$ and $\dot\Delta_a=-m\dot\omega$, so

$$
\frac{\dot b^{(0)}}{b^{(0)}}=\frac{\dot V}{V}-\frac{\dot\Delta}{\Delta}=\frac{m\dot\omega}{\Delta}-(l+1)\frac{\dot r}{r}
\;\;\Longrightarrow\;\;
\frac{b^{(1)}}{b^{(0)}}=\frac{i\,\kappa_a}{\Delta_a^{2}},\qquad
\kappa_a\equiv m_a\dot\omega-\Delta_a\,(l_a+1)\frac{\dot r}{r}
\tag{3.5}
$$

Rather than truncating, the geometric series is resummed, $b_{\rm eq}=b^{(0)}/\big(1-b^{(1)}/b^{(0)}\big)$:

$$
\boxed{\ b_a^{\rm(eq)}=\frac{\Delta_a\,\omega_{a0}\,V_a}{\Delta_a^{2}-i\,\kappa_a}\ }
\tag{3.6}
$$

(Paper 2 Eq. 18). It is finite at $\Delta_a=0$, where the orbital decay acts as an effective damping. Using (1.4),

$$
\kappa_a=\dot\omega\left[m_a+\tfrac23(l_a+1)\frac{\Delta_a}{\omega}\right]
\tag{3.7}
$$

which is positive for $\Delta_a>0$ and $m_a\ge0$ during an inspiral ($\dot\omega>0$), and also for $m_a=-2$ at the values here.

**Tidal lag.** The imaginary part of $b_a$ is the lag. To lowest order $\arg b^{\rm(eq)}_a=\arctan(\kappa_a/\Delta_a^2)$ (for
$\Delta_a>0$), and the harmonic's equatorial pattern is $|b|\cos(m\phi_s+\arg b)$, which peaks at $\phi_s=-\arg b/m$: **behind the
companion**. Without orbital decay $\kappa=0$ and (3.6) reduces to Paper 1's real $\omega_{a0}V_a/\Delta_a$.

### 3.3 The transfer function

Define $T_a\equiv C_a/V_a$ (with $C_a=b_a$). The equilibrium transfer function is

$$
T_a^{\rm(eq)}=\frac{z\,\omega_{a0}\,\Delta_a}{\Delta_a^{2}-i\kappa_a}
\tag{3.8}
$$

where $z$ is the GR redshift factor of §3.4 ($z=1$ in the Newtonian paper formulas). $T_a\to1$ in the adiabatic limit
($\omega\to0$, no spin). For $\Delta_a>0$, $|T_a|>1$ when the mode frequency is above the forcing (finite-frequency
enhancement $\omega_{a0}/\Delta_a$) and the response is ordered $|T_{+m}|>|T_0|>|T_{-m}|$ because $\Delta_{+m}<\Delta_0<\Delta_{-m}$.
This is the object that Paper 2 calls the frequency-dependent effective Love number (Eq. 42).

### 3.4 GR and nonlinear corrections, as implemented

The code multiplies and shifts (3.6) with three phenomenological corrections:

$$
\Delta_a=z\,(\omega_a+m_a\Omega_s)-m_a\omega+m_a\,\omega_{\rm FD}+A_{\rm nl}\Big(\frac{\omega}{\omega_{\rm ref}}\Big)^{2}\omega_{a0},\qquad
b_a^{\rm(eq)}=\frac{z\,\Delta_a\,\omega_{a0}V_a}{\Delta_a^{2}-i\kappa_a}
\tag{3.9}
$$

with, at second post-Newtonian order in $x$ (the expansions of "S+21", as named in the code comments),

$$
z=1+\tfrac{x}{2}X_2(X_1-3)+\tfrac{x^2}{24}X_2\!\left(5X_1^2X_2-9X_1-6X_1X_2-27\right),\quad
\frac{\omega_{\rm FD}}{\omega}=\tfrac{x}{2}X_2(X_1+3)-\tfrac{x^2}{24}X_2\!\left(X_1^2X_2-45X_1+30X_1X_2-27\right)
\tag{3.10}
$$

and $\omega_{\rm ref}=\pi\cdot10^3\ \mathrm{rad/s}$. These are **inputs adopted from the user's code, not derived here**. Their effect is a
modest detuning shift: for the heavy star's $(2,2)$ mode the Newtonian margin $\Delta/\omega_a=0.241$ becomes $0.218$ with
$z$, $\omega_{\rm FD}$ and $A_{\rm nl}$ (§6.1).

*Nonlinear shift.* For the $l=2$ modes the three-mode and four-mode couplings of Paper 1 lower the mode frequency by
$\Delta\omega_a/\omega_a=\sqrt{\pi/5}\,(J_2+4\kappa_2I)\,\frac{M'}{M}\frac{\Omega^2R^3}{GM_t}$ (Paper 1 Eq. 16). With the paper's Table 1
($J_2=-0.21$, $\kappa_2=-0.45$, $I=0.32$) the prefactor is $-0.623$. Since $\Omega^2R^3/GM_t=(R/r)^3\propto\omega^2$ this is of the form
$A_{\rm nl}(\omega/\omega_{\rm ref})^2$ with

$$
A_{\rm nl}=-0.623\,\frac{M'}{M}\,\frac{\omega_{\rm ref}^{2}R^{3}}{GM_t}
\tag{3.11}
$$

For equal masses and the paper's $1.3+1.3\,M_\odot$, $R=12$ km this evaluates to $-3.04\times10^{-2}$, which is where the generator's
default $-3\times10^{-2}$ comes from. For the sculpture stars: heavy $-0.02145$, light $-0.0431$. $A_{\rm nl}=0$ for $l=3$ (no
established value).

### 3.5 The dynamical tide

Write $b=b^{\rm(eq)}+b^{\rm(dyn)}$. Subtracting the equilibrium solution from (3.2) leaves a driven equation for the part that oscillates
at the mode frequency, with source $i\omega_{a0}V^{\rm dyn}=i\omega_{a0}V-\dot b_{\rm eq}-i\Delta_ab_{\rm eq}$ (Paper 2 Eq. 20). That source is
sharply peaked at resonance, so the response is an accumulated Fresnel-type integral. In terms of the time since resonance
$\tau$ and the dimensionless variable $u$ (Paper 2 Eqs. 24–26, using (1.5) at $\omega_r$ = resonance):

$$
u=\sqrt{\tfrac{m}{2}\dot\omega_r}\;\tau,\qquad
F(u)=\int_{-\infty}^{u}\frac{2i\,(1+a_1u')}{\left(1+a_2u'+2iu'^{2}\right)^{2}}\,e^{-iu'^{2}}\,du'
\tag{3.12}
$$

$$
a_1=\frac{2l-9}{6}\,\varepsilon,\qquad a_2=-\frac{2(l+1)}{3}\,\varepsilon,\qquad
\varepsilon=\sqrt{\frac{2}{m}\,\frac{\dot\omega_r}{\omega_r^{2}}}=\sqrt{\frac{2}{m\,\omega_rt_{\rm gw,r}}}
\tag{3.13}
$$

($a_1,a_2\ll1$ because $\omega_rt_{\rm gw,r}\gg1$). The amplitude and phase are (Paper 2 Eqs. 23, 30, with the resonance condition
$\omega_a+m\Omega_s\to m\omega_r$):

$$
b_a^{\rm(dyn)}=\frac{z_r\,\omega_{a0}\,V_{a,r}}{\sqrt{\tfrac m2\dot\omega_r}}\;F(u)\;
\exp\!\Big\{i\,m\big[\phi(\omega)-\phi_r-\omega_r\tau(\omega)\big]\Big\},\qquad V_{a,r}=V_a\big|_{\omega=\omega_r}
\tag{3.14}
$$

The resonance frequency $\omega_r$ is the root of $\Delta_a(\omega)=0$. `get_ba_l3` finds it by fixed-point iteration,
$\omega\leftarrow\big[z(\omega)(\omega_a+m\Omega_s)+A_{\rm nl}(\omega/\omega_{\rm ref})^2\omega_{a0}\big]/\big(m[1-\omega_{\rm FD}/\omega]\big)$;
`get_ba_l2` uses a 1PN closed form $\omega_r=\big(1+\eta x_r+A_{\rm nl}(\omega_N/\omega_{\rm ref})^2\big)\omega_N$ with
$\omega_N=\omega_a/2+\Omega_s$, which differs from the exact root by 0.5–2.3% with GR on (heavy $(2,2)$: 1636 vs 1616 Hz;
this shifts where the dynamical tide switches on, not its amplitude). Two implementation facts matter when reading the code:

- the code accumulates $F$ by a trapezoid rule from the first grid point, not from $-\infty$; the integrand decays like $u^{-4}$ so this is
  accurate when the grid starts well below resonance (`get_ba_l3` warns otherwise; `get_ba_l2` does not), with an error of $\sim10^{-4}$ of $|C|$
  at the snapshot (§8);
- only the prograde mode $m=l$ carries a dynamical tide ($m=+2$ for $l=2$, $m=+3$ for $l=3$); the others never resonate before merger.

*At the snapshot* the dynamical part is small but not zero: $u=-3.91$ (heavy $(2,2)$) and $u=-2.83$ (heavy $(3,3)$), giving
$|b^{\rm dyn}|/|C|=0.34\%$ and $0.92\%$ (light: $1.2\%$ and $0.66\%$). It becomes dominant only if the snapshot is placed past resonance,
which is why the regime gate (§7) refuses that.

### 3.6 Assembling the amplitude array

`ba1` $=$ `ba_eq1` $+$ the dynamical term on the prograde mode. The returned arrays are indexed by mode:

| generator | order of the first index | $I_a$ used |
|---|---|---|
| `get_ba_l2` | $[m=+2,\ 0,\ -2]$ | $\sqrt{5k_2/4\pi}$, **same sign for all $m$** ($I_{2,-2}=(-1)^2I_{2,2}$) |
| `get_ba_l3` | $[m=+3,\ +1,\ -1,\ -3]$ | $\sqrt{7k_3/4\pi}\,[+,+,-,-]$, **signed** by $I_{l,-m}=(-1)^mI_{l,m}$ (Paper 2, after Eq. 11) |

For odd $l$ both $W_{l,-m}=-W_{lm}$ and $I_{l,-m}=-I_{lm}$, so $V_{l,-m}=V_{lm}$ — the drive is symmetric between $\pm m$, as it must be
for a real potential. The shape uses $C_a/I_a$ mode by mode, so these signs must travel with the amplitudes
(`data/default/star_modes.toml` stores $I$ per mode).

---

## 4. From amplitudes to the deformed surface

### 4.1 Displacement field and frames

By (2.3) and (3.3), in the frame that corotates with the orbit,

$$
\boldsymbol\xi(\mathbf x)=\sum_{\omega_a>0}2\,\mathrm{Re}\left[C_a\,\boldsymbol\xi_a(\mathbf x;\phi_s)\right]
\tag{4.1}
$$

with the companion on the $+x$ axis ($\theta=\pi/2$, $\phi_s=0$), $+z$ along the orbital angular momentum, and the orbit turning in the
$+\phi$ sense.

### 4.2 Linear surface

The stellar surface is a material surface. To first order in $\boldsymbol\xi$ its radial position is $R+\xi_r(R,\theta,\phi_s)$; the tangential part
slides fluid along the unperturbed sphere and does not change the shape at this order. Hence $\xi_h(R)$ is never needed and

$$
\frac{\Delta R(\theta,\phi_s)}{R}=\sum_{\omega_a>0}\frac{2\,\xi_{r,a}(R)}{R}\;\mathrm{Re}\left[C_a\,Y_{lm_a}(\theta,\phi_s)\right]
\tag{4.2}
$$

### 4.3 Closing the eigenfunction problem: calibration to the equilibrium tide

(4.2) needs $\xi_{r,a}(R)$, a property of the stellar model that the amplitudes do not carry. All modes of one $l$ share the same radial
function, so there is **one unknown constant per $l$**, $\hat\xi_l\equiv\xi_{r,a}(R)/R$. It is fixed by demanding that the adiabatic limit of
(4.2) reproduces the equilibrium tide (2.9).

In the adiabatic limit $C_a\to V_a=\frac{M'}{M}W_{lm}I_a(R/r)^{l+1}$, and the addition theorem
$\sum_mW_{lm}Y_{lm}(\theta,\phi_s)=P_l(\cos\gamma)$ turns (4.2) into $2\hat\xi_l\,I_l\,\frac{M'}{M}(R/r)^{l+1}P_l(\cos\gamma)$. Equating with (2.9):

$$
\hat\xi_l=\frac{h_l}{2I_l}\qquad\Longrightarrow\qquad
\boxed{\ \frac{\Delta R}{R}=\frac{h_l}{I_l}\;\mathrm{Re}\sum_aC_a\,Y_{lm_a}(\theta,\phi_s)\ }
\tag{4.3}
$$

Consequences:

- $C_a\propto I_l$, so only the ratio $C_a/I_l$ enters. The arbitrary sign of the eigenfunction normalization (Paper 1, footnote 1) cancels.
- The adiabatic limit is exact by construction; all dynamics, lag and spin live in $C_a$.
- It is an *effective* $\hat\xi_l$: the f-mode alone has a true $\xi_r(R)/R$ that is 2% ($l=2$: 2.411 vs the calibrated 2.363) to 3.5% ($l=3$:
  2.581 vs 2.492) larger, because the p-modes (small $I$, large surface $\xi_r$) subtract $\sim2\text{–}3\%$ of the adiabatic surface sum
  (§5.3). The p-modes are far above the driving frequency and respond adiabatically, so folding them into the f-mode is the right place.

### 4.4 The surface as equilibrium tide times transfer functions

Insert $C_a=T_aV_a$ with $V_a$ from (2.6) into (4.3). The overlap cancels:

$$
\hat C_a\equiv\frac{C_a}{I_a}=T_a\,\frac{M'}{M}\Big(\frac{R}{r}\Big)^{l+1}W_{lm}
\quad\Longrightarrow\quad
\frac{\Delta R}{R}=\sum_lh_l\,\frac{M'}{M}\Big(\frac{R}{r}\Big)^{l+1}\mathrm{Re}\sum_mT_{lm}\,W_{lm}\,Y_{lm}(\theta,\phi_s)
\tag{4.4}
$$

Setting every $T=1$ gives back the equilibrium tide $h_l\frac{M'}{M}(R/r)^{l+1}P_l(\cos\gamma)$. The star's shape is therefore the
static tidal bulge with each azimuthal component amplified and rotated by its own $T_{lm}$.

### 4.5 The $l=2$ surface in closed form

Write $Y_{22}=N_2\sin^2\theta\,e^{2i\phi_s}$, $Y_{2,-2}=Y_{22}^*$, $Y_{20}=N_0(3\cos^2\theta-1)$ (real) and $\mathbf n=(n_x,n_y,n_z)$. Because
$\mathrm{Re}\,Y_{22}=\mathrm{Re}\,Y_{2,-2}$ and $\mathrm{Im}\,Y_{22}=-\mathrm{Im}\,Y_{2,-2}$,

$$
\mathrm{Re}\!\left[C_aY_{22}+C_bY_{22}^*\right]=N_2\sin^2\theta\,\big[S_c\cos2\phi_s-S_s\sin2\phi_s\big],\qquad
S_c=\mathrm{Re}(C_a+C_b),\ \ S_s=\mathrm{Im}(C_a-C_b)
\tag{4.5}
$$

with $C_a,C_b,C_c$ the $(2,+2),(2,-2),(2,0)$ amplitudes. The $m=\pm2$ modes enter only through $C_a+C_b$ and $C_a-C_b$. The imaginary part of
the $m=0$ amplitude drops out, because $Y_{20}$ is real; it encodes fluid velocity ($\dot{\boldsymbol\xi}\propto\omega\,\mathrm{Im}\,q$), not shape. With
$\sin^2\theta\cos2\phi=n_x^2-n_y^2$ and $\sin^2\theta\sin2\phi=2n_xn_y$, (4.3) becomes, for $K=h_2/I_2$,

$$
\frac{\Delta R}{R}=e_{xx}n_x^2+e_{yy}n_y^2+e_{zz}n_z^2+2e_{xy}n_xn_y,\qquad
\begin{aligned}
e_{xx}&=K\,(N_2S_c-N_0\mathrm{Re}\,C_c), &\quad e_{yy}&=K\,(-N_2S_c-N_0\mathrm{Re}\,C_c),\\
e_{zz}&=2KN_0\,\mathrm{Re}\,C_c, & e_{xy}&=-KN_2S_s
\end{aligned}
\tag{4.6}
$$

a **traceless quadratic form**: a triaxial shape with $z$ a principal axis and the $x$–$y$ axes rotated by the lag. Equivalently, in transfer
functions, with $\varepsilon_t\equiv\frac{M'}{M}(R/r)^3$ and $T_c+iT_s\equiv T_{+2}+T_{-2}^*$ (using $W_{22}N_2=\frac38$, $W_{20}N_0=-\frac14$):

$$
e_{xx}=h_2\varepsilon_t\big[\tfrac38T_c+\tfrac14\mathrm{Re}\,T_0\big],\quad
e_{yy}=h_2\varepsilon_t\big[-\tfrac38T_c+\tfrac14\mathrm{Re}\,T_0\big],\quad
e_{zz}=-\tfrac12h_2\varepsilon_t\,\mathrm{Re}\,T_0,\quad
e_{xy}=-\tfrac38h_2\varepsilon_tT_s
\tag{4.7}
$$

Check: $T=1$ gives $(e_{xx},e_{yy},e_{zz})=h_2\varepsilon_t(1,-\tfrac12,-\tfrac12)$, i.e. $P_2(n_x)$ times $h_2\varepsilon_t$.

**Principal axes and the lag.** In the equatorial plane $\theta=\pi/2$,

$$
\frac{\Delta R}{R}\Big|_{\rm eq}=K\big\{N_2\rho\cos(2\phi_s+\delta)-N_0\mathrm{Re}\,C_c\big\},\qquad
\rho\cos\delta=S_c,\ \rho\sin\delta=S_s,\qquad
\psi=\frac\delta2=\tfrac12\,\mathrm{atan2}(S_s,S_c)
\tag{4.8}
$$

The long in-plane axis sits at $\phi_s=-\psi$ ($\psi>0$ = behind the companion), at $\phi_s=-\psi+\pi/2$ for the short one, and $z$ is
principal. (If the overlap is negative, $S_c,S_s$ are replaced by $S_c/I_2,S_s/I_2$, i.e. $\psi$ is computed from $\hat C=C/I$ so that it does not depend on the sign convention.)

### 4.6 General $(l,m)$: combined amplitudes, the lag reference and $X_{\rm lag}$

For $m>0$ the pair $\pm m$ contributes $\mathrm{Re}[\hat C_{+m}Y_{lm}+\hat C_{-m}Y_{l,-m}]$. With $Y_{l,-m}=(-1)^mY_{lm}^*$ and
$\mathrm{Re}\,z^*=\mathrm{Re}\,z$ this is $\mathrm{Re}[D_mY_{lm}]$ with

$$
D_m=\hat C_{+m}+(-1)^m\,\hat C_{-m}^{*}\quad(m>0),\qquad D_0=\mathrm{Re}\,\hat C_0,\qquad
\frac{\Delta R_{\rm tide}}{R}=\sum_lh_l\sum_{m\ge0}\mathrm{Re}\left[D_{lm}\,Y_{lm}(\theta,\phi_s)\right]
\tag{4.9–4.10}
$$

For $l=2$, $m=2$: $D_2=(C_a+C_b^*)/I_2=(S_c+iS_s)/I_2$, which ties (4.10) to (4.5). The code is `combined_amplitudes` and
`delta_r_over_R`.

**Lag of a general harmonic.** $Y_{lm}(\pi/2,\phi)=Y_{lm}(\pi/2,0)\,e^{im\phi}$, so the equatorial pattern is $|D_m||Y_{lm}|\cos(m\phi_s+\arg\tilde D_m)$ with

$$
\tilde D_m=s_{lm}D_m,\qquad s_{lm}=\mathrm{sign}\,Y_{lm}(\tfrac\pi2,0),\qquad
\psi_m=\frac{\arg\tilde D_m}{m}\quad(\text{behind the companion for }\psi_m>0)
\tag{4.11}
$$

The reference sign is necessary: $s_{22}=s_{31}=+1$ but $s_{33}=-1$ (Condon–Shortley), so an *adiabatic* $(3,3)$ amplitude $D_3$ is real and
**negative**. Measuring $\arg D_m$ instead would report a spurious $60^\circ$ lag for it. (Test V10's adiabatic limit guards this.)

**Lag exaggeration.** The knob $X_{\rm lag}$ rotates each harmonic's pattern and leaves its amplitude alone:

$$
\tilde D_m\ \to\ |D_m|\,e^{iX_{\rm lag}\arg\tilde D_m},\qquad\psi_m\to X_{\rm lag}\,\psi_m,\qquad |D_m|\ \text{unchanged}
\tag{4.12}
$$

valid for $|X_{\rm lag}\psi_m|<180^\circ/m$ ($\arg$ is not unwrapped). $X_{\rm lag}=1$ is the physics.

### 4.7 Spin oblateness (an addition to the papers' physics)

The papers' amplitudes are first order in $\Omega_s/\omega_1$ and omit the star's rotational flattening. Rigid rotation adds the centrifugal
potential $\Phi_c=-\tfrac12\Omega_s^2r^2\sin^2\theta=-\tfrac13\Omega_s^2r^2\left[1-P_2(\cos\theta)\right]$; its $l=2$ part is $+\tfrac13\Omega_s^2r^2P_2(\cos\theta)$,
$\theta$ measured from the spin axis ($\pm z$). The response of the surface to *any* static $l=2$ potential is the same shape number $h_2$ as for the tide
(§2.3), so with $\Phi_{\rm tid}\to\Phi_c$ in the equipotential argument

$$
\frac{\Delta R_{\rm spin}}{R}=-\varepsilon_\Omega\,P_2(\cos\theta),\qquad
\varepsilon_\Omega=\frac{h_2}{3}\,\frac{\Omega_s^{2}R^{3}}{GM}
\tag{4.13}
$$

i.e. poles $-\varepsilon_\Omega$, equator $+\varepsilon_\Omega/2$, independent of the sign of the spin; in the quadratic form of (4.6) it adds
$e_{xx},e_{yy}\mathrel{+}=\varepsilon_\Omega/2$, $e_{zz}\mathrel{-}=\varepsilon_\Omega$. *Check:* a uniform-density star ($h_2=5/2$) gives
$R_{\rm eq}-R_{\rm pol}=\tfrac32\varepsilon_\Omega R=\tfrac54\,qR$ with $q=\Omega^2R^3/GM$, the classical Maclaurin small-rotation flattening. The $l=0$ part of
$\Phi_c$ (a uniform radius shift) is dropped. Neither $X$ nor $X_{\rm lag}$ scales this term.

### 4.8 The printed surface and its frames

$$
r(\mathbf n)=R\left[\,1+X\,\frac{\Delta R_{\rm tide}(\mathbf n)}{R}+s_{\rm obl}\,\frac{\Delta R_{\rm spin}(\mathbf n)}{R}\right],\qquad
\frac{\Delta V}{V}=\Big\langle\Big(\frac{r}{R}\Big)^{3}\Big\rangle_{S^2}-1
\tag{4.14}
$$

with $X=X_{\rm lag}=1$ and $s_{\rm obl}=1$ for the sculpture build. The mean of $\Delta R_{\rm tide}$ over the sphere is zero (no $l=0$ term), so
$\Delta V/V$ is second order; the code reports it and does not correct it.

*Frames.* All formulas above place the companion on $+x$. The heavy star (companion at $+x$ in the world) uses them as they stand. For the light star
(companion at $-x$) the star-frame result is rotated by $\pi$ about $z$, $(x,y,z)\to(-x,-y,z)$: the map is an involution, preserves the lag relative to the companion,
leaves $l=2$ unchanged (since $m$ is even) and flips odd-$l$ terms. *Mesh:* each unit vector $\mathbf n$ of an icosphere vertex is displaced to
$\mathbf n\,r(\mathbf n)$; the post socket is drilled along $-z$ at the pole $r(0,0,-1)$, which includes the oblateness.

---

## 5. The stellar-model inputs ($I_l$, $h_l$, $\omega_{a0}$)

The amplitudes of §3 need $I_l$, $\omega_{a0}$ and $k_l$; the surface needs $h_l$. For the $\Gamma=2$ ($n=1$) Newtonian polytrope used by Paper 1
(compactness 0.16 in the paper, treated as Newtonian here) all four are computed in `src/polytrope.py`, in units $G=M=R=1$.

### 5.1 Mode problem

For a neutrally stratified star with adiabatic index equal to the polytropic index, the Eulerian force balance
$-\omega^2\boldsymbol\xi=-\nabla(h'+\Phi')$ with $h'=2K\rho'$ (since $P=K\rho^2$) makes the displacement a potential flow,
$\boldsymbol\xi=\nabla\chi/\omega^2$, $\chi=2K\rho'+\Phi'$. With $\rho'=-\nabla\!\cdot(\rho\boldsymbol\xi)$ and Poisson's equation,

$$
\frac{\omega^2}{2K}(\chi-\Phi')=-\nabla\!\cdot(\rho\nabla\chi),\qquad
\nabla^2\Phi'=\frac{4\pi G}{2K}(\chi-\Phi'),\qquad 2K=\frac{4GR^2}{\pi}
\tag{5.1}
$$

Two radial ODEs for $(\chi,\Phi')$ at fixed $l$ (regular at the centre, exterior matching $\Phi'_r+(l+1)\Phi'/R=0$ at the surface) are solved by Chebyshev
collocation as a generalized eigenproblem in $\omega^2$. $\xi_r=\chi'/\omega^2$ and $\xi_h=\chi/(\omega^2r)$; normalizing by (2.2) gives
$I_l$ from (2.4) and $\xi_r(R)$.

### 5.2 Exact static tide

At $\omega\to0$ the surface is an equipotential and $h'+\Phi'+\Phi_{\rm tid}=\text{const}$, so $\psi\equiv\Phi'+\Phi_{\rm tid}$ obeys
$\nabla^2\psi=-(4\pi G/2K)\,\psi=-(\pi/R)^2\psi$, a Helmholtz equation with regular solution $\psi=A\,j_l(\pi r/R)Y_{lm}$. Matching $\Phi'$ and $\Phi'_r$
to the exterior ($\Phi'\propto r^{-(l+1)}$, $\Phi_{\rm tid}\propto r^l$; $\rho(R)=0$ so there is no surface layer) and defining
$y\equiv\pi j_l'(\pi)/j_l(\pi)$:

$$
k_l^{\rm Love}=\frac{l-y}{l+1+y},\qquad h_l=1+k_l^{\rm Love},\qquad k_l=\tfrac12k_l^{\rm Love}
\tag{5.2}
$$

For $l=2$ this is $k_2^{\rm Love}=(15-\pi^2)/\pi^2=0.51982$, i.e. $k_2=0.25991$.

### 5.3 Numbers

| | $l=2$ | $l=3$ | $l=4$ |
|---|---|---|---|
| f-mode $\omega/\omega_0$, $\omega_0=\sqrt{GM/R^3}$ | 1.2270 | 1.6983 | 2.0365 |
| f-mode $I_l$ (magnitude) | 0.32153 | 0.24341 | 0.20758 |
| $k_l$ exact static (5.2) | 0.25991 | 0.10645 | 0.06024 |
| $\frac{4\pi}{2l+1}\sum_nI_n^2$ over modes | 0.25991 | 0.10645 | 0.06024 |
| $h_l=1+2k_l$ | **1.51982** | **1.21291** | 1.12048 |
| $\sum_n2\,\xi_{r,n}(R)\,I_n/R$ over modes | 1.5198 | 1.2129 | 1.1205 |
| f-mode $\xi_r(R)/R$ (true) / calibrated $h_l/2I_l$ | 2.411 / 2.363 | 2.581 / 2.492 | 2.813 / 2.699 |

The two sum-rule rows are independent of the static solution (5.2): they involve only mode quantities, and they reproduce $k_l$ and $h_l$ to 4–5 digits.
That is the numerical proof of (2.7), (2.9) and of the factor 2 in (2.3). The first p-mode overlap, $-5.5\times10^{-3}$, matches the paper's Table 1.

### 5.4 Scaling to the sculpture stars

The generator takes $k_l$, $h_l$, $\omega/\omega_0$ from the table, sets $\omega_{a0}=(\omega/\omega_0)\sqrt{GM/R^3}$ for each star, and
$\Lambda_2=k_2/(1.5\,C^5)$ with $C=GM/Rc^2$ (the form `get_ba_l2` expects); `get_ba_l3` takes $k_3$ directly. The overlaps it uses are
$I_2=\sqrt{5k_2/4\pi}=0.32158$ and $I_3=\sqrt{7k_3/4\pi}=0.24352$ (what the generator computes internally; the f-mode-only values in the table differ by
0.02% and 0.05%, immaterial since only $C/I$ matters). Resulting $\omega_{a0}/2\pi$: heavy 2164.7 Hz ($l=2$), 2996.2 Hz ($l=3$); light 1763 Hz, 2440 Hz.
$A_{\rm nl}$ is from (3.11).

---

## 6. Worked example: the heavy star, step by step

Heavy star: $M=1.6\,M_\odot$, $R=12$ km, $M'=1.2\,M_\odot$, spin $-500$ Hz ($\Omega_s=-3141.59\ \mathrm{s^{-1}}$), snapshot $f_{\rm GW}=1023.58$ Hz.
Reproduce with `python verification/math_walkthrough.py --star heavy`.

### 6.1 Orbit and detuning

From §1.1: $\omega=3215.66$, $r/R=2.7500$, $\dot\omega=2.7018\times10^5$, $\dot r/r=-56.01$. GR factors (3.10) at $x=0.12529$, $X_1=0.5714$, $X_2=0.4286$:
$z=0.92557$, $\omega_{\rm FD}/\omega=0.10857$. Mode: $\omega_{a0}=13601.3\ \mathrm{s^{-1}}$ (2164.7 Hz), $\omega_a=\omega_{a0}-\Omega_s=16742.9$, $A_{\rm nl}=-0.02145$.

Detuning of the $(2,+2)$ mode, term by term in (3.9) [$\mathrm{s^{-1}}$]:

| $z(\omega_a+2\Omega_s)$ | $-2\omega$ | $+2\,\omega_{\rm FD}$ | $A_{\rm nl}(\omega/\omega_{\rm ref})^2\omega_{a0}$ | $\Delta_2$ |
|---|---|---|---|---|
| $+9681.18$ | $-6431.33$ | $+698.23$ | $-305.60$ | $\mathbf{+3642.49}$ |

(Newtonian $\Delta=4028.38$; margins $\Delta/\omega_a=0.218$ with GR/NL, $0.241$ Newtonian.)

### 6.2 Drive, decay term and equilibrium amplitude

$I_2=0.321581$, $W_{22}=+0.97081$, $M'/M=0.75$, $(R/r)^3=1/2.75^3$:
$V_{+2}=\tfrac34\cdot0.048084\cdot0.97081\cdot0.321581=\mathbf{1.12586\times10^{-2}}$.
Decay term (3.5): $m\dot\omega=5.4036\times10^5$ and $-\Delta(l+1)\dot r/r=+6.1208\times10^5$, so $\kappa=\mathbf{1.1524\times10^6}\ \mathrm{s^{-2}}$,
$\kappa/\Delta^2=0.0869$ ($\arctan=4.96^\circ$). Eq. (3.6):
$b^{\rm(eq)}_{+2}=\dfrac{0.92557\cdot3642.5\cdot13601.3\cdot1.12586\times10^{-2}}{3642.5^2-i\,1.1524\times10^{6}}=3.8620\times10^{-2}+3.355\times10^{-3}i$.

Dynamical tide (§3.5): resonance $f_{\rm GW}\approx1636$ Hz (code) / 1616 Hz (exact root), $u=-3.91$, $F=3.531\times10^{-4}+3.059\times10^{-4}i$, $|b^{\rm dyn}|=1.33\times10^{-4}$, so
$C_{+2}=\mathbf{3.853\times10^{-2}+3.456\times10^{-3}i}$ ($|C|=3.869\times10^{-2}$; the $|b^{\rm dyn}|/|C|=0.34\%$).

### 6.3 Transfer functions and the other modes

| mode | $V_a$ | $\Delta_a\ [\mathrm{s^{-1}}]$ | $\kappa_a\ [\mathrm{s^{-2}}]$ | $T_a=C_a/V_a$ | $\lvert T_a\rvert$ | $\arg T_a$ |
|---|---|---|---|---|---|---|
| $(2,+2)$ | $+1.1259\!\times\!10^{-2}$ | $+3642$ | $1.152\!\times\!10^{6}$ | $3.4226+0.3070i$ | 3.436 | $5.13^\circ$ |
| $(2,0)$ | $-9.1926\!\times\!10^{-3}$ | $+12283$ | $2.064\!\times\!10^{6}$ | $1.0247+0.0140i$ | 1.025 | $0.78^\circ$ |
| $(2,-2)$ | $+1.1259\!\times\!10^{-2}$ | $+20924$ | $2.976\!\times\!10^{6}$ | $0.6016+0.0041i$ | 0.602 | $0.39^\circ$ |
| $(3,+3)$ | $-2.3918\!\times\!10^{-3}$ | $+3009$ | $1.485\!\times\!10^{6}$ | $5.6708+0.9666i$ | 5.753 | $9.67^\circ$ |
| $(3,+1)$ | $+1.8527\!\times\!10^{-3}$ | $+12620$ | $3.098\!\times\!10^{6}$ | $1.3802+0.0269i$ | 1.381 | $1.11^\circ$ |
| $(3,-1)$ | $+1.8527\!\times\!10^{-3}$ | $+22230$ | $4.711\!\times\!10^{6}$ | $0.7838+0.0075i$ | 0.784 | $0.55^\circ$ |
| $(3,-3)$ | $-2.3918\!\times\!10^{-3}$ | $+31840$ | $6.323\!\times\!10^{6}$ | $0.5472+0.0034i$ | 0.547 | $0.36^\circ$ |

The ordering $|T_{+m}|>|T_0|>|T_{-m}|$ of §3.3 is visible, and the prograde $(2,+2)$, $(3,+3)$ modes carry most of the lag because they are the ones
nearest resonance.

### 6.4 $l=2$ surface

$\varepsilon_t=\frac{M'}{M}(R/r)^3=0.036065$, $h_2\varepsilon_t=0.054813$, $T_c+iT_s=T_{+2}+T_{-2}^*=4.02424+0.30288i$, $\mathrm{Re}\,T_0=1.02469$. From (4.7):

$$
e_{xx}=0.054813\,[\tfrac38\cdot4.02424+\tfrac14\cdot1.02469]=\mathbf{+0.09675},\quad
e_{yy}=\mathbf{-0.06867},\quad e_{zz}=\mathbf{-0.02808},\quad e_{xy}=-\tfrac38\cdot0.054813\cdot0.30288=\mathbf{-0.006225}
$$

The same values from (4.6) with $K=h_2/I_2=4.7261$, $S_c=4.5307\times10^{-2}$, $S_s=3.410\times10^{-3}$. The lag (4.8): $\delta=\mathrm{atan2}(S_s,S_c)=4.304^\circ$,
$\psi=\mathbf{2.152^\circ}$. Compared with the *adiabatic* tip $h_2\varepsilon_t=0.0548$, the dynamics raise the bulge toward the companion to $0.0968$ (×1.77).

### 6.5 $l=3$ and the full tidal surface

$\hat C_{3,+3}=C/I=-5.570\times10^{-2}-9.494\times10^{-3}i$, $\hat C_{3,-3}=+5.375\times10^{-3}+3.35\times10^{-5}i$, so $D_3=\hat C_{+3}-\hat C_{-3}^*=-6.107\times10^{-2}-9.460\times10^{-3}i$,
$\tilde D_3=-D_3$ ($s_{33}=-1$), $\arg\tilde D_3=8.80^\circ$, $\psi_3=\mathbf{2.935^\circ}$; $\psi_1=0.513^\circ$. The $l=3$ tip toward the companion is $+0.03736$
(adiabatic $h_3\frac{M'}{M}(R/r)^4=0.01591$, ×2.35). Along the axes the tidal surface is:

| direction | $+x$ | $+y$ | $\pm z$ | $-x$ |
|---|---|---|---|---|
| $\Delta R_{\rm tide}/R$ | $+0.1341$ | $-0.0639$ | $-0.0281$ | $+0.0594$ |

($+x$ and $-x$ differ because the $l=3$ part is odd.)

### 6.6 Spin oblateness and the final radius

$\varepsilon_\Omega=\frac{h_2}{3}\frac{\Omega_s^2R^3}{GM}=0.04069$; pole $-0.04069$, equator $+0.02034$. Adding (4.13) to the tidal surface:

| direction | $+x$ | $+y$ | $\pm z$ | $-x$ |
|---|---|---|---|---|
| $\Delta R/R$ (tide + spin) | $+0.1545$ | $-0.0436$ | $-0.0688$ | $+0.0797$ |

($+y$: $-0.0639+0.0203$; $z$: $-0.0281-0.0407$.) With $R=8.691$ mm on the sculpture ($S=0.72428$ mm/km) the tip toward the companion is
$8.691\times1.1545=10.03$ mm, the bottom pole $8.691\times0.9312=8.09$ mm.

### 6.7 The light star, same pipeline

$M=1.2\,M_\odot$, $R=12.5$ km, $M'=1.6\,M_\odot$, spin $-200$ Hz: $r/R=2.640$, $\varepsilon_t=0.07246$, $h_2\varepsilon_t=0.11013$;
$T_{+2}=3.5742+0.4990i$, $T_{-2}=0.6165+0.0052i$, $T_0=1.0527+0.0188i$, $T_c+iT_s=4.1907+0.4939i$; $(e_{xx},e_{yy},e_{zz},e_{xy})=(0.20206,-0.14409,-0.05797,-0.020397)$,
$\psi=3.361^\circ$, $\psi_3=2.518^\circ$; $\varepsilon_\Omega=0.00981$; tide tip $+0.2658$, with spin term $+0.2707$.

---

## 7. What the math does and does not contain

**Included.** Newtonian star with $l=2$, $l=3$ f-modes (p-modes folded in at the adiabatic level, §4.3); point-mass companion on a quasi-circular
inspiral; finite-frequency, orbital-decay and resonant (dynamical-tide) response of each mode; spin through $\omega_a$, $\Delta_a$ and the
optional centrifugal term; GR and nonlinear-frequency shifts as phenomenological terms (3.9)–(3.11).

**Approximations and their sizes** (all checked or measured in §8):

1. *Linear surface.* Second-order terms are dropped; the light star's largest tidal $\Delta R/R$ is 0.27, so they are $\sim(0.27)^2$ of the radius. The gate (V11 of
   `neutron_star_shape.md`) warns above 0.3 and above $|\Delta V/V|=5\%$.
2. *Truncation of Eq. 18.* Relative error $\sim(\omega_at_{\rm gw})^{-2}$ away from resonance and up to $\sim4\%$ for the prograde mode close to it:
   1.0% for $(2,2)$ and 3.7–4.2% for $(3,3)$ at the snapshot against the exact ODE.
3. *Resonance location.* The gate uses the **Newtonian** $\Delta_a$. GR and the nonlinear shift move the margin by $\le0.05$ for the resonant modes
   here ($\Delta_2/\omega_a$: heavy $0.241\to0.218$, light $0.275\to0.223$) and never change its sign in the final configuration; for the rejected
   light star at $-500$ Hz the exact $(3,3)$ resonance is at 997 Hz, still below the 1023.6 Hz snapshot.
4. *Calibration.* $\hat\xi_l=h_l/2I_l$ is 2–3.5% below the f-mode's true surface displacement (§4.3).
5. *Stellar structure* is the Newtonian $\Gamma=2$ polytrope (compactness enters only through $z$, $\omega_{\rm FD}$ and the unit conversions); a different EOS
   changes $k_l$, $h_l$, $I_l$, $\omega_{a0}$ but not the structure of the formulas.
6. *Spin.* Eigenfunctions are not corrected for rotation (first order in $\Omega_s/\omega_1$); the oblateness term is first order in $\Omega_s^2$.
7. *Snapshot.* $b^{\rm dyn}$ oscillates in the orbit frame, so past resonance the surface depends on exactly where the snapshot is taken; at the final snapshot
   it is $0.3$–$1.2\%$ of $|C|$.
8. *Not included:* $l=0$ and $l\ge4$ modes, g-modes (the polytrope has $N=0$), mode couplings in the surface (they enter only through $A_{\rm nl}$), eccentricity, and
   tidal dissipation beyond what orbital decay produces.

---

## 8. Verification map

| Equation / claim | Check | Result |
|---|---|---|
| (3.6) from the series (3.4)–(3.5) | `test_resummation_reproduces_eq18_in_the_point_mass_limit` (finite-difference $\dot b^{(0)}$) | agrees to $10^{-8}$; $\arg b=\arctan(\kappa/\Delta^2)$ |
| (3.2) solved by (3.6)+(3.14) | `verification/verify_amplitudes.py` (time-domain ODE, sculpture stars) | non-resonant modes $\le5\times10^{-4}$; prograde 1.0% / 3.7–4.2% |
| Independent chain §1–§3 vs the generator | `test_independent_chain_reproduces_the_generated_amplitudes` | non-resonant $\le10^{-12}$ (float64 round-off; measured $3\times10^{-15}$), prograde $\le4.4\times10^{-4}$ (trapezoid accumulation of $F$) |
| (3.14) is small at the snapshot | `test_dynamical_tide_is_small_but_nonzero_at_the_snapshot` | 0.3–1.2% of $\lvert C\rvert$ |
| resonance root | `test_resonance_root_and_closed_form` | $\Delta_a(\omega_r)=0$ to $10^{-12}$; 1PN form within 3% |
| (2.7), (2.9), (5.2) | `test_polytrope_shape_numbers`; `polytrope.py` $l=2,3,4$ | $k_l$, $h_l$ from mode sums and static solution agree to $2\times10^{-5}$ |
| (4.3) adiabatic limit | `test_v1_adiabatic_limit` ($l=2,3$, both overlap signs) | $10^{-12}$ |
| (4.4) $I$ cancels | `test_surface_equals_equilibrium_tide_times_transfer_function` | $10^{-13}$ |
| (4.6)–(4.7) closed form | `test_v5_quadratic_form_matches_harmonics`; `test_l2_closed_form_in_transfer_functions`; `test_W_N_products` | $10^{-13}$; $10^{-12}$ |
| (4.8) axes and lag | `test_v6_principal_axes_follow_lag_angle`; `test_v8_lag_sign_and_adiabatic_limit` | exact to $10^{-8}$ deg; $\psi>0$ and $\to0$ as $\dot\omega\to0$ |
| (4.11)–(4.12) | `test_v10_*` (incl. zero lag for adiabatic $(3,3)$) | exact |
| (4.13) | `test_v12_*` (axisymmetry, spin-sign independence, Maclaurin $5q/4$) | exact |
| (4.14) mesh | `test_v7_*` (bare-mesh volume vs $\langle(r/R)^3\rangle$; gap; socket at the pole) | $2\times10^{-3}$ (mesh discretization) |

Run: `python -m pytest tests` (48 tests), `python verification/math_walkthrough.py`, `python verification/verify_amplitudes.py`.

---

## Appendix A: pitfalls

1. **Mode order differs between the generators.** `get_ba_l2` returns $[+2,0,-2]$, `get_ba_l3` $[+3,+1,-1,-3]$. The tool `make_star_modes.py` maps both to explicit $(l,m)$.
2. **The overlap $I$ is signed for odd $l$** ($I_{l,-m}=(-1)^mI_{l,m}$), and the shape divides by it mode by mode. Using $|I|$ throughout flips the $-m$ contributions.
3. **$\psi$ and $D_m$ must be built from $C/I$**, not $C$: with $I<0$ the raw $C$ flips sign and $\psi$ jumps by $90^\circ$ for $l=2$.
4. **The lag reference is $s_{lm}D_m$**, not $D_m$ ($s_{33}=-1$).
5. **The factor 2 of $\mathrm{Re}$** (2.3) is absorbed in $k_l=4\pi I^2/(2l+1)$ and in $\hat\xi_l=h_l/2I_l$; it is not an extra factor in (4.3).
6. **$h_l=1+2k_l$, not $1+k_l$**, because the paper's $k_l$ is half the classical Love number.
7. **$\mathrm{Im}\,C_c$ for $m=0$ is not a shape change** (velocity only); only $\mathrm{Re}\,C_c$ enters.
8. **Companion on $-x$:** compute with the companion on $+x$, then rotate by $\pi$ about $z$ (odd $l$ flips).
9. **Gate margins are Newtonian**; the amplitudes themselves use the GR/NL-corrected $\Delta_a$.
