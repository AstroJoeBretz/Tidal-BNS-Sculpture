# Binary Neutron Star Gravitational-Wave Sculpture — Design

Two 3D-printed, tidally deformed neutron stars on posts above a CNC-milled hardwood
panel carrying the gravitational-wave field of their own orbit.

This document describes the system as built: the physical model, the geometry, the
manufacturing constraints and the checks. The neutron-star surfaces are specified in
[`neutron_star_shape.md`](neutron_star_shape.md) and derived in [`math.md`](math.md).
Section numbers are referenced from the code (`docs/design.md sec N`).

---

## 0. Design intent, and what is not physical

Everything in this piece is at a single true scale $S$ (mm per km) except the following:

- the **wave speed used for retardation** is $0.35\,c$ (§3.5), which winds the spiral more
  tightly than the real speed of light would;
- the **field is evaluated off the orbital plane** (`Z_EVAL_MM`, §3.4) and its amplitude is
  normalized to a chosen relief depth (§4.3), so the relief depth carries no strain value;
- the **neutron-star spins** (−500 Hz heavy, −200 Hz light, anti-aligned) are design choices, picked to
  make the tidal lag and the flattening legible (`neutron_star_shape.md` §3.6).

Not exaggerated: the orbital frequency, the separation, the stellar radii, the mass
ratio, the barycentre offset and the handedness of the spiral. The tidal deformation of
the stars is **not** exaggerated either ($X=1$): the linear f-mode response of a
$\Gamma=2$ polytrope at this separation is a visible $1$–$2.5$ mm bulge (§6.2).

**Snapshot semantics.** The panel is the gravitational-wave field frozen at $t_0 = 0$.
The stars are at their $t_0$ positions. Because the field at radius $r$ was emitted at
$t_0 - r/c$, the outer arms show the binary as it was several orbits earlier — that is
the point of carving retarded time into wood.

---

## 1. Master parameters

### 1.1 Physical system

| Symbol | Value | Notes |
|---|---|---|
| $M_1$ | 1.6 M☉ | heavier star |
| $M_2$ | 1.2 M☉ | lighter star, $q = 0.75$ |
| $M$ | 2.8 M☉ | total |
| $\mu$ | 0.6857 M☉ | reduced mass |
| $R_1$ | 12.0 km | NS radii are nearly mass-independent |
| $R_2$ | 12.5 km | lighter star slightly larger |
| $a$ | **33.0 km** | orbital separation — see §1.3 |
| $\Omega$ | 3215.7 rad/s | $\sqrt{GM/a^3}$, Keplerian |
| $f_{\rm orb}$ | 511.8 Hz | |
| $f_{\rm GW}$ | 1023.6 Hz | $=2 f_{\rm orb}$ |
| $\lambda_{\rm GW}$ | 292.89 km | $=c / f_{\rm GW}$ |
| $v_{\rm orb}/c$ | 0.354 | $\Omega a/c$; mildly relativistic, the Newtonian orbit is ~10% off (§10) |
| $e$ | 0 | circular orbit assumed |

Constants: $G = 6.67430\times10^{-11}$, $c = 2.99792458\times10^{8}$, $M_\odot = 1.98847\times10^{30}$ (SI).

### 1.2 Master scale

$$
S=\frac{\Lambda}{\lambda_{\rm GW}}=0.72428\ \mathrm{mm/km},\qquad
\Lambda=\texttt{SCALE\_LAMBDA\_MM}=150\sqrt2=212.13\ \mathrm{mm}
$$

One GW wavelength is $212.13$ mm on the finished piece. The reference length $\Lambda$ is a
parameter independent of the panel size, so resizing the panel does not rescale the stars, rods or
any printed part; a smaller panel simply windows the same field. The 254 mm panel's half-diagonal
(179.6 mm) is $0.85\,\lambda_{\rm GW}$.

Derived geometry (mm on the finished piece):

| Quantity | Value |
|---|---|
| separation $a$ | 23.901 |
| star 1 radius $R_1$ | 8.691 (Ø17.4) |
| star 2 radius $R_2$ | 9.053 (Ø18.1) |
| star 1 centre | $x = -10.243$ |
| star 2 centre | $x = +13.658$ |
| $\lambda_{\rm GW}$ | 212.13 |
| near-zone radius $\lambda/2\pi$ | 33.76 |

The barycentre offset is real: the heavier star sits closer to the origin. This asymmetry
is visible and should not be "cleaned up."

### 1.3 Why $a = 33$ km

$\lambda_{\rm GW}/a=\pi c/v_{\rm orb}\gtrsim 8$ for any binary neutron star, so a single global scale
cannot make both the binary and many wavelengths large simultaneously. Tightening the
spiral shrinks the stars. $a = 33$ km is the best compromise:

- near contact ($R_1 + R_2 = 24.5$ km), so tidal deformation is large — the subject of the piece;
- stars land at Ø17–18 mm, large enough to print the $l=2,3$ structure legibly;
- a 2.36 mm surface gap between the deformed stars (§6.2), so they read as two objects;
- one full two-armed wrap per 212.13 mm of radius (the 254 mm panel reaches 0.85 of it at
  the corners and 0.60 at the edge midpoints).

To retune, change `binary.a_km` and `panel.scale_lambda_mm` in a configuration file (§9) and re-run
`build_all.py`, which regenerates the star amplitudes and whose validation step re-checks the
parameters that depend on them.

---

## 2. Coordinate frames

**World frame** — right-handed, origin at the binary barycentre.

- $x, y$ — the orbital plane, and also the plane of the wood panel;
- $z$ — orbital angular momentum, pointing up out of the panel;
- $z = 0$ — the **top surface of the wood stock** (CAM datum; all cuts are $z \le 0$);
- at $t_0 = 0$ both stars lie on the $x$ axis (heavy at $-x$, light at $+x$); the binary rotates in $+\hat\phi$.

**Star frame** — per star, origin at the star's centre, companion on $+x$, $+z$ the orbital
angular momentum (`neutron_star_shape.md` §2). The light star's frame is mapped to the world by a
rotation of $\pi$ about $z$. The deformed surface is triaxial and slightly rotated by the tidal
lag, so the full $(l,m)$ harmonic content is used, not only $m=0$.

**Panel plane.** The panel is the orbital plane $z=0$. The field is a scalar potential
difference, not a strain tensor, so it has no cross-polarization component to worry about;
evaluating it off the plane at $z = $ `Z_EVAL_MM` is the mechanism that controls the dynamic range
(§3.4).

---

## 3. Wave field model

### 3.1 Why not a point quadrupole

The exact point-quadrupole outgoing solution (a spherical Hankel function $h_2(kr)$) is exact
only for a point source at the origin. The source here is two masses separated by
$a = 24$ mm, and the visible plate begins at only 2–3 $a$. In that range no point-multipole
treatment of an extended two-body source is trustworthy, whatever radial function is attached to
it, and its $1/r^3$ growth toward the origin would need an arbitrary blanked disk. The field below
keeps each star as what it is — a point mass on its own orbit.

### 3.2 Retarded per-star Newtonian potential

Each star is a point mass on its own circular orbit, and its potential is evaluated at its own
**self-consistently retarded position**:

$$
\Phi(\mathbf x,t)=-\frac{G\,m_1}{\lvert\mathbf x-\mathbf x_1(t_{{\rm ret},1})\rvert}
-\frac{G\,m_2}{\lvert\mathbf x-\mathbf x_2(t_{{\rm ret},2})\rvert},\qquad
t-t_{{\rm ret},i}=\frac{\lvert\mathbf x-\mathbf x_i(t_{{\rm ret},i})\rvert}{c_{\rm eff}}
$$

The implicit equation is solved per star by fixed-point iteration on the retardation phase
$\varphi_i=k_{1/2}\,r_i$, with $k_{1/2}=\Omega/c_{\rm eff}$ at scale (half the GW wavenumber: that one is
defined for the $\omega_{\rm GW}=2\Omega$ radiation, but each star's own orbital retardation runs at
$\Omega$). The star's retarded angular position is $\theta_{i0}-\varphi_i$. The iteration contracts
geometrically ($d(k_{1/2}r)/d\varphi\approx k_{1/2}d_{\rm star}\approx0.18$ here), so a fixed count
(`RETARDED_ITERS = 10`) converges to double precision.

The potential is finite everywhere except at the two stellar surfaces, so no near-zone window is
needed. Near each star the light-crossing time across a near-zone distance is negligible against the
orbital period, so $t_{{\rm ret},i}\approx t$ and it reduces to ordinary two-body Newtonian gravity —
the physics that also shapes the stars (§6). Far away, expanding in multipoles is the standard
derivation of the quadrupole radiation formula, and because the orbit is exactly periodic, exact
retardation resums the whole derivative series rather than truncating at quadrupole order.

### 3.3 What gets carved: retarded minus instantaneous

$$
\delta\Phi(\mathbf x)=\Phi_{\rm retarded}(\mathbf x,t_0)-\Phi_{\rm instantaneous}(\mathbf x,t_0)
$$

$\Phi_{\rm instantaneous}$ is the same two-body potential with **no** retardation — both stars
evaluated at their $t_0=0$ positions directly. This isolates the propagation-delay content:

- **Near each star**, retarded $\approx$ instantaneous (§3.2), so $\delta\Phi\to0$ smoothly on its
  own — this is what removes the divergence, not a window.
- **Far away**, the static quadrupole-like piece of $\Phi_{\rm instantaneous}$ ($\sim1/r^3$) decays
  faster than the radiative tail of $\Phi_{\rm retarded}$ ($\sim1/r$), so $\delta\Phi\to\Phi_{\rm retarded}$
  and the trailing spiral emerges.

The carved scalar is $u_{\rm raw}=\delta\Phi/c^2$ — dimensionless, in the spirit of the weak-field
relation $h\sim\Phi/c^2$.

**Honesty check.** This is an ansatz, not a solved wave equation: its numerical prefactor does not
match the textbook quadrupole formula's $4G/c^4\,\ddot Q_{ij}$. That is acceptable because the
relief is normalized to a depth budget (§4.3) whichever field model feeds it, and the piece makes no
claim about real strain magnitudes (§0). What matters is the *structure*: the correct near-zone limit,
the correct qualitative far-zone radiation character, and no artificial divergence.

### 3.4 Dynamic range and the z-offset

Even finite, $\delta\Phi$ has a dynamic-range problem: right next to a star the retardation *phase*
is small, but it multiplies the star's own large potential, so the near-star region dwarfs the
far-field radiative content. Evaluated in-plane ($z=0$) the near-star peaks are about 4000× the mean edge value on the 0.5 mm grid (the value is grid-dependent, since the potential is singular at each star).

**Fix: evaluate the field off-plane**, at $(x,y,\texttt{Z\_EVAL\_MM})$. This acts as a softening
length: it bounds the minimum source distance to $\lvert\texttt{Z\_EVAL\_MM}\rvert$, which is what
controls the dynamic range. Larger magnitudes suppress the near-star peaks *and* relax milling
slope/curvature together, but push the retardation phase toward a near-constant offset dominated by
`Z_EVAL_MM` itself, which washes out the in-plane spiral structure if taken too far.

Measured on the current configuration (254 mm panel, $A_{\rm pp}=20$ mm, $c_{\rm eff}=0.35c$):

| `Z_EVAL_MM` | edge mean $\lvert u\rvert/N$ | max slope | min concave $R_c$ | look |
|---|---|---|---|---|
| −30 | 0.065 | 43.1° | 12.4 mm (7.8×) | tight bowtie, background nearly flat |
| −50 | 0.187 | 34.8° | 17.5 mm (11.0×) | same shape, softer |
| **−70** | **0.251** | **36.9°** | **16.5 mm (10.4×)** | **clear trailing two-armed pinwheel, visible background gradient — chosen** |
| −100 | 0.435 | 35.2° | 19.5 mm (12.3×) | pattern reads as more diffuse/rounded |

($N$ is the normalization of §4.3; $R_c$ in parentheses is in units of the ball radius.)

`Z_EVAL_MM = −70` was chosen by inspection of the rendered heightmap, not a formula: it is the
largest softening depth that still reads unambiguously as a two-armed spiral rather than a bowtie
or a diffuse blob. It is an artistic tuning parameter, not a physically derived softening length.

### 3.5 Hypothetical wave speed — `C_EFF_FACTOR`

`binary.retardation_k_half(b, c_eff)` decouples the retardation wavenumber from the real speed of
light. It changes only how tightly the pattern winds, not the plate/star geometry, which stays fixed
at the scale of §1.2. `params.C_EFF_FACTOR = 0.35` multiplies the retardation phase at every radius by
$1/0.35=2.86$ relative to the real speed of light. Set `C_EFF_FACTOR = 1.0` to recover the real
speed and the roughly single-wrap pattern. Re-verify `A_PP_MM` against `validate.py` whenever it
changes: both slope and curvature depend on it.

---

## 4. Panel geometry

### 4.1 Stock

| | |
|---|---|
| Footprint | **254 × 254 mm (10″ × 10″)** |
| Thickness | **36 mm** (≈ 1.42″) |
| Species | **Hard maple** (recommended) — tight grain, holds crisp detail, pale enough that shadows carry the pattern. Cherry acceptable. |
| Avoid | Baltic birch ply (contour rings fight the spiral), MDF (perfect detail, dead surface), open-pore oak/ash |

### 4.2 Relief

| | |
|---|---|
| Peak-to-valley $A_{\rm pp}$ | **20.0 mm** |
| Datum | $z = 0$ at top of stock; relief occupies $z\in[-20,0]$ (the field as built reaches −19.29 mm) |
| Mean plane | $z = -10.0$ mm |
| Material below deepest valley | 16.0 mm budget, 16.7 mm as built (36 mm stock) |

Slope and curvature scale linearly with $A_{\rm pp}$ (§5.1). **Post-hole clearance** (§7) is
checked explicitly (`check_post_hole_clearance`, `MIN_HOLE_CLEARANCE_MM = 3.0` mm): at the stars the
relief floor is −9.45 / −10.33 mm, leaving 18.5 / 17.7 mm of wood below the holes in the 36 mm stock,
so the post holes do not limit $A_{\rm pp}$. The remaining limits are the slope guard (36.9° of 45°) and
the curvature margin (10.4× the ball radius), §5.1; the real constraint on going deeper is surface
quality on steep faces, not a hard failure.

### 4.3 Height field

$$
\begin{aligned}
z_{\rm eval}&=\texttt{Z\_EVAL\_MM}=-70\ \mathrm{mm}&&\text{softening depth, §3.4}\\
c_{\rm eff}&=0.35\,c&&\text{§3.5}\\
u_{\rm raw}(x,y)&=\big[\Phi_{\rm ret}(x,y,z_{\rm eval};c_{\rm eff})-\Phi_{\rm inst}(x,y,z_{\rm eval})\big]/c^2&&\text{§3.2–3.3}\\
u(x,y)&=\mathrm{clamp}\!\big(u_{\rm raw}/N,\,-1,\,+1\big)\\
z(x,y)&=-\tfrac12A_{\rm pp}\,\big(1-u\big)&&\Rightarrow z\in[-A_{\rm pp},0]
\end{aligned}
$$

**No edge window:** the field is carved all the way to the panel's cut edge rather than tapering
flat over the outer margin; a taper's own curvature would be the tightest feature on the plate.

**Normalization constant — computed, not chosen:** $N=\max\lvert u_{\rm raw}\rvert$ over the plate sampled at
the 0.5 mm analysis step (`ANALYSIS_STEP_MM`) $=5.403\times10^{-4}$, at $(x,y)=(-88.5,-23.0)$ mm. $N$, the
milling checks (§5) and the star-height solve (§7) all use that analysis grid, so the carved field does not
depend on the export resolution `panel.grid_step_mm` (§8). The field at the stars is small: at
`Z_EVAL_MM = −70`, $u_{\rm raw}/N=+0.055$ at the heavy star and $-0.033$ at the light star, so both
stars sit near the mid-plane of the relief (the post holes start at −9.45 / −10.33 mm). The extremes lie on a
ring of radius ≈ 90–100 mm. Local maxima of $\lvert u\rvert$:

| location $(x,y)$ mm | $u_{\rm raw}/N$ |
|---|---|
| (−88.5, −23.0) | **+1.000 (global peak)** |
| (−27.5, +96.0) | −0.929 |
| (+92.0, +22.5) | +0.924 |
| (+28.5, −91.0) | −0.907 |
| (+40.0, +16.5) | −0.512 |

**Rim values**, $c_{\rm eff}=0.35c$: $\lvert u/N\rvert$ averages **0.25** along the cut edge (range
0.00–0.46) and is 0.19–0.35 at the four corners ($\lvert x\rvert,\lvert y\rvert>110$ mm: mean 0.36). At
$A_{\rm pp}=20$ mm the edge swings $\pm2.5$ mm on average ($\pm4.6$ mm at its strongest) and the corners
$\pm1.9$ to $\pm3.5$ mm.

---

## 5. Milling constraints (1/8″ ball, $R_b = 1.5875$ mm)

### 5.1 Derived limits

**Scallop height** $h_s = R_b - \sqrt{R_b^2 - (s/2)^2}$:

| stepover $s$ | scallop |
|---|---|
| 0.50 mm | 0.020 mm |
| **0.75 mm** | **0.045 mm** — use this |
| 1.00 mm | 0.081 mm |

**Slope and concave radius of curvature** at $A_{\rm pp}=20.0$ mm, $c_{\rm eff}=0.35c$,
$z_{\rm eval}=-70$ mm:

| | value | location |
|---|---|---|
| max slope | **36.9°** (guard 45°) | (−81.5, +58.0) mm |
| min concave $R_c$ | **16.5 mm (10.4× the ball radius)** | (−96.5, +31.5) mm |

Both extremes lie on the ring at $r\approx100$ mm (§4.3). A 2-D sweep of $A_{\rm pp}$ on the 254 mm
plate gives, for slope / $R_c$ / $R_c$ over ball radius: 16 mm 31.0° / 20.6 / 13.0×; 18 mm 34.0° / 18.3 / 11.5×;
**20 mm 36.9° / 16.5 / 10.4×**; 24 mm 42.0° / 13.7 / 8.6×; 28 mm 46.4° / 11.8 / 7.4×. Curvature is
never the limit (floor $2R_b=3.175$ mm); the 45° slope guard is reached at ≈ 26.7 mm.

**Cost of steep faces.** The stepover is 0.75 mm in the horizontal; along a face at slope $\theta$ the
pass spacing is $s/\cos\theta$, so at 36.9° it is 0.94 mm and the scallop rises from **0.045 mm to
0.071 mm** (≈ the 1.0 mm row above). Use a 0.5 mm stepover on the steepest faces if that matters
(0.020 mm flat, ≈ 0.031 mm at 37°), at ~1.5× the finishing time. Tearout across the grain also gets
worse as faces steepen, so the sharp-bit, climb-mill and 45°-to-grain advice in §5.2 matters. A 1/8″
ball and a 1/8″ shank clear a 37° face, but check the tool's flute length against the 19.3 mm depth
range (the roughing pass should remove most of it).

Single 3-axis finishing pass, no rest machining, no steep-wall strategy (those only become relevant
above ~45–60°).

**Chordal sag** of the exported mesh: $h^2/(8R_c)$ is 0.0019 mm at 0.5 mm spacing and ≈ 0.03 mm at the default
2 mm spacing ($R_c=16.5$ mm). Measured against the field itself, the 2 mm plate surface deviates by at most
0.030 mm (rms 0.009 mm), below the 0.045 mm scallop of the finishing pass.

### 5.2 Toolpath

1. **Rough** — 1/4″ flat endmill, adaptive or horizontal, 0.5 mm stock-to-leave.
2. **Finish** — 1/8″ ball, parallel raster, 0.75 mm stepover, ~2500 mm/min, **oriented 45° to the grain**
   to split tearout risk between the two bad quadrants.
3. **Border — no roundover.** The relief runs a live, non-trivial value right to the cut line
   (§4.3: $\lvert u\rvert$ averages 0.25 of peak along the edge, up to 0.46). A roundover bit would cut
   into the pattern itself, so the top edge is left sharp, square to the perimeter cut.
4. **Post holes** — see §7.

Estimated finishing time: 339 passes × 254 mm ≈ 86 m of path ≈ 35 min at 2500 mm/min.

**Grain.** The spiral crosses grain at every angle; the against-grain quadrants will fuzz. Sharp
bit, climb-mill the finish pass, and hand-sand with a foam block 180→320 following the ridges (never
across them).

---

## 6. Neutron star geometry

### 6.1 Surface

The surface $r(\theta,\phi_s)$ of each star is built from the complex f-mode amplitudes of Yu et al.
(arXiv:2211.07002, 2404.00147) at the snapshot orbital frequency, by Love-number-calibrated linear
displacement (`neutron_star_shape.md`, derivation in `math.md`):

$$
r(\theta,\phi_s)=R\left[1+X\,\frac{\Delta R_{\rm tide}(\theta,\phi_s)}{R}
+s_{\rm obl}\,\frac{\Delta R_{\rm spin}(\theta)}{R}\right],\qquad
\frac{\Delta R_{\rm tide}}{R}=\sum_l\frac{h_l}{I_l}\,\mathrm{Re}\sum_a C_a\,Y_{lm_a}
$$

with $l=2,3$, the complex amplitudes $C_a$ stored in the generated `star_modes.toml`, the shape numbers $h_l$ of a
$\Gamma=2$ polytrope, $X=1$ (no exaggeration) and spin oblateness on. The surface is triaxial, its
long axis trails the companion by the tidal lag $\psi$, and the $l=3$ term makes the near-side bulge
larger than the far-side one.

### 6.2 What is printed

| | heavy | light |
|---|---|---|
| Mass, radius, companion | 1.6 M☉, 12.0 km, 1.2 M☉ | 1.2 M☉, 12.5 km, 1.6 M☉ |
| Spin (anti-aligned) | −500 Hz | −200 Hz |
| Separation at snapshot $r/R$ | 2.75 | 2.64 |
| Tidal lag $\psi$ ($l=2$) | 2.15° | 3.37° |
| $\Delta R/R$ range (tide + oblateness) | [−0.070, +0.155] | [−0.154, +0.272] |
| Deviation from the sphere, mm | −0.61 … +1.35 | −1.39 … +2.46 |
| $\min r/\max r$ | 0.805 | 0.665 |
| $\Delta V/V$ (with oblateness) | +0.92% | +3.0% |
| Polar radius $r_{\rm pole}$ | 8.094 mm | 8.440 mm |

The line-of-centres gap between the deformed stars is **2.363 mm**, equal to the minimum 3-D vertex
distance, against 6.16 mm between ideal spheres. The light star is the more deformed one (tide
$\propto M'/M\cdot(R/r)^3$), and the one that limits the spin (`neutron_star_shape.md` §8).

### 6.3 Mesh and printing

| | |
|---|---|
| Base mesh | icosphere, subdivision level 6 (81,920 faces before the socket cut; 81,478 heavy / 81,526 light after) |
| Facet edge at $R\approx9$ mm | ~0.17 mm — well below any print resolution |
| Process | **SLA/resin strongly preferred** at this size. FDM acceptable at 0.1 mm layers / 0.4 mm nozzle, but layer lines will read across the smooth bulges. |
| Orientation | free — SLA, one piece, on supports (not socket-down: there is no bed flat). Put support contacts on the underside, the half that faces the plate in the finished piece, and keep the bore clear of supports. |

Each star is a single print; splitting it into glued hemispheres would save little on SLA and would put
a glue seam on the widest circumference. The icosphere is displaced radially by $r(\theta,\phi)$, then the
post socket is a boolean difference. Use `trimesh` with the `manifold3d` backend — hand-editing the mesh
around the socket will break watertightness.

### 6.4 Post socket

- **Axis:** world $-\hat z$ (straight down), entering the underside of the star at the pole on the spin
  axis. The pole radius is the actual deformed surface (`star.polar_radius`), including the spin
  oblateness and the tidal compression along $z$: 7% below $R$ for the heavy star.
- **Socket:** **Ø3.375 mm × 6.0 mm deep** (`POST_DIAMETER_MM + SOCKET_CLEARANCE_MM`: 0.2 mm clearance
  on the 1/8″ = 3.175 mm rod, which also covers the epoxy gap).
- **Depth is measured from the pole point** $z=-r_{\rm pole}$. The bore removes that point, so the
  lowest material is the bore rim, within about 0.1 mm of it.
- **The socket stays in the lower half of the star; the rod never reaches the equator.** The blind end
  sits 2.09 mm (heavy) / 2.44 mm (light) below the equatorial plane. `star.build_star` refuses a socket
  that would reach the equator.
- **No flat.** The star is not seated on a post shoulder (the rod is a plain cylinder); its height is
  set by the rod bottoming out in the socket and by the fixture while the epoxy cures. The bore opens
  onto the curved, nearly flat pole. Test-fit a printed socket on a real rod before printing both
  stars: resin holes often print 0.1 mm undersize.

---

## 7. Posts and assembly

| | |
|---|---|
| Material | **1/8″ (3.175 mm) brass rod** |
| Panel holes | Ø3.175 mm (a 1/8″ bit) × 8.0 mm deep, at $(-10.243, 0)$ and $(+13.658, 0)$; a CAM operation, not in the plate STL |
| Hole location | at each star's $(x,y)$ panel position, drilled from the *local* relief surface (there is no flat zone under the stars, §4.3) |
| Star centre height | $z=+3.85$ mm, both stars level: derived by `build_all.py` as the smallest multiple of 0.05 mm that makes the closest approach of either star to the plate $\ge$ **5 mm** (`posts.star_min_gap_mm`) |

**What "5 mm" means.** The height is the *closest approach between a star's surface and the plate
surface* (3-D, measured on the placed meshes against the plate grid, `validate.check_star_height`), not
the star-centre height and not the clearance above the $z=0$ datum. The plate is carved, so those
differ: the relief under the stars is ≈ −10 mm (−9.45 / −10.33), not 0. Heavy star 5.05 mm, light star
5.49 mm at the derived centre height of 3.85 mm (the heavy star governs; a level pair means the light one cannot
also be exactly 5.0 — it sits 0.5 mm higher above its floor). For comparison, the same stars with their
*centre* at 5 mm would sit at 6.2 / 6.6 mm. The stars' centres are below the high ridges of the relief
($z=0$), a deliberate consequence of the deep valley they sit in.

The relief at the two hole positions is read from the 0.5 mm plate grid (interpolated; `plate.relief_at`):

| | local relief $z$ at hole | panel hole bottom (8 mm deeper) | rod cut length | material below hole (36 mm stock) | rod visible, relief → star pole |
|---|---|---|---|---|---|
| heavy star | −9.45 mm | −17.45 mm | **19.2 mm** | 18.5 mm | 5.2 mm |
| light star | −10.33 mm | −18.33 mm | **19.7 mm** | 17.7 mm | 5.7 mm |

The two rods differ by 0.5 mm (the stars are level, but the light star's floor is 0.9 mm deeper), so
cut each to its row. The rod length is the 8 mm in the panel + the visible length + the 6 mm in the
socket; the visible length is *pole above relief*, not centre height. Bond with slow-cure epoxy; the
sockets are the only alignment feature, so **fixture the bulge axis along $\pm\hat x$ while curing** — a
rotated star is the most likely assembly error and the most visible. The tidal-lag rotation of each
star's long axis ($\psi\approx2.2^\circ$ heavy, $3.4^\circ$ light) is already built into the mesh, so
the fixture is still the $\pm\hat x$ line of centres.

**Where the numbers come from.** `plate.relief_at` interpolates the plate grid. Do not call
`wavefield.height_field` on a single point: it normalizes by the maximum over the points it is given,
so a lone point clamps to exactly $0$ or $-A_{\rm pp}$.

---

## 8. Deliverables

All STL exports: **binary format, millimetre units, watertight, outward normals, no
self-intersections.** `validate.py` verifies `is_watertight` and `is_winding_consistent` before they are
written.

| File | Contents | Approx. size |
|---|---|---|
| `ns_heavy_1p6.stl` | 1.6 M☉ star, body frame, socket at $-z$, centred at origin | ~4 MB |
| `ns_light_1p2.stl` | 1.2 M☉ star, same convention | ~4 MB |
| `ns_heavy_1p6_sphere.stl`, `ns_light_1p2_sphere.stl` | **Comparison stars, no deformation**: ideal spheres of radius $R_0$ (8.691 / 9.053 mm), no tide, no spin, no lag; same Ø3.375 × 6.0 mm socket from the pole point, no flat. Not placed in `assembly_preview.stl`, and not covered by `check_star_height` — see below. | ~4 MB each |
| `plate_relief.stl` | **closed solid**: relief top + 4 side walls + flat bottom, 254 × 254 mm, bottom at $z=-36$ mm. **No post holes** — a separate CAM operation (§7). | ~3 MB (2 mm grid) |
| `plate_heightmap.png` | 16-bit greyscale, 128×128 at the 2 mm grid, with a `.json` sidecar giving the $z$-range — alternate input for CAM that prefers depth maps | ~30 kB |
| `assembly_preview.stl` | everything positioned, posts included — **render only, not for manufacture** | ~11 MB |
| `star_modes.toml` | the generated per-star mode amplitudes the star meshes were built from | ~3 kB |

**Sphere stars on the same posts.** With the same centre height ($z=3.85$ mm) a sphere's pole is about
0.6 mm lower than the deformed star's (heavy 0.60, light 0.61), so their closest approach to the plate
is 4.51 / 4.89 mm, not 5. For $\ge5$ mm use centre $z=4.35$ (heavy) and $4.00$ (light) mm, giving rods of
19.1 / 19.3 mm; at 3.85 mm the rods are 18.6 / 19.1 mm. The sphere–sphere surface gap is 6.16 mm.

**Plate sampling grid: 2 mm** (`panel.grid_step_mm`; 128×128 → 32,258 top triangles, 65,532 in the closed
solid). This is coarser than the 0.75 mm finishing stepover; the surface stays within 0.03 mm of the field
(§5.1), and CAM interpolates linearly between samples. Set `grid_step_mm = 0.5` for a 1.04M-face, ~52 MB
plate that oversamples the stepover. The field and every check are independent of this setting (§4.3). The plate STL must be a closed solid, not an open surface —
Fusion, VCarve and most CAM will reject or mis-offset an open mesh.

---

## 9. Implementation

**Stack:** Python 3.11+ (uses `tomllib`), `numpy`, `scipy` ≥ 1.15, `trimesh`, `manifold3d`, `pillow`.

```
build_all.py              generate every deliverable into out/<config name>/ and run the validation
configs/default.toml      the reference design; every key a user can change, with its meaning
configs/equal_mass.toml   example alternative: 1.4 + 1.4 Msun, nonspinning
src/params.py             reads the configuration (--config / TIDAL_NS_CONFIG) into constants
src/cli_config.py         selects the configuration file from the command line
src/binary.py             Kepler solve, Omega, lambda_GW, scale S, star positions
src/wavefield.py          retarded per-star potential, height field z(x, y)
src/plate.py              sample grid -> watertight solid -> STL; heightmap PNG
src/tidal_waveform.py     f-mode amplitudes b_a (l=2, l=3) along the inspiral (Yu et al. 2022, 2024)
src/polytrope.py          Gamma=2 polytrope: mode overlaps I_l, Love numbers k_l, shape numbers h_l
src/tidal_shape.py        amplitudes -> surface radius r(theta, phi); regime gate (V11); TOML loader
src/star.py               icosphere, radial displacement, socket boolean -> STL
src/assembly.py           preview scene
src/validate.py           post-generation checks
tools/make_star_modes.py  amplitudes at the snapshot -> <out>/star_modes.toml (build_all.py calls it)
data/<config>/star_modes.toml  reference amplitudes of the default and equal-mass configurations, for the tests and verification
tests/                    pytest suite
verification/             independent re-derivations and a time-domain check of the amplitudes
```

`validate.py` hard-fails on: a non-watertight mesh; $\max\lvert\nabla z\rvert$ above `MAX_SLOPE_DEG`
(45°, a margin above genuine tool capability, not a milling limit itself); minimum concave
$R_c<2R_b$; either panel post hole leaving less than `MIN_HOLE_CLEARANCE_MM` (3.0 mm, §7) of wood below
it; either star closer to the plate than `STAR_MIN_GAP_MM` (5.0 mm, §7); a star–star surface gap
$\le0$; relief exceeding the depth budget; and a linear-regime error from the star shapes
(`neutron_star_shape.md` V11). These are the failure modes that only show up after the wood is cut.
A check is only as good as its inputs: it must read the relief from the full plate grid (§7).

---

## 10. Build notes and known limitations

1. **Test cut first.** Cut a 100 mm square of the outer field before committing 254 mm of maple. This
   shows whether the rim ($\pm2.5$ mm average edge swing, §4.3) reads clearly enough or whether
   `Z_EVAL_MM` should move toward −100 (§3.4).
2. **Panel edge treatment** (roundover, chamfer or square, and whether to bevel the underside) is
   left to the builder; the design leaves the top edge square (§5.2).
3. **Wood movement** across 254 mm of maple is ~1–2 mm seasonally. It does not affect the relief, but
   the post holes shift slightly relative to each other. Bond the posts after the panel has acclimated.
4. **`Z_EVAL_MM` is tuned by eye** (§3.4). If $A_{\rm pp}$, the plate size or `C_EFF_FACTOR` change,
   re-render a few candidate values and re-check `validate.py`.
5. **Newtonian orbit at $v/c=0.35$.** The Keplerian $\Omega(a)$ is ~10% off. A post-Newtonian correction
   would shift $f_{\rm GW}$, $\lambda_{\rm GW}$, $S$ and the star trajectories fed into the field model
   (both use the same `binary.solve()` orbit). Cosmetically irrelevant, but worth noting if the piece
   is ever labelled with numbers.
6. **Missing Liénard–Wiechert velocity factor.** The retarded potential (§3.2) uses the delayed
   *position* of each star but not the $1/(1-\hat n\cdot\mathbf v/c)$ beaming factor of an exact retarded
   scalar-wave solution. At $v/c=0.35$ this factor differs from 1 by up to ~35% depending on viewing
   angle, which would make the two trailing arms visibly unequal in strength. Deferred as a refinement.
7. **Amplitude normalization is not GR-calibrated** (§3.3). Do not label the piece with a claimed strain
   value.
8. **Star shapes are first order.** The linear surface is an approximation at the printed amplitudes
   (the light star's tip is +0.27 R; second-order terms are $\sim0.27^2$), and the stellar model is a
   Newtonian $\Gamma=2$ polytrope. `neutron_star_shape.md` §8 states the limits.
