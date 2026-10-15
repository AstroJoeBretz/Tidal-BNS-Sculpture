# Tidal NS Sculpture

Generate the files for a gravitational-wave sculpture of a binary neutron star: a CNC-milled
hardwood panel carrying the retarded field of the orbit, with two 3D-printed, **tidally deformed**
neutron stars on brass posts above it.

![Relief heightmap of the panel](docs/images/plate_heightmap.png)

*The milling heightmap produced by `build_all.py` (light = high, dark = deep; 254 × 254 mm panel,
20 mm relief).*

The default design is a 1.6 + 1.2 M☉ binary at a 33 km separation, a $\lambda_{\rm GW}$ = 212 mm
scale, and stars that are shaped from the f-mode amplitudes of Yu et al.
([arXiv:2211.07002](https://arxiv.org/abs/2211.07002),
[arXiv:2404.00147](https://arxiv.org/abs/2404.00147)): triaxial, with the bulge lagging the
companion and, optionally, an anti-aligned spin. Everything is Python and runs in seconds.

## What you get

`python build_all.py` writes into `out/` (not committed; ~125 MB):

| File | Contents |
|---|---|
| `plate_relief.stl` | closed solid of the carved panel, 254 × 254 × 36 mm, for CAM (no post holes) |
| `plate_heightmap.png` + `.json` | the same relief as a 16-bit depth map, with its $z$-range |
| `ns_heavy_1p6.stl`, `ns_light_1p2.stl` | the two stars, watertight, with a Ø3.375 mm post socket |
| `ns_*_sphere.stl` | undeformed comparison stars of the same size |
| `assembly_preview.stl` | everything positioned, for rendering only |

It also runs a validation pass (milling slope and curvature, wood left under the post holes, star
height over the carved relief, star–star gap, linear-regime check) and fails loudly on a bad design.

## Quick start

Python 3.11 or newer.

```bash
pip install -r requirements.txt
python build_all.py        # ~20 s
python -m pytest           # ~10 s
```

`data/star_modes.toml` (the per-star tidal amplitudes) is committed, so building needs nothing else.

## Make your own

The system lives in [`src/params.py`](src/params.py): masses, radii, separation, panel size, relief
depth, milling and post parameters. To change the binary:

1. Edit `M1_MSUN`, `M2_MSUN`, `R1_KM`, `R2_KM`, `A_KM` in `src/params.py`.
2. Regenerate the star amplitudes (spins in Hz, negative = anti-aligned):

   ```bash
   python tools/make_star_modes.py --spin-heavy -500 --spin-light -200
   ```

   The generator refuses configurations outside the linear regime (for example a spin that puts the
   snapshot past a mode's resonance) and writes nothing.
3. `python build_all.py`, and read the validation output.

Things to know:

- The star height `STAR_CENTER_Z_MM` is a constant: it is the level height that gives the required
  5 mm clearance over the carved relief. If the validation reports the clearance is violated after a
  change, re-derive it ([`docs/design.md`](docs/design.md) §7).
- `Z_EVAL_MM`, `C_EFF_FACTOR` and `A_PP_MM` set the *look* of the panel and its milling margins; they
  interact (§3.4–§5.1 of the design document). Re-run the validation after changing any of them.
- Output file names say `1p6` and `1p2` for the default masses; they are not derived from `params.py`.
- Tidal exaggeration `X`, lag exaggeration `X_lag` and the spin-oblateness switch are per-star fields
  of `data/star_modes.toml`, written by `tools/make_star_modes.py` (`star_block`).

## How it works

```
params.py ──► binary.py  (orbit, scale S, star positions)
                 │
                 ├──► wavefield.py ──► plate.py ──► plate_relief.stl, heightmap
                 │      retarded per-star potential, carved as retarded − instantaneous
                 │
tools/make_star_modes.py ──► data/star_modes.toml ──► tidal_shape.py ──► star.py ──► ns_*.stl
   polytrope.py (stellar model)                          amplitudes → r(θ, φ)
   tidal_waveform.py (f-mode amplitudes)
                 │
                 └──► validate.py   checks everything against the milling and assembly constraints
```

## Documentation

- [`docs/design.md`](docs/design.md) — the system as built: field model, panel geometry, milling,
  star geometry, posts and assembly, deliverables, known limitations.
- [`docs/neutron_star_shape.md`](docs/neutron_star_shape.md) — specification of the deformed stellar
  surface, with its validation checks (V1–V12).
- [`docs/math.md`](docs/math.md) — derivation of every formula from the orbit to the surface, with
  a worked example for the sculpture's heavy star.

## Verification

- `python -m pytest` — the shape math, the amplitude pipeline, the regime gate and an independent
  float64 re-derivation of the worked example.
- `python verification/math_walkthrough.py` — prints the worked example of `docs/math.md` step by step.
- `python verification/verify_amplitudes.py` — compares the amplitude generator with a time-domain
  integration of the mode equation (a few minutes).

## Limitations

This is an artwork built on real physics, not a simulation of a binary neutron star merger:

- The panel is not a strain map. The field is a retarded Newtonian potential difference, evaluated
  off the orbital plane, with a wave speed of $0.35\,c$, and normalized to the relief depth.
- The stars are a linear, Newtonian $\Gamma=2$ polytrope response at one orbital frequency; the orbit
  is Keplerian even though $v/c\approx0.35$.
- Do not quote strain values from the piece. See `docs/design.md` §10 for the full list.

## Credits and licence

The amplitude model implements Yu, Weinberg, Arras, Kwon & Venumadhav (arXiv:2211.07002) and Yu,
Arras & Weinberg (arXiv:2404.00147); please cite those papers if you use the stellar shapes.

[MIT](LICENSE).
