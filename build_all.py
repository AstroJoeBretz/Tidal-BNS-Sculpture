"""Generate all sculpture deliverables for one configuration.

Run:  python build_all.py [--config configs/default.toml] [--out-root out]

Writes plate, stars, heightmap, preview and the star amplitudes into <out-root>/<config name>/,
then validates the design. See docs/design.md sec 8-9 for the file list and module layout.
"""
import argparse
import sys
from pathlib import Path

ROOT = Path(__file__).parent
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT / "tools"))

import cli_config  # noqa: E402

cli_config.select_from_argv()

import numpy as np                # noqa: E402

import assembly as ASM            # noqa: E402
import binary                     # noqa: E402
import make_star_modes as MSM     # noqa: E402
import params as P                # noqa: E402
import plate as PLATE             # noqa: E402
import star as STAR               # noqa: E402
import tidal_shape as TS          # noqa: E402
import validate as VALIDATE       # noqa: E402


def mass_tag(m_msun):
    """1.6 -> '1p6' (used in the star file names)."""
    return f"{m_msun:g}".replace(".", "p")


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--config", type=Path, default=None,
                    help="configuration TOML (default configs/default.toml)")
    ap.add_argument("--out-root", type=Path, default=ROOT / "out",
                    help="parent directory of the per-configuration output directory")
    args = ap.parse_args(argv)

    out = args.out_root / P.CONFIG_NAME
    out.mkdir(parents=True, exist_ok=True)
    heavy_name = f"ns_heavy_{mass_tag(P.M1_MSUN)}"
    light_name = f"ns_light_{mass_tag(P.M2_MSUN)}"
    print(f"configuration: {P.CONFIG_PATH}\noutput:        {out}\n")

    b = binary.solve()
    print(f"S = {b.S_mm_per_km:.5f} mm/km | a = {b.a_mm:.3f} mm | "
          f"R1 = {b.R1_mm:.3f} mm | R2 = {b.R2_mm:.3f} mm | "
          f"x1 = {b.x1_mm:.3f} mm | x2 = {b.x2_mm:.3f} mm")
    print(f"f_orb = {b.f_orb:.1f} Hz | f_GW = {b.f_gw:.1f} Hz | "
          f"lambda_GW = {b.lambda_gw_km:.2f} km\n")

    # --- plate: field normalized on the analysis grid, exported on GRID_STEP_MM ---
    N = PLATE.field_norm(b)
    Xa, Ya, _, _ = PLATE.sample_grid(P.ANALYSIS_STEP_MM)
    Za = PLATE.height_grid(b, Xa, Ya, norm=N)
    if P.GRID_STEP_MM == P.ANALYSIS_STEP_MM:
        X, Y, Z = Xa, Ya, Za
    else:
        X, Y, _, _ = PLATE.sample_grid(P.GRID_STEP_MM)
        Z = PLATE.height_grid(b, X, Y, norm=N)
    plate_mesh = PLATE.grid_to_solid(X, Y, Z, bottom_z=-P.STOCK_THICKNESS_MM)
    plate_mesh.export(out / "plate_relief.stl")
    PLATE.write_heightmap(Z, out / "plate_heightmap.png", out / "plate_heightmap.json")
    print(f"plate_relief.stl: {len(plate_mesh.faces)} faces ({P.GRID_STEP_MM:g} mm grid), "
          f"watertight={plate_mesh.is_watertight}\n")

    # --- stars: shape from the f-mode amplitudes (docs/neutron_star_shape.md) ---
    modes_file = out / "star_modes.toml"
    if MSM.generate(modes_file) != 0:
        raise SystemExit("star amplitudes failed the linear-regime gate; see above")
    print()
    models = TS.load_star_modes(modes_file)
    heavy_model, light_model = models["heavy"], models["light"]
    heavy_fn = heavy_model.radius_fn(b.R1_mm)
    light_fn = light_model.radius_fn(b.R2_mm)
    heavy = STAR.build_star(heavy_fn)
    light = STAR.build_star(light_fn)
    heavy.export(out / f"{heavy_name}.stl")
    light.export(out / f"{light_name}.stl")
    for name, model, mesh in ((heavy_name, heavy_model, heavy), (light_name, light_model, light)):
        psi = TS.lag_angles(model.modes, model.X_lag).get((2, 2), 0.0)
        print(f"{name}.stl: spin={model.spin_hz:.0f} Hz  lag psi={psi:.2f} deg  X={model.X:g} X_lag={model.X_lag:g} "
              f"oblateness={'on' if model.spin_oblateness else 'off'}  "
              f"faces={len(mesh.faces)} watertight={mesh.is_watertight}")
    print()

    # --- ideal-sphere comparison stars: radius R0, no tide, no spin, same socket (docs/design.md sec 8) ---
    for R_mm, fname in ((b.R1_mm, f"{heavy_name}_sphere"), (b.R2_mm, f"{light_name}_sphere")):
        sphere = STAR.build_star(lambda n, R=R_mm: np.full(len(n), R))
        VALIDATE.check_mesh(sphere, fname)
        sphere.export(out / f"{fname}.stl")
    print()

    # --- star height: the level centre height that gives STAR_MIN_GAP_MM over the relief ---
    z_center = ASM.star_center_z(Xa, Ya, Za, [(heavy, b.x1_mm), (light, b.x2_mm)])
    print(f"star centre height z = {z_center:.2f} mm (closest approach >= {P.STAR_MIN_GAP_MM:g} mm, "
          f"rounded up to {ASM.STAR_HEIGHT_QUANTUM_MM:g} mm)\n")

    # --- assembly preview ---
    zp_h = -STAR.polar_radius(heavy_fn)
    zp_l = -STAR.polar_radius(light_fn)
    zl_h = PLATE.relief_at(Xa, Ya, Za, b.x1_mm, 0.0)
    zl_l = PLATE.relief_at(Xa, Ya, Za, b.x2_mm, 0.0)
    preview = ASM.build(plate_mesh, heavy, light, b, zp_h, zp_l, zl_h, zl_l, z_center)
    preview.export(out / "assembly_preview.stl")
    print(f"assembly_preview.stl: {len(preview.faces)} faces")
    for name, zp, zl in (("heavy", zp_h, zl_h), ("light", zp_l, zl_l)):
        z_bot, z_top = ASM.post_span(zp, zl, z_center)
        print(f"  {name} post (O{P.POST_DIAMETER_MM:.3f}): relief at hole z = {zl:.2f} mm, "
              f"rod {z_top - z_bot:.1f} mm long, {z_top - zl - P.SOCKET_DEPTH_MM:.1f} mm "
              f"visible from the relief to the star's pole")
    print()

    # --- validate (milling checks on the analysis grid, not the export mesh) ---
    VALIDATE.run(b, Xa, Ya, Za, heavy, light, plate_mesh, heavy_model, light_model, z_center,
                 heavy_name=heavy_name, light_name=light_name)


if __name__ == "__main__":
    main()
