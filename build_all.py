"""Generate all sculpture deliverables into out/. Run: python build_all.py

See docs/design.md sec 8-9 for the file list and module layout.
"""
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0,str(Path(__file__).parent / "src"))

import assembly as ASM     # noqa: E402
import binary               # noqa: E402
import params as P          # noqa: E402
import plate as PLATE       # noqa: E402
import star as STAR         # noqa: E402
import tidal_shape as TS    # noqa: E402
import validate as VALIDATE  # noqa: E402


def main():
    out = Path(__file__).parent / "out"
    out.mkdir(exist_ok=True)

    b = binary.solve()
    print(f"S = {b.S_mm_per_km:.5f} mm/km | a = {b.a_mm:.3f} mm | "
          f"R1 = {b.R1_mm:.3f} mm | R2 = {b.R2_mm:.3f} mm | "
          f"x1 = {b.x1_mm:.3f} mm | x2 = {b.x2_mm:.3f} mm")
    print(f"f_orb = {b.f_orb:.1f} Hz | f_GW = {b.f_gw:.1f} Hz | "
          f"lambda_GW = {b.lambda_gw_km:.2f} km\n")

    # --- plate ---
    X, Y, xs, ys = PLATE.sample_grid()
    Z = PLATE.height_grid(b, X, Y)
    plate_mesh = PLATE.grid_to_solid(X, Y, Z, bottom_z=-P.STOCK_THICKNESS_MM)
    plate_mesh.export(out / "plate_relief.stl")
    PLATE.write_heightmap(Z, out / "plate_heightmap.png", out / "plate_heightmap.json")
    print(f"plate_relief.stl: {len(plate_mesh.faces)} faces, "
          f"watertight={plate_mesh.is_watertight}\n")

    # --- stars: shape from the f-mode amplitudes (docs/neutron_star_shape.md) ---
    modes_file = Path(__file__).parent / "data" / "star_modes.toml"
    if not modes_file.exists():
        raise SystemExit(f"{modes_file} not found: run  python tools/make_star_modes.py  first")
    models = TS.load_star_modes(modes_file)
    heavy_model, light_model = models["heavy"], models["light"]
    heavy_fn = heavy_model.radius_fn(b.R1_mm)
    light_fn = light_model.radius_fn(b.R2_mm)
    heavy = STAR.build_star(heavy_fn)
    light = STAR.build_star(light_fn)
    heavy.export(out / "ns_heavy_1p6.stl")
    light.export(out / "ns_light_1p2.stl")
    for name, model, mesh in (("ns_heavy_1p6.stl", heavy_model, heavy), ("ns_light_1p2.stl", light_model, light)):
        psi = TS.lag_angles(model.modes, model.X_lag).get((2, 2), 0.0)
        print(f"{name}: spin={model.spin_hz:.0f} Hz  lag psi={psi:.2f} deg  X={model.X:g} X_lag={model.X_lag:g} "
              f"oblateness={'on' if model.spin_oblateness else 'off'}  "
              f"faces={len(mesh.faces)} watertight={mesh.is_watertight}")
    print()

    # --- ideal-sphere comparison stars: radius R0, no tide, no spin, same socket (docs/design.md sec 8) ---
    for name, R_mm, fname in (("heavy", b.R1_mm, "ns_heavy_1p6_sphere.stl"),
                              ("light", b.R2_mm, "ns_light_1p2_sphere.stl")):
        sphere = STAR.build_star(lambda n, R=R_mm: np.full(len(n), R))
        VALIDATE.check_mesh(sphere, fname[:-4])
        sphere.export(out / fname)
    print()

    # --- assembly preview ---
    zp_h = -STAR.polar_radius(heavy_fn)
    zp_l = -STAR.polar_radius(light_fn)
    zl_h = PLATE.relief_at(X, Y, Z, b.x1_mm, 0.0)
    zl_l = PLATE.relief_at(X, Y, Z, b.x2_mm, 0.0)
    preview = ASM.build(plate_mesh, heavy, light, b, zp_h, zp_l, zl_h, zl_l)
    preview.export(out / "assembly_preview.stl")
    print(f"assembly_preview.stl: {len(preview.faces)} faces")
    for name, zp, zl in (("heavy", zp_h, zl_h), ("light", zp_l, zl_l)):
        z_bot, z_top = ASM.post_span(zp, zl)
        print(f"  {name} post (O{P.POST_DIAMETER_MM:.3f}): relief at hole z = {zl:.2f} mm, "
              f"rod {z_top - z_bot:.1f} mm long, {z_top - zl - P.SOCKET_DEPTH_MM:.1f} mm "
              f"visible from the relief to the star's pole")
    print()

    # --- validate ---
    VALIDATE.run(b, X, Y, Z, heavy, light, plate_mesh, heavy_model, light_model)


if __name__ == "__main__":
    main()
