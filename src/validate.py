"""Post-generation checks. Hard-fails on the mistakes that only show up after
you've cut the wood or printed the resin. docs/design.md sec 9."""
import numpy as np
from scipy.spatial import cKDTree

import params as P
import plate as PLATE
import tidal_shape as TS

# Not a tool limit -- single-pass 3-axis ball-nose finishing is fine well
# beyond this (docs/design.md sec 5.1 notes steep-wall/rest-machining only becomes a
# concern above ~45-60 deg). Set to 45 deg, the low end of that range; the current
# relief measures 36.9 deg. Steeper faces cost finish (scallop ~2x at 45 deg for
# the same stepover) and raise tearout risk across the grain -- that is what this
# guard is for.
# The genuine hard floor is the concave-curvature check below (ball radius
# vs surface curvature), not this one.
MAX_SLOPE_DEG = P.MAX_SLOPE_DEG

# Minimum wood left behind a panel post hole (docs/design.md sec 7).
# There is no flat zone under the stars in the current field model, so each
# panel hole is drilled from wherever the local relief actually is, read from
# the full plate grid (plate.relief_at). Slightly more conservative than the
# ~2.2 mm precedent already accepted for the star's own socket wall (sec
# 6.4), since this location sees a wider, thinner span of surrounding stock.
MIN_HOLE_CLEARANCE_MM = P.MIN_HOLE_CLEARANCE_MM


def check_mesh(mesh, name):
    ok = mesh.is_watertight and mesh.is_winding_consistent
    print(f"[{'PASS' if ok else 'FAIL'}] {name}: watertight={mesh.is_watertight} "
          f"winding_consistent={mesh.is_winding_consistent} faces={len(mesh.faces)}")
    if not ok:
        raise AssertionError(f"{name} is not a clean manifold solid")


def check_slope_and_curvature(X, Y, Z, ball_radius_mm=P.BALL_RADIUS_MM):
    dx = X[0, 1] - X[0, 0]
    dy = Y[1, 0] - Y[0, 0]
    dzdy, dzdx = np.gradient(Z, dy, dx)
    slope = np.hypot(dzdx, dzdy)
    max_slope_deg = np.degrees(np.arctan(slope.max()))

    d2zdx2 = np.gradient(dzdx, dx, axis=1)
    d2zdy2 = np.gradient(dzdy, dy, axis=0)
    curvature = np.maximum(np.abs(d2zdx2), np.abs(d2zdy2))
    min_Rc = 1.0 / max(curvature.max(), 1e-12)

    slope_ok = max_slope_deg <= MAX_SLOPE_DEG
    curv_ok = min_Rc >= 2 * ball_radius_mm

    print(f"[{'PASS' if slope_ok else 'FAIL'}] max slope = {max_slope_deg:.2f} deg "
          f"(limit {MAX_SLOPE_DEG:.1f} deg)")
    print(f"[{'PASS' if curv_ok else 'FAIL'}] min concave R_c = {min_Rc:.3f} mm "
          f"(need >= {2*ball_radius_mm:.3f} mm, margin {min_Rc/ball_radius_mm:.1f}x ball radius)")

    if not slope_ok:
        raise AssertionError(f"max slope {max_slope_deg:.2f} deg exceeds {MAX_SLOPE_DEG} deg")
    if not curv_ok:
        raise AssertionError(f"min concave R_c {min_Rc:.3f} mm below 2*ball_radius")


def check_depth_budget(Z, tol=1e-6):
    z_min, z_max = Z.min(), Z.max()
    ok = (z_min >= -P.A_PP_MM - tol) and (z_max <= 0.0 + tol)
    print(f"[{'PASS' if ok else 'FAIL'}] relief depth range = [{z_min:.4f}, {z_max:.4f}] mm "
          f"(budget [-{P.A_PP_MM:.1f}, 0.0] mm)")
    if not ok:
        raise AssertionError("relief exceeds the depth budget")


def check_post_hole_clearance(b, X, Y, Z):
    """Material remaining below each panel post hole, drilled from the local
    relief surface (not a flat reference -- docs/design.md sec 7)."""
    ok = True
    for name, x_mm in [("heavy", b.x1_mm), ("light", b.x2_mm)]:
        z_local = PLATE.relief_at(X, Y, Z, x_mm, 0.0)
        hole_bottom_z = z_local - P.PANEL_HOLE_DEPTH_MM
        material = P.STOCK_THICKNESS_MM + hole_bottom_z
        this_ok = material >= MIN_HOLE_CLEARANCE_MM
        ok = ok and this_ok
        print(f"[{'PASS' if this_ok else 'FAIL'}] {name} star post hole: local relief z = "
              f"{z_local:.3f} mm, material below hole = {material:.3f} mm "
              f"(need >= {MIN_HOLE_CLEARANCE_MM:.1f} mm)")
    if not ok:
        raise AssertionError("a panel post hole leaves too little material below it "
                              "-- reduce A_pp or increase stock thickness")


def check_star_height(b, X, Y, Z, heavy_mesh, light_mesh, z_center):
    """Closest approach of each placed star (flat, socket and all) to the plate surface must be
    at least STAR_MIN_GAP_MM (docs/design.md sec 7).

    Minimum vertex-to-vertex distance on the analysis grid (X, Y, Z at ANALYSIS_STEP_MM); with
    that 0.5 mm spacing and a ~5 mm gap this overestimates the surface distance by < 0.02 mm.
    """
    tree = cKDTree(np.column_stack([X.ravel(), Y.ravel(), Z.ravel()]))
    ok = True
    for name, mesh, x_mm in (("heavy", heavy_mesh, b.x1_mm), ("light", light_mesh, b.x2_mm)):
        placed = mesh.vertices + np.array([x_mm, 0.0, z_center])
        gap = float(tree.query(placed)[0].min())
        this_ok = gap >= P.STAR_MIN_GAP_MM
        ok = ok and this_ok
        print(f"[{'PASS' if this_ok else 'FAIL'}] {name} star height: centre z = {z_center:.2f} mm, "
              f"closest approach to plate = {gap:.3f} mm (need >= {P.STAR_MIN_GAP_MM:.1f} mm)")
    if not ok:
        raise AssertionError("a star is closer to the plate than STAR_MIN_GAP_MM "
                             "-- raise the star height")


def check_star_gap(b, heavy_mesh, light_mesh, heavy_fn, light_fn):
    """docs/neutron_star_shape.md V7: the deformed stars must not touch.

    Reports the line-of-centres gap (tips toward the companion) and the true 3-D gap, the
    minimum vertex-to-vertex distance between the two placed meshes (an upper bound on the
    surface distance; the mesh spacing is ~0.17 mm).
    """
    tip_h = float(heavy_fn(np.array([[1.0, 0.0, 0.0]]))[0])
    tip_l = float(light_fn(np.array([[-1.0, 0.0, 0.0]]))[0])
    gap_line = b.a_mm - tip_h - tip_l

    vh = heavy_mesh.vertices + np.array([b.x1_mm, 0.0, 0.0])
    vl = light_mesh.vertices + np.array([b.x2_mm, 0.0, 0.0])
    gap_3d, _ = cKDTree(vh).query(vl)
    gap_3d = float(gap_3d.min())

    ok = gap_line > 0 and gap_3d > 0
    print(f"[{'PASS' if ok else 'FAIL'}] star gap: line of centres = {gap_line:.3f} mm "
          f"(tip_heavy={tip_h:.3f} tip_light={tip_l:.3f} a={b.a_mm:.3f}), minimum 3-D = {gap_3d:.3f} mm")
    if not ok:
        raise AssertionError("tidally deformed stars overlap -- reduce the tidal amplitudes, "
                             "the exaggeration X, or the spin, or increase a_km")


def check_star_regimes(heavy_model, light_model):
    """docs/neutron_star_shape.md V11: linear-regime gate for each star; errors fail, warnings print."""
    for model in (heavy_model, light_model):
        rep = TS.check_regime(model, quiet=True)      # raises RegimeError for past-resonance modes
        margins = ", ".join(f"(l={r['l']},m={r['m']:+d}) {r['delta_over_omega']:.2f}"
                            for r in rep["modes"] if r["m"] > 0)
        status = "WARN" if rep["warnings"] else "PASS"
        print(f"[{status}] {model.name} linear regime: Delta/omega {margins}; "
              f"tide dR/R in [{rep['tide_range'][0]:+.3f}, {rep['tide_range'][1]:+.3f}], "
              f"dV/V = {rep['dV_tide']:+.4f} (tide) / {rep['dV_total']:+.4f} (with spin oblateness), "
              f"min r/max r = {rep['min_over_max_r']:.3f}")
        for msg in rep["warnings"]:
            print(f"        warning: {msg}")


def run(b, X, Y, Z, heavy_mesh, light_mesh, plate_mesh, heavy_model, light_model, z_center,
        heavy_name="ns_heavy", light_name="ns_light"):
    """X, Y, Z: the plate relief on the analysis grid (ANALYSIS_STEP_MM), not the export grid."""
    check_mesh(plate_mesh, "plate_relief")
    check_mesh(heavy_mesh, heavy_name)
    check_mesh(light_mesh, light_name)
    check_depth_budget(Z)
    check_slope_and_curvature(X, Y, Z)
    check_post_hole_clearance(b, X, Y, Z)
    check_star_height(b, X, Y, Z, heavy_mesh, light_mesh, z_center)
    check_star_regimes(heavy_model, light_model)
    check_star_gap(b, heavy_mesh, light_mesh,
                   heavy_model.radius_fn(b.R1_mm), light_model.radius_fn(b.R2_mm))
    print("\nAll checks passed.")
