"""Build a tidally-deformed neutron star mesh with post socket. docs/design.md sec 6.3-6.4.

The surface shape comes from tidal_shape.py (docs/neutron_star_shape.md): this module only turns a
*radius callable* r(n_world) [mm] into a watertight icosphere solid with a blind socket at
the bottom pole. The callable is expressed in the world frame, i.e. the companion-direction
mapping (docs/neutron_star_shape.md sec 5) is already applied -- see tidal_shape.StarModel.radius_fn.

World frame: +z up, the post socket is drilled along world -z at the star's bottom pole
(the spin axis), straight into the curved surface (no flat). It is printed on SLA supports in
any orientation; lock rotation about the vertical axis during assembly so the long in-plane
axis lands along the line of centres (docs/design.md sec 7); with a nonzero tidal lag that axis is
rotated by psi from the line of centres, and is already built into the mesh.

The socket depth is measured from the pole point z = -r_pole (the surface on the spin axis).
The bore removes that point: the lowest material is the bore rim, rho^2 / 2 r_pole ~ 0.2 mm
higher. The socket must stay in the lower half of the star (it never reaches the equator).
"""
import numpy as np
import trimesh

import params as P


def polar_radius(radius_fn):
    """Surface radius [mm] at the bottom pole, n = (0, 0, -1)."""
    return float(np.asarray(radius_fn(np.array([[0.0, 0.0, -1.0]])))[0])


def socket_bottom_top_z(r_pole):
    """World-frame z [mm] of the pole point and of the socket's blind end (both star-frame)."""
    return -r_pole, -r_pole + P.SOCKET_DEPTH_MM


def build_star(radius_fn, subdivisions=P.ICOSPHERE_SUBDIVISIONS, socket=True):
    """radius_fn: unit vectors (N, 3) in the world frame -> radius [mm] (N,).

    socket=False returns the bare displaced icosphere (used to compare mesh volume with
    tidal_shape.volume_drift).
    """
    sphere = trimesh.creation.icosphere(subdivisions=subdivisions, radius=1.0)
    n = sphere.vertices / np.linalg.norm(sphere.vertices, axis=1, keepdims=True)
    r = np.asarray(radius_fn(n), dtype=float)
    if np.min(r) <= 0:
        raise ValueError("star radius is not positive")
    sphere.vertices = n * r[:, None]
    if not socket:
        sphere.fix_normals()
        return sphere

    z_pole, z_top = socket_bottom_top_z(polar_radius(radius_fn))
    if z_top >= 0.0:
        raise ValueError(f"socket ({P.SOCKET_DEPTH_MM} mm from the pole) would reach the equator; "
                         "the post must stay in the lower half of the star")

    # drill the blind post socket, overshooting below the pole for a clean cut
    overshoot = 1.0
    z_bottom = z_pole - overshoot
    height = z_top - z_bottom
    cut = trimesh.creation.cylinder(radius=P.SOCKET_DIAMETER_MM / 2, height=height)
    cut.apply_translation([0, 0, (z_top + z_bottom) / 2])
    star = trimesh.boolean.difference([sphere, cut], engine="manifold")

    star.fix_normals()
    return star


if __name__ == "__main__":
    from pathlib import Path

    import binary
    import tidal_shape as TS

    b = binary.solve()
    models = TS.load_star_modes(Path(__file__).resolve().parent.parent / "data" / "star_modes.toml")
    for name, R_mm in (("heavy", b.R1_mm), ("light", b.R2_mm)):
        fn = models[name].radius_fn(R_mm)
        mesh = build_star(fn)
        print(f"{name}: R0={R_mm:.3f} mm  r_pole={polar_radius(fn):.3f} mm  "
              f"faces={len(mesh.faces)} watertight={mesh.is_watertight}")
