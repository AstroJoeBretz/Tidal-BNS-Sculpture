"""Positioned preview of the full sculpture: plate + stars + posts.

Render only -- not for manufacture (docs/design.md sec 8). Simple concatenation, no
booleans between plate and stars; they don't intersect (the stars float
above the panel on posts).
"""
import numpy as np
import trimesh
from scipy.optimize import brentq
from scipy.spatial import cKDTree

import binary
import params as P

STAR_HEIGHT_QUANTUM_MM = 0.05


def star_center_z(X, Y, Z, placements, quantum=STAR_HEIGHT_QUANTUM_MM, gap=None):
    """Level star-centre height [mm] that puts every star's closest approach to the plate at
    least `gap` (default STAR_MIN_GAP_MM) above the relief, rounded up to `quantum`.

    placements: [(mesh, x_mm), ...]. X, Y, Z: the relief on the analysis grid. The closest
    approach is the minimum vertex-to-vertex distance (an overestimate of the surface distance
    by < 0.02 mm on the 0.5 mm grid); it is monotone in the height, so each star is a bracketed
    root; the tallest requirement governs (docs/design.md sec 7).
    """
    gap = P.STAR_MIN_GAP_MM if gap is None else gap
    pts = np.column_stack([X.ravel(), Y.ravel(), Z.ravel()])
    needed = []
    for mesh, x_mm in placements:
        r_max = float(np.linalg.norm(mesh.vertices, axis=1).max())
        near = np.abs(pts[:, 0] - x_mm) < r_max + gap + 3.0
        near &= np.abs(pts[:, 1]) < r_max + gap + 3.0
        tree = cKDTree(pts[near])
        lower = mesh.vertices[mesh.vertices[:, 2] < 0.0]      # the plate is always below the star
        verts = lower + np.array([x_mm, 0.0, 0.0])

        def f(zc, verts=verts, tree=tree):
            return float(tree.query(verts + np.array([0.0, 0.0, zc]))[0].min()) - gap

        # walk down from a height at which the star clears the whole relief, to the first height
        # that violates the gap, then refine; above the plate the gap is monotone in the height
        hi = float(pts[near, 2].max()) + r_max + gap + 1.0
        lo = hi - 1.0
        while f(lo) >= 0.0:
            hi, lo = lo, lo - 1.0
            if lo < pts[near, 2].min() - r_max - gap:
                raise ValueError("no star height satisfies the clearance over the plate")
        needed.append(brentq(f, lo, hi, xtol=1e-5))
    z = max(needed)
    return float(np.ceil(z / quantum - 1e-9) * quantum)


def post_span(z_pole, z_local, z_center):
    """(z_bottom, z_top) of the rod for a star centred at z_center whose bottom pole sits at
    z_pole in the star frame (= -r_pole), in a panel hole drilled PANEL_HOLE_DEPTH_MM down from
    the local relief z_local. The rod bottoms out in the socket, SOCKET_DEPTH_MM above the
    pole point."""
    return z_local - P.PANEL_HOLE_DEPTH_MM, z_center + z_pole + P.SOCKET_DEPTH_MM


def post_mesh(x_mm, z_pole, z_local, z_center):
    z_bottom, z_top = post_span(z_pole, z_local, z_center)
    height = z_top - z_bottom
    post = trimesh.creation.cylinder(radius=P.POST_DIAMETER_MM / 2, height=height)
    post.apply_translation([x_mm, 0, (z_top + z_bottom) / 2])
    return post


def build(plate_mesh, heavy_mesh, light_mesh, b: binary.Binary, z_pole_heavy, z_pole_light,
          z_local_heavy, z_local_light, z_center):
    heavy = heavy_mesh.copy()
    heavy.apply_translation([b.x1_mm, 0, z_center])

    light = light_mesh.copy()
    light.apply_translation([b.x2_mm, 0, z_center])

    post1 = post_mesh(b.x1_mm, z_pole_heavy, z_local_heavy, z_center)
    post2 = post_mesh(b.x2_mm, z_pole_light, z_local_light, z_center)

    return trimesh.util.concatenate([plate_mesh, heavy, light, post1, post2])
