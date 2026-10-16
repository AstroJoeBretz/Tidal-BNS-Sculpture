"""Build the milled-panel solid and heightmap. docs/design.md sec 4, 8."""
import json

import numpy as np
import trimesh
from PIL import Image

import binary
import params as P
import wavefield as W


def sample_grid(step=P.GRID_STEP_MM):
    half = P.PLATE_SIZE_MM / 2
    n = int(round(P.PLATE_SIZE_MM / step)) + 1
    xs = np.linspace(-half, half, n)
    ys = np.linspace(-half, half, n)
    X, Y = np.meshgrid(xs, ys, indexing="xy")
    return X, Y, xs, ys


def field_norm(b: binary.Binary):
    """Normalization N of the carved field (docs/design.md sec 4.3): the maximum |u_raw| over the
    plate sampled at ANALYSIS_STEP_MM, independent of the export resolution GRID_STEP_MM."""
    X, Y, _, _ = sample_grid(P.ANALYSIS_STEP_MM)
    k_half = binary.retardation_k_half(b, P.C_EFF)
    return float(np.max(np.abs(W.delta_phi(X, Y, P.Z_EVAL_MM, b, k_half))) / P.C ** 2)


def height_grid(b: binary.Binary, X, Y, norm=None):
    """Relief z(x, y) in mm on the given grid, normalized by `norm` (default: field_norm(b))."""
    return W.height_field(X, Y, b, norm=field_norm(b) if norm is None else norm)


def _wall_faces(top_line, bottom_line, flip):
    """Two triangles per segment along a polyline, param increasing with j.

    With param increasing along the line and z(top) > z(bottom), the
    unflipped triangle order (a, b, d) / (a, d, c) yields outward normal
    z_hat x t_hat (derived by hand, sec 9 -- see grid_to_solid). Pass
    flip=True on the two sides where that direction is the inward one.
    """
    faces = []
    for j in range(len(top_line) - 1):
        a, b = top_line[j], top_line[j + 1]
        c, d = bottom_line[j], bottom_line[j + 1]
        if flip:
            faces.append((a, d, b))
            faces.append((a, c, d))
        else:
            faces.append((a, b, d))
            faces.append((a, d, c))
    return faces


def grid_to_solid(X, Y, Z, bottom_z):
    """Top relief surface + 4 side walls + flat bottom -> one watertight solid.

    Winding is fixed by hand (not via trimesh.fix_normals(), which is
    correct but ~O(minutes) on a ~1.4M-face mesh at the spec's 0.5 mm grid):

    Top surface (i0, i1, i2) / (i0, i2, i3) gives +z (outward, correct).
    Bottom is the same x,y layout at a flat z, so it needs the mirror-image
    winding, (j0, j2, j1) / (j0, j3, j2), to face -z (outward).
    Walls: with idx rows/cols built so the tangential parameter always
    increases with j, the unflipped order gives outward normal z_hat x
    t_hat. That is already outward on the +y edge (t_hat=+x -> +y) and the
    -x edge (t_hat=+y -> -x); it is inward -- needs flip=True -- on the
    -y edge (t_hat=+x -> +y, wrong sign) and the +x edge (t_hat=+y -> -x,
    wrong sign).
    """
    ny, nx = Z.shape
    top = np.column_stack([X.ravel(), Y.ravel(), Z.ravel()])
    idx = np.arange(nx * ny).reshape(ny, nx)
    n_top = nx * ny

    bottom = top.copy()
    bottom[:, 2] = bottom_z
    bidx = idx + n_top

    faces = []

    # top surface, two triangles per grid cell, wound for +z outward
    i0 = idx[:-1, :-1].ravel()
    i1 = idx[:-1, 1:].ravel()
    i2 = idx[1:, 1:].ravel()
    i3 = idx[1:, :-1].ravel()
    faces.append(np.column_stack([i0, i1, i2]))
    faces.append(np.column_stack([i0, i2, i3]))

    # side walls around the 4 edges of the grid
    faces += _wall_faces(idx[0, :], bidx[0, :], flip=True)    # y = -half
    faces += _wall_faces(idx[-1, :], bidx[-1, :], flip=False)  # y = +half
    faces += _wall_faces(idx[:, 0], bidx[:, 0], flip=False)    # x = -half
    faces += _wall_faces(idx[:, -1], bidx[:, -1], flip=True)   # x = +half

    # flat bottom, mirroring the top grid's connectivity (so every
    # wall-bottom edge is shared with a bottom face) but wound for -z outward
    j0 = bidx[:-1, :-1].ravel()
    j1 = bidx[:-1, 1:].ravel()
    j2 = bidx[1:, 1:].ravel()
    j3 = bidx[1:, :-1].ravel()
    faces.append(np.column_stack([j0, j2, j1]))
    faces.append(np.column_stack([j0, j3, j2]))

    faces = np.vstack(faces)
    vertices = np.vstack([top, bottom])

    return trimesh.Trimesh(vertices=vertices, faces=faces, process=True)


def relief_at(X, Y, Z, x_mm, y_mm):
    """Relief z [mm] at (x, y), bilinearly interpolated from the plate grid (Z on X, Y).

    wavefield.height_field cannot be used for this: it normalizes by the maximum of the
    points it is handed, so a lone point (or a small patch) comes back clamped to 0 or -A_pp.
    """
    from scipy.interpolate import RegularGridInterpolator

    f = RegularGridInterpolator((Y[:, 0], X[0, :]), Z)
    return float(f([[y_mm, x_mm]])[0])


def write_heightmap(Z, png_path, json_path):
    """16-bit greyscale PNG + JSON sidecar with the z-range. docs/design.md sec 8."""
    z_min, z_max = -P.A_PP_MM, 0.0
    norm = np.clip((Z - z_min) / (z_max - z_min), 0.0, 1.0)
    img16 = np.round(norm * 65535).astype(np.uint16)
    # row 0 = +y (north/top edge of the plate in a top-down view)
    img16 = np.flipud(img16)
    Image.fromarray(img16, mode="I;16").save(png_path)

    meta = {
        "z_min_mm": z_min,
        "z_max_mm": z_max,
        "plate_size_mm": P.PLATE_SIZE_MM,
        "grid_step_mm": P.GRID_STEP_MM,
        "width": img16.shape[1],
        "height": img16.shape[0],
        "row0": "y = +plate_size_mm/2 (north edge)",
        "col0": "x = -plate_size_mm/2 (west edge)",
        "value_to_z_mm": "z = z_min_mm + (value/65535)*(z_max_mm - z_min_mm)",
    }
    with open(json_path, "w") as f:
        json.dump(meta, f, indent=2)


if __name__ == "__main__":
    b = binary.solve()
    X, Y, xs, ys = sample_grid()
    Z = height_grid(b, X, Y)
    print(f"grid: {X.shape}, z range [{Z.min():.4f}, {Z.max():.4f}] mm")
    mesh = grid_to_solid(X, Y, Z, bottom_z=-P.STOCK_THICKNESS_MM)
    print(f"faces={len(mesh.faces)} watertight={mesh.is_watertight} "
          f"winding_consistent={mesh.is_winding_consistent}")
