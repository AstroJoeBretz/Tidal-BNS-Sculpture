"""Positioned preview of the full sculpture: plate + stars + posts.

Render only -- not for manufacture (docs/design.md sec 8). Simple concatenation, no
booleans between plate and stars; they don't intersect (the stars float
above the panel on posts).
"""
import trimesh

import binary
import params as P


def post_span(z_pole, z_local):
    """(z_bottom, z_top) of the rod for a star whose bottom pole sits at z_pole in the star
    frame (= -r_pole), in a panel hole drilled PANEL_HOLE_DEPTH_MM down from the local relief
    z_local. The rod bottoms out in the socket, SOCKET_DEPTH_MM above the pole point."""
    return z_local - P.PANEL_HOLE_DEPTH_MM, P.STAR_CENTER_Z_MM + z_pole + P.SOCKET_DEPTH_MM


def post_mesh(x_mm, z_pole, z_local):
    z_bottom, z_top = post_span(z_pole, z_local)
    height = z_top - z_bottom
    post = trimesh.creation.cylinder(radius=P.POST_DIAMETER_MM / 2, height=height)
    post.apply_translation([x_mm, 0, (z_top + z_bottom) / 2])
    return post


def build(plate_mesh, heavy_mesh, light_mesh, b: binary.Binary, z_pole_heavy, z_pole_light,
          z_local_heavy, z_local_light):
    heavy = heavy_mesh.copy()
    heavy.apply_translation([b.x1_mm, 0, P.STAR_CENTER_Z_MM])

    light = light_mesh.copy()
    light.apply_translation([b.x2_mm, 0, P.STAR_CENTER_Z_MM])

    post1 = post_mesh(b.x1_mm, z_pole_heavy, z_local_heavy)
    post2 = post_mesh(b.x2_mm, z_pole_light, z_local_light)

    return trimesh.util.concatenate([plate_mesh, heavy, light, post1, post2])
