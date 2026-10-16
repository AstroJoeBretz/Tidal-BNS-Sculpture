"""Parameters for the binary neutron star GW sculpture, read from a TOML configuration.

The configuration file is chosen by the environment variable TIDAL_NS_CONFIG (the command-line
tools set it from `--config`), defaulting to configs/default.toml. Every key is optional; omitted
keys take the defaults below, and unknown keys are an error. Names derived from the config:
CONFIG_PATH, CONFIG_NAME (the file stem, used for the output directory).

Single source of input constants. Everything derivable (scale factor, star positions, wavenumber,
normalization, tidal amplitudes, star height) is computed in binary.py / wavefield.py / star.py /
build_all.py, not repeated here. See docs/design.md.
"""
import os
import tomllib
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parent.parent

# --- physical constants (SI) ---
G = 6.67430e-11
C = 2.99792458e8
M_SUN = 1.98847e30

# --- configuration keys and defaults (see configs/default.toml for the meaning of each) ---
_DEFAULTS = {
    "binary.m1_msun": 1.6, "binary.m2_msun": 1.2, "binary.r1_km": 12.0, "binary.r2_km": 12.5,
    "binary.a_km": 33.0,
    "stars.tidal_exaggeration": 1.0, "stars.lag_exaggeration": 1.0, "stars.spin_oblateness": True,
    "stars.heavy.spin_hz": -500.0, "stars.light.spin_hz": -200.0,
    "field.z_eval_mm": -70.0, "field.c_eff_factor": 0.35, "field.retarded_iters": 10,
    "panel.plate_size_mm": 254.0, "panel.stock_thickness_mm": 36.0, "panel.relief_pp_mm": 20.0,
    "panel.grid_step_mm": 2.0, "panel.scale_lambda_mm": 150.0 * np.sqrt(2),
    "milling.ball_radius_mm": 3.175 / 2, "milling.stepover_mm": 0.75, "milling.max_slope_deg": 45.0,
    "posts.post_diameter_mm": 3.175, "posts.socket_clearance_mm": 0.2, "posts.socket_depth_mm": 6.0,
    "posts.panel_hole_depth_mm": 8.0, "posts.star_min_gap_mm": 5.0, "posts.min_hole_clearance_mm": 3.0,
    "mesh.icosphere_subdivisions": 6,
}


def _flatten(d, prefix=""):
    for k, v in d.items():
        if isinstance(v, dict):
            yield from _flatten(v, f"{prefix}{k}.")
        else:
            yield f"{prefix}{k}", v


def load_config(path):
    """Read a configuration TOML into a flat {"section.key": value} dict with defaults filled in."""
    with open(path, "rb") as f:
        given = dict(_flatten(tomllib.load(f)))
    unknown = sorted(set(given) - set(_DEFAULTS))
    if unknown:
        raise ValueError(f"{path}: unknown configuration key(s): {', '.join(unknown)}")
    cfg = dict(_DEFAULTS)
    for k, v in given.items():
        want = type(_DEFAULTS[k])
        if want is float and isinstance(v, int) and not isinstance(v, bool):
            v = float(v)
        if want is np.float64:
            want, v = float, float(v)
        if not isinstance(v, want):
            raise ValueError(f"{path}: {k} must be {want.__name__}, got {v!r}")
        cfg[k] = v
    if cfg["binary.m1_msun"] < cfg["binary.m2_msun"]:
        raise ValueError(f"{path}: binary.m1_msun must be >= binary.m2_msun (star 1 is the heavier one)")
    return cfg


CONFIG_PATH = Path(os.environ.get("TIDAL_NS_CONFIG", ROOT / "configs" / "default.toml")).resolve()
CONFIG_NAME = CONFIG_PATH.stem
_cfg = load_config(CONFIG_PATH)

# --- binary system (docs/design.md sec 1.1) ---
M1_MSUN = _cfg["binary.m1_msun"]   # heavier star
M2_MSUN = _cfg["binary.m2_msun"]   # lighter star
R1_KM = _cfg["binary.r1_km"]       # heavier star radius
R2_KM = _cfg["binary.r2_km"]       # lighter star radius
A_KM = _cfg["binary.a_km"]         # orbital separation (docs/design.md sec 1.3)

# --- tidal deformation (docs/neutron_star_shape.md) ---
# These are inputs to tools/make_star_modes.py, which writes the per-star amplitudes to
# <out>/star_modes.toml.
STAR_X = _cfg["stars.tidal_exaggeration"]
STAR_X_LAG = _cfg["stars.lag_exaggeration"]
SPIN_OBLATENESS = _cfg["stars.spin_oblateness"]
SPIN_HEAVY_HZ = _cfg["stars.heavy.spin_hz"]
SPIN_LIGHT_HZ = _cfg["stars.light.spin_hz"]

# --- panel (docs/design.md sec 4) ---
PLATE_SIZE_MM = _cfg["panel.plate_size_mm"]
STOCK_THICKNESS_MM = _cfg["panel.stock_thickness_mm"]
# Master scale reference (sec 1.2): one lambda_GW is this many mm on the finished piece, which
# fixes S = SCALE_LAMBDA_MM / lambda_GW[km]. It is independent of the plate size, so the stars,
# separation and every printed part keep their size when the panel is resized.
SCALE_LAMBDA_MM = _cfg["panel.scale_lambda_mm"]
# Relief peak-to-valley, sec 4.2. Slope and curvature scale linearly with it. Re-verify (the
# build's validation does) whenever A_pp, Z_EVAL_MM, C_EFF_FACTOR or STOCK_THICKNESS_MM change.
A_PP_MM = _cfg["panel.relief_pp_mm"]
GRID_STEP_MM = _cfg["panel.grid_step_mm"]   # sampling of the exported plate STL and heightmap
ANALYSIS_STEP_MM = 0.5                      # sampling used for the field normalization and the milling checks;
                                            # independent of GRID_STEP_MM, so the carved field does not
                                            # depend on the export resolution (sec 4.3)

# The field is carved all the way to the panel's cut edge: no edge window or border taper.

# --- field model: retarded per-star Newtonian potential (docs/design.md sec 3) ---
# Z_EVAL_MM: the field is evaluated off the orbital/panel plane, at (x, y, Z_EVAL_MM). This acts
# as a softening length that bounds the minimum source distance and so controls the dynamic
# range (sec 3.4). Chosen by inspection, not derived.
Z_EVAL_MM = _cfg["field.z_eval_mm"]
RETARDED_ITERS = _cfg["field.retarded_iters"]   # fixed-point iterations for the retarded-time solve;
                                                 # contraction factor is ~(k/2)*d_star ~ 0.18/iteration

# Hypothetical wave speed, applied only to each star's own retardation phase
# (binary.retardation_k_half) -- NOT to the plate/star geometry, which stays fixed at whatever
# scale the real C produced (sec 1.2). Smaller values wind the pattern more tightly.
# Set to 1.0 to recover the real speed of light.
C_EFF_FACTOR = _cfg["field.c_eff_factor"]
C_EFF = C * C_EFF_FACTOR

# --- milling (docs/design.md sec 5) ---
BALL_RADIUS_MM = _cfg["milling.ball_radius_mm"]
STEPOVER_MM = _cfg["milling.stepover_mm"]
MAX_SLOPE_DEG = _cfg["milling.max_slope_deg"]

# --- posts / assembly (docs/design.md sec 7) ---
POST_DIAMETER_MM = _cfg["posts.post_diameter_mm"]
SOCKET_CLEARANCE_MM = _cfg["posts.socket_clearance_mm"]   # diametral clearance + epoxy gap in the star
SOCKET_DIAMETER_MM = POST_DIAMETER_MM + SOCKET_CLEARANCE_MM
SOCKET_DEPTH_MM = _cfg["posts.socket_depth_mm"]           # from the pole point; lower half only (star.build_star)
PANEL_HOLE_DEPTH_MM = _cfg["posts.panel_hole_depth_mm"]   # drilled from the local relief surface, not z = 0
# Panel hole diameter = POST_DIAMETER_MM (a 1/8 in bit); it is a CAM operation, not modelled in the plate STL.
STAR_MIN_GAP_MM = _cfg["posts.star_min_gap_mm"]           # required closest approach of each star's surface
                                                           # to the plate; the level star centre height that
                                                           # delivers it is derived by build_all (sec 7)
MIN_HOLE_CLEARANCE_MM = _cfg["posts.min_hole_clearance_mm"]

# --- star mesh (docs/design.md sec 6.3) ---
ICOSPHERE_SUBDIVISIONS = _cfg["mesh.icosphere_subdivisions"]
