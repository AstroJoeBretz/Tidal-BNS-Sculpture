"""Configuration loading, regeneration of the reference amplitudes, and the derived star height."""
import sys
import tomllib
from pathlib import Path

import numpy as np
import pytest
import trimesh

import assembly as ASM
import params as P

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "tools"))


def write(tmp_path, text):
    p = tmp_path / "c.toml"
    p.write_text(text, encoding="utf-8")
    return p


def test_default_config_file_matches_the_built_in_defaults():
    cfg = P.load_config(ROOT / "configs" / "default.toml")
    for k, v in P._DEFAULTS.items():
        assert cfg[k] == pytest.approx(v, rel=0, abs=0), k


def test_partial_config_fills_defaults(tmp_path):
    cfg = P.load_config(write(tmp_path, "[binary]\na_km = 35.0\n"))
    assert cfg["binary.a_km"] == 35.0 and cfg["binary.m1_msun"] == P._DEFAULTS["binary.m1_msun"]


@pytest.mark.parametrize("text, match", [
    ("[binary]\nm1_mass = 1.6\n", "unknown"),
    ("[binary]\nm1_msun = 1.0\nm2_msun = 1.2\n", ">="),
    ("[panel]\ngrid_step_mm = \"fine\"\n", "must be"),
])
def test_bad_config_is_rejected(tmp_path, text, match):
    with pytest.raises(ValueError, match=match):
        P.load_config(write(tmp_path, text))


def test_regenerated_amplitudes_match_the_committed_reference(tmp_path):
    import make_star_modes as MSM
    out = tmp_path / "modes.toml"
    assert MSM.generate(out, verbose=False) == 0
    new = tomllib.loads(out.read_text())
    ref = tomllib.load(open(ROOT / "data" / "star_modes.toml", "rb"))
    for sn, sr in zip(new["star"], ref["star"]):
        assert sn["name"] == sr["name"] and sn["spin_hz"] == sr["spin_hz"]
        for mn, mr in zip(sn["mode"], sr["mode"]):
            assert mn["C_re"] == pytest.approx(mr["C_re"], rel=1e-9, abs=1e-15)
            assert mn["C_im"] == pytest.approx(mr["C_im"], rel=1e-9, abs=1e-15)


def test_star_center_height_over_a_flat_plate():
    xs = np.linspace(-40, 40, 161)
    X, Y = np.meshgrid(xs, xs)
    Z = np.full_like(X, -10.0)
    ball = trimesh.creation.icosphere(subdivisions=4, radius=8.0)
    # lowest point 8 mm below the centre; the 5 mm gap puts the centre at -10 + 8 + 5 = 3 mm
    z = ASM.star_center_z(X, Y, Z, [(ball, 0.0)], gap=5.0)
    assert 3.0 <= z <= 3.1
