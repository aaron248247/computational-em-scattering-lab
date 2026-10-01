import json
import numpy as np
import pytest
from em_rf_lab.config import ConfigError, SimulationConfig, load_config, save_config
from em_rf_lab import geometry as G


def test_default_config_valid():
    SimulationConfig()


@pytest.mark.parametrize("change", [
    dict(polarization="TEM"), dict(radius=-1.0), dict(eps_r=0.5), dict(resolution=2),
    dict(ka_min=3.0, ka_max=2.0), dict(n_freq=0), dict(pml=0.0), dict(pml=40.0),
    dict(n2f_half_width=0.5), dict(n2f_half_width=30.0), dict(cell_x=8.0, pml=3.0, n2f_half_width=2.5),
    dict(beam_waist=-1.0), dict(courant=0.9), dict(decay_tol=0.0), dict(farfield_distance=10.0),
    dict(min_run_time=5000.0), dict(n_angles=2), dict(decay_probe_radius=1.2), dict(decay_probe_radius=0.0),
])
def test_invalid_configs_rejected(change):
    with pytest.raises(ConfigError):
        SimulationConfig(**change)


def test_hash_deterministic_and_name_independent():
    a, b = SimulationConfig(name="a"), SimulationConfig(name="b")
    assert a.config_hash() == b.config_hash() == SimulationConfig().config_hash()
    assert a.config_hash() != a.replace(resolution=31.0).config_hash()


def test_json_roundtrip(tmp_path):
    cfg = SimulationConfig(name="rt", resolution=25.0, polarization="TE")
    save_config(cfg, tmp_path / "c.json")
    assert load_config(tmp_path / "c.json") == cfg


def test_unknown_key_rejected(tmp_path):
    (tmp_path / "c.json").write_text(json.dumps({"resolutionn": 3}))
    with pytest.raises(ConfigError):
        load_config(tmp_path / "c.json")


def test_ka_freq_roundtrip_and_values():
    assert G.ka_to_freq(2 * np.pi) == pytest.approx(1.0)
    ka = np.linspace(1, 6, 7)
    np.testing.assert_allclose(G.freq_to_ka(G.ka_to_freq(ka, 2.0), 2.0), ka)


def test_grids_deterministic_and_inclusive():
    cfg = SimulationConfig()
    ka = G.ka_grid(cfg)
    assert ka[0] == cfg.ka_min and ka[-1] == cfg.ka_max and len(ka) == cfg.n_freq
    np.testing.assert_array_equal(ka, G.ka_grid(cfg))
    assert np.any(np.isclose(ka, 3.0)) and np.any(np.isclose(ka, 4.5))


def test_angle_grid():
    a = G.angle_grid(91)
    assert a[0] == 0 and a[-1] == pytest.approx(np.pi) and len(a) == 91
    with pytest.raises(ValueError):
        G.angle_grid(2)


def test_points_per_wavelength_and_timestep():
    assert G.points_per_wavelength(20, 0.5, 2.0) == pytest.approx(20.0)
    assert G.courant_timestep(SimulationConfig(resolution=40.0)) == pytest.approx(0.5 / 40)
    with pytest.raises(ValueError):
        G.vacuum_wavelength(0.0)


def test_near2far_segments_closed_outward():
    segs = G.near2far_segments(2.5)
    assert len(segs) == 4
    for s in segs:
        c = np.array(s["center"])
        axis = 0 if s["direction"] == "x" else 1
        assert np.sign(c[axis]) == s["weight"]          # outward normal
        assert s["size"][axis] == 0 and s["size"][1 - axis] == 5.0
