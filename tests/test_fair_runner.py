import importlib.util
import os
import sys
import types

import pytest

_fair_dir = os.path.join(
    os.path.dirname(__file__), "..", "dynamic_characterization", "fair"
)
_dc = sys.modules.setdefault(
    "dynamic_characterization", types.ModuleType("dynamic_characterization")
)
_pkg = types.ModuleType("dynamic_characterization.fair")
_pkg.__path__ = [_fair_dir]
_pkg.FAIR_IMPORT_ERROR_MSG = "install dynamic_characterization[fair]"
sys.modules["dynamic_characterization.fair"] = _pkg
_spec = importlib.util.spec_from_file_location(
    "dynamic_characterization.fair.runner",
    os.path.join(_fair_dir, "runner.py"),
)
runner = importlib.util.module_from_spec(_spec)
sys.modules["dynamic_characterization.fair.runner"] = runner
_spec.loader.exec_module(runner)


def test_require_fair_raises_when_absent():
    if importlib.util.find_spec("fair") is not None:
        pytest.skip("fair installed")
    with pytest.raises(ImportError, match=r"dynamic_characterization\[fair\]"):
        runner.require_fair()


def test_find_calibration_files_returns_tuple():
    pytest.importorskip("fair")
    result = runner.find_calibration_files()
    # Always a 2-tuple; each entry is a path or None depending on the install.
    assert isinstance(result, tuple)
    assert len(result) == 2
    for entry in result:
        assert entry is None or isinstance(entry, str)


def test_get_calibration_files_returns_existing_paths():
    pytest.importorskip("fair")
    # Downloads the calibration 1.4.1 pair on first use (cached afterwards).
    params, props = runner.get_calibration_files()
    assert os.path.exists(params)
    assert os.path.exists(props)


def test_unit_factor_co2_vs_other():
    assert runner._unit_factor("CO2 FFI") == 1e-12
    assert runner._unit_factor("CO2 AFOLU") == 1e-12
    assert runner._unit_factor("CH4") == 1e-9
    assert runner._unit_factor("Sulfur") == 1e-9


def test_run_fair_rejects_bad_output():
    import pytest

    with pytest.raises(ValueError, match="output must be"):
        runner.run_fair(
            "ssp245", None, __import__("numpy").array([2030, 2031]), output="bogus"
        )


def test_batch_size_is_at_least_one_and_shrinks_with_size():
    import numpy as np  # noqa: F401  (kept local to the test)

    small = runner._batch_size(n_timebounds=80, n_configs=841, n_species=61)
    large = runner._batch_size(n_timebounds=350, n_configs=841, n_species=61)
    huge = runner._batch_size(n_timebounds=10**6, n_configs=841, n_species=61)
    assert small >= large >= 1
    assert huge == 1


def test_cache_key_distinguishes_perturbations():
    import numpy as np

    years = np.array([2030, 2031])
    baseline = runner._cache_key("ssp245", years, None)
    assert baseline == runner._cache_key("ssp245", years, {})
    co2 = runner._cache_key("ssp245", years, {"CO2 FFI": np.array([1.0, 0.0])})
    ch4 = runner._cache_key("ssp245", years, {"CH4": np.array([1.0, 0.0])})
    other_amount = runner._cache_key(
        "ssp245", years, {"CO2 FFI": np.array([2.0, 0.0])}
    )
    other_marker = runner._cache_key("ssp126", years, {"CO2 FFI": np.array([1.0, 0.0])})
    assert len({baseline, co2, ch4, other_amount, other_marker}) == 5
    # Same perturbation, same key: that is what makes the memoization work.
    assert co2 == runner._cache_key("ssp245", years, {"CO2 FFI": np.array([1.0, 0.0])})


def _fake_template(background_gtco2=40.0, n_years=10):
    """Template stub with a constant CO2 background, for scale calculations."""
    import numpy as np

    return runner._Template(
        species=[],
        properties={},
        configs=[],
        specie_axis=["CO2 FFI"],
        emissions=np.full((n_years, 1), background_gtco2),
        forcing=None,
        species_configs=None,
        climate_configs=None,
        ebms=None,
    )


def test_perturbation_scale_lifts_lca_sized_pulses():
    import numpy as np

    years = np.arange(1750, 1760)
    tmpl = _fake_template()
    pulse = np.zeros(len(years))
    pulse[1] = 1e4  # 10 t CO2, a typical LCA magnitude

    scale = runner._perturbation_scale(tmpl, {"CO2 FFI": pulse}, years)
    # Scaled peak emission is the target fraction of the background.
    scaled_peak = pulse.max() * runner._unit_factor("CO2 FFI") * scale
    assert scaled_peak == pytest.approx(runner._LINEARIZATION_TARGET * 40.0)


def test_perturbation_scale_leaves_large_inventories_alone():
    import numpy as np

    years = np.arange(1750, 1760)
    tmpl = _fake_template()
    huge = np.zeros(len(years))
    huge[0] = 1e16  # 10 Gt CO2: already far above the noise floor

    assert runner._perturbation_scale(tmpl, {"CO2 FFI": huge}, years) == 1.0
    assert runner._perturbation_scale(tmpl, None, years) == 1.0
    assert runner._perturbation_scale(tmpl, {}, years) == 1.0


def test_ch4_method_matches_the_calibration():
    """The 1.4.1 ensemble is calibrated for Thornhill CH4 chemistry.

    Under FaIR's default ("leach2021") the calibration's per-species
    ch4_lifetime_chemical_sensitivity values are ignored, which silently
    removes the CH4-lifetime response to NOx, VOC and N2O.
    """
    assert runner._CH4_METHOD == "thornhill2021"


def test_nox_shortens_methane_lifetime():
    """A NOx perturbation must move CH4 - the check that caught the wrong mode."""
    import numpy as np
    import pytest

    pytest.importorskip("fair")
    years = np.arange(2030, 2051)
    perturbation = {"NOx": np.full(len(years), 10e9)}  # +10 Mt NOx/yr, unmistakable
    tmpl = runner._template("ssp245", int(years[-1]))
    f = runner._new_run(tmpl, 2, int(years[0]), int(years[-1]))
    runner._apply_state(f, runner._spinup_state("ssp245", int(years[-1]), int(years[0])))
    runner._apply_perturbation(f, tmpl, 1, perturbation, years, 1.0)
    with runner._RUN_LOCK:
        f.run(progress=False)

    ch4 = tmpl.specie_axis.index("CH4")
    concentration = f.concentration.values
    delta = np.median(concentration[-1, 1, :, ch4] - concentration[-1, 0, :, ch4])
    assert np.isfinite(delta)
    assert delta < 0  # more NOx -> shorter CH4 lifetime -> less CH4
