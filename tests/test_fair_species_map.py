import importlib.util
import os
import sys
import types

_fair_dir = os.path.join(
    os.path.dirname(__file__), "..", "dynamic_characterization", "fair"
)
sys.modules.setdefault(
    "dynamic_characterization", types.ModuleType("dynamic_characterization")
)
_pkg = types.ModuleType("dynamic_characterization.fair")
_pkg.__path__ = [_fair_dir]
sys.modules["dynamic_characterization.fair"] = _pkg
_spec = importlib.util.spec_from_file_location(
    "dynamic_characterization.fair.species_map",
    os.path.join(_fair_dir, "species_map.py"),
)
species_map = importlib.util.module_from_spec(_spec)
sys.modules["dynamic_characterization.fair.species_map"] = species_map
_spec.loader.exec_module(species_map)


def test_co2_resolves():
    sp, sign = species_map.resolve_species("Carbon dioxide, fossil")
    assert sp == "CO2 FFI"
    assert sign == 1


def test_co2_uptake_negative_sign():
    sp, sign = species_map.resolve_species("Carbon dioxide, in air")
    assert sp in {"CO2 FFI", "CO2 AFOLU"}
    assert sign == -1


def test_methane_resolves():
    sp, _ = species_map.resolve_species("Methane, fossil")
    assert sp == "CH4"


def test_precursor_maps_to_response_channel():
    m = species_map.load_species_map()
    assert "NOx" in m["precursors"]
    assert m["precursors"]["NOx"] in {"Ozone", "Aerosol-radiation interactions",
                                      "Aerosol-cloud interactions"}


def test_unmappable_flow_returns_none():
    sp, _ = species_map.resolve_species("Occupation, forest")
    assert sp is None
