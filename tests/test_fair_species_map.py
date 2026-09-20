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
    sp, sign = species_map.resolve_species("Methane, fossil")
    assert sp == "CH4"
    assert sign == 1


def test_precursor_maps_to_response_channel():
    m = species_map.load_species_map()
    assert "NOx" in m["precursors"]
    assert m["precursors"]["NOx"] in {"Ozone", "Aerosol-radiation interactions",
                                      "Aerosol-cloud interactions"}


def test_unmappable_flow_returns_none():
    sp, sign = species_map.resolve_species("Occupation, forest")
    assert sp is None
    assert sign == 1


def test_halocarbons_are_not_swallowed_by_the_methane_rule():
    """Substring matching would route every *methane to CH4."""
    for name, expected in (
        ("Tetrafluoromethane", "CF4"),
        ("Trifluoromethane", "HFC-23"),
        ("Dichloromethane", "CH2Cl2"),
        ("Dichlorodifluoromethane", "CFC-12"),
        ("Bromotrifluoromethane", "Halon-1301"),
        ("NMVOC, non-methane volatile organic compounds", "VOC"),
    ):
        assert species_map.resolve_species(name)[0] == expected


def test_fossil_and_biogenic_co2_are_separated():
    assert species_map.resolve_species("Carbon dioxide, fossil")[0] == "CO2 FFI"
    for name in (
        "Carbon dioxide, non-fossil",
        "Carbon dioxide, in air",
        "Carbon dioxide, non-fossil, resource correction",
        "Carbon dioxide, to soil or biomass stock",
    ):
        assert species_map.resolve_species(name)[0] == "CO2 AFOLU"
    # Uptake and resource-correction flows remove CO2 from the atmosphere.
    for name in (
        "Carbon dioxide, in air",
        "Carbon dioxide, non-fossil, resource correction",
        "Carbon dioxide, to soil or biomass stock",
    ):
        assert species_map.resolve_species(name)[1] == -1
    assert species_map.resolve_species("Carbon dioxide, non-fossil")[1] == 1


def test_cas_fallback_resolves_unknown_naming():
    # Zero-padded (ecoinvent) and plain spellings both work.
    assert species_map.resolve_species("Ethane, hexafluoro-", "000076-16-4")[0] == "C2F6"
    assert species_map.resolve_species("whatever", "76-16-4")[0] == "C2F6"


def test_mass_basis_converts_to_the_unit_fair_expects():
    _, sign, factor = species_map.resolve_flow("Nitric oxide")
    assert (sign, round(factor, 3)) == (1.0, 1.534)  # NO reported, NO2 expected
    _, _, factor = species_map.resolve_flow("Nitrogen oxides")
    assert factor == 1.0


def test_every_mapped_species_exists_in_fair():
    import csv
    import importlib.util

    if importlib.util.find_spec("fair") is None:
        import pytest

        pytest.skip("fair not installed")
    from fair.structure.units import desired_emissions_units

    mapping = species_map.load_species_map()
    targets = set(mapping["by_name"].values()) | set(mapping["by_cas"].values())
    unknown = targets - set(desired_emissions_units)
    assert not unknown, f"not FAIR species: {sorted(unknown)}"
