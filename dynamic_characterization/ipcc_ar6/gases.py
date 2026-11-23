"""
Pulse-response specifications (:class:`~.ar6_core.GasSpec`) of individual gases for the
``"ar6"`` parameter set.

The values come from ``data/ar6_parameters.json`` and, for halocarbons, ``data/hodnebrog20.csv``
(the input file of the AR6 metric code). See ``ar6_core`` for the formulation.
"""

import csv
import json
import os
from functools import lru_cache
from typing import Dict, Optional, Tuple

from dynamic_characterization.flows import KNOWN_CAS_TYPOS, normalize_cas

from .ar6_core import (
    CCFParameters,
    GasSpec,
    co2_from_decay,
    co2_spec,
    exponential_spec,
    kg_per_ppb,
    load_parameters,
    precursor_component,
)

_DATA_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "data")

BIOGENIC_OXIDATION_MODES = ("ar6", "neutral", "explicit")

# CAS numbers that are blank in hodnebrog20.csv, as given in AR6 Table 7.SM.7 (by acronym)
_CAS_FROM_AR6_TABLE = {
    "HFO-1336mzz(E)": "66711-86-2",
    "HFE-569sf2 (HFE-7200)": "163702-05-4",
}


def co2() -> GasSpec:
    """1 kg CO2 (no climate-carbon feedback: the impulse response already contains the carbon cycle)."""
    return co2_spec()


def co2_uptake() -> GasSpec:
    """1 kg CO2 taken up from the atmosphere."""
    return co2_spec().scaled(-1.0, name="CO2 uptake")


def _ch4_radiative_efficiency_per_kg() -> float:
    p = load_parameters()
    ch4 = p["ch4"]
    re = (
        ch4["radiative_efficiency_direct_W_m2_ppb"]
        + ch4["radiative_efficiency_ozone_W_m2_ppb"]
        + ch4["radiative_efficiency_stratospheric_h2o_W_m2_ppb"]
    )
    return re / kg_per_ppb(p["constants"]["M_ch4"])


def ch4(
    origin: str = "fossil",
    biogenic_oxidation_co2: str = "ar6",
    oxidation_yield: Optional[float] = None,
    ccf: Optional[CCFParameters] = None,
) -> GasSpec:
    """
    1 kg CH4, including ozone and stratospheric water vapour (AR6 Table 7.SM.7).

    Parameters
    ----------
    origin : {"fossil", "non-fossil"}
        Fossil CH4 adds the CO2 from its oxidation (AR6 Eq. 7.SM.5.8; Table 7.15 "fossil").
    biogenic_oxidation_co2 : {"ar6", "neutral", "explicit"}
        Only for non-fossil CH4:

        - ``"ar6"``: oxidation CO2 is added and the carbon that was taken up by the biomass
          (M_CO2/M_CH4 kg CO2 per kg CH4) is subtracted at the time of emission, as in AR6
          Table 7.15 "non-fossil" (GWP100 27.0).
        - ``"neutral"``: neither term, i.e. the CH4 pulse alone (GWP100 27.9).
        - ``"explicit"``: oxidation CO2 is added, uptake is assumed to be modelled in the
          inventory (GWP100 29.8, as fossil).
    oxidation_yield : float, optional
        Fraction of the CH4 that is oxidised to CO2 (AR6: 0.75, range 0.5-1).
    ccf : CCFParameters, optional
        Climate-carbon feedback added to the CH4 forcing (AR6 adds it), or None. It is not
        applied to the oxidation CO2.
    """
    if origin not in ("fossil", "non-fossil"):
        raise ValueError(f"origin must be 'fossil' or 'non-fossil', not {origin!r}")
    if biogenic_oxidation_co2 not in BIOGENIC_OXIDATION_MODES:
        raise ValueError(
            f"biogenic_oxidation_co2 must be one of {BIOGENIC_OXIDATION_MODES}, not {biogenic_oxidation_co2!r}"
        )
    p = load_parameters()
    constants = p["constants"]
    yield_ = p["ch4"]["oxidation_yield"] if oxidation_yield is None else oxidation_yield
    if not 0 <= yield_ <= 1:
        raise ValueError(f"oxidation_yield must be between 0 and 1, not {yield_}")
    co2_per_ch4 = constants["M_co2"] / constants["M_ch4"]

    spec = exponential_spec("CH4", _ch4_radiative_efficiency_per_kg(), p["ch4"]["lifetime"], ccf=ccf)
    add_oxidation = origin == "fossil" or biogenic_oxidation_co2 in ("ar6", "explicit")
    if add_oxidation and yield_ > 0:
        spec = spec + co2_from_decay(yield_ * co2_per_ch4, p["ch4"]["oxidation_lifetime"])
    if origin == "non-fossil" and biogenic_oxidation_co2 == "ar6":
        spec = spec + co2_spec().scaled(-co2_per_ch4)
    return spec.renamed(f"CH4, {origin}")


def n2o(ccf: Optional[CCFParameters] = None) -> GasSpec:
    """1 kg N2O, including its effects on ozone and on the CH4 lifetime (AR6 Table 7.SM.7)."""
    p = load_parameters()
    n = p["n2o"]
    ch4_re_ppb = _ch4_radiative_efficiency_per_kg() * kg_per_ppb(p["constants"]["M_ch4"])
    re = (
        n["radiative_efficiency_direct_W_m2_ppb"]
        + n["radiative_efficiency_ozone_W_m2_ppb"]
        + n["ch4_per_n2o"] * ch4_re_ppb
    )
    return exponential_spec("N2O", re / kg_per_ppb(p["constants"]["M_n2o"]), n["lifetime"], ccf=ccf)


def co_oxidation(lifetime: float = 2.0 / 12.0) -> GasSpec:
    """
    1 kg CO as the CO2 from its oxidation only (as the legacy ``characterize_co``, but with a
    finite CO lifetime of about two months, AR6 Table 6.1). Its ozone and CH4 effects are not
    included.
    """
    constants = load_parameters()["constants"]
    co2_per_co = constants["M_co2"] / 28.01e-3
    return co2_from_decay(co2_per_co, lifetime, name="CO (oxidation only)")


@lru_cache(maxsize=None)
def halocarbon_table() -> Dict[str, dict]:
    """
    Rows of ``hodnebrog20.csv`` with a lifetime, keyed by canonical CAS number.

    Applies the CFC-11/12 tropospheric adjustment and the lifetime overrides of
    ``ar6_parameters.json`` (SF6: 1000 yr as in AR6 Table 7.15).
    """
    adjustment = load_parameters()["halocarbons"]["tropospheric_adjustment"]
    lifetime_overrides = load_parameters()["halocarbons"]["lifetime_overrides"]
    table = {}
    with open(os.path.join(_DATA_DIR, "hodnebrog20.csv"), newline="", encoding="utf-8-sig") as f:
        for row in csv.DictReader(f):
            acronym = (row.get("Acronym") or "").strip()
            cas = normalize_cas(row.get("CASRN")) or normalize_cas(_CAS_FROM_AR6_TABLE.get(acronym))
            try:
                lifetime = float(row["Lifetime (yr)"])
                re = float(row["RE (W m-2 ppb-1)"])
                mass = float(row["Molar mass"])
            except (TypeError, ValueError):
                continue
            if cas is None or cas in table:
                continue
            if cas in lifetime_overrides:
                lifetime = lifetime_overrides[cas]["lifetime"]
            table[cas] = {
                "name": row["Name"].strip(),
                "acronym": acronym,
                "lifetime": lifetime,
                "radiative_efficiency": re * (1 + adjustment.get(acronym, 0.0)),
                "molar_mass": mass,
            }
    for correct, typo in KNOWN_CAS_TYPOS.items():
        if typo in table:
            table.setdefault(correct, table[typo])
    return table


def halocarbon(cas: str, ccf: Optional[CCFParameters] = None) -> Optional[GasSpec]:
    """1 kg of the halocarbon with this CAS number (AR6 Table 7.SM.7), or None if not listed."""
    row = halocarbon_table().get(normalize_cas(cas))
    if row is None:
        return None
    return exponential_spec(
        row["acronym"] or row["name"],
        row["radiative_efficiency"] / kg_per_ppb(row["molar_mass"]),
        row["lifetime"],
        ccf=ccf,
    )


# --- indirect (precursor) gases -------------------------------------------------------------

PRECURSOR_SPECIES = ("CO", "VOC", "H2", "NOx", "nitrate")


@lru_cache(maxsize=None)
def load_precursor_parameters() -> dict:
    """The published precursor parameter sets in ``data/precursor_parameters.json``."""
    with open(os.path.join(_DATA_DIR, "precursor_parameters.json")) as f:
        return json.load(f)


def precursor_parameter_set(species: str, name: Optional[str] = None) -> dict:
    """The parameter set ``name`` of ``species`` (default set if None)."""
    data = load_precursor_parameters()
    if species not in data["sets"]:
        raise ValueError(f"Unknown precursor species {species!r}; available: {sorted(data['sets'])}")
    name = name or data["defaults"][species]
    sets = data["sets"][species]
    if name not in sets:
        raise ValueError(f"Unknown parameter set {name!r} for {species}; available: {sorted(sets)}")
    return sets[name]


def precursor(species: str, name: Optional[str] = None, ccf: Optional[CCFParameters] = None) -> GasSpec:
    """
    Indirect forcing of 1 kg of a precursor, per kg of the set's ``basis`` species (e.g. kg C
    for the Fry et al. 2012 NMVOC set, kg N for NOx). See ``data/precursor_parameters.json``.

    ``ccf`` adds the climate-carbon feedback (the published sets do not include it).
    """
    p = precursor_parameter_set(species, name)
    components = tuple(
        precursor_component(
            c["steady_state_forcing"],
            c["response_time"],
            p["precursor_lifetime"],
            ccf=ccf,
        )
        for c in p["components"]
    )
    return GasSpec(name=f"{species} ({name or load_precursor_parameters()['defaults'][species]})", components=components)


def precursor_component_specs(species: str, name: Optional[str] = None) -> Tuple[Tuple[str, GasSpec], ...]:
    """The components of a parameter set as separate specs ``(component name, spec)``, for checks."""
    p = precursor_parameter_set(species, name)
    return tuple(
        (
            c["name"],
            GasSpec(
                name=c["name"],
                components=(precursor_component(c["steady_state_forcing"], c["response_time"], p["precursor_lifetime"]),),
            ),
        )
        for c in p["components"]
    )


def mass_basis_factor(key: str) -> float:
    """kg of the basis species per kg of flow, e.g. ``"kg N per kg NO"``."""
    return load_precursor_parameters()["mass_basis"][key]


def co(
    origin: str = "fossil",
    parameter_set: Optional[str] = None,
    biogenic_oxidation_co2: str = "neutral",
    ccf: Optional[CCFParameters] = None,
    lifetime: float = 2.0 / 12.0,
) -> GasSpec:
    """
    1 kg CO: indirect forcing (O3, CH4, sulfate) plus, for fossil CO, the CO2 from its oxidation.

    Non-fossil CO adds the oxidation CO2 only with ``biogenic_oxidation_co2="explicit"`` (the
    carbon uptake is then assumed to be in the inventory); otherwise it is treated as carbon
    neutral, as in the static IPCC 2021 factors (2.49 non-fossil, 4.06 fossil).
    """
    spec = precursor("CO", parameter_set, ccf=ccf)
    if origin == "fossil" or biogenic_oxidation_co2 == "explicit":
        spec = spec + co_oxidation(lifetime=lifetime)
    return spec.renamed(f"CO, {origin} ({parameter_set or load_precursor_parameters()['defaults']['CO']})")


def voc(parameter_set: Optional[str] = None, ccf: Optional[CCFParameters] = None) -> GasSpec:
    """1 kg VOC; converted to the basis of the set (kg C for Fry et al. 2012)."""
    p = precursor_parameter_set("VOC", parameter_set)
    factor = mass_basis_factor("kg C per kg VOC") if p["basis"] == "kg C" else 1.0
    return precursor("VOC", parameter_set, ccf=ccf).scaled(factor, name=f"VOC ({parameter_set or 'default'})")


def h2(parameter_set: Optional[str] = None, ccf: Optional[CCFParameters] = None) -> GasSpec:
    """1 kg H2 (no direct forcing; CH4, O3 and stratospheric H2O effects)."""
    return precursor("H2", parameter_set, ccf=ccf)


def nox(
    emitted_as: str = "NO",
    parameter_set: Optional[str] = None,
    nitrate: bool = False,
    nitrate_parameter_set: Optional[str] = None,
    ccf: Optional[CCFParameters] = None,
) -> GasSpec:
    """
    1 kg of NO or NO2 (``emitted_as``), from a per-kg-N parameter set. Optionally adds the
    direct forcing of nitrate aerosol. Sign and size of the NOx metric are very uncertain.
    """
    if emitted_as not in ("NO", "NO2"):
        raise ValueError("emitted_as must be 'NO' or 'NO2'")
    spec = precursor("NOx", parameter_set, ccf=ccf)
    if nitrate:
        spec = spec + precursor("nitrate", nitrate_parameter_set, ccf=ccf)
    return spec.scaled(mass_basis_factor(f"kg N per kg {emitted_as}"), name=f"NOx as {emitted_as}")
