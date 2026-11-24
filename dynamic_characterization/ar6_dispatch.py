"""Default characterization functions of the ``"ar6"`` parameter set, matched to biosphere flows."""

from typing import Callable, Dict, Iterable, NamedTuple, Optional

from loguru import logger

from dynamic_characterization import ar6_prospective
from dynamic_characterization.flows import is_co2_uptake, normalize_cas
from dynamic_characterization.ipcc_ar6 import ar6_core, gases
from dynamic_characterization.ipcc_ar6.radiative_forcing import (
    characterize_co2_ar6,
    characterize_co2_uptake_ar6,
    make_characterization_function,
)
from dynamic_characterization.options import CharacterizationOptions

SPECIES_BY_CAS = {
    "124-38-9": "CO2",
    "74-82-8": "CH4",
    "10024-97-2": "N2O",
    "630-08-0": "CO",
    "1333-74-0": "H2",
    "10102-43-9": "NO",
    "10102-44-0": "NO2",
    "11104-93-1": "NOx",
}

# exact (lower-case) names, for flows without a CAS number. Never match by substring:
# "methane" is also part of "Methane, dichloro-", "Tetrafluoromethane", ...
SPECIES_BY_NAME = {
    "carbon dioxide": "CO2",
    "carbon dioxide, fossil": "CO2",
    "carbon dioxide, non-fossil": "CO2",
    "carbon dioxide, from soil or biomass stock": "CO2",
    "carbon dioxide, to soil or biomass stock": "CO2",
    "carbon dioxide, in air": "CO2",
    "methane": "CH4",
    "methane, fossil": "CH4",
    "methane, non-fossil": "CH4",
    "methane, from soil or biomass stock": "CH4",
    "dinitrogen monoxide": "N2O",
    "carbon monoxide": "CO",
    "carbon monoxide, fossil": "CO",
    "carbon monoxide, non-fossil": "CO",
    "carbon monoxide, from soil or biomass stock": "CO",
    "hydrogen": "H2",
    "voc, volatile organic compounds": "VOC",
    "nmvoc, non-methane volatile organic compounds": "VOC",
    "nitric oxide": "NO",
    "nitrogen dioxide": "NO2",
    "nitrogen oxides": "NOx",
}

_NAMED_SPECIES = ("CO2", "CO2 uptake", "CH4", "N2O", "CO", "H2", "VOC", "NO", "NO2", "NOx")


class FlowClass(NamedTuple):
    species: Optional[str]  # one of _NAMED_SPECIES, a canonical CAS, or None
    origin: Optional[str]  # "fossil", "non-fossil" or None (not stated in the flow name)


def origin_from_name(name: str) -> Optional[str]:
    """Fossil or non-fossil origin as stated in an ecoinvent-style flow name, else None.

    "from soil or biomass stock" counts as fossil, as in the ecoinvent IPCC 2021 implementation
    (carbon released from a stock, not recently taken up).
    """
    name = name.lower()
    if ", non-fossil" in name:
        return "non-fossil"
    if ", fossil" in name or "from soil or biomass stock" in name:
        return "fossil"
    return None


def classify_flow(node, options: Optional[CharacterizationOptions] = None) -> FlowClass:
    """Species and origin of a biosphere flow, for the ``"ar6"`` default functions."""
    options = options or CharacterizationOptions()
    name = (node.get("name") or "").strip()
    origin = origin_from_name(name)
    cas = normalize_cas(node.get("CAS number"))

    if cas is not None:
        species = SPECIES_BY_CAS.get(cas, cas)
    else:
        species = SPECIES_BY_NAME.get(name.lower())

    if species == "CO2" and is_co2_uptake(node):
        species = "CO2 uptake"

    override = options.override_for(node)
    if override:
        species = override.get("species", species)
        if species is not None and species not in _NAMED_SPECIES:
            species = normalize_cas(species) or species
        origin = override.get("origin", origin)
    return FlowClass(species, origin)


def function_for(
    flow_class: FlowClass,
    options: CharacterizationOptions,
    flow_label: str = "",
    prospective: bool = False,
) -> Optional[Callable]:
    """
    The ``"ar6"`` characterization function for a classified flow, or None.

    With ``prospective=True`` the function follows the scenario of ``prospective.set_scenario``
    (see ``ar6_prospective``): CO2 terms use the scenario CO2 response, CH4 and N2O the
    scenario radiative efficiencies.
    """
    ccf = ar6_core.ccf_parameters(options.ccf_parameters) if options.climate_carbon_feedback else None
    precursor_ccf = ccf if options.precursor_ccf else None
    species, origin = flow_class

    if species is None:
        return None
    method = options.prospective_method
    if species == "CO2":
        return ar6_prospective.co2_function(method) if prospective else characterize_co2_ar6
    if species == "CO2 uptake":
        return ar6_prospective.co2_function(method, uptake=True) if prospective else characterize_co2_uptake_ar6
    if prospective and method == "watanabe2026" and species in ("CH4", "N2O"):
        if species == "CH4" and origin is None:
            logger.warning(
                "Flow {!r} does not state a fossil or non-fossil origin; the Watanabe & Cherubini "
                "(2026) CH4 metric does not distinguish them.",
                flow_label,
            )
        return ar6_prospective.paper_function(species)
    if species == "CH4":
        if origin is None:
            origin = options.unspecified_origin
            logger.warning(
                "Flow {!r} does not state a fossil or non-fossil origin; treated as {} CH4 "
                "(set 'unspecified_origin' or 'flow_overrides' to change this).",
                flow_label,
                origin,
            )
        ch4_options = dict(
            origin=origin,
            biogenic_oxidation_co2=options.biogenic_oxidation_co2,
            oxidation_yield=options.ch4_oxidation_yield,
            ccf=ccf,
        )
        if prospective:
            return ar6_prospective.make_prospective_function(
                ar6_prospective.Recipe(kind="CH4", method="ar6", options=tuple(sorted(ch4_options.items())))
            )
        spec = gases.ch4(**ch4_options)
    elif species == "N2O":
        if prospective:
            return ar6_prospective.make_prospective_function(
                ar6_prospective.Recipe(kind="N2O", method="ar6", options=(("ccf", ccf),))
            )
        spec = gases.n2o(ccf=ccf)
    elif species == "CO":
        if "CO" in options.include_precursors:
            if origin is None:
                origin = options.unspecified_origin
            spec = gases.co(
                origin=origin,
                parameter_set=options.precursor_parameter_set("CO"),
                biogenic_oxidation_co2="explicit" if options.biogenic_oxidation_co2 == "explicit" else "neutral",
                ccf=precursor_ccf,
            )
        else:
            spec = gases.co_oxidation()
    elif species == "VOC":
        if "VOC" not in options.include_precursors:
            return None
        spec = gases.voc(options.precursor_parameter_set("VOC"), ccf=precursor_ccf)
    elif species == "H2":
        if "H2" not in options.include_precursors:
            return None
        spec = gases.h2(options.precursor_parameter_set("H2"), ccf=precursor_ccf)
    elif species in ("NO", "NO2", "NOx"):
        if not options.include_nox:
            return None
        spec = gases.nox(
            emitted_as="NO" if species == "NO" else "NO2",
            parameter_set=options.precursor_parameter_set("NOx"),
            nitrate=options.nox_nitrate,
            nitrate_parameter_set=options.precursor_parameter_set("nitrate"),
            ccf=precursor_ccf,
        )
    else:
        spec = gases.halocarbon(species, ccf=ccf)
        if spec is None:
            return None
    if prospective:
        return ar6_prospective.fixed(spec, method=method)
    return make_characterization_function(spec)


SCENARIO_SPECIES = ("CO2", "CO2 uptake", "CH4", "N2O")


def build_functions(
    nodes: Iterable,
    options: CharacterizationOptions,
    characterize_uptake: bool = True,
    prospective: bool = False,
    fallback_to_ipcc: bool = True,
) -> Dict[int, Callable]:
    """
    ``{node.id: function}`` for all nodes that have an ``"ar6"`` characterization function.

    With ``prospective=True`` the functions follow the prospective scenario; with
    ``fallback_to_ipcc=False`` only the gases with scenario data (CO2, CH4, N2O) are included.
    """
    functions = {}
    for node in nodes:
        flow_class = classify_flow(node, options)
        if flow_class.species == "CO2 uptake" and not characterize_uptake:
            continue
        if prospective and not fallback_to_ipcc and flow_class.species not in SCENARIO_SPECIES:
            continue
        label = f"{node.get('name')} {tuple(node.get('categories') or ())}"
        func = function_for(flow_class, options, flow_label=label, prospective=prospective)
        if func is not None:
            functions[node.id] = func
    return functions
