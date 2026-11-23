"""Default characterization functions of the ``"ar6"`` parameter set, matched to biosphere flows."""

from typing import Callable, Dict, Iterable, NamedTuple, Optional

from loguru import logger

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
}


class FlowClass(NamedTuple):
    species: Optional[str]  # "CO2", "CO2 uptake", "CH4", "N2O", "CO", a canonical CAS, or None
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
    options = options or CharacterizationOptions(parameter_set="ar6")
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
        if species not in (None, "CO2", "CO2 uptake", "CH4", "N2O", "CO"):
            species = normalize_cas(species) or species
        origin = override.get("origin", origin)
    return FlowClass(species, origin)


def function_for(flow_class: FlowClass, options: CharacterizationOptions, flow_label: str = "") -> Optional[Callable]:
    """The ``"ar6"`` characterization function for a classified flow, or None."""
    ccf = ar6_core.ccf_parameters(options.ccf_parameters) if options.climate_carbon_feedback else None
    species, origin = flow_class

    if species is None:
        return None
    if species == "CO2":
        return characterize_co2_ar6
    if species == "CO2 uptake":
        return characterize_co2_uptake_ar6
    if species == "CH4":
        if origin is None:
            origin = options.unspecified_origin
            logger.warning(
                "Flow {!r} does not state a fossil or non-fossil origin; treated as {} CH4 "
                "(set 'unspecified_origin' or 'flow_overrides' to change this).",
                flow_label,
                origin,
            )
        spec = gases.ch4(
            origin=origin,
            biogenic_oxidation_co2=options.biogenic_oxidation_co2,
            oxidation_yield=options.ch4_oxidation_yield,
            ccf=ccf,
        )
    elif species == "N2O":
        spec = gases.n2o(ccf=ccf)
    elif species == "CO":
        spec = gases.co_oxidation()
    else:
        spec = gases.halocarbon(species, ccf=ccf)
        if spec is None:
            return None
    return make_characterization_function(spec)


def build_functions(
    nodes: Iterable,
    options: CharacterizationOptions,
    characterize_uptake: bool = True,
) -> Dict[int, Callable]:
    """``{node.id: function}`` for all nodes that have an ``"ar6"`` characterization function."""
    functions = {}
    for node in nodes:
        flow_class = classify_flow(node, options)
        if flow_class.species == "CO2 uptake" and not characterize_uptake:
            continue
        label = f"{node.get('name')} {tuple(node.get('categories') or ())}"
        func = function_for(flow_class, options, flow_label=label)
        if func is not None:
            functions[node.id] = func
    return functions
