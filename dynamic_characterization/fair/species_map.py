"""Resolve dynamic-inventory flows to FAIR species, signs and mass bases."""

import copy
import os
from functools import lru_cache
from typing import Dict, Optional, Tuple

import yaml


def _data_path() -> str:
    return os.path.join(os.path.dirname(__file__), "data", "species_map.yaml")


@lru_cache(maxsize=1)
def _load_species_map_cached() -> Dict:
    with open(_data_path()) as fh:
        raw = yaml.safe_load(fh)
    raw.setdefault("by_name", {})
    raw.setdefault("by_cas", {})
    raw.setdefault("signs", {})
    raw.setdefault("precursors", {})
    raw.setdefault("mass_basis", {})
    raw["mass_basis"].setdefault("by_name", {})
    raw["mass_basis"].setdefault("by_cas", {})
    raw["by_cas"] = {_normalize_cas(k): v for k, v in raw["by_cas"].items()}
    raw["mass_basis"]["by_cas"] = {
        _normalize_cas(k): v for k, v in raw["mass_basis"]["by_cas"].items()
    }
    # convenient alias kept for callers that only care about the flow -> species table
    raw["species"] = raw["by_name"]
    return raw


def load_species_map() -> Dict:
    """
    Load the flow mapping from YAML, including precursor metadata.

    The returned dict includes a 'precursors' key with informational metadata
    documenting which forcing channel each precursor drives in FAIR; these are
    not used for routing (FAIR applies precursor responses natively).
    """
    return copy.deepcopy(_load_species_map_cached())


def _normalize_cas(cas: Optional[str]) -> str:
    """CAS numbers are written zero-padded in some databases ('000074-82-8')."""
    if not cas:
        return ""
    cas = str(cas).strip()
    head, _, tail = cas.partition("-")
    return f"{head.lstrip('0') or '0'}-{tail}" if tail else cas


@lru_cache(maxsize=8192)
def resolve_flow(
    flow_name: str, cas: Optional[str] = None
) -> Tuple[Optional[str], float, float]:
    """
    Map a biosphere flow to a FAIR species, an emission sign and a mass factor.

    The flow name is matched exactly (case-insensitively); if that fails, the
    CAS number is tried, which keeps other naming conventions working. Uptake
    flows (CO2 in air, CO2 to soil or biomass stock, resource corrections) get
    sign -1. The mass factor converts the reported mass to the basis FAIR
    expects (e.g. NO -> NO2).

    Returns ``(species_or_None, sign, mass_factor)``.
    """
    mapping = _load_species_map_cached()
    name = (flow_name or "").strip().lower()
    cas_key = _normalize_cas(cas)

    species = mapping["by_name"].get(name)
    if species is None and cas_key:
        species = mapping["by_cas"].get(cas_key)
    if species is None:
        return None, 1.0, 1.0

    sign = mapping["signs"].get(name)
    if sign is None:
        # Fallback for naming conventions the table does not cover.
        sign = -1.0 if ("in air" in name or "uptake" in name) else 1.0
    factor = mapping["mass_basis"]["by_name"].get(
        name, mapping["mass_basis"]["by_cas"].get(cas_key, 1.0)
    )
    return species, float(sign), float(factor)


def resolve_species(
    flow_name: str, cas: Optional[str] = None
) -> Tuple[Optional[str], int]:
    """
    Map a flow name to a FAIR species and emission sign.

    Returns ``(species_or_None, sign)``. See :func:`resolve_flow` for the
    mass-basis factor that some flows also carry.
    """
    species, sign, _ = resolve_flow(flow_name, cas)
    return species, int(sign)
