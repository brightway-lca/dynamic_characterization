"""Resolve dynamic-inventory flow names to FAIR species and signs."""

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
    raw.setdefault("species", {})
    raw.setdefault("precursors", {})
    raw.setdefault("negative_sign_species", [])
    # convenient alias used by callers/tests
    raw["signs"] = {sp: -1 for sp in raw["negative_sign_species"]}
    return raw


def load_species_map() -> Dict:
    """
    Load species mapping from YAML, including precursor metadata.

    The returned dict includes a 'precursors' key with informational metadata
    documenting which forcing channel each precursor drives in FAIR; these are
    not used for routing (FAIR applies precursor responses natively).
    """
    return copy.deepcopy(_load_species_map_cached())


def resolve_species(
    flow_name: str, cas: Optional[str] = None
) -> Tuple[Optional[str], int]:
    """
    Map a flow name to a FAIR species and emission sign.

    Returns (species_or_None, sign). Uptake flows ("in air" CO2 resources,
    explicit uptake) get sign -1. Unmappable flows return (None, sign).
    """
    name = flow_name.lower()
    mapping = load_species_map()

    sign = 1
    if "in air" in name or "uptake" in name:
        sign = -1

    for needle, species in mapping["species"].items():
        if needle in name:
            if species in mapping["signs"]:
                sign = mapping["signs"][species]
            return species, sign
    return None, sign
