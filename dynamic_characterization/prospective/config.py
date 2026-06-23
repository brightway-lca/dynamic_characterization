"""Scenario configuration for prospective characterization factors."""

from typing import Dict, List, Optional, Tuple

_P = frozenset({"prospective"})
_PF = frozenset({"prospective", "fair"})
_F = frozenset({"fair"})

# (iam, ssp, rcp) -> {"metrics": frozenset, "fair_marker": str | None}
SCENARIO_REGISTRY: Dict[Tuple[str, str, str], Dict] = {
    # --- prospective (Watanabe IAM-SSP-RCP) scenarios ---
    ("AIM", "SSP3", "4.5"): {"metrics": _P, "fair_marker": None},
    ("AIM", "SSP3", "6.0"): {"metrics": _P, "fair_marker": None},
    ("AIM", "SSP3", "8.5"): {"metrics": _P, "fair_marker": None},
    ("GCAM4", "SSP4", "2.6"): {"metrics": _P, "fair_marker": None},
    ("GCAM4", "SSP4", "4.5"): {"metrics": _P, "fair_marker": None},
    ("GCAM4", "SSP4", "6.0"): {"metrics": _PF, "fair_marker": "ssp460"},
    ("GCAM4", "SSP4", "8.5"): {"metrics": _P, "fair_marker": None},
    ("IMAGE", "SSP1", "2.6"): {"metrics": _PF, "fair_marker": "ssp126"},
    ("IMAGE", "SSP1", "4.5"): {"metrics": _P, "fair_marker": None},
    ("IMAGE", "SSP1", "8.5"): {"metrics": _P, "fair_marker": None},
    ("MESSAGE", "SSP2", "2.6"): {"metrics": _P, "fair_marker": None},
    ("MESSAGE", "SSP2", "4.5"): {"metrics": _PF, "fair_marker": "ssp245"},
    ("MESSAGE", "SSP2", "6.0"): {"metrics": _P, "fair_marker": None},
    ("MESSAGE", "SSP2", "8.5"): {"metrics": _P, "fair_marker": None},
    ("REMIND", "SSP5", "2.6"): {"metrics": _P, "fair_marker": None},
    ("REMIND", "SSP5", "4.5"): {"metrics": _P, "fair_marker": None},
    ("REMIND", "SSP5", "6.0"): {"metrics": _P, "fair_marker": None},
    ("REMIND", "SSP5", "8.5"): {"metrics": _PF, "fair_marker": "ssp585"},
    # --- FAIR-native marker scenarios (IAM-agnostic) ---
    ("FAIR", "SSP1", "1.9"): {"metrics": _F, "fair_marker": "ssp119"},
    ("FAIR", "SSP1", "2.6"): {"metrics": _F, "fair_marker": "ssp126"},
    ("FAIR", "SSP2", "4.5"): {"metrics": _F, "fair_marker": "ssp245"},
    ("FAIR", "SSP3", "7.0"): {"metrics": _F, "fair_marker": "ssp370"},
    ("FAIR", "SSP4", "3.4"): {"metrics": _F, "fair_marker": "ssp434"},
    ("FAIR", "SSP4", "6.0"): {"metrics": _F, "fair_marker": "ssp460"},
    ("FAIR", "SSP5", "3.4-over"): {"metrics": _F, "fair_marker": "ssp534-over"},
    ("FAIR", "SSP5", "8.5"): {"metrics": _F, "fair_marker": "ssp585"},
}

# Backwards-compatible alias: the set of valid (iam, ssp, rcp) keys.
VALID_SCENARIOS = set(SCENARIO_REGISTRY.keys())

_METRIC_FAMILY = {
    "pGWP": "prospective",
    "pGTP": "prospective",
    "prospective_radiative_forcing": "prospective",
    "fair_radiative_forcing": "fair",
    "fair_temperature": "fair",
}


def metric_family(metric: str) -> str:
    """Return 'prospective' or 'fair' for a scenario-dependent metric."""
    try:
        return _METRIC_FAMILY[metric]
    except KeyError as exc:
        raise ValueError(f"Metric {metric!r} is not scenario-dependent.") from exc


def available_scenarios(metric: Optional[str] = None) -> List[Tuple[str, str, str]]:
    """
    List scenarios, optionally filtered by metric or metric family.

    Parameters
    ----------
    metric : str, optional
        A metric string ('pGWP', 'fair_temperature', ...) or a family
        ('prospective' / 'fair'). If None, returns all scenario keys.
    """
    if metric is None:
        return sorted(SCENARIO_REGISTRY)
    family = metric if metric in {"prospective", "fair"} else metric_family(metric)
    return sorted(
        key
        for key, entry in SCENARIO_REGISTRY.items()
        if family in entry["metrics"]
    )

# Module-level state
_current_scenario: Optional[Dict[str, str]] = None


def set_scenario(iam: str, ssp: str, rcp: str) -> None:
    """
    Set the current (iam, ssp, rcp) scenario.

    Use ``available_scenarios(metric)`` to see which scenarios support a
    given metric. Prospective metrics (pGWP, pGTP,
    prospective_radiative_forcing) and FAIR metrics (fair_radiative_forcing,
    fair_temperature) are available for different scenario subsets; FAIR-only
    scenarios use ``iam="FAIR"``.
    """
    global _current_scenario
    if (iam, ssp, rcp) not in SCENARIO_REGISTRY:
        raise ValueError(
            f"Invalid scenario combination: ({iam}, {ssp}, {rcp}). "
            f"Valid scenarios: {sorted(SCENARIO_REGISTRY)}"
        )
    _current_scenario = {"iam": iam, "ssp": ssp, "rcp": rcp}


def get_scenario() -> Dict[str, str]:
    """
    Get the current scenario configuration.

    Returns
    -------
    dict
        Current scenario with keys: iam, ssp, rcp

    Raises
    ------
    RuntimeError
        If no scenario has been set
    """
    if _current_scenario is None:
        raise RuntimeError(
            "No scenario set. Call prospective.set_scenario(iam, ssp, rcp) first."
        )
    return _current_scenario.copy()


def reset_scenario() -> None:
    """Reset scenario to None (for testing)."""
    global _current_scenario
    _current_scenario = None


def scenario_supports(metric_family_name: str) -> bool:
    """Whether the current scenario supports the given metric family."""
    s = get_scenario()
    entry = SCENARIO_REGISTRY[(s["iam"], s["ssp"], s["rcp"])]
    return metric_family_name in entry["metrics"]


def current_fair_marker() -> str:
    """Return the FAIR marker for the current scenario, or raise."""
    s = get_scenario()
    entry = SCENARIO_REGISTRY[(s["iam"], s["ssp"], s["rcp"])]
    marker = entry["fair_marker"]
    if marker is None:
        raise ValueError(
            f"Scenario ({s['iam']}, {s['ssp']}, {s['rcp']}) does not support "
            f"fair metrics. Fair-capable scenarios: "
            f"{available_scenarios('fair')}"
        )
    return marker
