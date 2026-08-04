"""Scenario configuration for prospective and FAIR characterization."""

from typing import Dict, List, Optional, Tuple

_P = frozenset({"prospective"})
_PF = frozenset({"prospective", "fair"})

# Prospective (Watanabe IAM-SSP-RCP) scenarios.
# (iam, ssp, rcp) -> {"metrics": frozenset, "fair_marker": str | None}
# A few of these line up with an exact FAIR SSP marker and therefore also
# support the FAIR metrics.
SCENARIO_REGISTRY: Dict[Tuple[str, str, str], Dict] = {
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
}

# FAIR-native SSP marker scenarios. These are not tied to any IAM; FAIR runs
# the marker's own emission pathway. Selected with ``set_fair_scenario``.
# (ssp, rcp) -> FAIR marker string
FAIR_NATIVE_SCENARIOS: Dict[Tuple[str, str], str] = {
    ("SSP1", "1.9"): "ssp119",
    ("SSP1", "2.6"): "ssp126",
    ("SSP2", "4.5"): "ssp245",
    ("SSP3", "7.0"): "ssp370",
    ("SSP4", "3.4"): "ssp434",
    ("SSP4", "6.0"): "ssp460",
    ("SSP5", "3.4-over"): "ssp534-over",
    ("SSP5", "8.5"): "ssp585",
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
    List the prospective ``(iam, ssp, rcp)`` scenarios, optionally filtered.

    Parameters
    ----------
    metric : str, optional
        A metric string ('pGWP', 'fair_temperature', ...) or a family
        ('prospective' / 'fair'). If None, returns all prospective scenarios.

    Notes
    -----
    For ``metric="fair"`` this returns only the IAM scenarios that also line up
    with a FAIR marker. The IAM-agnostic FAIR-native scenarios (set with
    :func:`set_fair_scenario`) are listed by :func:`available_fair_scenarios`.
    """
    if metric is None:
        return sorted(SCENARIO_REGISTRY)
    family = metric if metric in {"prospective", "fair"} else metric_family(metric)
    return sorted(
        key for key, entry in SCENARIO_REGISTRY.items() if family in entry["metrics"]
    )


def available_fair_scenarios() -> List[Tuple[str, str]]:
    """List the FAIR-native ``(ssp, rcp)`` scenarios for :func:`set_fair_scenario`."""
    return sorted(FAIR_NATIVE_SCENARIOS)


# Module-level state. ``iam`` is None for FAIR-native scenarios.
_current_scenario: Optional[Dict[str, Optional[str]]] = None


def set_scenario(iam: str, ssp: str, rcp: str) -> None:
    """
    Set the current prospective ``(iam, ssp, rcp)`` scenario.

    Supports the prospective metrics (pGWP, pGTP,
    prospective_radiative_forcing); four of these scenarios additionally
    support the FAIR metrics. See ``available_scenarios(metric)``. To select a
    FAIR background that is not tied to an IAM, use :func:`set_fair_scenario`.
    """
    global _current_scenario
    if (iam, ssp, rcp) not in SCENARIO_REGISTRY:
        raise ValueError(
            f"Invalid scenario combination: ({iam}, {ssp}, {rcp}). "
            f"Valid scenarios: {sorted(SCENARIO_REGISTRY)}"
        )
    _current_scenario = {"iam": iam, "ssp": ssp, "rcp": rcp}


def set_fair_scenario(ssp: str, rcp: str) -> None:
    """
    Set the current FAIR-native ``(ssp, rcp)`` background scenario.

    FAIR runs the SSP marker's own emission pathway, so no IAM is involved.
    Supports only the FAIR metrics (fair_radiative_forcing, fair_temperature).
    See :func:`available_fair_scenarios` for the valid pairs.
    """
    global _current_scenario
    if (ssp, rcp) not in FAIR_NATIVE_SCENARIOS:
        raise ValueError(
            f"Invalid FAIR scenario: ({ssp}, {rcp}). "
            f"Valid FAIR scenarios: {available_fair_scenarios()}"
        )
    _current_scenario = {"iam": None, "ssp": ssp, "rcp": rcp}


def get_scenario() -> Dict[str, Optional[str]]:
    """
    Get the current scenario configuration.

    Returns
    -------
    dict
        Current scenario with keys: iam, ssp, rcp. ``iam`` is ``None`` for a
        FAIR-native scenario set via :func:`set_fair_scenario`.

    Raises
    ------
    RuntimeError
        If no scenario has been set.
    """
    if _current_scenario is None:
        raise RuntimeError(
            "No scenario set. Call prospective.set_scenario(iam, ssp, rcp) or "
            "prospective.set_fair_scenario(ssp, rcp) first."
        )
    return _current_scenario.copy()


def reset_scenario() -> None:
    """Reset scenario to None (for testing)."""
    global _current_scenario
    _current_scenario = None


def scenario_supports(metric_family_name: str) -> bool:
    """Whether the current scenario supports the given metric family."""
    s = get_scenario()
    if s["iam"] is None:  # FAIR-native scenario
        return metric_family_name == "fair"
    entry = SCENARIO_REGISTRY[(s["iam"], s["ssp"], s["rcp"])]
    return metric_family_name in entry["metrics"]


def current_fair_marker() -> str:
    """Return the FAIR marker for the current scenario, or raise."""
    s = get_scenario()
    if s["iam"] is None:  # FAIR-native scenario
        return FAIR_NATIVE_SCENARIOS[(s["ssp"], s["rcp"])]
    marker = SCENARIO_REGISTRY[(s["iam"], s["ssp"], s["rcp"])]["fair_marker"]
    if marker is None:
        raise ValueError(
            f"Scenario ({s['iam']}, {s['ssp']}, {s['rcp']}) does not support "
            f"fair metrics. Use set_fair_scenario(ssp, rcp) with one of "
            f"{available_fair_scenarios()}, or a fair-capable IAM scenario "
            f"from {available_scenarios('fair')}."
        )
    return marker
