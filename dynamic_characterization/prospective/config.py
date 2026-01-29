"""Scenario configuration for prospective characterization factors."""

from typing import Dict, Optional, Set, Tuple

# Valid IAM-SSP-RCP combinations from the paper
# Each IAM is associated with one SSP but multiple RCPs
VALID_SCENARIOS: Set[Tuple[str, str, str]] = {
    # AIM - SSP3
    ("AIM", "SSP3", "4.5"),
    ("AIM", "SSP3", "6.0"),
    ("AIM", "SSP3", "8.5"),
    # GCAM4 - SSP4
    ("GCAM4", "SSP4", "2.6"),
    ("GCAM4", "SSP4", "4.5"),
    ("GCAM4", "SSP4", "6.0"),
    ("GCAM4", "SSP4", "8.5"),
    # IMAGE - SSP1
    ("IMAGE", "SSP1", "2.6"),
    ("IMAGE", "SSP1", "4.5"),
    ("IMAGE", "SSP1", "8.5"),
    # MESSAGE - SSP2
    ("MESSAGE", "SSP2", "2.6"),
    ("MESSAGE", "SSP2", "4.5"),
    ("MESSAGE", "SSP2", "6.0"),
    ("MESSAGE", "SSP2", "8.5"),
    # REMIND - SSP5
    ("REMIND", "SSP5", "2.6"),
    ("REMIND", "SSP5", "4.5"),
    ("REMIND", "SSP5", "6.0"),
    ("REMIND", "SSP5", "8.5"),
}

# Module-level state
_current_scenario: Optional[Dict[str, str]] = None


def set_scenario(iam: str, ssp: str, rcp: str) -> None:
    """
    Set the current IAM-SSP-RCP scenario for characterization.

    Parameters
    ----------
    iam : str
        Integrated Assessment Model: AIM, GCAM4, IMAGE, MESSAGE, or REMIND
    ssp : str
        Shared Socioeconomic Pathway: SSP1, SSP2, SSP3, SSP4, or SSP5
    rcp : str
        Representative Concentration Pathway: 2.6, 4.5, 6.0, or 8.5

    Raises
    ------
    ValueError
        If the combination is not valid per Watanabe et al. (2026).
    """
    global _current_scenario

    if (iam, ssp, rcp) not in VALID_SCENARIOS:
        raise ValueError(
            f"Invalid scenario combination: ({iam}, {ssp}, {rcp}). "
            f"Each IAM is associated with a specific SSP. "
            f"Valid combinations: IMAGE-SSP1, MESSAGE-SSP2, AIM-SSP3, GCAM4-SSP4, REMIND-SSP5"
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
