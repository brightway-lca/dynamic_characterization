"""Data loading functions for prospective characterization SI Excel files."""

import os
from functools import lru_cache
from typing import Dict

import numpy as np
import pandas as pd


def _get_data_dir() -> str:
    """Return path to prospective data directory."""
    return os.path.join(os.path.dirname(__file__), "data")


@lru_cache(maxsize=1)
def load_irf_ch4() -> np.ndarray:
    """
    Load impulse response function for CH4.

    Returns 101-element array (years 0-100) of IRF values.
    All RCP scenarios use the same CH4 lifetime (11.8 years).
    """
    filepath = os.path.join(_get_data_dir(), "es5c12391_si_009.xlsx")
    df = pd.read_excel(filepath)
    # All columns are identical, use first IRF column
    return df.iloc[:, 1].values


@lru_cache(maxsize=1)
def load_irf_co2() -> Dict[str, np.ndarray]:
    """
    Load impulse response functions for CO2 by RCP scenario.

    Returns dict mapping RCP name to 101-element IRF array.
    CO2 IRF varies by RCP due to carbon cycle feedbacks.
    """
    filepath = os.path.join(_get_data_dir(), "es5c12391_si_010.xlsx")
    df = pd.read_excel(filepath)

    return {
        "RCP26": df.iloc[:, 1].values,
        "RCP45": df.iloc[:, 2].values,
        "RCP60": df.iloc[:, 3].values,
        "RCP85": df.iloc[:, 4].values,
    }


@lru_cache(maxsize=1)
def load_irf_n2o() -> np.ndarray:
    """
    Load impulse response function for N2O.

    Returns 100-element array (years 0-99) of IRF values.
    All scenarios use the same N2O lifetime (109 years).
    """
    filepath = os.path.join(_get_data_dir(), "es5c12391_si_011.xlsx")
    df = pd.read_excel(filepath)
    # All columns are identical, use first IRF column
    return df.iloc[:, 1].values


# Mapping of SI file numbers to IAM names
_RE_CH4_FILES = {
    "AIM": "es5c12391_si_012.xlsx",
    "GCAM4": "es5c12391_si_013.xlsx",
    "IMAGE": "es5c12391_si_014.xlsx",
    "MESSAGE": "es5c12391_si_015.xlsx",
    "REMIND": "es5c12391_si_016.xlsx",
}

_RE_CO2_FILES = {
    "AIM": "es5c12391_si_017.xlsx",
    "GCAM4": "es5c12391_si_018.xlsx",
    "IMAGE": "es5c12391_si_019.xlsx",
    "MESSAGE": "es5c12391_si_020.xlsx",
    "REMIND": "es5c12391_si_021.xlsx",
}

_RE_N2O_FILES = {
    "AIM": "es5c12391_si_022.xlsx",
    "GCAM4": "es5c12391_si_023.xlsx",
    "IMAGE": "es5c12391_si_024.xlsx",
    "MESSAGE": "es5c12391_si_025.xlsx",
    "REMIND": "es5c12391_si_026.xlsx",
}


def _parse_re_column_name(col: str) -> tuple:
    """
    Parse RE column name to extract SSP and RCP.

    Examples:
        "AIM - SSP3 - 4.5" -> ("SSP3", "4.5")
        "GCAM4 -SSP4 - 2.6" -> ("SSP4", "2.6")
        "MESSAGE-GLOBIOM -SSP2 - 4.5" -> ("SSP2", "4.5")
        "REMIND-MAGPIE -SSP5 - 8.5" -> ("SSP5", "8.5")
    """
    import re

    # Find SSP pattern (SSP followed by digit)
    ssp_match = re.search(r"(SSP\d)", col)
    if not ssp_match:
        return None
    ssp = ssp_match.group(1)

    # Find RCP pattern (number like 2.6, 4.5, 6.0, 8.5 at end)
    rcp_match = re.search(r"(\d+\.\d+)\s*$", col)
    if not rcp_match:
        return None
    rcp = rcp_match.group(1)

    return (ssp, rcp)


def _load_re_file(filepath: str) -> Dict[str, Dict[str, np.ndarray]]:
    """
    Load RE file and return nested dict: {ssp: {rcp: array}}.

    Years in file run from 2020-2150 (131 rows).
    """
    df = pd.read_excel(filepath)
    result = {}

    for col in df.columns[1:]:  # Skip 'Year' column
        parsed = _parse_re_column_name(col)
        if parsed:
            ssp, rcp = parsed
            if ssp not in result:
                result[ssp] = {}
            result[ssp][rcp] = df[col].values

    # Store years as well
    result["_years"] = df["Year"].values

    return result


@lru_cache(maxsize=5)
def load_re_ch4(iam: str = "IMAGE") -> Dict[str, Dict[str, np.ndarray]]:
    """
    Load radiative efficiency data for CH4 by IAM.

    Parameters
    ----------
    iam : str
        IAM name: AIM, GCAM4, IMAGE, MESSAGE, or REMIND

    Returns
    -------
    Dict with structure {ssp: {rcp: array}} plus "_years" key
    """
    if iam not in _RE_CH4_FILES:
        raise ValueError(f"Unknown IAM: {iam}. Valid: {list(_RE_CH4_FILES.keys())}")
    filepath = os.path.join(_get_data_dir(), _RE_CH4_FILES[iam])
    return _load_re_file(filepath)


@lru_cache(maxsize=5)
def load_re_co2(iam: str = "IMAGE") -> Dict[str, Dict[str, np.ndarray]]:
    """
    Load radiative efficiency data for CO2 by IAM.

    Parameters
    ----------
    iam : str
        IAM name: AIM, GCAM4, IMAGE, MESSAGE, or REMIND

    Returns
    -------
    Dict with structure {ssp: {rcp: array}} plus "_years" key
    """
    if iam not in _RE_CO2_FILES:
        raise ValueError(f"Unknown IAM: {iam}. Valid: {list(_RE_CO2_FILES.keys())}")
    filepath = os.path.join(_get_data_dir(), _RE_CO2_FILES[iam])
    return _load_re_file(filepath)


@lru_cache(maxsize=5)
def load_re_n2o(iam: str = "IMAGE") -> Dict[str, Dict[str, np.ndarray]]:
    """
    Load radiative efficiency data for N2O by IAM.

    Parameters
    ----------
    iam : str
        IAM name: AIM, GCAM4, IMAGE, MESSAGE, or REMIND

    Returns
    -------
    Dict with structure {ssp: {rcp: array}} plus "_years" key
    """
    if iam not in _RE_N2O_FILES:
        raise ValueError(f"Unknown IAM: {iam}. Valid: {list(_RE_N2O_FILES.keys())}")
    filepath = os.path.join(_get_data_dir(), _RE_N2O_FILES[iam])
    return _load_re_file(filepath)
