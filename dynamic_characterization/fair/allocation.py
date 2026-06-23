"""Signed per-kg allocation of FAIR responses back to inventory flows."""

import numpy as np


def per_kg_response(
    response_series: np.ndarray, signed_cumulative_emissions: np.ndarray
) -> np.ndarray:
    """Per-kg response = response / total signed cumulative emission.

    Element-wise; 0 where the denominator is 0.
    """
    response_series = np.asarray(response_series, dtype="float64")
    denom = np.asarray(signed_cumulative_emissions, dtype="float64")
    out = np.zeros_like(response_series)
    nz = denom != 0
    out[nz] = response_series[nz] / denom[nz]
    return out


def allocate_to_flows(
    flow_signed_emissions: np.ndarray, per_kg: np.ndarray
) -> np.ndarray:
    """Allocate a species response to flows by signed cumulative emission.

    flow_signed_emissions: (n_flows, n_years) cumulative signed emissions.
    per_kg: (n_years,). Returns (n_flows, n_years).
    """
    flow_signed_emissions = np.asarray(flow_signed_emissions, dtype="float64")
    per_kg = np.asarray(per_kg, dtype="float64")
    return flow_signed_emissions * per_kg[None, :]


def safe_nanpercentile(stack: np.ndarray, quantiles, axis: int) -> np.ndarray:
    """np.nanpercentile with all-NaN slices coerced to 0."""
    stack = np.asarray(stack, dtype="float64")
    with np.errstate(all="ignore"):
        out = np.nanpercentile(stack, quantiles, axis=axis)
    return np.nan_to_num(out, nan=0.0)
