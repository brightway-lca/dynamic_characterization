# dynamic_characterization/fair/core.py
"""Orchestrate a FAIR run from a dynamic inventory into a long DataFrame."""

from concurrent.futures import ThreadPoolExecutor
from typing import Dict, Optional, Tuple

import numpy as np
import pandas as pd

from ..prospective import config
from . import allocation, runner, species_map


def _emission_years(df: pd.DataFrame, time_horizon: Optional[int]) -> np.ndarray:
    start = int(df["date"].min().year)
    if time_horizon is not None:
        end = int(df["date"].max().year) + int(time_horizon)
    else:
        end = 2100
    end = max(end, start + 1)
    return np.arange(start, end + 1)


def _inventory_emissions_by_species(
    df: pd.DataFrame, years: np.ndarray
) -> Tuple[Dict[str, np.ndarray], list]:
    """Sum signed inventory emissions per FAIR species per year.

    Returns (by_species, flow_records) where flow_records is a list of
    dicts {flow, activity, species, signed_yearly (n_years,)}.
    """
    year_index = {int(y): i for i, y in enumerate(years)}
    by_species: Dict[str, np.ndarray] = {}
    flow_records = []

    # resolve species once per (flow, name) — here flow id used as name proxy
    for (flow, activity), grp in df.groupby(["flow", "activity"], sort=False):
        # flow name: callers pass a real name; fall back to str(flow)
        name = str(grp.get("flow_name", pd.Series([str(flow)])).iloc[0])
        sp, sign = species_map.resolve_species(name)
        if sp is None:
            continue
        signed = np.zeros(len(years), dtype="float64")
        for _, row in grp.iterrows():
            yi = year_index.get(int(pd.Timestamp(row["date"]).year))
            if yi is None:
                continue
            signed[yi] += sign * float(row["amount"])
        by_species.setdefault(sp, np.zeros(len(years), dtype="float64"))
        by_species[sp] += signed
        flow_records.append(
            {
                "flow": flow,
                "activity": activity,
                "species": sp,
                "signed": signed,
            }
        )
    return by_species, flow_records


def characterize_with_fair(
    dynamic_inventory_df: pd.DataFrame,
    output: str = "radiative_forcing",
    quantiles=(2.5, 25, 50, 75, 97.5),
    time_horizon: Optional[int] = None,
    workers: Optional[int] = None,
) -> pd.DataFrame:
    """Run FAIR and return ΔRF/ΔT per (year, flow, activity, quantile)."""
    runner.require_fair()  # fail fast with a clear message
    marker = config.current_fair_marker()

    years = _emission_years(dynamic_inventory_df, time_horizon)
    by_species, flow_records = _inventory_emissions_by_species(
        dynamic_inventory_df, years
    )
    if not by_species:
        return pd.DataFrame(
            columns=["date", "amount", "flow", "activity", "quantile"]
        )

    baseline = runner.run_fair(marker, None, years, output)  # (cfg, yr)

    # One perturbed run per species, split by sign, in parallel.
    def _species_response(sp: str) -> np.ndarray:
        delta = by_species[sp]
        pos = np.maximum(delta, 0.0)
        neg = np.minimum(delta, 0.0)
        resp = np.zeros_like(baseline)
        for part in (pos, neg):
            if not np.any(part):
                continue
            run = runner.run_fair(marker, {sp: part}, years, output)
            resp = resp + (run - baseline)
        return resp  # (cfg, yr)

    species_list = list(by_species)
    with ThreadPoolExecutor(max_workers=workers) as ex:
        responses = dict(zip(species_list, ex.map(_species_response, species_list)))

    rows = []
    qs = list(quantiles)
    for rec in flow_records:
        sp = rec["species"]
        species_cum = np.cumsum(by_species[sp])
        flow_cum = np.cumsum(rec["signed"])
        resp = responses[sp]  # (cfg, yr)
        # per-config per-kg, allocate this flow's share, then quantiles over cfg
        per_kg = resp / np.where(species_cum == 0, np.nan, species_cum)[None, :]
        flow_resp = np.nan_to_num(per_kg * flow_cum[None, :], nan=0.0)  # (cfg, yr)
        q_arr = allocation.safe_nanpercentile(flow_resp, qs, axis=0)  # (nq, yr)
        for qi, q in enumerate(qs):
            for yi, year in enumerate(years):
                rows.append(
                    {
                        "date": np.datetime64(f"{int(year)}-01-01"),
                        "amount": q_arr[qi, yi],
                        "flow": rec["flow"],
                        "activity": rec["activity"],
                        "quantile": q,
                    }
                )

    out = (
        pd.DataFrame(rows)
        .astype({"date": "datetime64[s]", "amount": "float64"})
        .query("amount != 0")[["date", "amount", "flow", "activity", "quantile"]]
        .sort_values(by=["quantile", "date", "flow"])
        .reset_index(drop=True)
    )
    return out
