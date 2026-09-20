# dynamic_characterization/fair/core.py
"""Orchestrate a FAIR run from a dynamic inventory into a long DataFrame.

Everything here is sized for real dynamic inventories (hundreds of thousands
of rows): the inventory is aggregated with vectorized pandas/numpy, the FAIR
runs are batched and cached in :mod:`.runner`, and the per-flow attribution
never touches a Python loop over rows.
"""

from typing import Dict, List, Optional, Tuple

import numpy as np
import pandas as pd
from loguru import logger

from ..prospective import config
from . import allocation, runner, species_map

_OUT_COLUMNS = ["date", "amount", "flow", "activity", "quantile"]

# Species below this total mass get no FAIR run: solving the inventory leaves
# numerical dust (picograms of an F-gas) that costs a model run and cannot
# move any result, even for the most potent species.
_NEGLIGIBLE_KG = 1e-12


def _emission_years(df: pd.DataFrame, time_horizon: Optional[int]) -> np.ndarray:
    start = int(df["date"].min().year)
    if time_horizon is not None:
        end = int(df["date"].max().year) + int(time_horizon)
    else:
        end = 2100
    end = max(end, start + 1)
    return np.arange(start, end + 1)


def _flow_names(df: pd.DataFrame) -> pd.Series:
    """Flow names to map to FAIR species, falling back to the flow id."""
    if "flow_name" in df.columns:
        names = df["flow_name"]
        if names.isna().any():
            names = names.fillna(df["flow"].astype(str))
        return names
    return df["flow"].astype(str)


def _flow_cas(df: pd.DataFrame, names: pd.Series) -> Dict[str, Optional[str]]:
    """CAS number per distinct flow name, when the inventory carries one."""
    if "flow_cas" not in df.columns:
        return {}
    cas = (
        pd.DataFrame({"name": names, "cas": df["flow_cas"]})
        .dropna()
        .drop_duplicates("name")
    )
    return dict(zip(cas["name"], cas["cas"]))


def _warn_about_unmapped(df: pd.DataFrame, names: pd.Series, mapped: np.ndarray):
    """Report the largest flows that no FAIR species covers.

    Silence here is how under-coverage hides: an unmapped greenhouse gas simply
    does not contribute, and the result still looks plausible.
    """
    unmapped = ~mapped
    if not unmapped.any():
        return
    amounts = df["amount"].to_numpy(dtype="float64")
    by_name = (
        pd.Series(np.abs(amounts[unmapped]), index=names[unmapped].to_numpy())
        .groupby(level=0)
        .sum()
        .sort_values(ascending=False)
    )
    top = ", ".join(f"{name} ({mass:.4g} kg)" for name, mass in by_name.head(5).items())
    logger.info(
        "FAIR characterization skipped {} of {} inventory flows that no species "
        "in the map covers, {:.4g} kg in total. Largest: {}.",
        int(by_name.size),
        int(names.nunique()),
        float(by_name.sum()),
        top,
    )


def _aggregate_inventory(
    df: pd.DataFrame, years: np.ndarray
) -> Tuple[List[str], np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Aggregate a dynamic inventory to per-(flow, activity) yearly emissions.

    Returns ``(species_names, species_of_pair, flow_ids, activity_ids,
    yearly)`` where ``yearly`` is an ``(n_pairs, n_years)`` array of signed
    emissions in kg and ``species_of_pair`` indexes ``species_names``.
    """
    n_years = len(years)

    # Species resolution happens once per distinct flow name, not per row.
    names = _flow_names(df)
    cas_by_name = _flow_cas(df, names)
    name_codes, unique_names = pd.factorize(names, sort=False)
    resolved = [
        species_map.resolve_flow(str(name), cas_by_name.get(name))
        for name in unique_names
    ]
    species_names = sorted({sp for sp, _, _ in resolved if sp is not None})
    species_index = {sp: i for i, sp in enumerate(species_names)}
    name_species = np.array(
        [species_index.get(sp, -1) if sp is not None else -1 for sp, _, _ in resolved],
        dtype="int64",
    )
    # The sign carries the mass-basis conversion (e.g. NO reported as NO2).
    name_sign = np.array([sign * factor for _, sign, factor in resolved], dtype="float64")

    row_species = name_species[name_codes]
    row_sign = name_sign[name_codes]
    _warn_about_unmapped(df, names, row_species >= 0)
    year_of_row = df["date"].dt.year.to_numpy()
    year_index = year_of_row - int(years[0])
    keep = (row_species >= 0) & (year_index >= 0) & (year_index < n_years)
    if not keep.any():
        return species_names, np.empty(0, "int64"), np.empty(0), np.empty(0), None

    flow_codes, flow_uniques = pd.factorize(df["flow"], sort=False)
    activity_codes, activity_uniques = pd.factorize(df["activity"], sort=False)
    n_activities = len(activity_uniques)
    pair_keys = flow_codes.astype("int64") * n_activities + activity_codes
    pair_codes, pair_uniques = pd.factorize(pair_keys[keep], sort=False)
    n_pairs = len(pair_uniques)

    signed = row_sign[keep] * df["amount"].to_numpy(dtype="float64")[keep]
    yearly = np.bincount(
        pair_codes * n_years + year_index[keep],
        weights=signed,
        minlength=n_pairs * n_years,
    ).reshape(n_pairs, n_years)

    species_of_pair = np.zeros(n_pairs, dtype="int64")
    species_of_pair[pair_codes] = row_species[keep]  # constant within a pair
    flow_ids = np.asarray(flow_uniques)[pair_uniques // n_activities]
    activity_ids = np.asarray(activity_uniques)[pair_uniques % n_activities]
    return species_names, species_of_pair, flow_ids, activity_ids, yearly


def _inventory_emissions_by_species(
    df: pd.DataFrame, years: np.ndarray
) -> Tuple[Dict[str, np.ndarray], list]:
    """Sum signed inventory emissions per FAIR species per year.

    Returns ``(by_species, flow_records)`` where ``flow_records`` is a list of
    dicts ``{flow, activity, species, signed}``. Kept as a readable view on
    :func:`_aggregate_inventory`, which is what the characterization uses.
    """
    species_names, species_of_pair, flow_ids, activity_ids, yearly = (
        _aggregate_inventory(df, years)
    )
    if yearly is None:
        return {}, []
    by_species = np.zeros((len(species_names), len(years)), dtype="float64")
    np.add.at(by_species, species_of_pair, yearly)
    records = [
        {
            "flow": flow_ids[i],
            "activity": activity_ids[i],
            "species": species_names[species_of_pair[i]],
            "signed": yearly[i],
        }
        for i in range(len(flow_ids))
    ]
    return {name: by_species[i] for i, name in enumerate(species_names)}, records


def _species_responses(
    marker: str,
    species_names: List[str],
    by_species: np.ndarray,
    years: np.ndarray,
    output: str,
) -> Dict[str, np.ndarray]:
    """ΔRF/ΔT per species, as ``(n_configs, n_years)`` ensemble arrays.

    Positive and negative parts of a species' emissions are perturbed
    separately (uptake is not the mirror image of a release), and all
    perturbations plus the baseline go into one batched, cached FAIR call.
    The runner scales small perturbations up to keep the difference out of the
    float64 noise floor, so each response is divided by that factor again.
    """
    perturbations: List[Optional[Dict[str, np.ndarray]]] = [None]  # baseline first
    owners: List[Tuple[int, str]] = []
    for i, name in enumerate(species_names):
        delta = by_species[i]
        if np.abs(delta).max() < _NEGLIGIBLE_KG:
            continue
        for part in (np.maximum(delta, 0.0), np.minimum(delta, 0.0)):
            if not np.any(part):
                continue
            perturbations.append({name: part})
            owners.append((i, name))

    runs = runner.run_perturbations(marker, years, perturbations)
    baseline = runs[0][output]

    responses = {
        name: np.zeros_like(baseline) for name in species_names
    }  # (cfg, n_years)
    for (_, name), run in zip(owners, runs[1:]):
        responses[name] = responses[name] + (run[output] - baseline) / run["scale"]
    return responses


def _attribute_to_flows(
    responses: Dict[str, np.ndarray],
    species_names: List[str],
    species_of_pair: np.ndarray,
    by_species: np.ndarray,
    yearly: np.ndarray,
    quantiles,
) -> np.ndarray:
    """Allocate species responses to flows and take ensemble quantiles.

    Returns an ``(n_pairs, n_quantiles, n_years)`` array. A flow's share is a
    scalar per year, so the quantiles are taken once per species on the per-kg
    response and then scaled - with the quantile order flipped where a flow's
    cumulative emission is negative (uptake), since ``percentile(c * x)`` is
    ``c * percentile_{100-q}(x)`` for ``c < 0``.
    """
    qs = list(quantiles)
    flipped = [100.0 - q for q in qs]
    out = np.zeros((yearly.shape[0], len(qs), yearly.shape[1]), dtype="float64")
    flow_cum = np.cumsum(yearly, axis=1)

    for i, name in enumerate(species_names):
        rows = np.flatnonzero(species_of_pair == i)
        if rows.size == 0:
            continue
        species_cum = np.cumsum(by_species[i])
        denominator = np.where(species_cum == 0, np.nan, species_cum)
        per_kg = responses[name] / denominator[None, :]  # (cfg, n_years)
        q_pos = allocation.safe_nanpercentile(per_kg, qs, axis=0)  # (nq, n_years)
        q_neg = allocation.safe_nanpercentile(per_kg, flipped, axis=0)
        cum = flow_cum[rows][:, None, :]  # (n_rows, 1, n_years)
        out[rows] = np.where(cum >= 0, q_pos[None, :, :], q_neg[None, :, :]) * cum
    return np.nan_to_num(out, nan=0.0)


def _to_long_dataframe(
    values: np.ndarray,
    flow_ids: np.ndarray,
    activity_ids: np.ndarray,
    years: np.ndarray,
    quantiles,
) -> pd.DataFrame:
    """Long frame of (date, amount, flow, activity, quantile), zeros dropped."""
    n_pairs, n_quantiles, n_years = values.shape
    nonzero = np.flatnonzero(values.reshape(-1) != 0)
    if nonzero.size == 0:
        return pd.DataFrame(columns=_OUT_COLUMNS)

    pair_idx, quantile_idx, year_idx = np.unravel_index(
        nonzero, (n_pairs, n_quantiles, n_years)
    )
    dates = np.asarray(
        [np.datetime64(f"{int(year)}-01-01", "s") for year in years]
    )
    out = pd.DataFrame(
        {
            "date": dates[year_idx],
            "amount": values.reshape(-1)[nonzero],
            "flow": flow_ids[pair_idx],
            "activity": activity_ids[pair_idx],
            "quantile": np.asarray(list(quantiles), dtype="float64")[quantile_idx],
        }
    )
    return (
        out.astype({"date": "datetime64[s]", "amount": "float64"})
        .sort_values(by=["quantile", "date", "flow"])
        .reset_index(drop=True)
    )


def characterize_with_fair(
    dynamic_inventory_df: pd.DataFrame,
    output: str = "radiative_forcing",
    quantiles=(2.5, 25, 50, 75, 97.5),
    time_horizon: Optional[int] = None,
    workers: Optional[int] = None,  # kept for API compatibility; runs are batched
) -> pd.DataFrame:
    """Run FAIR and return ΔRF/ΔT per (year, flow, activity, quantile)."""
    # Guard: empty inventory must return early, before require_fair()
    EMPTY = pd.DataFrame(columns=_OUT_COLUMNS)
    if dynamic_inventory_df.empty:
        return EMPTY

    years = _emission_years(dynamic_inventory_df, time_horizon)
    (
        species_names,
        species_of_pair,
        flow_ids,
        activity_ids,
        yearly,
    ) = _aggregate_inventory(dynamic_inventory_df, years)
    if yearly is None or not species_names:
        logger.warning(
            "No inventory flows could be mapped to FAIR species. Ensure the "
            "inventory carries flow names (a 'flow_name' column) and that the "
            "species map covers them."
        )
        return EMPTY

    runner.require_fair()  # fail fast with a clear message
    marker = config.current_fair_marker()

    by_species = np.zeros((len(species_names), len(years)), dtype="float64")
    np.add.at(by_species, species_of_pair, yearly)

    responses = _species_responses(marker, species_names, by_species, years, output)
    values = _attribute_to_flows(
        responses, species_names, species_of_pair, by_species, yearly, quantiles
    )
    return _to_long_dataframe(values, flow_ids, activity_ids, years, quantiles)
