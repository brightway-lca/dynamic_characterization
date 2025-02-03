import json
import os
import warnings
from collections.abc import Collection
from datetime import datetime
from typing import Callable, Dict, Tuple

import bw2data as bd
import numpy as np
import pandas as pd
from bw2data.utils import UnknownObject

from dynamic_characterization.classes import CharacterizedRow
from dynamic_characterization.ipcc_ar6.radiative_forcing import (
    characterize_ch4,
    characterize_co,
    characterize_co2,
    characterize_co2_uptake,
    characterize_n2o,
    create_generic_characterization_function,
)

def characterize(
    dynamic_inventory_df: pd.DataFrame,
    metric: str = "radiative_forcing",
    characterization_functions: Dict[int, Callable] = None,
    base_lcia_method: Tuple[str, ...] = None,
    time_horizon: int = 100,
    fixed_time_horizon: bool = False,
    time_horizon_start: datetime = datetime.now(),
    characterization_function_co2: Callable = None,
) -> pd.DataFrame:
    """
    Fully vectorized characterization function that applies characterization functions efficiently.
    """

    if metric not in {"radiative_forcing", "GWP"}:
        raise ValueError(f"Metric must be 'radiative_forcing' or 'GWP', not {metric}")

    if characterization_functions is None:
        warnings.warn(
            "No custom dynamic characterization functions provided. Using default functions."
        )
        if base_lcia_method is None:
            raise ValueError("Please provide an LCIA method for default characterization functions.")
        characterization_functions = create_characterization_functions_from_method(base_lcia_method)

    if metric == "GWP" and characterization_function_co2 is None:
        characterization_function_co2 = characterize_co2

    # Convert DataFrame to structured NumPy array
    structured_array = dynamic_inventory_df.to_records(index=False)

    # Prepare storage
    all_dates, all_amounts, all_flows, all_activities = [], [], [], []

    # Process each flow separately
    for flow_id, char_func in characterization_functions.items():
        # Select relevant rows (batch processing)
        mask = structured_array.flow == flow_id
        if not np.any(mask):
            continue

        selected_rows = structured_array[mask]  # Keep all relevant rows

        # Compute dynamic time horizons for all selected emissions at once
        dynamic_time_horizons = _calculate_dynamic_time_horizon_vectorized(
            emission_dates=selected_rows.date,
            time_horizon_start=time_horizon_start,
            time_horizon=time_horizon,
            fixed_time_horizon=fixed_time_horizon,
        )

        # Apply characterization function to all rows in batch
        if metric == "radiative_forcing":
            characterized_results = [
                _characterize_radiative_forcing(
                    characterization_functions, date, amount, flow, activity, period
                ) 
                for date, amount, flow, activity, period in zip(
                    selected_rows.date, selected_rows.amount, selected_rows.flow,
                    selected_rows.activity, dynamic_time_horizons
                )
            ]

        elif metric == "GWP":
            characterized_results = [
                _characterize_gwp(
                    characterization_functions=characterization_functions,
                    date=date,
                    amount=amount,
                    flow=flow,
                    activity=activity,
                    original_time_horizon=time_horizon,
                    dynamic_time_horizon=period,
                    characterization_function_co2=characterization_function_co2,
                ) 
                for date, amount, flow, activity, period in zip(
                    selected_rows.date, selected_rows.amount, selected_rows.flow,
                    selected_rows.activity, dynamic_time_horizons
                )
            ]

        # Store results efficiently (unpacking tuples correctly)
        all_dates.append(np.concatenate([res[0] for res in characterized_results]))
        all_amounts.append(np.concatenate([res[1] for res in characterized_results]))
        all_flows.append(np.concatenate([np.full(len(res[0]), res[2]) for res in characterized_results]))
        all_activities.append(np.concatenate([np.full(len(res[0]), res[3]) for res in characterized_results]))

    if not all_dates:
        raise ValueError("No flows to characterize. Check time horizon and available characterization functions.")

    # Convert results to a DataFrame
    characterized_inventory = pd.DataFrame({
        "date": np.concatenate(all_dates),
        "amount": np.concatenate(all_amounts),
        "flow": np.concatenate(all_flows),
        "activity": np.concatenate(all_activities),
    })

    return (
        characterized_inventory
        .astype({"date": "datetime64[s]", "amount": "float64"})
        .query("amount != 0")
        .sort_values(by=["date", "amount"])
        .reset_index(drop=True)
    )


def _calculate_dynamic_time_horizon_vectorized(emission_dates, time_horizon_start, time_horizon, fixed_time_horizon):
    """
    Fully vectorized computation of dynamic time horizons.
    """
    if fixed_time_horizon:
        # Compute Levasseur approach in batch
        end_time_horizon = time_horizon_start + pd.DateOffset(years=time_horizon)
        end_time_horizon = np.datetime64(end_time_horizon)

        # Compute horizon lengths in years (vectorized)
        delta_years = (end_time_horizon - emission_dates) / np.timedelta64(1, "Y")
        return np.maximum(0, np.round(delta_years)).astype(int)

    else:
        # Conventional approach: All emissions get the same time_horizon
        return np.full_like(emission_dates, time_horizon, dtype=int)


def create_characterization_functions_from_method(
    base_lcia_method: Tuple[str, ...]
) -> dict:
    """
    Add default dynamic characterization functions for CO2, CH4, N2O and other GHGs, based on IPCC AR6 Chapter 7 decay curves.

    Please note: Currently, only CO2, CH4 and N2O include climate-carbon feedbacks.

    This has not yet been added for other GHGs. Refer to https://esd.copernicus.org/articles/8/235/2017/esd-8-235-2017.html"
    "Methane, non-fossil" is currently also excluded from the default characterization functions, as it has a different static CF than fossil methane and we need to check the correct value (#TODO)

    Parameters
    ----------
    base_lcia_method : tuple
        Tuple of the selected the LCIA method, e.g. `("EF v3.1", "climate change", "global warming potential (GWP100)")`.

    Returns
    -------
    None but adds default dynamic characterization functions to the `characterization_functions` attribute of the DynamicCharacterization object.

    """

    characterization_functions = dict()

    filepath = os.path.join(
        os.path.dirname(os.path.abspath(__file__)),
        "ipcc_ar6",
        "data",
        "decay_multipliers.json",
    )

    with open(filepath) as json_file:
        decay_multipliers = json.load(json_file)

    # look up which GHGs are characterized in the selected static LCA method
    method_data = bd.Method(base_lcia_method).load()

    biosphere_db = bd.Database(bd.config.biosphere)

    # the bioflow-identifier stored in the method data can be the database id or the tuple (database, code)
    def get_bioflow_node(identifier):
        if (
            isinstance(identifier, Collection) and len(identifier) == 2
        ):  # is (probably) tuple of (database, code)
            try:
                biosphere_node = biosphere_db.get(
                    database=identifier[0], code=identifier[1]
                )
            except UnknownObject as e:
                raise UnknownObject(
                    f"Failed to set up the default characterization functions because a biosphere node was not found. Make sure the biosphere is set up correctly or provide a characterization_functions. Original error: {e}"
                )
            return biosphere_node

        elif isinstance(identifier, int):  # id is an int
            return biosphere_db.get(id=identifier)
        else:
            raise ValueError(
                "The flow-identifier stored in the selected method is neither an id nor the tuple (database, code). No automatic matching possible."
            )

    bioflow_nodes = set(get_bioflow_node(identifier) for identifier, _ in method_data)

    for node in bioflow_nodes:
        if "carbon dioxide" in node["name"].lower():
            if "soil" in node.get("categories", []):
                characterization_functions[node.id] = (
                    characterize_co2_uptake  # negative emission because uptake by soil
                )

            else:
                characterization_functions[node.id] = characterize_co2

        elif (
            "methane, fossil" in node["name"].lower()
            or "methane, from soil or biomass stock" in node["name"].lower()
        ):
            # TODO Check why "methane, non-fossil" has a CF of 27 instead of 29.8, currently excluded
            characterization_functions[node.id] = characterize_ch4

        elif "dinitrogen monoxide" in node["name"].lower():
            characterization_functions[node.id] = characterize_n2o

        elif "carbon monoxide" in node["name"].lower():
            characterization_functions[node.id] = characterize_co

        else:
            cas_number = node.get("CAS number")
            if cas_number:
                decay_series = decay_multipliers.get(cas_number)
                if decay_series is not None:
                    characterization_functions[node.id] = (
                        create_generic_characterization_function(np.array(decay_series))
                    )
    return characterization_functions


def _characterize_radiative_forcing(
    characterization_functions, date, amount, flow, activity, time_horizon
) -> CharacterizedRow:
    return characterization_functions[flow](date, amount, flow, activity, time_horizon)


def _characterize_gwp(
    characterization_functions,
    date, amount, flow, activity,
    original_time_horizon,
    dynamic_time_horizon,
    characterization_function_co2,
) -> CharacterizedRow:
    _, radiative_forcing_ghg, _, _ = characterization_functions[flow](
        date, amount, flow, activity,
        dynamic_time_horizon,
    )

    # calculate reference radiative forcing for 1 kg of CO2
    _, radiative_forcing_co2, _, _ = characterization_function_co2(
        date, 1, flow, activity, original_time_horizon
    )

    ghg_integral = radiative_forcing_ghg.sum()
    co2_integral = radiative_forcing_co2.sum()
    co2_equiv = ghg_integral / co2_integral

    return date, co2_equiv, flow, activity