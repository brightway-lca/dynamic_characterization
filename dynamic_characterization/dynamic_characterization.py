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
    characterization_function_dict: Dict[int, Callable] = None,
    base_lcia_method: Tuple[str, ...] = None,
    time_horizon: int = 100,
    fixed_time_horizon: bool = False,
    time_horizon_start: datetime = datetime.now(),
    characterization_function_co2: Callable = None,
) -> pd.DataFrame:
    """
    Characterizes the dynamic inventory, formatted as a Dataframe, by evaluating each emission (row in DataFrame) using given dynamic characterization functions.

    Available metrics are radiative forcing [W/m2] and GWP [kg CO2eq], defaulting to `radiative_forcing`.

    In case users don't provide own dynamic characterization functions, it adds dynamic characterization functions from the timex submodule
    for the GHGs mentioned in the IPCC AR6 Chapter 7, if these GHG are also characterized in the selected static LCA method.

    This method applies the dynamic characterization functions to each row of the dynamic_inventory_df for the duration of `time_horizon`, defaulting to 100 years.
    The `fixed_time_horizon` parameter determines whether the evaluation time horizon for all emissions is calculated from the
    functional unit (`fixed_time_horizon=True`), regardless of when the actual emission occurs, or from the time of the emission itself(`fixed_time_horizon=False`).
    The former is the implementation of the Levasseur approach (https://doi.org/10.1021/es9030003), while the latter is how conventional LCA is done.
    The Levasseur approach means that earlier emissions are characterized for a longer time period than later emissions.

    Parameters
    ----------
    dynamic_inventory_df : pd.DataFrame
        Dynamic inventory, formatted as a DataFrame, which contains the timing, id and amount of emissions and the emitting activity.
    metric : str, optional
        The metric for which the dynamic LCIA should be calculated. Default is "GWP". Available: "GWP" and "radiative_forcing". Default is "radiative_forcing".
    characterization_function_dict : dict, optional
        A dictionary of the form {biosphere_flow_id: dynamic_characterization_function} allowing users to specify their own functions and what flows to apply them to.
        Default is none, in which case a set of default functions are added based on the base_lcia_method.
    base_lcia_method : tuple, optional
        Tuple of the selcted the LCIA method, e.g. `("EF v3.1", "climate change", "global warming potential (GWP100)")`. This is
        required for adding the default characterization functions and can be kept empty if custom ones are provided.
    time_horizon: int, optional
        Length of the time horizon for the dynamic characterization. Default is 100 years.
    fixed_time_horizon: bool, optional
        If True, the time horizon is calculated from the time of the functional unit (FU) instead of the time of emission. Default is False.
    time_horizon_start: pd.Timestamp, optional
        The starting timestamp of the time horizon for the dynamic characterization. Only needed for fixed time horizons. Default is datetime.now().
    characterization_function_co2: Callable, optional
        Characterization function for CO2. This is required for the GWP calculation. If None is given, we try using timex' default CO2 function.

    Returns
    -------
    pd.DataFrame
        characterized dynamic inventory
    """

    if metric not in {"radiative_forcing", "GWP"}:
        raise ValueError(
            f"Metric must be either 'radiative_forcing' or 'GWP', not {metric}"
        )

    if not characterization_function_dict:
        warnings.warn(
            "No custom dynamic characterization functions provided. Using default dynamic characterization functions.\
                The flows that are characterized are based on the selection of the initially chosen impact category.\
                You can look up the mapping in the bw_timex.dynamic_characterizer.characterization_function_dict."
        )
        if not base_lcia_method:
            raise ValueError(
                "Please provide an LCIA method to base the default dynamic characterization functions on."
            )
        characterization_function_dict = (
            create_characterization_function_dict_from_method(base_lcia_method)
        )

    if metric == "GWP" and not characterization_function_co2:
        characterization_function_co2 = characterize_co2

    characterized_inventory_data = []

    for row in dynamic_inventory_df.itertuples(index=False):

        # skip uncharacterized biosphere flows
        if row.flow not in characterization_function_dict.keys():
            continue

        dynamic_time_horizon = _calculate_dynamic_time_horizon(
            emission_date=row.date,
            time_horizon_start=time_horizon_start,
            time_horizon=time_horizon,
            fixed_time_horizon=fixed_time_horizon,
        )

        if metric == "radiative_forcing":  # radiative forcing in W/m2
            characterized_inventory_data.append(
                _characterize_radiative_forcing(
                    characterization_function_dict, row, dynamic_time_horizon
                )
            )

        if metric == "GWP":  # scale radiative forcing to GWP [kg CO2 equivalent]
            characterized_inventory_data.append(
                _characterize_gwp(
                    characterization_function_dict=characterization_function_dict,
                    row=row,
                    original_time_horizon=time_horizon,
                    dynamic_time_horizon=dynamic_time_horizon,
                    characterization_function_co2=characterization_function_co2,
                )
            )

    if not characterized_inventory_data:
        raise ValueError(
            "There are no flows to characterize. Please make sure your time horizon matches the timing of emissions and make sure there are characterization functions for the flows in the dynamic inventories."
        )

    characterized_inventory = (
        pd.DataFrame(characterized_inventory_data)
        .explode(["amount", "date"])
        .astype({"date": "datetime64[s]", "amount": "float64"})
        .query("amount != 0")[["date", "amount", "flow", "activity"]]
        .sort_values(by=["date", "amount"])
        .reset_index(drop=True)
    )

    return characterized_inventory


def create_characterization_function_dict_from_method(
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
    None but adds default dynamic characterization functions to the `characterization_function_dict` attribute of the DynamicCharacterization object.

    """

    characterization_function_dict = dict()

    filepath = os.path.join(
        os.path.dirname(os.path.abspath(__file__)),
        "timex",
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
                    f"Failed to set up the default characterization functions because a biosphere node was not found. Make sure the biosphere is set up correctly or provide a characterization_function_dict. Original error: {e}"
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
                characterization_function_dict[node.id] = (
                    characterize_co2_uptake  # negative emission because uptake by soil
                )

            else:
                characterization_function_dict[node.id] = characterize_co2

        elif (
            "methane, fossil" in node["name"].lower()
            or "methane, from soil or biomass stock" in node["name"].lower()
        ):
            # TODO Check why "methane, non-fossil" has a CF of 27 instead of 29.8, currently excluded
            characterization_function_dict[node.id] = characterize_ch4

        elif "dinitrogen monoxide" in node["name"].lower():
            characterization_function_dict[node.id] = characterize_n2o

        elif "carbon monoxide" in node["name"].lower():
            characterization_function_dict[node.id] = characterize_co

        else:
            cas_number = node.get("CAS number")
            if cas_number:
                decay_series = decay_multipliers.get(cas_number)
                if decay_series is not None:
                    characterization_function_dict[node.id] = (
                        create_generic_characterization_function(np.array(decay_series))
                    )
    return characterization_function_dict


def _calculate_dynamic_time_horizon(
    emission_date: pd.Timestamp,
    time_horizon_start: pd.Timestamp,
    time_horizon: int,
    fixed_time_horizon: bool,
) -> datetime:
    """
    Calculate the dynamic time horizon for the dynamic characterization of an emission.
    Distinguishes between the Levasseur approach (fixed_time_horizon = True) and the conventional approach (fixed_time_horizon = False).

    Parameters
    ----------
    emission_date: pd.Timestamp
        The date of the emission
    time_horizon_start: pd.Timestamp
        Start timestamp of the time horizon
    time_horizon: int
        Length of the time horizon in years
    fixed_time_horizon: bool
        If True, the time horizon is calculated from the time of the functional unit (FU) instead of the time of emission

    Returns
    -------
    datetime.datetime
        dynamic time horizon for the specific emission
    """
    if fixed_time_horizon:
        # Levasseur approach: time_horizon for all emissions starts at timing of FU + time_horizon
        # e.g. an emission occuring n years before FU is characterized for time_horizon+n years
        end_time_horizon = time_horizon_start + pd.DateOffset(years=time_horizon)
        emission_datetime = emission_date.to_pydatetime()

        return max(0, round((end_time_horizon - emission_datetime).days / 365.25))

    else:
        # conventional approach, emission is calculated from t emission for the length of time horizon
        return time_horizon


def _characterize_radiative_forcing(
    characterization_function_dict, row, time_horizon
) -> CharacterizedRow:
    return characterization_function_dict[row.flow](row, time_horizon)


def _characterize_gwp(
    characterization_function_dict,
    row,
    original_time_horizon,
    dynamic_time_horizon,
    characterization_function_co2,
) -> CharacterizedRow:
    radiative_forcing_ghg = characterization_function_dict[row.flow](
        row,
        dynamic_time_horizon,
    )

    # calculate reference radiative forcing for 1 kg of CO2
    radiative_forcing_co2 = characterization_function_co2(
        row._replace(amount=1), original_time_horizon
    )

    ghg_integral = radiative_forcing_ghg.amount.sum()
    co2_integral = radiative_forcing_co2.amount.sum()
    co2_equiv = ghg_integral / co2_integral

    return CharacterizedRow(
        date=row.date,
        amount=co2_equiv,
        flow=row.flow,
        activity=row.activity,
    )
