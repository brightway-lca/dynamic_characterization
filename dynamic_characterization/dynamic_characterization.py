import json
import os
import warnings
from collections.abc import Collection
from datetime import datetime
from typing import Callable, Dict, Tuple
from loguru import logger

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
from dynamic_characterization.prospective import agwp, agtp
from dynamic_characterization.prospective.radiative_forcing import (
    characterize_ch4 as prospective_characterize_ch4,
    characterize_co2 as prospective_characterize_co2,
    characterize_co2_uptake as prospective_characterize_co2_uptake,
    characterize_n2o as prospective_characterize_n2o,
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
    time_varying_re: bool = False,
    fallback_to_ipcc: bool = True,
) -> pd.DataFrame:
    """
    Characterizes the dynamic inventory, formatted as a Dataframe, by evaluating each emission (row in DataFrame) using given dynamic characterization functions.

    Available metrics are radiative forcing [W/m2] and GWP [kg CO2eq], defaulting to `radiative_forcing`.
    Additional prospective metrics pGWP, pGTP, and prospective_radiative_forcing use Watanabe et al. (2026)
    scenario-based characterization factors. For these, set the scenario first via prospective.set_scenario().

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
        The metric for which the dynamic LCIA should be calculated. Available: "radiative_forcing", "GWP", "pGWP", "pGTP", "prospective_radiative_forcing". Default is "radiative_forcing".
        - "radiative_forcing": Returns W/m2 time series using IPCC AR6 functions
        - "GWP": Returns kg CO2eq using IPCC AR6 functions
        - "pGWP": Returns kg CO2eq using Watanabe et al. (2026) prospective GWP (requires scenario)
        - "pGTP": Returns kg CO2eq using Watanabe et al. (2026) prospective GTP (requires scenario)
        - "prospective_radiative_forcing": Returns W/m2 time series using Watanabe et al. (2026) functions (requires scenario)
    characterization_functions : dict, optional
        A dictionary of the form {biosphere_flow_id: dynamic_characterization_function} allowing users to specify their own functions and what flows to apply them to.
        Default is none, in which case a set of default functions are added based on the base_lcia_method.
    base_lcia_method : tuple, optional
        Tuple of the selected the LCIA method, e.g. `("EF v3.1", "climate change", "global warming potential (GWP100)")`. This is
        required for adding the default characterization functions and can be kept empty if custom ones are provided.
    time_horizon: int, optional
        Length of the time horizon for the dynamic characterization. Default is 100 years.
    fixed_time_horizon: bool, optional
        If True, the time horizon is calculated from the time of the functional unit (FU) instead of the time of emission. Default is False.
    time_horizon_start: pd.Timestamp, optional
        The starting timestamp of the time horizon for the dynamic characterization. Only needed for fixed time horizons. Default is datetime.now().
    characterization_function_co2: Callable, optional
        Characterization function for CO2. This is required for the GWP calculation. If None is given, we try using timex' default CO2 function.
    time_varying_re: bool, optional
        Only for prospective metrics (pGWP/pGTP/prospective_radiative_forcing). If True, use radiative efficiency that evolves over the decay period.
        If False (default), use fixed RE from emission year (IPCC standard approach).
    fallback_to_ipcc: bool, optional
        Only for prospective metrics. If True (default), use IPCC AR6 characterization functions for GHGs not available
        in the Watanabe module (e.g., CO and other GHGs). If False, only GHGs available in Watanabe (CO2, CH4, N2O) are characterized.

    Returns
    -------
    pd.DataFrame
        characterized dynamic inventory
    """

    valid_metrics = {"radiative_forcing", "GWP", "pGWP", "pGTP", "prospective_radiative_forcing"}
    if metric not in valid_metrics:
        raise ValueError(
            f"Metric must be one of {valid_metrics}, not {metric}"
        )

    # For prospective metrics, use Watanabe characterization functions
    use_prospective = metric in {"pGWP", "pGTP", "prospective_radiative_forcing"}

    if not characterization_functions:
        logger.info(
            "No custom dynamic characterization functions provided. Using default dynamic \
            characterization functions. The flows that are characterized are based on the selection\
                of the initially chosen impact category."
        )
        if not base_lcia_method:
            raise ValueError(
                "Please provide an LCIA method to base the default dynamic characterization \
                functions on."
            )
        characterization_functions = create_characterization_functions_from_method(
            base_lcia_method, use_prospective=use_prospective, fallback_to_ipcc=fallback_to_ipcc
        )

    if metric == "GWP" and not characterization_function_co2:
        characterization_function_co2 = characterize_co2

    characterized_inventory_data = []

    for row in dynamic_inventory_df.itertuples(index=False):

        # skip uncharacterized biosphere flows
        if row.flow not in characterization_functions.keys():
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
                    characterization_functions, row, dynamic_time_horizon
                )
            )

        elif metric == "GWP":  # scale radiative forcing to GWP [kg CO2 equivalent]
            characterized_inventory_data.append(
                _characterize_gwp(
                    characterization_functions=characterization_functions,
                    row=row,
                    original_time_horizon=time_horizon,
                    dynamic_time_horizon=dynamic_time_horizon,
                    characterization_function_co2=characterization_function_co2,
                )
            )

        elif metric == "pGWP":  # prospective GWP using Watanabe
            characterized_inventory_data.append(
                _characterize_pgwp(
                    characterization_functions=characterization_functions,
                    row=row,
                    original_time_horizon=time_horizon,
                    dynamic_time_horizon=dynamic_time_horizon,
                    time_varying_re=time_varying_re,
                )
            )

        elif metric == "pGTP":  # prospective GTP using Watanabe
            characterized_inventory_data.append(
                _characterize_pgtp(
                    characterization_functions=characterization_functions,
                    row=row,
                    original_time_horizon=time_horizon,
                    dynamic_time_horizon=dynamic_time_horizon,
                    time_varying_re=time_varying_re,
                )
            )

        elif metric == "prospective_radiative_forcing":  # prospective radiative forcing using Watanabe
            characterized_inventory_data.append(
                _characterize_prospective_radiative_forcing(
                    characterization_functions=characterization_functions,
                    row=row,
                    time_horizon=dynamic_time_horizon,
                    time_varying_re=time_varying_re,
                )
            )

    if not characterized_inventory_data:
        logger.warning(
            "There are no flows to characterize. Please make sure your time horizon matches the "
            "timing of emissions and make sure there are characterization functions for the flows "
            "in the dynamic inventories."
        )
        return pd.DataFrame(columns=["date", "amount", "flow", "activity"])

    characterized_inventory = (
        pd.DataFrame(characterized_inventory_data)
        .explode(["amount", "date"])
        .astype({"date": "datetime64[s]", "amount": "float64"})
        .query("amount != 0")[["date", "amount", "flow", "activity"]]
        .sort_values(by=["date", "amount"])
        .reset_index(drop=True)
    )

    return characterized_inventory


def create_characterization_functions_from_method(
    base_lcia_method: Tuple[str, ...],
    use_prospective: bool = False,
    fallback_to_ipcc: bool = True,
) -> dict:
    """
    Add default dynamic characterization functions for CO2, CH4, N2O and other GHGs, based on IPCC
    AR6 Chapter 7 decay curves. If use_prospective=True, uses Watanabe et al. (2026) functions
    for CO2, CH4, and N2O (requires scenario to be set first).

    Please note: Currently, only CO2, CH4 and N2O include climate-carbon feedbacks.
    This has not yet been added for other GHGs.
    Refer to https://esd.copernicus.org/articles/8/235/2017/esd-8-235-2017.html for more info.
    "Methane, non-fossil" is currently also excluded from the default characterization functions,
    as it has a different static CF than fossil methane and we need to check the correct value.

    Parameters
    ----------
    base_lcia_method : tuple
        Tuple of the selected the LCIA method, e.g. `("EF v3.1", "climate change", "global warming potential (GWP100)")`.
    use_prospective : bool, optional
        If True, use Watanabe et al. (2026) prospective characterization functions for CO2, CH4, N2O.
        Default is False.
    fallback_to_ipcc : bool, optional
        Only relevant when use_prospective=True. If True (default), use IPCC AR6 characterization functions
        for GHGs not available in the Watanabe module (e.g., CO and other GHGs). If False, only GHGs
        available in Watanabe (CO2, CH4, N2O) are characterized.

    Returns
    -------
    dict
        Dictionary mapping biosphere flow IDs to characterization functions.

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
                    f"Failed to set up the default characterization functions because a biosphere \
                    node was not found. Make sure the biosphere is set up correctly or provide \
                    a correct mapping in 'characterization_functions'. Original error: {e}"
                ) from e
            return biosphere_node

        elif isinstance(identifier, int):  # id is an int
            return biosphere_db.get(id=identifier)
        else:
            raise ValueError(
                "The flow-identifier stored in the selected method is neither an id nor the tuple \
                (database, code). No automatic matching possible."
            )

    bioflow_nodes = set(get_bioflow_node(identifier) for identifier, _ in method_data)

    # Select characterization functions based on whether prospective metrics are used
    if use_prospective:
        co2_func = prospective_characterize_co2
        co2_uptake_func = prospective_characterize_co2_uptake
        ch4_func = prospective_characterize_ch4
        n2o_func = prospective_characterize_n2o
    else:
        co2_func = characterize_co2
        co2_uptake_func = characterize_co2_uptake
        ch4_func = characterize_ch4
        n2o_func = characterize_n2o

    for node in bioflow_nodes:
        if "carbon dioxide" in node["name"].lower():
            if "soil" in node.get("categories", []):
                characterization_functions[node.id] = (
                    co2_uptake_func  # negative emission because uptake by soil
                )
            elif "in air" in node.get("categories", []) and node.get("type", []) == 'natural resource':
                # CO2 as a natural resource in air is assumed to be used for uptake in CDR processes
                characterization_functions[node.id] = (
                    co2_uptake_func  # negative emission because uptake by CDR processes
                )
            else:
                characterization_functions[node.id] = co2_func

        elif (
            "methane, fossil" in node["name"].lower()
            or "methane, non-fossil" in node["name"].lower()
            or "methane, from soil or biomass stock" in node["name"].lower()
        ):
            characterization_functions[node.id] = ch4_func

        elif "dinitrogen monoxide" in node["name"].lower():
            characterization_functions[node.id] = n2o_func

        elif "carbon monoxide" in node["name"].lower():
            # CO is not available in Watanabe module, use IPCC if fallback is enabled
            if not use_prospective or fallback_to_ipcc:
                characterization_functions[node.id] = characterize_co

        else:
            # Other GHGs from decay_multipliers are not available in Watanabe module
            if not use_prospective or fallback_to_ipcc:
                cas_number = node.get("CAS number")
                if cas_number:
                    decay_series = decay_multipliers.get(cas_number)
                    if decay_series is not None:
                        characterization_functions[node.id] = (
                            create_generic_characterization_function(np.array(decay_series))
                        )
    return characterization_functions


def _calculate_dynamic_time_horizon(
    emission_date: pd.Timestamp,
    time_horizon_start: pd.Timestamp,
    time_horizon: int,
    fixed_time_horizon: bool,
) -> datetime:
    """
    Calculate the dynamic time horizon for the dynamic characterization of an emission.
    Distinguishes between the Levasseur approach (fixed_time_horizon = True) and the conventional
    approach (fixed_time_horizon = False).

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
    characterization_functions, row, time_horizon
) -> CharacterizedRow:
    return characterization_functions[row.flow](row, time_horizon)


def _characterize_gwp(
    characterization_functions,
    row,
    original_time_horizon,
    dynamic_time_horizon,
    characterization_function_co2,
) -> CharacterizedRow:
    radiative_forcing_ghg = characterization_functions[row.flow](
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


def _characterize_pgwp(
    characterization_functions,
    row,
    original_time_horizon,
    dynamic_time_horizon,
    time_varying_re: bool = False,
) -> CharacterizedRow:
    """
    Calculate prospective GWP using Watanabe et al. (2026) characterization.

    Uses AGWP_gas / AGWP_CO2 to calculate kg CO2 equivalent.
    Emission year is extracted from the row's date.
    """
    # Get emission year from the row date
    emission_date = row.date
    emission_year = int(str(emission_date.to_numpy())[:4])

    # Calculate AGWP for the gas using its characterization function
    radiative_forcing_ghg = characterization_functions[row.flow](
        row,
        dynamic_time_horizon,
        time_varying_re=time_varying_re,
    )

    # Calculate reference AGWP for 1 kg of CO2
    radiative_forcing_co2 = prospective_characterize_co2(
        row._replace(amount=1),
        original_time_horizon,
        time_varying_re=time_varying_re,
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


def _characterize_pgtp(
    characterization_functions,
    row,
    original_time_horizon,
    dynamic_time_horizon,
    time_varying_re: bool = False,
) -> CharacterizedRow:
    """
    Calculate prospective GTP using Watanabe et al. (2026) characterization.

    Uses AGTP_gas / AGTP_CO2 to calculate kg CO2 equivalent.
    Emission year is extracted from the row's date.
    """
    # Get emission year from the row date
    emission_date = row.date
    emission_year = int(str(emission_date.to_numpy())[:4])

    # For GTP, we need to use the AGTP functions directly since they
    # compute temperature change potential, not just integrated radiative forcing
    # First, determine which gas we're dealing with based on the characterization function
    char_func = characterization_functions[row.flow]

    # Get the AGTP for this gas
    if char_func in (prospective_characterize_co2, prospective_characterize_co2_uptake):
        agtp_gas = agtp.agtp_co2(
            emission_year=emission_year,
            time_horizon=dynamic_time_horizon,
            time_varying_re=time_varying_re,
        )
        # For CO2, pGTP = 1.0 by definition
        if char_func == prospective_characterize_co2_uptake:
            agtp_gas = -agtp_gas
    elif char_func == prospective_characterize_ch4:
        agtp_gas = agtp.agtp_ch4(
            emission_year=emission_year,
            time_horizon=dynamic_time_horizon,
            time_varying_re=time_varying_re,
        )
    elif char_func == prospective_characterize_n2o:
        agtp_gas = agtp.agtp_n2o(
            emission_year=emission_year,
            time_horizon=dynamic_time_horizon,
            time_varying_re=time_varying_re,
        )
    else:
        # For other GHGs, fall back to using integrated RF as proxy
        # This may not be as accurate but provides a fallback
        radiative_forcing_ghg = char_func(
            row,
            dynamic_time_horizon,
        )
        agtp_gas = radiative_forcing_ghg.amount.sum()

    # Calculate reference AGTP for 1 kg of CO2
    agtp_co2 = agtp.agtp_co2(
        emission_year=emission_year,
        time_horizon=original_time_horizon,
        time_varying_re=time_varying_re,
    )

    co2_equiv = row.amount * agtp_gas / agtp_co2

    return CharacterizedRow(
        date=row.date,
        amount=co2_equiv,
        flow=row.flow,
        activity=row.activity,
    )


def _characterize_prospective_radiative_forcing(
    characterization_functions,
    row,
    time_horizon,
    time_varying_re: bool = False,
) -> CharacterizedRow:
    """
    Calculate prospective radiative forcing using Watanabe et al. (2026) characterization.

    For GHGs available in Watanabe (CO2, CH4, N2O), uses scenario-based radiative efficiencies.
    For other GHGs (when fallback_to_ipcc=True), uses standard IPCC AR6 functions.
    """
    char_func = characterization_functions[row.flow]

    # Check if this is a Watanabe function (they accept time_varying_re parameter)
    prospective_functions = (
        prospective_characterize_co2,
        prospective_characterize_co2_uptake,
        prospective_characterize_ch4,
        prospective_characterize_n2o,
    )

    if char_func in prospective_functions:
        # Watanabe functions support time_varying_re
        return char_func(row, time_horizon, time_varying_re=time_varying_re)
    else:
        # IPCC fallback functions don't support time_varying_re
        return char_func(row, time_horizon)
