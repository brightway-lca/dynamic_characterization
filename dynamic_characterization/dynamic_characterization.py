import inspect
from collections import namedtuple
from collections.abc import Collection
from datetime import datetime
from functools import lru_cache
from typing import Callable, Dict, Tuple
from loguru import logger

import bw2data as bd
import numpy as np
import pandas as pd
from bw2data.utils import UnknownObject

from dynamic_characterization import ar6_prospective
from dynamic_characterization.ar6_dispatch import build_functions as build_ar6_functions
from dynamic_characterization.classes import CharacterizedRow
from dynamic_characterization.ipcc_ar6.radiative_forcing import characterize_co2
from dynamic_characterization.options import CharacterizationOptions
from dynamic_characterization.prospective import watanabe2026

# Building the default characterization functions issues one biosphere DB query per method flow.
# The result only depends on the BW project, the method, and the options, so it is memoized
# across calls.
_CHARACTERIZATION_FUNCTION_CACHE = {}


def clear_characterization_function_cache() -> None:
    """Clear the memoized default characterization functions (e.g. after editing a method)."""
    _CHARACTERIZATION_FUNCTION_CACHE.clear()


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
    characterize_biogenic_uptake: bool = True,
    characterization_options: dict = None,
) -> pd.DataFrame:
    """
    Characterizes the dynamic inventory, formatted as a Dataframe, by evaluating each emission (row in DataFrame) using given dynamic characterization functions.

    Available metrics are radiative forcing [W/m2] and GWP [kg CO2eq], defaulting to `radiative_forcing`.
    The prospective metrics pGWP, pGTP and prospective_radiative_forcing follow a future scenario
    (Barbosa Watanabe & Cherubini 2026); set it first with ``prospective.set_scenario()``.

    If no characterization functions are given, default ones are built for the flows of
    ``base_lcia_method``: the IPCC AR6 emission metrics (WG1 Ch. 7) for CO2, CH4, N2O and the
    halocarbons, incl. the climate-carbon feedback and the CO2 from fossil-CH4 oxidation, plus the
    indirect effects of CO, VOC and H2 (see ``characterization_options`` and the documentation).

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
        - "pGWP": Returns kg CO2eq, prospective GWP (requires scenario)
        - "pGTP": Returns kg CO2eq, prospective GTP (requires scenario)
        - "prospective_radiative_forcing": Returns W/m2 time series for the scenario (requires scenario)
        Prospective metrics cover at most 100 years; longer horizons are truncated with a warning.
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
        Only for prospective metrics (pGWP/pGTP/prospective_radiative_forcing). If True, the radiative efficiencies
        of CO2, CH4 and N2O follow the scenario over the decay period. If False (default), those of the emission
        year are used throughout (as in Watanabe & Cherubini 2026).
    fallback_to_ipcc: bool, optional
        Only for prospective metrics. If True (default), GHGs without scenario data (halocarbons, CO, VOC, H2, ...)
        are characterized with their IPCC AR6 / literature functions. If False, only CO2, CH4 and N2O are characterized.
    characterize_biogenic_uptake: bool, optional
        Whether to characterize biogenic uptake flows (e.g. CO2 uptake by soil or biomass, or CDR processes).
        If True, the default characterization functions include uptake flows and these are characterized with a negative emission profile. Default is True.
        Set to False if you want to explicitly model uptake outside of the dynamic characterization (e.g. from a forest model)
    characterization_options: dict, optional
        Options for the default characterization functions, see
        :class:`dynamic_characterization.options.CharacterizationOptions`, e.g.
        ``{"include_nox": True}``, ``{"climate_carbon_feedback": False}`` or, for the prospective
        metrics, ``{"prospective_method": "ar6"}``.

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

    options = CharacterizationOptions.from_value(characterization_options)

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
            base_lcia_method,
            use_prospective=use_prospective,
            fallback_to_ipcc=fallback_to_ipcc,
            characterize_uptake=characterize_biogenic_uptake,
            characterization_options=characterization_options,
        )

    if metric == "GWP" and not characterization_function_co2:
        characterization_function_co2 = characterize_co2

    # Fast path: for the (default) radiative_forcing metric without a fixed time
    # horizon, every row of a given flow uses the same characterization function
    # with the same period. The per-emission response is linear in the emission
    # amount and its shape is independent of the emission date, so we can evaluate
    # each characterization function once per flow and broadcast over all its rows
    # instead of iterating the inventory row by row.
    if metric == "radiative_forcing" and not fixed_time_horizon:
        vectorized = _vectorized_radiative_forcing(
            dynamic_inventory_df, characterization_functions, time_horizon
        )
        if vectorized is not None:
            return vectorized

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
                    options=options,
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
                    options=options,
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


_VectorizationRow = namedtuple(
    "_VectorizationRow", ["date", "amount", "flow", "activity"]
)


def _vectorized_radiative_forcing(
    dynamic_inventory_df: pd.DataFrame,
    characterization_functions: Dict[int, Callable],
    period: int,
):
    """
    Vectorized evaluation of the ``radiative_forcing`` metric without a fixed
    time horizon.

    Each characterization function is evaluated once per flow (in cumulative
    mode, giving the signed per-unit decay array) and then broadcast over all
    inventory rows of that flow. ``np.diff`` is applied to ``amount * decay``
    exactly as in the per-row implementation, so the result is bit-for-bit
    identical to the row-by-row loop.

    Returns the characterized inventory DataFrame, or ``None`` to signal that
    the caller should fall back to the generic row-by-row loop (nothing to
    characterize, a function without a ``cumulative`` keyword, or an unexpected
    response length).
    """
    flow_values = dynamic_inventory_df["flow"].to_numpy()
    known = np.fromiter(
        (f in characterization_functions for f in flow_values),
        dtype=bool,
        count=len(flow_values),
    )
    if not known.any():
        return None  # let the row loop emit the "nothing to characterize" warning

    df = dynamic_inventory_df.loc[known]
    dates_all = df["date"].to_numpy()
    amounts_all = df["amount"].to_numpy().astype("float64")
    flows_all = df["flow"].to_numpy()
    activities_all = df["activity"].to_numpy()

    # timedelta offsets identical to those built inside every IPCC AR6 / generic
    # characterization function
    offsets = np.arange(start=0, stop=period, dtype="timedelta64[Y]").astype(
        "timedelta64[s]"
    )

    date_blocks = []
    amount_blocks = []
    flow_blocks = []
    activity_blocks = []

    for flow_id in pd.unique(flows_all):
        func = characterization_functions[flow_id]

        # The fast path calls the function in cumulative mode to obtain the
        # signed per-unit decay array. If it does not accept `cumulative`, or
        # returns an unexpected length, fall back to the safe row loop.
        try:
            params = inspect.signature(func).parameters
        except (TypeError, ValueError):
            return None
        if "cumulative" not in params:
            return None

        mask = flows_all == flow_id
        sub_dates = dates_all[mask]
        sub_amounts = amounts_all[mask]
        sub_activities = activities_all[mask]

        sample = _VectorizationRow(
            date=pd.Timestamp(sub_dates[0]),
            amount=1.0,
            flow=flow_id,
            activity=sub_activities[0],
        )
        unit_cumulative = np.asarray(
            func(sample, period, cumulative=True).amount, dtype="float64"
        )
        if unit_cumulative.shape != (period,):
            return None

        # forcing[i, k] = amount_i * decay_k, then marginalized with np.diff,
        # exactly reproducing the per-row `np.diff(amount * decay, prepend=0)`
        forcing = sub_amounts[:, None] * unit_cumulative[None, :]
        forcing = np.diff(forcing, prepend=0, axis=1)

        n = sub_amounts.shape[0]
        date_blocks.append((sub_dates[:, None] + offsets[None, :]).reshape(-1))
        amount_blocks.append(forcing.reshape(-1))
        flow_blocks.append(np.full(n * period, flow_id))
        activity_blocks.append(np.repeat(sub_activities, period))

    characterized_inventory = (
        pd.DataFrame(
            {
                "date": np.concatenate(date_blocks),
                "amount": np.concatenate(amount_blocks),
                "flow": np.concatenate(flow_blocks),
                "activity": np.concatenate(activity_blocks),
            }
        )
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
    characterize_uptake: bool = True,
    characterization_options: dict = None,
) -> dict:
    """
    Default dynamic characterization functions for the flows of a static LCIA method.

    The functions reproduce the IPCC AR6 emission metrics (WG1 Ch. 7, Tables 7.15 and 7.SM.7):
    CO2, fossil and non-fossil CH4 (incl. oxidation CO2), N2O and all halocarbons of Table 7.SM.7,
    with the climate-carbon feedback; CO, VOC and H2 include their indirect effects from published
    parameter sets, NO/NOx optionally. Flows are matched by CAS number, then by exact name. With
    ``use_prospective=True`` they follow the scenario of ``prospective.set_scenario()`` (CO2, CH4
    and N2O after Barbosa Watanabe & Cherubini 2026). See the documentation for details.

    Parameters
    ----------
    base_lcia_method : tuple
        Tuple of the selected the LCIA method, e.g. `("EF v3.1", "climate change", "global warming potential (GWP100)")`.
    use_prospective : bool, optional
        If True, the functions follow the prospective scenario (metrics pGWP, pGTP,
        prospective_radiative_forcing). Default is False.
    fallback_to_ipcc : bool, optional
        Only relevant when use_prospective=True. If True (default), GHGs without scenario data are
        characterized with their IPCC AR6 / literature functions. If False, only CO2, CH4 and N2O are.
    characterize_uptake: bool, optional
        Whether to include characterization functions for biogenic uptake flows (e.g. CO2 uptake by soil or biomass, or CDR processes).
        If True (default), the default characterization functions include uptake. Set to False
        if you want to explicitly model uptake outside of the dynamic characterization (e.g. from a forest model)
    characterization_options: dict, optional
        Options for the default characterization functions, see
        :class:`dynamic_characterization.options.CharacterizationOptions`.

    Returns
    -------
    dict
        Dictionary mapping biosphere flow IDs to characterization functions.

    """

    options = CharacterizationOptions.from_value(characterization_options)

    cache_key = (
        bd.projects.current,
        tuple(base_lcia_method),
        use_prospective,
        fallback_to_ipcc,
        characterize_uptake,
        options,
    )
    cached = _CHARACTERIZATION_FUNCTION_CACHE.get(cache_key)
    if cached is not None:
        return cached

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

    characterization_functions = build_ar6_functions(
        bioflow_nodes,
        options,
        characterize_uptake=characterize_uptake,
        prospective=use_prospective,
        fallback_to_ipcc=fallback_to_ipcc,
    )
    _CHARACTERIZATION_FUNCTION_CACHE[cache_key] = characterization_functions
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


@lru_cache(maxsize=None)
def _accepts_time_varying_re(char_func) -> bool:
    try:
        return "time_varying_re" in inspect.signature(char_func).parameters
    except (TypeError, ValueError):
        return False


def _apply_characterization_function(
    char_func,
    row,
    time_horizon,
    time_varying_re: bool = False,
) -> CharacterizedRow:
    """Call a characterization function, passing ``time_varying_re`` to those that accept it."""
    if time_varying_re and _accepts_time_varying_re(char_func):
        return char_func(row, time_horizon, time_varying_re=time_varying_re)
    return char_func(row, time_horizon)


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
    options: CharacterizationOptions = None,
) -> CharacterizedRow:
    """
    Prospective GWP: the forcing of the emission summed over its (dynamic) time horizon, divided
    by that of 1 kg CO2 in the same scenario and year over the original time horizon.
    """
    options = options or CharacterizationOptions()
    radiative_forcing_ghg = _apply_characterization_function(
        characterization_functions[row.flow],
        row,
        dynamic_time_horizon,
        time_varying_re=time_varying_re,
    )
    radiative_forcing_co2 = ar6_prospective.co2_function(options.prospective_method)(
        row._replace(amount=1), original_time_horizon, time_varying_re=time_varying_re
    )
    co2_equiv = radiative_forcing_ghg.amount.sum() / radiative_forcing_co2.amount.sum()
    return CharacterizedRow(date=row.date, amount=co2_equiv, flow=row.flow, activity=row.activity)


def _characterize_pgtp(
    characterization_functions,
    row,
    original_time_horizon,
    dynamic_time_horizon,
    time_varying_re: bool = False,
    options: CharacterizationOptions = None,
) -> CharacterizedRow:
    """
    Prospective GTP: AGTP of the emission at its (dynamic) time horizon, divided by the AGTP of
    1 kg CO2 in the same scenario and year at the original time horizon.
    """
    options = options or CharacterizationOptions()
    char_func = characterization_functions[row.flow]
    reference = ar6_prospective.co2_function(options.prospective_method)
    if hasattr(char_func, "agtp_for"):
        agtp_gas = char_func.agtp_for(row.date, dynamic_time_horizon, time_varying_re=time_varying_re)
    else:
        agtp_gas = _agtp_from_radiative_forcing(char_func, row, dynamic_time_horizon, time_varying_re=time_varying_re)
    agtp_co2 = reference.agtp_for(row.date, original_time_horizon, time_varying_re=time_varying_re)
    co2_equiv = row.amount * agtp_gas / agtp_co2
    return CharacterizedRow(date=row.date, amount=co2_equiv, flow=row.flow, activity=row.activity)


def _agtp_from_radiative_forcing(
    char_func,
    row,
    time_horizon: int,
    time_varying_re: bool = False,
) -> float:
    """
    AGTP at ``time_horizon`` [K per unit emission] of a custom characterization function, for
    pGTP. The function is evaluated for a unit emission and its annual forcing is passed through
    the annual temperature model of ``prospective.watanabe2026`` (as for the CO2 reference). A
    function that provides its own ``agtp(time_horizon)`` uses that instead.
    """
    if time_horizon <= 0:
        return 0.0
    own_agtp = getattr(char_func, "agtp", None)
    if callable(own_agtp):
        return float(own_agtp(time_horizon))
    marginal = np.asarray(
        _apply_characterization_function(
            char_func, row._replace(amount=1), time_horizon, time_varying_re=time_varying_re
        ).amount,
        dtype="float64",
    )
    return float(watanabe2026.agtp_from_rf(marginal)[time_horizon - 1])


def _characterize_prospective_radiative_forcing(
    characterization_functions,
    row,
    time_horizon,
    time_varying_re: bool = False,
) -> CharacterizedRow:
    """Prospective radiative forcing of an emission (functions built with use_prospective=True)."""
    return _apply_characterization_function(
        characterization_functions[row.flow],
        row,
        time_horizon,
        time_varying_re=time_varying_re,
    )


def compare_with_static(
    base_lcia_method: Tuple[str, ...],
    time_horizon: int = 100,
    characterization_options: dict = None,
    characterize_uptake: bool = True,
) -> pd.DataFrame:
    """
    Compare the static characterization factors of a method with the dynamic GWP of its
    default characterization functions, flow by flow.

    The dynamic GWP of a flow is the forcing of a 1 kg pulse summed over ``time_horizon`` years,
    divided by that of 1 kg CO2 with the same parameter set (as the "GWP" metric of
    :func:`characterize` does).

    Parameters
    ----------
    base_lcia_method : tuple
        The static LCIA method, e.g. ``("ecoinvent-3.12", "IPCC 2021", "climate change",
        "global warming potential (GWP100)")``.
    time_horizon : int, optional
        Time horizon of the dynamic GWP in years. Default 100.
    characterization_options : dict, optional
        As for :func:`characterize`, e.g. ``{"include_nox": True}``.
    characterize_uptake : bool, optional
        As for :func:`create_characterization_functions_from_method`.

    Returns
    -------
    pd.DataFrame
        One row per method flow: ``flow`` (id), ``name``, ``categories``, ``static_cf``,
        ``dynamic_gwp`` (NaN if the flow has no default function) and ``ratio``
        (dynamic / static), sorted by name.
    """
    options = CharacterizationOptions.from_value(characterization_options)
    functions = create_characterization_functions_from_method(
        base_lcia_method,
        characterize_uptake=characterize_uptake,
        characterization_options=options,
    )
    co2_function = characterize_co2

    unit = _VectorizationRow(date=pd.Timestamp("2020-01-01"), amount=1.0, flow=None, activity=None)
    co2_integral = co2_function(unit, time_horizon).amount.sum()

    biosphere_db = bd.Database(bd.config.biosphere)
    rows = []
    for identifier, static_cf in bd.Method(base_lcia_method).load():
        if isinstance(identifier, int):
            node = biosphere_db.get(id=identifier)
        else:
            node = biosphere_db.get(database=identifier[0], code=identifier[1])
        func = functions.get(node.id)
        dynamic = (
            func(unit._replace(flow=node.id), time_horizon).amount.sum() / co2_integral
            if func is not None
            else np.nan
        )
        rows.append(
            {
                "flow": node.id,
                "name": node["name"],
                "categories": tuple(node.get("categories") or ()),
                "static_cf": float(static_cf),
                "dynamic_gwp": float(dynamic),
                "ratio": float(dynamic) / float(static_cf) if static_cf else np.nan,
            }
        )
    return pd.DataFrame(rows).sort_values(["name", "categories"]).reset_index(drop=True)
