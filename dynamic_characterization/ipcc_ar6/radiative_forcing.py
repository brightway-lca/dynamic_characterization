from functools import lru_cache

import numpy as np
import pandas as pd

from dynamic_characterization.classes import CharacterizedRow

from . import ar6_core, gases


def IRF_co2(year) -> callable:
    """
    Impulse Resonse Function (IRF) of CO2

    Parameters
    ----------
    year : int
        The year after emission for which the IRF is calculated.

    Returns
    -------
    float
        The IRF value for the given year.

    """
    alpha_0, alpha_1, alpha_2, alpha_3 = 0.2173, 0.2240, 0.2824, 0.2763
    tau_1, tau_2, tau_3 = 394.4, 36.54, 4.304
    exponentials = lambda year, alpha, tau: alpha * tau * (1 - np.exp(-year / tau))
    return (
        alpha_0 * year
        + exponentials(year, alpha_1, tau_1)
        + exponentials(year, alpha_2, tau_2)
        + exponentials(year, alpha_3, tau_3)
    )


# The decay-multiplier arrays below depend only on `period`, not on the emission
# row, so they were previously recomputed for every row of the inventory. They are
# memoized here and returned read-only; callers only ever multiply by them.
@lru_cache(maxsize=None)
def _co2_decay_multipliers(period: int) -> np.ndarray:
    radiative_efficiency_ppb = 1.33e-5
    M_co2 = 44.01
    M_air = 28.97
    m_atmosphere = 5.135e18
    radiative_efficiency_kg = (
        radiative_efficiency_ppb * M_air / M_co2 * 1e9 / m_atmosphere
    )
    arr = np.array(
        [radiative_efficiency_kg * IRF_co2(year) for year in range(period)]
    )
    arr.setflags(write=False)
    return arr


@lru_cache(maxsize=None)
def _co_decay_multipliers(period: int) -> np.ndarray:
    radiative_efficiency_ppb = 1.33e-5
    M_co2 = 44.01
    M_co = 28.01
    M_air = 28.97
    m_atmosphere = 5.135e18
    radiative_efficiency_kg = (
        radiative_efficiency_ppb * M_air / M_co2 * 1e9 / m_atmosphere
    )
    arr = np.array(
        [
            M_co2 / M_co * radiative_efficiency_kg * IRF_co2(year)
            for year in range(period)
        ]
    )
    arr.setflags(write=False)
    return arr


@lru_cache(maxsize=None)
def _ch4_decay_multipliers(period: int) -> np.ndarray:
    radiative_efficiency_ppb = 5.7e-4
    M_ch4 = 16.04
    M_air = 28.97
    m_atmosphere = 5.135e18
    radiative_efficiency_kg = (
        radiative_efficiency_ppb * M_air / M_ch4 * 1e9 / m_atmosphere
    )
    tau = 11.8
    arr = np.array(
        [
            radiative_efficiency_kg * tau * (1 - np.exp(-year / tau))
            for year in range(period)
        ]
    )
    arr.setflags(write=False)
    return arr


@lru_cache(maxsize=None)
def _n2o_decay_multipliers(period: int) -> np.ndarray:
    radiative_efficiency_ppb = 2.8e-3
    M_n2o = 44.01
    M_air = 28.97
    m_atmosphere = 5.135e18
    radiative_efficiency_kg = (
        radiative_efficiency_ppb * M_air / M_n2o * 1e9 / m_atmosphere
    )
    tau = 109
    arr = np.array(
        [
            radiative_efficiency_kg * tau * (1 - np.exp(-year / tau))
            for year in range(period)
        ]
    )
    arr.setflags(write=False)
    return arr


def characterize_co2(
    series,
    period: int | None = 100,
    cumulative: bool | None = False,
) -> CharacterizedRow:
    """
    Calculate the cumulative or marginal radiative forcing (CRF) from CO2 for each year in a given period.

    Based on characterize_co2 from bw_temporalis, but updated numerical values from IPCC AR6 Ch7 & SM.

    If `cumulative` is True, the cumulative CRF is calculated. If `cumulative` is False, the marginal CRF is calculated.
    Takes a single row of the TimeSeries Pandas DataFrame (corresponding to a set of (`date`/`amount`/`flow`/`activity`).
    For each year in the given period, the CRF is calculated.
    Units are watts/square meter/kilogram of CO2.

    Returns
    -------
    A CharacterizedRow object (namedtuple) with the following fields:
    - date: datetime64[s]
    - amount: float
    - flow: str
    - activity: str

    See also
    --------
    Joos2013: Relevant scientific publication on CRF: https://doi.org/10.5194/acp-13-2793-2013
    Schivley2015: Relevant scientific publication on the numerical calculation of CRF: https://doi.org/10.1021/acs.est.5b01118
    Forster2023: Updated numerical values from IPCC AR6 Chapter 7 (Table 7.15): https://doi.org/10.1017/9781009157896.009
    """

    date_beginning: np.datetime64 = series.date.to_numpy()
    dates_characterized: np.ndarray = date_beginning + np.arange(
        start=0, stop=period, dtype="timedelta64[Y]"
    ).astype("timedelta64[s]")

    decay_multipliers = _co2_decay_multipliers(period)

    forcing = np.array(series.amount * decay_multipliers, dtype="float64")

    if not cumulative:
        forcing = np.diff(forcing, prepend=0)

    return CharacterizedRow(
        date=np.array(dates_characterized, dtype="datetime64[s]"),
        amount=forcing,
        flow=series.flow,
        activity=series.activity,
    )


def characterize_co2_uptake(
    series,
    period: int | None = 100,
    cumulative: bool | None = False,
) -> CharacterizedRow:
    """
    The same as characterize_co2, but with a negative sign for uptake of CO2.

    Based on characterize_co2 from bw_temporalis, but updated numerical values from IPCC AR6 Ch7 & SM.

    Calculate the negative cumulative or marginal radiative forcing (CRF) from CO2-uptake for each year in a given period.

    If `cumulative` is True, the cumulative CRF is calculated. If `cumulative` is False, the marginal CRF is calculated.
    Takes a single row of the TimeSeries Pandas DataFrame (corresponding to a set of (`date`/`amount`/`flow`/`activity`).
    For each year in the given period, the CRF is calculated.
    Units are watts/square meter/kilogram of CO2.

    Returns
    -------
    A CharacterizedRow object (namedtuple) with the following fields:
    - date: datetime64[s]
    - amount: float
    - flow: str
    - activity: str

    See also
    --------
    Joos2013: Relevant scientific publication on CRF: https://doi.org/10.5194/acp-13-2793-2013
    Schivley2015: Relevant scientific publication on the numerical calculation of CRF: https://doi.org/10.1021/acs.est.5b01118
    Forster2023: Updated numerical values from IPCC AR6 Chapter 7 (Table 7.15): https://doi.org/10.1017/9781009157896.009
    """

    date_beginning: np.datetime64 = series.date.to_numpy()
    dates_characterized: np.ndarray = date_beginning + np.arange(
        start=0, stop=period, dtype="timedelta64[Y]"
    ).astype("timedelta64[s]")

    decay_multipliers = _co2_decay_multipliers(period)

    forcing = np.array(series.amount * decay_multipliers, dtype="float64")

    # flip the sign of the characterization function for CO2 uptake and not release
    forcing = -forcing

    if not cumulative:
        forcing = np.diff(forcing, prepend=0)

    return CharacterizedRow(
        date=np.array(dates_characterized, dtype="datetime64[s]"),
        amount=forcing,
        flow=series.flow,
        activity=series.activity,
    )


def characterize_co(
    series,
    period: int | None = 100,
    cumulative: bool | None = False,
) -> CharacterizedRow:
    """
    Calculate the cumulative or marginal radiative forcing (CRF) from CO for each year in a given period.

    This is exactly the same function as for CO2, it's just scaled by the ratio of molar masses of CO and CO2. This is because CO is very short-lived (lifetime ~2 months) and we assume that it completely reacts to CO2 within the first year.

    Only this oxidation to CO2 is modelled. The indirect effects of CO on ozone and methane, which
    make up most of its GWP in the literature (AR5 Table 8.A.4), are not included.

    Based on characterize_co2 from bw_temporalis, but updated numerical values from IPCC AR6 Ch7 & SM.

    Calculate the cumulative or marginal radiative forcing (CRF) from CO2 for each year in a given period.

    If `cumulative` is True, the cumulative CRF is calculated. If `cumulative` is False, the marginal CRF is calculated.
    Takes a single row of the TimeSeries Pandas DataFrame (corresponding to a set of (`date`/`amount`/`flow`/`activity`).
    For each year in the given period, the CRF is calculated.
    Units are watts/square meter/kilogram of CO2.

    Returns
    -------
    A CharacterizedRow object (namedtuple) with the following fields:
    - date: datetime64[s]
    - amount: float
    - flow: str
    - activity: str

    See also
    --------
    Joos2013: Relevant scientific publication on CRF: https://doi.org/10.5194/acp-13-2793-2013
    Schivley2015: Relevant scientific publication on the numerical calculation of CRF: https://doi.org/10.1021/acs.est.5b01118
    Forster2023: Updated numerical values from IPCC AR6 Chapter 7 (Table 7.15): https://doi.org/10.1017/9781009157896.009
    """

    date_beginning: np.datetime64 = series.date.to_numpy()
    dates_characterized: np.ndarray = date_beginning + np.arange(
        start=0, stop=period, dtype="timedelta64[Y]"
    ).astype("timedelta64[s]")

    decay_multipliers = _co_decay_multipliers(period)

    forcing = np.array(series.amount * decay_multipliers, dtype="float64")

    if not cumulative:
        forcing = np.diff(forcing, prepend=0)

    return CharacterizedRow(
        date=np.array(dates_characterized, dtype="datetime64[s]"),
        amount=forcing,
        flow=series.flow,
        activity=series.activity,
    )


def characterize_ch4(
    series,
    period: int = 100,
    cumulative=False,
) -> CharacterizedRow:
    """
    Calculate the cumulative or marginal radiative forcing from CH4 for each year in a given period.

    Based on characterize_methane from bw_temporalis, but updated numerical values from IPCC AR6 
    Ch7 & SM.

    This DOES include indirect effects of CH4 on ozone and water vapor, but DOES NOT include the
    decay to CO2. Therefore, this also treats fossil and biogenic methane the same. It also DOES
    NOT include the climate-carbon feedback that AR6 adds to its metrics. For more info,
    see IPCC AR6 Ch7 under "Carbon Cycle Responses and Other Indirect Contributions" on page 1014.

    If `cumulative` is True, the cumulative CRF is calculated. If `cumulative` is False, the 
    marginal CRF is calculated.

    Takes a single row of the TimeSeries Pandas DataFrame (corresponding to a set of 
    (`date`/`amount`/`flow`/`activity`). For each year in the given period, the CRF is calculated.
    Units are watts/square meter/kilogram of CH4.

    Parameters
    ----------
    series : array-like
        A single row of the TimeSeries dataframe.
    period : int, optional
        Time period for calculation (number of years), by default 100
    cumulative : bool, optional
        Should the RF amounts be summed over time?

    Returns
    -------
    A CharacterizedRow object (namedtuple) with the following fields:
    - date: datetime64[s]
    - amount: float
    - flow: str
    - activity: str

    See also
    --------
    Joos2013: Relevant scientific publication on CRF: https://doi.org/10.5194/acp-13-2793-2013
    Schivley2015: Relevant scientific publication on the numerical calculation of CRF: https://doi.org/10.1021/acs.est.5b01118
    Forster2023: Updated numerical values from IPCC AR6 Chapter 7 (Table 7.15): https://doi.org/10.1017/9781009157896.009
    """

    date_beginning: np.datetime64 = series.date.to_numpy()
    dates_characterized: np.ndarray = date_beginning + np.arange(
        start=0, stop=period, dtype="timedelta64[Y]"
    ).astype("timedelta64[s]")

    decay_multipliers = _ch4_decay_multipliers(period)

    forcing = np.array(series.amount * decay_multipliers, dtype="float64")

    if not cumulative:
        forcing = np.diff(forcing, prepend=0)

    return CharacterizedRow(
        date=np.array(dates_characterized, dtype="datetime64[s]"),
        amount=forcing,
        flow=series.flow,
        activity=series.activity,
    )


def characterize_n2o(
    series,
    period: int = 100,
    cumulative=False,
) -> CharacterizedRow:
    """
    Calculate the cumulative or marginal radiative forcing (CRF) from N2O for each year in a given period.

    Based on characterize_methane from bw_temporalis, but updated numerical values from IPCC AR6 Ch7 & SM.

    If `cumulative` is True, the cumulative CRF is calculated. If `cumulative` is False, the marginal CRF is calculated.
    Takes a single row of the TimeSeries Pandas DataFrame (corresponding to a set of (`date`/`amount`/`flow`/`activity`).
    For earch year in the given period, the CRF is calculated.
    Units are watts/square meter/kilogram of N2O.

    Parameters
    ----------
    series : array-like
        A single row of the TimeSeries dataframe.
    period : int, optional
        Time period for calculation (number of years), by default 100
    cumulative : bool, optional
        Should the RF amounts be summed over time?

    Returns
    -------
    A CharacterizedRow object (namedtuple) with the following fields:
    - date: datetime64[s]
    - amount: float
    - flow: str
    - activity: str

    See also
    --------
    Joos2013: Relevant scientific publication on CRF: https://doi.org/10.5194/acp-13-2793-2013
    Schivley2015: Relevant scientific publication on the numerical calculation of CRF: https://doi.org/10.1021/acs.est.5b01118
    Forster2023: Updated numerical values from IPCC AR6 Chapter 7 (Table 7.15): https://doi.org/10.1017/9781009157896.009
    """

    date_beginning: np.datetime64 = series.date.to_numpy()
    dates_characterized: np.ndarray = date_beginning + np.arange(
        start=0, stop=period, dtype="timedelta64[Y]"
    ).astype("timedelta64[s]")

    decay_multipliers = _n2o_decay_multipliers(period)

    forcing = np.array(series.amount * decay_multipliers, dtype="float64")
    if not cumulative:
        forcing = np.diff(forcing, prepend=0)

    return CharacterizedRow(
        date=np.array(dates_characterized, dtype="datetime64[s]"),
        amount=forcing,
        flow=series.flow,
        activity=series.activity,
    )


def create_generic_characterization_function(decay_series) -> CharacterizedRow:
    """
    Creates a characterization function for a GHG based on a decay series, by calling the nested method `characterize_generic()`.

    Parameters
    ----------
    decay_series : np.ndarray
        A decay series for a specific GHG. This is retrieved from `../data/decay_multipliers.pkl`

    Returns
    -------
    A function called `characterize_generic`, which in turn returns a TimeSeries dataframe that contains the forcing of the emission of the row over the given period based on the decay series of that biosphere flow.

    """

    def characterize_generic(
        series,
        period: int = 100,
        cumulative=False,
    ) -> CharacterizedRow:
        """
        Uses lookup generated in /dev/calculate_metrics.ipynb
        Data originates from https://doi.org/10.1029/2019RG000691

        Parameters
        ----------
        series : array-like
            A single row of the dynamic inventory dataframe.
        period : int, optional
            Time period for calculation (number of years), by default 100
        cumulative : bool,
            cumulative impact

        Returns
        -------
        A TimeSeries dataframe that contains the forcing of the point emission from the row for each year in the given period.
          date: datetime64[s]
          amount: float (forcing at this timestep)
          flow: str
          activity: str

        See also
        --------
        Joos2013: Relevant scientific publication on CRF: https://doi.org/10.5194/acp-13-2793-2013
        Schivley2015: Relevant scientific publication on the numerical calculation of CRF: https://doi.org/10.1021/acs.est.5b01118
        Forster2023: Updated numerical values from IPCC AR6 Chapter 7 (Table 7.15): https://doi.org/10.1017/9781009157896.009

        """

        date_beginning: np.datetime64 = series.date.to_numpy()

        dates_characterized: np.ndarray = date_beginning + np.arange(
            start=0, stop=period, dtype="timedelta64[Y]"
        ).astype("timedelta64[s]")

        decay_multipliers = decay_series[:period]

        forcing = np.array(series.amount * decay_multipliers, dtype="float64")

        if not cumulative:
            forcing = np.diff(forcing, prepend=0)

        return CharacterizedRow(
            date=np.array(dates_characterized, dtype="datetime64[s]"),
            amount=forcing,
            flow=series.flow,
            activity=series.activity,
        )

    return characterize_generic


# --- "ar6" parameter set ---------------------------------------------------------------------


@lru_cache(maxsize=None)
def make_characterization_function(spec):
    """
    Characterization function for a pulse response from ``ar6_core`` (a ``GasSpec``).

    The returned function has the same signature and output as the legacy functions in this
    module, ``(series, period=100, cumulative=False) -> CharacterizedRow``, with one
    difference in the time convention: the value stamped at ``date + k`` years is the forcing
    integrated over the year *starting* then, so the annual values add up to ``AGWP(period)``
    (the legacy functions add up to ``AGWP(period - 1)``).

    The function also exposes ``.spec``, ``.agwp(time_horizon)`` and ``.agtp(time_horizon)``
    (per kg). It is memoized per spec, so the same spec always gives the same function object.
    """
    def characterize_ar6(series, period: int = 100, cumulative: bool = False) -> CharacterizedRow:
        date_beginning: np.datetime64 = series.date.to_numpy()
        dates_characterized: np.ndarray = date_beginning + np.arange(
            start=0, stop=period, dtype="timedelta64[Y]"
        ).astype("timedelta64[s]")

        forcing = np.array(series.amount * ar6_core.annual_cumulative(spec, period), dtype="float64")
        if not cumulative:
            forcing = np.diff(forcing, prepend=0)

        return CharacterizedRow(
            date=np.array(dates_characterized, dtype="datetime64[s]"),
            amount=forcing,
            flow=series.flow,
            activity=series.activity,
        )

    characterize_ar6.__name__ = f"characterize_ar6[{spec.name}]"
    characterize_ar6.__qualname__ = characterize_ar6.__name__
    characterize_ar6.spec = spec
    characterize_ar6.agwp = lambda time_horizon: float(ar6_core.agwp(spec, float(time_horizon)))
    characterize_ar6.agtp = lambda time_horizon: float(ar6_core.agtp(spec, float(time_horizon)))
    return characterize_ar6


def _co2_ar6_function():
    return make_characterization_function(gases.co2())


def characterize_co2_ar6(series, period: int = 100, cumulative: bool = False) -> CharacterizedRow:
    """
    Radiative forcing of CO2 with the ``"ar6"`` parameter set (AR6 Table 7.SM.7 CO2 response).

    Differs from :func:`characterize_co2` only in the radiative efficiency (Meinshausen et al.
    2020 at 2019 concentrations, +5 % adjustment) and in the time convention: the annual values
    add up to ``AGWP(period)``. Used as the GWP reference of the ``"ar6"`` parameter set.
    """
    return _co2_ar6_function()(series, period, cumulative=cumulative)


def characterize_co2_uptake_ar6(series, period: int = 100, cumulative: bool = False) -> CharacterizedRow:
    """As :func:`characterize_co2_ar6`, with a negative sign for CO2 uptake."""
    return make_characterization_function(gases.co2_uptake())(series, period, cumulative=cumulative)
