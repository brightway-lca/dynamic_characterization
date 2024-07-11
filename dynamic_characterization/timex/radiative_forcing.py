from functools import cache
import numpy as np
import pandas as pd

M_air = 28.97  # g/mol, dry air
m_atmosphere = 5.135e18  # kg [Trenberth and Smith, 2005]

# Cache dictionary to store precomputed decay multipliers
_decay_multipliers_co2_cache = {}
_decay_multipliers_ch4_cache = {}
_decay_multipliers_n2o_cache = {}

@cache
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

@cache
def radiative_efficiency_per_kg_co2():
    radiative_efficiency_ppb = 1.33e-5  # W/m2/ppb; 2019 background CO2 concentration; IPCC AR6 Table 7.15
    M_co2 = 44.01  # g/mol
    return radiative_efficiency_ppb * M_air / M_co2 * 1e9 / m_atmosphere  # W/m2/kg-CO2

@cache
def radiative_efficiency_per_kg_ch4():
    radiative_efficiency_ppb = 5.7e-4  # # W/m2/ppb; 2019 background cch4 concentration; IPCC AR6 Table 7.15. This number includes indirect effects.
    M_ch4 = 16.04  # g/mol
    return radiative_efficiency_ppb * M_air / M_ch4 * 1e9 / m_atmosphere # W/m2/kg-CH4
    
@cache
def radiative_efficiency_per_kg_n2o():
    radiative_efficiency_ppb = 2.8e-3  # # W/m2/ppb; 2019 background cch4 concentration; IPCC AR6 Table 7.15. This number includes indirect effects.
    M_n2o = 44.01  # g/mol
    return radiative_efficiency_ppb * M_air / M_n2o * 1e9 / m_atmosphere # W/m2/kg-N2O

@cache
def calculate_decay_multipliers_co2(period):
    global _decay_multipliers_co2_cache

    if period in _decay_multipliers_co2_cache:
        return _decay_multipliers_co2_cache[period]

    # Find the longest precomputed period that is shorter than the requested period
    max_cached_period = max((p for p in _decay_multipliers_co2_cache if p < period), default=0)
    multipliers = list(_decay_multipliers_co2_cache.get(max_cached_period, []))

    radiative_efficiency_kg = radiative_efficiency_per_kg_co2()

    # Compute any additional multipliers needed
    for year in range(max_cached_period, period):
        multipliers.append(radiative_efficiency_kg * IRF_co2(year))

    # Cache the result for future use
    _decay_multipliers_co2_cache[period] = np.array(multipliers)
    
    return _decay_multipliers_co2_cache[period]

@cache
def calculate_decay_multipliers_ch4(period):
    global _decay_multipliers_ch4_cache
    tau = 11.8  # Lifetime (years)

    if period in _decay_multipliers_ch4_cache:
        return _decay_multipliers_ch4_cache[period]

    # Find the longest precomputed period that is shorter than the requested period
    max_cached_period = max((p for p in _decay_multipliers_ch4_cache if p < period), default=0)
    multipliers = list(_decay_multipliers_ch4_cache.get(max_cached_period, []))

    radiative_efficiency_kg = radiative_efficiency_per_kg_ch4()

    # Compute any additional multipliers needed
    for year in range(max_cached_period, period):
        multipliers.append(radiative_efficiency_kg * tau * (1 - np.exp(-year / tau)))

    # Cache the result for future use
    _decay_multipliers_ch4_cache[period] = np.array(multipliers)
    
    return _decay_multipliers_ch4_cache[period]

@cache
def calculate_decay_multipliers_n2o(period):
    global _decay_multipliers_n2o_cache
    tau = 109  # Lifetime (years)

    if period in _decay_multipliers_n2o_cache:
        return _decay_multipliers_n2o_cache[period]

    # Find the longest precomputed period that is shorter than the requested period
    max_cached_period = max((p for p in _decay_multipliers_n2o_cache if p < period), default=0)
    multipliers = list(_decay_multipliers_n2o_cache.get(max_cached_period, []))

    radiative_efficiency_kg = radiative_efficiency_per_kg_n2o()

    # Compute any additional multipliers needed
    for year in range(max_cached_period, period):
        multipliers.append(radiative_efficiency_kg * tau * (1 - np.exp(-year / tau)))

    # Cache the result for future use
    _decay_multipliers_n2o_cache[period] = np.array(multipliers)
    
    return _decay_multipliers_n2o_cache[period]

@cache
def characterize_co2(
    row: tuple,
    period: int | None = 100,
    cumulative: bool | None = False,
) -> pd.DataFrame:
    """
    Calculate the cumulative or marginal radiative forcing (CRF) from CO2 for each year in a given period.
    
    Based on characterize_co2 from bw_temporalis, but updated numerical values from IPCC AR6 Ch7 & SM.

    If `cumulative` is True, the cumulative CRF is calculated. If `cumulative` is False, the marginal CRF is calculated.
    Takes a single row of the TimeSeries Pandas DataFrame (corresponding to a set of (`date`/`amount`/`flow`/`activity`).
    For each year in the given period, the CRF is calculated.
    Units are watts/square meter/kilogram of CO2.

    Returns
    -------
    A TimeSeries dataframe with the following columns:
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
    date, amount, flow, activity = row

    date_beginning: np.datetime64 = date.to_numpy()
    date_characterized: np.ndarray = date_beginning + np.arange(
        start=0, stop=period, dtype="timedelta64[Y]"
    ).astype("timedelta64[s]")
    
    decay_multipliers = calculate_decay_multipliers_co2(period)

    forcing = pd.Series(data=amount * decay_multipliers, dtype="float64")

    if not cumulative:
        forcing = forcing.diff(periods=1).fillna(0)

    return pd.DataFrame(
        {
            "date": pd.Series(data=date_characterized, dtype="datetime64[s]"),
            "amount": forcing,
            "flow": flow,
            "activity": activity,
        }
    )

@cache
def characterize_co2_uptake(
    row,
    period: int | None = 100,
    cumulative: bool | None = False,
) -> pd.DataFrame:
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
    A TimeSeries dataframe with the following columns:
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
    date, amount, flow, activity = row

    date_beginning: np.datetime64 = date.to_numpy()
    date_characterized: np.ndarray = date_beginning + np.arange(
        start=0, stop=period, dtype="timedelta64[Y]"
    ).astype("timedelta64[s]")

    decay_multipliers: np.ndarray = calculate_decay_multipliers_co2(period)

    forcing = pd.Series(data=amount * decay_multipliers, dtype="float64")

    forcing = (
        -forcing
    )  # flip the sign of the characterization function for CO2 uptake and not release

    if not cumulative:
        forcing = forcing.diff(periods=1).fillna(0)

    return pd.DataFrame(
        {
            "date": pd.Series(data=date_characterized, dtype="datetime64[s]"),
            "amount": forcing,
            "flow": flow,
            "activity": activity,
        }
    )

@cache
def characterize_co(
    row,
    period: int | None = 100,
    cumulative: bool | None = False,
) -> pd.DataFrame:
    """
    Calculate the cumulative or marginal radiative forcing (CRF) from CO for each year in a given period.

    This is exactly the same function as for CO2, it's just scaled by the ratio of molar masses of CO and CO2. This is because CO is very short-lived (lifetime ~2 months) and we assume that it completely reacts to CO2 within the first year.

    Based on characterize_co2 from bw_temporalis, but updated numerical values from IPCC AR6 Ch7 & SM.

    Calculate the cumulative or marginal radiative forcing (CRF) from CO2 for each year in a given period.

    If `cumulative` is True, the cumulative CRF is calculated. If `cumulative` is False, the marginal CRF is calculated.
    Takes a single row of the TimeSeries Pandas DataFrame (corresponding to a set of (`date`/`amount`/`flow`/`activity`).
    For each year in the given period, the CRF is calculated.
    Units are watts/square meter/kilogram of CO2.

    Returns
    -------
    A TimeSeries dataframe with the following columns:
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
    date, amount, flow, activity = row

    M_co2 = 44.01  # g/mol
    M_co = 28.01   # g/mol

    date_beginning: np.datetime64 = date.to_numpy()
    date_characterized: np.ndarray = date_beginning + np.arange(
        start=0, stop=period, dtype="timedelta64[Y]"
    ).astype("timedelta64[s]")
    
    decay_multipliers: np.ndarray = M_co2 / M_co * calculate_decay_multipliers_co2(period)

    forcing = pd.Series(data=amount * decay_multipliers, dtype="float64")

    if not cumulative:
        forcing = forcing.diff(periods=1).fillna(0)

    return pd.DataFrame(
        {
            "date": pd.Series(data=date_characterized, dtype="datetime64[s]"),
            "amount": forcing,
            "flow": flow,
            "activity": activity,
        }
    )

@cache
def characterize_ch4(
    row,
    period: int = 100,
    cumulative=False,
) -> pd.DataFrame:
    """
    Calculate the cumulative or marginal radiative forcing (CRF) from CH4 for each year in a given period. 

    Based on characterize_methane from bw_temporalis, but updated numerical values from IPCC AR6 Ch7 & SM.
    
    This DOES include indirect effects of CH4 on ozone and water vapor, but DOES NOT include the decay to CO2. 
    For more info on that, see the deprecated version of bw_temporalis.

    If `cumulative` is True, the cumulative CRF is calculated. If `cumulative` is False, the marginal CRF is calculated.
    Takes a single row of the TimeSeries Pandas DataFrame (corresponding to a set of (`date`/`amount`/`flow`/`activity`).
    For earch year in the given period, the CRF is calculated.
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
    A TimeSeries dataframe with the following columns:
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
    date, amount, flow, activity = row

    date_beginning: np.datetime64 = date.to_numpy()
    date_characterized: np.ndarray = date_beginning + np.arange(
        start=0, stop=period, dtype="timedelta64[Y]"
    ).astype("timedelta64[s]")

    decay_multipliers = calculate_decay_multipliers_ch4(period)

    forcing = pd.Series(data=amount * decay_multipliers, dtype="float64")
    
    if not cumulative:
        forcing = forcing.diff(periods=1).fillna(0)

    return pd.DataFrame(
        {
            "date": pd.Series(data=date_characterized, dtype="datetime64[s]"),
            "amount": forcing,
            "flow": flow,
            "activity": activity,
        }
    )

@cache
def characterize_n2o(
    row,
    period: int = 100,
    cumulative=False,
) -> pd.DataFrame:
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
    A TimeSeries dataframe with the following columns:
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
    date, amount, flow, activity = row

    date_beginning: np.datetime64 = date.to_numpy()
    date_characterized: np.ndarray = date_beginning + np.arange(
        start=0, stop=period, dtype="timedelta64[Y]"
    ).astype("timedelta64[s]")

    decay_multipliers = calculate_decay_multipliers_n2o(period)

    forcing = pd.Series(data=amount * decay_multipliers, dtype="float64")
    
    if not cumulative:
        forcing = forcing.diff(periods=1).fillna(0)

    return pd.DataFrame(
        {
            "date": pd.Series(data=date_characterized, dtype="datetime64[s]"),
            "amount": forcing,
            "flow": flow,
            "activity": activity,
        }
    )


def create_generic_characterization_function(decay_series) -> pd.DataFrame:
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
    @cache
    def characterize_generic(
        row,
        period: int = 100,
        cumulative=False,
    ) -> pd.DataFrame:
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
        date, amount, flow, activity = row

        date_beginning: np.datetime64 = date.to_numpy()

        dates_characterized: np.ndarray = date_beginning + np.arange(
            start=0, stop=period, dtype="timedelta64[Y]"
        ).astype("timedelta64[s]")

        decay_multipliers = decay_series[:period]

        forcing = pd.Series(data=amount * decay_multipliers, dtype="float64")

        if not cumulative:
            forcing = forcing.diff(periods=1).fillna(0)

        return pd.DataFrame(
            {
                "date": pd.Series(data=dates_characterized, dtype="datetime64[s]"),
                "amount": forcing,
                "flow": flow,
                "activity": activity,
            }
        )

    return characterize_generic
