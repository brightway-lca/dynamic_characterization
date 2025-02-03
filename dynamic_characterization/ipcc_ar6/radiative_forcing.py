import numpy as np
import pandas as pd

from dynamic_characterization.classes import CharacterizedRow


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


def characterize_co2(
    date,
    amount,
    flow,
    activity,
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

    # functional variables and units (from publications listed in docstring)
    radiative_efficiency_ppb = (
        1.33e-5  # W/m2/ppb; 2019 background co2 concentration; IPCC AR6 Table 7.15
    )

    # for conversion from ppb to kg-CO2
    M_co2 = 44.01  # g/mol
    M_air = 28.97  # g/mol, dry air
    m_atmosphere = 5.135e18  # kg [Trenberth and Smith, 2005]

    radiative_efficiency_kg = (
        radiative_efficiency_ppb * M_air / M_co2 * 1e9 / m_atmosphere
    )  # W/m2/kg-CO2

    dates_characterized: np.ndarray = date + np.arange(
        start=0, stop=period, dtype="timedelta64[Y]"
    ).astype("timedelta64[s]")

    decay_multipliers: np.ndarray = np.array(
        [radiative_efficiency_kg * IRF_co2(year) for year in range(period)]
    )

    forcing = np.array(amount * decay_multipliers, dtype="float64")

    if not cumulative:
        forcing = np.diff(forcing, prepend=0)

    return dates_characterized, forcing, flow, activity


def characterize_co2_uptake(
    date,
    amount,
    flow,
    activity,
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

    # functional variables and units (from publications listed in docstring)
    radiative_efficiency_ppb = (
        1.33e-5  # W/m2/ppb; 2019 background co2 concentration; IPCC AR6 Table 7.15
    )

    # for conversion from ppb to kg-CO2
    M_co2 = 44.01  # g/mol
    M_air = 28.97  # g/mol, dry air
    m_atmosphere = 5.135e18  # kg [Trenberth and Smith, 2005]

    radiative_efficiency_kg = (
        radiative_efficiency_ppb * M_air / M_co2 * 1e9 / m_atmosphere
    )  # W/m2/kg-CO2

    dates_characterized: np.ndarray = date + np.arange(
        start=0, stop=period, dtype="timedelta64[Y]"
    ).astype("timedelta64[s]")

    decay_multipliers: np.ndarray = np.array(
        [radiative_efficiency_kg * IRF_co2(year) for year in range(period)]
    )

    forcing = np.array(amount * decay_multipliers, dtype="float64")

    # flip the sign of the characterization function for CO2 uptake and not release
    forcing = -forcing

    if not cumulative:
        forcing = np.diff(forcing, prepend=0)

    return dates_characterized, forcing, flow, activity


def characterize_co(
    date,
    amount,
    flow,
    activity,
    period: int | None = 100,
    cumulative: bool | None = False,
) -> CharacterizedRow:
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

    # functional variables and units (from publications listed in docstring)
    radiative_efficiency_ppb = (
        1.33e-5  # W/m2/ppb; 2019 background co2 concentration; IPCC AR6 Table 7.15
    )

    # for conversion from ppb to kg-CO2
    M_co2 = 44.01  # g/mol
    M_co = 28.01  # g/mol
    M_air = 28.97  # g/mol, dry air
    m_atmosphere = 5.135e18  # kg [Trenberth and Smith, 2005]

    radiative_efficiency_kg = (
        radiative_efficiency_ppb * M_air / M_co2 * 1e9 / m_atmosphere
    )  # W/m2/kg-CO2

    dates_characterized: np.ndarray = date + np.arange(
        start=0, stop=period, dtype="timedelta64[Y]"
    ).astype("timedelta64[s]")

    decay_multipliers: np.ndarray = np.array(
        [
            M_co2 / M_co * radiative_efficiency_kg * IRF_co2(year)
            for year in range(period)
        ]  # <-- Scaling from co2 to co is done here
    )

    forcing = np.array(amount * decay_multipliers, dtype="float64")

    if not cumulative:
        forcing = np.diff(forcing, prepend=0)

    return dates_characterized, forcing, flow, activity


def characterize_ch4(
    date,
    amount,
    flow,
    activity,
    period: int = 100,
    cumulative=False,
) -> CharacterizedRow:
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

    # functional variables and units (from publications listed in docstring)
    radiative_efficiency_ppb = 5.7e-4  # # W/m2/ppb; 2019 background cch4 concentration; IPCC AR6 Table 7.15. This number includes indirect effects.

    # for conversion from ppb to kg-CH4
    M_ch4 = 16.04  # g/mol
    M_air = 28.97  # g/mol, dry air
    m_atmosphere = 5.135e18  # kg [Trenberth and Smith, 2005]

    radiative_efficiency_kg = (
        radiative_efficiency_ppb * M_air / M_ch4 * 1e9 / m_atmosphere
    )  # W/m2/kg-CH4

    tau = 11.8  # Lifetime (years)

    dates_characterized: np.ndarray = date + np.arange(
        start=0, stop=period, dtype="timedelta64[Y]"
    ).astype("timedelta64[s]")

    decay_multipliers: list = np.array(
        [
            radiative_efficiency_kg * tau * (1 - np.exp(-year / tau))
            for year in range(period)
        ]
    )

    forcing = np.array(amount * decay_multipliers, dtype="float64")

    if not cumulative:
        forcing = np.diff(forcing, prepend=0)

    return dates_characterized, forcing, flow, activity


def characterize_n2o(
    date,
    amount,
    flow,
    activity,
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

    # functional variables and units (from publications listed in docstring)
    radiative_efficiency_ppb = 2.8e-3  # # W/m2/ppb; 2019 background cch4 concentration; IPCC AR6 Table 7.15. This number includes indirect effects.

    # for conversion from ppb to kg-CH4
    M_n2o = 44.01  # g/mol
    M_air = 28.97  # g/mol, dry air
    m_atmosphere = 5.135e18  # kg [Trenberth and Smith, 2005]

    radiative_efficiency_kg = (
        radiative_efficiency_ppb * M_air / M_n2o * 1e9 / m_atmosphere
    )  # W/m2/kg-N2O

    tau = 109  # Lifetime (years)

    dates_characterized: np.ndarray = date + np.arange(
        start=0, stop=period, dtype="timedelta64[Y]"
    ).astype("timedelta64[s]")

    decay_multipliers: list = np.array(
        [
            radiative_efficiency_kg * tau * (1 - np.exp(-year / tau))
            for year in range(period)
        ]
    )

    forcing = np.array(amount * decay_multipliers, dtype="float64")
    if not cumulative:
        forcing = np.diff(forcing, prepend=0)

    return dates_characterized, forcing, flow, activity


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
        date,
        amount,
        flow,
        activity,
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



        dates_characterized: np.ndarray = date + np.arange(
            start=0, stop=period, dtype="timedelta64[Y]"
        ).astype("timedelta64[s]")

        decay_multipliers = decay_series[:period]

        forcing = np.array(amount * decay_multipliers, dtype="float64")

        if not cumulative:
            forcing = np.diff(forcing, prepend=0)

        return dates_characterized, forcing, flow, activity

    return characterize_generic
