"""
Prospective absolute global temperature change potentials (AGTP) of CO2, CH4 and N2O, in K kg-1.

Since version 2.0 these are computed with :mod:`.watanabe2026` (Barbosa Watanabe & Cherubini
2026, annual two-layer temperature model, with two errors of the published code corrected), for
the scenario of ``prospective.set_scenario()``. The signatures are unchanged; horizons are
limited to 100 years.
"""

from .agwp import _value


def agtp_co2(emission_year: int, time_horizon: int = 100, time_varying_re: bool = False) -> float:
    """AGTP of 1 kg CO2 emitted in ``emission_year`` (2020-2150, clamped outside)."""
    return _value("CO2", "agtp", emission_year, time_horizon, time_varying_re)


def agtp_ch4(emission_year: int, time_horizon: int = 100, time_varying_re: bool = False) -> float:
    """AGTP of 1 kg CH4, incl. ozone, stratospheric water vapour and climate-carbon feedback."""
    return _value("CH4", "agtp", emission_year, time_horizon, time_varying_re)


def agtp_n2o(emission_year: int, time_horizon: int = 100, time_varying_re: bool = False) -> float:
    """AGTP of 1 kg N2O, incl. ozone, its CH4 effect and climate-carbon feedback."""
    return _value("N2O", "agtp", emission_year, time_horizon, time_varying_re)
