"""
The pGWP and pGTP calculation of Barbosa Watanabe & Cherubini (2026),
Environ. Sci. Technol. 60, 3202-3215, doi:10.1021/acs.est.5c12391.

This is a port of the computational functions of the paper's code in Supporting Information
es5c12391_si_028.pdf (Section 3, "public release version, November 2025"):

- Code #3  AGWP of CO2          -> :func:`agwp_co2_series`
- Code #4  AGWP of CH4          -> :func:`agwp_ch4_series`
- Code #5  AGWP of N2O          -> :func:`agwp_n2o_series`
- Code #6  GWP20 and GWP100     -> :func:`pgwp`
- Code #7  AGTP of CO2          -> :func:`agtp_co2_series`
- Code #8  AGTP of CH4          -> :func:`agtp_ch4_series`
- Code #9  AGTP of N2O          -> :func:`agtp_n2o_series`
- Code #10 GTP50 and GTP100     -> :func:`pgtp`

The arithmetic, constants, time indexing and end-point treatment follow that code, apart from
the two corrections listed below. In short:

- RF of CH4 = (RE_CH4 + 1.4e-4 + 0.4e-4) x IRF (RE of the pulse year, fixed, as in the paper;
  optionally time-varying, see :func:`scenario_series`; the two constants
  are the AR6 ozone and stratospheric water vapour terms, no tropospheric adjustment of RE_CH4);
  RF of N2O = (RE_N2O + 5.5e-4 - 1.7 x (RE_CH4 + 1.8e-4)) x IRF; RF of CO2 = RE_CO2 x IRF.
- Temperature: annual explicit two-layer model (C_s = 7.7, C_d = 147 W yr m-2 K-1,
  lambda = 1.31 W m-2 K-1, gamma = 0.88 W m-2 K-1, efficacy 1.03).
- Climate-carbon feedback for CH4 and N2O: three boxes driven by the surface temperature
  (Gasser et al. 2017: gamma = 11.06 Gt CO2 yr-1 K-1, a = 0.6368 / 0.3322 / 0.031,
  tau = 2.376 / 30.14 / 490.1 yr), whose CO2 forces with RE_CO2 of the pulse year.
- AGWP by the trapezoidal rule, AGTP as the mean of consecutive annual temperatures.

Two features of the published code are corrected here, so the results differ from the SI tables:

- In Code #4, #5 and #8 the second and third feedback boxes decay from the previous value of the
  *first* box (``C2[t] = ... + C1[t-1] * exp(-1/30.14)``). Here, as in Code #9, every box decays
  from its own previous value.
- Code #3 and #7 take the column of the CO2 impulse-response file (RCP2.6, 4.5, 6.0, 8.5) at the
  position of the scenario among the IAM's scenarios, which picks the wrong RCP for AIM (RCP4.5,
  6.0, 8.5) and IMAGE (RCP2.6, 4.5, 8.5); e.g. IMAGE-SSP1-8.5 used the RCP6.0 response. Here the
  response of the scenario's own RCP is used.

Everything else follows the published code, which, with both features restored, reproduces every
cell of es5c12391_si_001-004.xlsx and SI Tables S.2-S.9 within their rounding to one decimal
(checked when this module was written). With the corrections, N2O pGTP of GCAM4, MESSAGE and
REMIND, which neither feature affects, still equals the tables; CH4 pGTP100 is 16-34 % higher
in the scenarios checked (2030).

GWP20/GWP100 and GTP50/GTP100 are the gas value divided by the CO2 value at row H - 1 of the
annual series.

Units: the code works in W m-2 per kg x 1e12 (pW m-2 kg-1) and K kg-1 x 1e12; the ``*_series``
functions return these values divided by 1e12, i.e. W m-2 yr kg-1 and K kg-1.
"""

from functools import lru_cache
from typing import Tuple

import numpy as np

from .config import get_scenario
from .data_loader import load_irf_ch4, load_irf_co2, load_irf_n2o, load_re_ch4, load_re_co2, load_re_n2o
from ._common import _get_year_index

_SCALE = 1e12  # the paper's code works in pW m-2 kg-1 and pK kg-1

# two-layer temperature model (Code #4-#9)
_C_S, _C_D, _LAMBDA, _GAMMA, _EFFICACY = 7.7, 147.0, 1.31, 0.88, 1.03
# climate-carbon feedback (Gasser et al. 2017, as in Code #4, #5, #8, #9)
_GAMMA_CO2 = 11.06  # Gt CO2 yr-1 K-1; times T in pK gives kg CO2
_CCF_A = (0.6368, 0.3322, 0.031)
_CCF_TAU = (2.376, 30.14, 490.1)
# indirect radiative efficiencies [W m-2 ppb-1] added to CH4 and N2O
_CH4_INDIRECT = 0.00014 + 0.00004
_N2O_OZONE = 0.00055
_N2O_CH4 = -1.7


def _per_kg(molar_mass: float) -> float:
    """The code's ppb -> kg conversion: 1e9 * 1e-6 * M_air / M / 5.13252."""
    return 10**9 * 10**-6 * 28.97 / molar_mass / 5.13252


def _at(re, t: int) -> float:
    """Radiative efficiency at step ``t``: a constant, or element ``t`` of an annual series
    (clamped to its last element)."""
    if np.ndim(re) == 0:
        return float(re)
    return float(re[min(t, len(re) - 1)])


def _temperature_step(T_s, T_d, forcing):
    T_s_new = T_s + forcing / _C_S + (-_LAMBDA) * T_s / _C_S - _EFFICACY * (_GAMMA / _EFFICACY) * (T_s - T_d) / _C_S
    T_d_new = T_d + (_GAMMA / _EFFICACY) * (T_s - T_d) / _C_D
    return T_s_new, T_d_new


def agtp_from_rf(rf_values) -> np.ndarray:
    """
    AGTP series [K kg-1] from an annual RF series [W m-2 kg-1] at t = 0, 1, ..., with the
    temperature scheme of Code #7 (as for CO2); used for gases outside the paper.
    """
    rf = [float(v) for v in rf_values]
    T_s, T_d = [0.0], [0.0]
    for t in range(1, len(rf)):
        T_s_new, T_d_new = _temperature_step(T_s[t - 1], T_d[t - 1], rf[t - 1])
        T_s.append(T_s_new)
        T_d.append(T_d_new)
    agtp = [(T_s[t] + T_s[t + 1]) / 2 for t in range(len(T_s) - 1)]
    agtp.append(agtp[-1] + (T_s[-1] - T_s[-2]) / 2)
    return np.asarray(agtp)


def _ccf_boxes(T_s, C):
    """One annual step of the three feedback boxes, each decaying from its own previous value."""
    return tuple(
        T_s * 10**-9 * _GAMMA_CO2 * 10**9 * a + previous * np.exp(-1 / tau)
        for a, tau, previous in zip(_CCF_A, _CCF_TAU, C)
    )


# --- Code #3: AGWP of CO2 ---------------------------------------------------------------------


def agwp_co2_series(re_co2: float, irf_co2, integration_years: int = 125) -> np.ndarray:
    """AGWP of 1 kg CO2 [W m-2 yr kg-1]; element i belongs to calendar year pulse + i."""
    irf = np.asarray(irf_co2, dtype="float64")
    rf = [_at(re_co2, i) * irf[i] * _per_kg(44.01) for i in range(min(integration_years, len(irf)))]
    agwp = [(rf[0] + rf[1]) / 2]
    for t in range(1, len(rf) - 1):
        agwp.append(agwp[t - 1] + (rf[t] + rf[t + 1]) / 2)
    agwp.append(agwp[-1] + rf[-1])
    return np.asarray(agwp) / _SCALE


# --- Code #4 and #8: CH4 ----------------------------------------------------------------------


def _ch4_chain(re_ch4: float, re_co2: float, irf_ch4) -> Tuple[list, list, list, list]:
    """factor1 (CH4 RF), factor2 (CCF RF), T_s and T_s' (temperature from the CCF CO2), t = 1..101."""
    irf = np.asarray(irf_ch4, dtype="float64")
    constant = 10**-6 * 28.97 / 16.04 * 1 / 5.13252
    co2_constant = 10**-6 * 28.97 / 44.01 * 1 / 5.13252
    T_s = T_d = T_s_p = T_d_p = 0.0
    C = (0.0, 0.0, 0.0)
    factor1, factor2, T_s_values, T_s_p_values = [], [], [], []
    for t in range(1, 102):
        f1 = (_at(re_ch4, t - 1) + _CH4_INDIRECT) * 10**9 * (irf[t - 1] * constant)
        T_s_new, T_d = _temperature_step(T_s, T_d, f1)
        T_s = T_s_new
        C = _ccf_boxes(T_s, C)
        f2 = sum(C) * co2_constant * _at(re_co2, t) * 10**9
        T_s_p, T_d_p = _temperature_step(T_s_p, T_d_p, f2)
        factor1.append(f1)
        factor2.append(f2)
        T_s_values.append(T_s)
        T_s_p_values.append(T_s_p)
    return factor1, factor2, T_s_values, T_s_p_values


def agwp_ch4_series(re_ch4: float, re_co2: float, irf_ch4) -> np.ndarray:
    """AGWP of 1 kg CH4 incl. ozone, stratospheric H2O and CCF [W m-2 yr kg-1], 100 values."""
    factor1, factor2, _, _ = _ch4_chain(re_ch4, re_co2, irf_ch4)
    agwp = [0.0]
    for t in range(1, 101):
        agwp.append(agwp[t - 1] + (factor1[t - 1] + factor1[t]) / 2 + (factor2[t - 1] + factor2[t]) / 2)
    return np.asarray(agwp[1:]) / _SCALE


def agtp_ch4_series(re_ch4: float, re_co2: float, irf_ch4) -> np.ndarray:
    """AGTP of 1 kg CH4 incl. ozone, stratospheric H2O and CCF [K kg-1], 100 values."""
    _, _, T_s, T_s_p = _ch4_chain(re_ch4, re_co2, irf_ch4)
    agtp = [(T_s[t] + T_s[t + 1]) / 2 + (T_s_p[t] + T_s_p[t + 1]) / 2 for t in range(len(T_s) - 1)]
    slope, intercept = np.polyfit(np.array([99, 100]), np.array(agtp[-2:]), 1)
    agtp.append(slope * 101 + intercept)
    return np.asarray(agtp[:-1]) / _SCALE


# --- Code #5 and #9: N2O ----------------------------------------------------------------------


def agwp_n2o_series(re_n2o: float, re_ch4: float, re_co2: float, irf_n2o) -> np.ndarray:
    """AGWP of 1 kg N2O incl. ozone, CH4 effect and CCF [W m-2 yr kg-1] (Code #5), 100 values."""
    irf = np.asarray(irf_n2o, dtype="float64")
    constant = 10**-6 * 28.97 / 44.013 * 1 / 5.13252
    T_s = T_d = 0.0
    C = (0.0, 0.0, 0.0)
    factor1, factor2 = [], []
    for t in range(1, 101):
        f1 = (_at(re_n2o, t - 1) + _N2O_OZONE + _N2O_CH4 * (_at(re_ch4, t - 1) + _CH4_INDIRECT)) * 10**9 * (irf[t - 1] * constant)
        T_s_new, T_d = _temperature_step(T_s, T_d, f1)
        T_s = T_s_new
        C = _ccf_boxes(T_s, C)
        factor1.append(f1)
        factor2.append(sum(C) * constant * _at(re_co2, t) * 10**9)
    cumulative, agwp = 0.0, []
    for t in range(len(factor1) - 1):
        cumulative += (factor1[t] + factor1[t + 1]) / 2 + (factor2[t] + factor2[t + 1]) / 2
        agwp.append(cumulative)
    agwp.append(agwp[-1] + factor1[-1] + factor2[-1])
    return np.asarray(agwp[:100]) / _SCALE


def agtp_n2o_series(re_n2o: float, re_ch4: float, re_co2: float, irf_n2o) -> np.ndarray:
    """AGTP of 1 kg N2O incl. ozone, CH4 effect and CCF [K kg-1] (Code #9), 100 values."""
    irf = np.asarray(irf_n2o, dtype="float64")
    n = irf.shape[0]
    constant = 10**-6 * 28.97 / 44.013 * 1 / 5.13252
    T_s, T_d, T_s_p, T_d_p = [0.0] * n, [0.0] * n, [0.0] * n, [0.0] * n
    C = (0.0, 0.0, 0.0)
    for t in range(n):
        f1 = (_at(re_n2o, t) + _N2O_OZONE + _N2O_CH4 * (_at(re_ch4, t) + _CH4_INDIRECT)) * 10**9 * (irf[t] * constant)
        f2 = 0.0
        if t > 0:
            T_s[t], T_d[t] = _temperature_step(T_s[t - 1], T_d[t - 1], f1)
            C = _ccf_boxes(T_s[t], C)
            f2 = sum(C) * constant * _at(re_co2, t) * 10**9
        # T_s_prime[t - 1] is T_s_prime[-1] (= 0) for t = 0, as in the code
        T_s_p[t], T_d_p[t] = _temperature_step(T_s_p[t - 1], T_d_p[t - 1], f2)
    agtp = [(T_s[t] + T_s[t + 1]) / 2 + (T_s_p[t] + T_s_p[t + 1]) / 2 for t in range(n - 1)]
    while len(agtp) < 100:
        agtp.append(agtp[-1] + (agtp[-1] - agtp[-2]))
    return np.asarray(agtp[:100]) / _SCALE


# --- Code #7: AGTP of CO2 ---------------------------------------------------------------------


def agtp_co2_series(re_co2: float, irf_co2, n_years: int = 100) -> np.ndarray:
    """AGTP of 1 kg CO2 [K kg-1] (Code #7), 100 values."""
    irf = np.asarray(irf_co2, dtype="float64")
    rf = [_at(re_co2, i) * irf[i] * _per_kg(44.01) for i in range(min(len(irf), n_years))]
    T_s, T_d = [0.0], [0.0]
    for t in range(1, len(rf)):
        rf_value = rf[t - 1] if t > 1 else rf[0]
        T_s_new, T_d_new = _temperature_step(T_s[t - 1], T_d[t - 1], rf_value)
        T_s.append(T_s_new)
        T_d.append(T_d_new)
    agtp = [(T_s[t] + T_s[t + 1]) / 2 for t in range(len(T_s) - 1)]
    agtp.append(agtp[-1] + (T_s[-1] - T_s[-2]) / 2)
    return np.asarray(agtp) / _SCALE


# --- scenario data and the metrics (Code #6, #10) ---------------------------------------------


def _forward_filled(series: np.ndarray) -> np.ndarray:
    """Replace missing values by the last available one."""
    if not np.isnan(series).any():
        return series
    filled = series.copy()
    for i in range(1, filled.size):
        if np.isnan(filled[i]):
            filled[i] = filled[i - 1]
    return filled


def _re(loader, iam: str, ssp: str, rcp: str, year: int, time_varying_re: bool = False):
    """
    RE of the pulse year, or with ``time_varying_re`` the annual series from the pulse year on
    (held at its last value, 2150, beyond the data). Pulse years outside the data are clamped
    with a warning, as elsewhere in the module. Years without data (REMIND CH4 after 2100) take
    the last available value.
    """
    data = loader(iam)
    series = _forward_filled(np.asarray(data[ssp][rcp], dtype="float64"))
    index = min(_get_year_index(year, np.asarray(data["_years"])), len(series) - 1)
    if not time_varying_re:
        return float(series[index])
    return series[index:]


def co2_irf_key(rcp: str) -> str:
    """Key of the CO2 impulse response of an RCP, e.g. "4.5" -> "RCP45"."""
    return f"RCP{rcp.replace('.', '')}"


@lru_cache(maxsize=None)
def scenario_series(
    gas: str, metric: str, iam: str, ssp: str, rcp: str, year: int, time_varying_re: bool = False
) -> np.ndarray:
    """
    AGWP or AGTP series (``metric`` "agwp" / "agtp") of CO2, CH4 or N2O for a scenario and pulse
    year. With ``time_varying_re`` the radiative efficiencies follow the scenario over the
    response (the RE of year ``pulse + t`` at time ``t``); otherwise, as in the paper, those of
    the pulse year are used throughout.
    """
    re_co2 = _re(load_re_co2, iam, ssp, rcp, year, time_varying_re)
    irf_co2 = load_irf_co2()[co2_irf_key(rcp)]
    if gas == "CO2":
        series = agwp_co2_series(re_co2, irf_co2) if metric == "agwp" else agtp_co2_series(re_co2, irf_co2)
    elif gas == "CH4":
        re_ch4 = _re(load_re_ch4, iam, ssp, rcp, year, time_varying_re)
        fn = agwp_ch4_series if metric == "agwp" else agtp_ch4_series
        series = fn(re_ch4, re_co2, load_irf_ch4())
    elif gas == "N2O":
        re_n2o = _re(load_re_n2o, iam, ssp, rcp, year, time_varying_re)
        re_ch4 = _re(load_re_ch4, iam, ssp, rcp, year, time_varying_re)
        if metric == "agwp":
            series = agwp_n2o_series(re_n2o, re_ch4, re_co2, load_irf_n2o())
        else:
            series = agtp_n2o_series(re_n2o, re_ch4, re_co2, load_irf_n2o())
    else:
        raise ValueError(f"gas must be 'CO2', 'CH4' or 'N2O', not {gas!r}")
    series = np.asarray(series)
    series.setflags(write=False)
    return series


def _scenario(iam, ssp, rcp):
    if iam is None:
        scenario = get_scenario()
        return scenario["iam"], scenario["ssp"], scenario["rcp"]
    return iam, ssp, rcp


def pgwp(gas: str, year: int, time_horizon: int = 100, iam=None, ssp=None, rcp=None, time_varying_re: bool = False) -> float:
    """
    pGWP of CH4 or N2O as in Code #6: gas AGWP / CO2 AGWP at row ``time_horizon - 1`` of the
    annual series (row 99 for GWP100). The scenario defaults to ``prospective.get_scenario()``.
    """
    if not 1 <= time_horizon <= 100:
        raise ValueError("time_horizon must be between 1 and 100 years")
    iam, ssp, rcp = _scenario(iam, ssp, rcp)
    gas_series = scenario_series(gas, "agwp", iam, ssp, rcp, year, time_varying_re)
    co2_series = scenario_series("CO2", "agwp", iam, ssp, rcp, year, time_varying_re)
    return float(gas_series[time_horizon - 1] / co2_series[time_horizon - 1])


def pgtp(gas: str, year: int, time_horizon: int = 100, iam=None, ssp=None, rcp=None, time_varying_re: bool = False) -> float:
    """pGTP of CH4 or N2O as in Code #10: gas AGTP / CO2 AGTP at row ``time_horizon - 1``."""
    if not 1 <= time_horizon <= 100:
        raise ValueError("time_horizon must be between 1 and 100 years")
    iam, ssp, rcp = _scenario(iam, ssp, rcp)
    gas_series = scenario_series(gas, "agtp", iam, ssp, rcp, year, time_varying_re)
    co2_series = scenario_series("CO2", "agtp", iam, ssp, rcp, year, time_varying_re)
    return float(gas_series[time_horizon - 1] / co2_series[time_horizon - 1])
