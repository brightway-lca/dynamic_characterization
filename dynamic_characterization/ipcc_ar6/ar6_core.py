"""
Response functions behind the IPCC AR6 emission metrics (WG1 Ch. 7, Section 7.6 and 7.SM.5).

This module is pure numpy and has no Brightway dependency. It follows the AR6 metric code
(https://github.com/chrisroadmap/ar6, ``src/ar6/metrics`` and notebook
``335_chapter7_generate_metrics.ipynb``, MIT licence), which produced AR6 Table 7.SM.7:

- the radiative forcing (RF) of a 1 kg pulse is a constant plus a sum of exponentials. This
  covers CO2 (Joos et al. 2013 impulse response, AR6 parameters), gases with a single
  lifetime, and the CO2 released by the oxidation of another gas (an exponential source
  convolved with the CO2 response is again a sum of exponentials). AGWP and AGTP of such a
  sum are evaluated analytically, as in the AR6 code;
- the climate-carbon feedback (CCF, Eqs. 7.SM.5.5-7.SM.5.6, Gasser et al. 2017) is added
  numerically on the AR6 time grid (0.1 yr), with the same discrete convolution as the AR6
  code.

All RFs are in W m-2 per kg emitted, AGWPs in W m-2 yr kg-1 and AGTPs in K kg-1.
"""

import json
import os
from dataclasses import dataclass
from functools import lru_cache
from typing import Optional, Tuple

import numpy as np

_PARAMETER_FILE = os.path.join(os.path.dirname(os.path.abspath(__file__)), "data", "ar6_parameters.json")

# time step of the numerical CCF convolution, as in the AR6 code (ts_per_year = 10)
STEPS_PER_YEAR = 10


@lru_cache(maxsize=None)
def load_parameters() -> dict:
    """The AR6 parameters bundled in ``ipcc_ar6/data/ar6_parameters.json`` (see its ``_source``)."""
    with open(_PARAMETER_FILE, encoding="utf-8") as f:
        return json.load(f)


def kg_per_ppb(molar_mass: float) -> float:
    """Mass of a gas [kg] that raises its mixing ratio by 1 ppb (molar mass in kg/mol)."""
    constants = load_parameters()["constants"]
    return 1e-9 * molar_mass / constants["M_air"] * constants["m_atmosphere"]


@dataclass(frozen=True)
class CCFParameters:
    """Climate-carbon feedback: CO2 flux response to warming (Gasser et al. 2017, Table 1)."""

    gamma: float  # kg C yr-1 K-1
    a: Tuple[float, ...]
    tau: Tuple[float, ...]  # yr


def ccf_parameters(name: str = "gasser2017_mean") -> CCFParameters:
    """CCF parameter set by name: ``gasser2017_mean`` (AR6), ``gasser2017_upper`` or ``gasser2017_lower``."""
    sets = load_parameters()["climate_carbon_feedback"]
    if name not in sets:
        raise ValueError(f"Unknown climate-carbon feedback parameters {name!r}; available: {sorted(sets)}")
    p = sets[name]
    return CCFParameters(gamma=p["gamma_kgC_per_yr_K"], a=tuple(p["a"]), tau=tuple(p["tau"]))


@dataclass(frozen=True)
class Component:
    """
    One part of a pulse response: RF = ``constant + sum(c * exp(-t / tau) for c, tau in terms)``.

    If ``ccf`` is set, the climate-carbon feedback of this part's warming is added on top. It
    is never set for CO2 (emitted or from oxidation), whose impulse response already contains
    the carbon cycle.
    """

    terms: Tuple[Tuple[float, float], ...] = ()
    constant: float = 0.0
    ccf: Optional[CCFParameters] = None

    def scaled(self, factor: float) -> "Component":
        return Component(
            terms=tuple((c * factor, tau) for c, tau in self.terms),
            constant=self.constant * factor,
            ccf=self.ccf,
        )


@dataclass(frozen=True)
class GasSpec:
    """Response of a 1 kg pulse emission, as a sum of :class:`Component` s."""

    name: str
    components: Tuple[Component, ...] = ()

    def __add__(self, other: "GasSpec") -> "GasSpec":
        return GasSpec(name=self.name, components=self.components + other.components)

    def scaled(self, factor: float, name: Optional[str] = None) -> "GasSpec":
        return GasSpec(name=name or self.name, components=tuple(c.scaled(factor) for c in self.components))

    def renamed(self, name: str) -> "GasSpec":
        return GasSpec(name=name, components=self.components)


# --- temperature and CO2 responses ---------------------------------------------------------


def _thermal() -> Tuple[np.ndarray, np.ndarray]:
    p = load_parameters()["temperature_response"]
    return np.asarray(p["q"]), np.asarray(p["d"])


def co2_spec() -> GasSpec:
    """RF of a 1 kg CO2 pulse (AR6 Table 7.SM.7 CO2 row; no CCF by construction)."""
    p = load_parameters()["co2"]
    A = p["radiative_efficiency_W_m2_ppm"] / (1e3 * kg_per_ppb(load_parameters()["constants"]["M_co2"]))
    a0, *a = p["a"]
    return GasSpec(
        name="CO2",
        components=(Component(constant=A * a0, terms=tuple((A * ai, taui) for ai, taui in zip(a, p["tau"]))),),
    )


def exponential_spec(name: str, radiative_efficiency_per_kg: float, lifetime: float, ccf=None) -> GasSpec:
    """RF of a gas with a single perturbation lifetime."""
    return GasSpec(name=name, components=(Component(terms=((radiative_efficiency_per_kg, lifetime),), ccf=ccf),))


def co2_from_decay(kg_co2: float, lifetime: float, name: str = "oxidation CO2") -> GasSpec:
    """
    RF of ``kg_co2`` of CO2 released at the rate ``exp(-t / lifetime) / lifetime`` (oxidation of
    a gas with that lifetime; AR6 Eq. 7.SM.5.8).

    Convolving that source with each CO2 box gives
    ``a_i * tau_i / (tau_i - lifetime) * (exp(-t/tau_i) - exp(-t/lifetime))`` and, for the
    constant box, ``a_0 * (1 - exp(-t/lifetime))``.
    """
    (co2,) = co2_spec().components
    terms = [(-co2.constant * kg_co2, lifetime)]
    constant = co2.constant * kg_co2
    for c, tau in co2.terms:
        factor = c * kg_co2 * tau / (tau - lifetime)
        terms += [(factor, tau), (-factor, lifetime)]
    return GasSpec(name=name, components=(Component(terms=tuple(terms), constant=constant),))


def precursor_component(
    steady_state_forcing: float,
    response_time: float,
    precursor_lifetime: float = 0.0,
    ccf: Optional[CCFParameters] = None,
) -> Component:
    """
    One forcing component of an indirect (precursor) emission, e.g. the ozone or CH4 change
    caused by 1 kg of CO, VOC, H2 or NOx (Fuglestvedt et al. 2010, Eqs. 10-14; Aamaas et al.
    2013; with a finite precursor lifetime as in Warwick et al. 2023, Eq. 5):

        RF(t) = F_ss * (exp(-t / alpha_M) - exp(-t / alpha_p)) / (alpha_M - alpha_p)

    ``steady_state_forcing`` (F_ss) is the forcing per unit sustained emission
    [W m-2 per kg yr-1], which equals the AGWP of a 1 kg pulse at an infinite horizon.
    ``response_time`` (alpha_M) is the time constant of the component (e.g. the CH4 perturbation
    lifetime), ``precursor_lifetime`` (alpha_p) that of the precursor; either may be 0, in which
    case the RF follows the other one alone.
    """
    if response_time < 0 or precursor_lifetime < 0:
        raise ValueError("time constants must not be negative")
    if response_time == 0 and precursor_lifetime == 0:
        raise ValueError("response_time and precursor_lifetime cannot both be zero")
    if np.isclose(response_time, precursor_lifetime, rtol=1e-9):
        raise ValueError("response_time and precursor_lifetime must differ")
    if precursor_lifetime == 0:
        terms = ((steady_state_forcing / response_time, response_time),)
    elif response_time == 0:
        terms = ((steady_state_forcing / precursor_lifetime, precursor_lifetime),)
    else:
        c = steady_state_forcing / (response_time - precursor_lifetime)
        terms = ((c, response_time), (-c, precursor_lifetime))
    return Component(terms=terms, ccf=ccf)


# --- analytic metrics of a sum of exponentials ---------------------------------------------


def _rf_direct(spec: Component, t) -> np.ndarray:
    t = np.asarray(t, dtype="float64")
    rf = np.full_like(t, spec.constant)
    for c, tau in spec.terms:
        rf = rf + c * np.exp(-t / tau)
    return rf


def _agwp_direct(spec: Component, t) -> np.ndarray:
    t = np.asarray(t, dtype="float64")
    agwp = spec.constant * t
    for c, tau in spec.terms:
        agwp = agwp + c * tau * -np.expm1(-t / tau)
    return agwp


def _agtp_direct(spec: Component, t) -> np.ndarray:
    """Convolution of the RF with R_T(t) = sum_j q_j / d_j exp(-t / d_j) (AR6 7.SM.5.2)."""
    t = np.asarray(t, dtype="float64")
    q, d = _thermal()
    agtp = np.zeros_like(t)
    for qj, dj in zip(q, d):
        agtp = agtp + spec.constant * qj * -np.expm1(-t / dj)
        for c, tau in spec.terms:
            if np.isclose(tau, dj, rtol=1e-12):
                agtp = agtp + c * qj * t / dj * np.exp(-t / dj)
            else:
                agtp = agtp + c * tau * qj * (np.exp(-t / tau) - np.exp(-t / dj)) / (tau - dj)
    return agtp


# --- climate-carbon feedback, numerically as in the AR6 code (gasser.py) ------------------


def _convolve(x: np.ndarray, y: np.ndarray) -> np.ndarray:
    """out[j] = sum_{i<=j} x[i] * y[j-i] (same length as x), via FFT."""
    n = x.size
    size = 1 << int(np.ceil(np.log2(2 * n - 1)))
    return np.fft.irfft(np.fft.rfft(x, size) * np.fft.rfft(y, size), size)[:n]


@lru_cache(maxsize=None)
def _ccf_on_grid(spec: Component, years: int) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """(rf_cc, agwp_cc, agtp_cc) on the grid t = 0, dt, ..., years (AR6 carbon_cycle_adjustment)."""
    t = np.linspace(0, years, years * STEPS_PER_YEAR + 1)
    dt = 1.0 / STEPS_PER_YEAR
    ccf = spec.ccf
    agtp = _agtp_direct(spec, t)

    r_f = np.zeros_like(t)
    r_f[0] = np.sum(ccf.a) / dt  # the delta term of the flux response
    for a, tau in zip(ccf.a, ccf.tau):
        r_f = r_f - (a / tau) * np.exp(-t / tau)
    flux_c = ccf.gamma * _convolve(agtp, r_f) * dt  # kg C yr-1

    constants = load_parameters()["constants"]
    flux_co2 = flux_c * constants["M_co2"] / constants["M_c"]  # kg CO2 yr-1
    (co2,) = co2_spec().components
    rf_cc = _convolve(flux_co2, _rf_direct(co2, t)) * dt
    agwp_cc = _convolve(flux_co2, _agwp_direct(co2, t)) * dt
    agtp_cc = _convolve(flux_co2, _agtp_direct(co2, t)) * dt
    for arr in (rf_cc, agwp_cc, agtp_cc):
        arr.setflags(write=False)
    return rf_cc, agwp_cc, agtp_cc


def _ccf_at(spec: Component, horizons, which: int) -> np.ndarray:
    horizons = np.atleast_1d(np.asarray(horizons, dtype="float64"))
    if spec.ccf is None:
        return np.zeros_like(horizons)
    idx = np.rint(horizons * STEPS_PER_YEAR).astype(int)
    if not np.allclose(idx, horizons * STEPS_PER_YEAR):
        raise ValueError("With the climate-carbon feedback, time horizons must be multiples of 0.1 yr.")
    years = int(np.ceil(horizons.max())) if horizons.size else 0
    return _ccf_on_grid(spec, max(years, 1))[which][idx]


# --- public metric functions ---------------------------------------------------------------


def _sum_over_components(spec: GasSpec, x, direct, which: int) -> np.ndarray:
    total = np.zeros(np.shape(x), dtype="float64")
    for component in spec.components:
        total = total + direct(component, x) + _ccf_at(component, x, which).reshape(np.shape(x))
    return total


def radiative_forcing(spec: GasSpec, t) -> np.ndarray:
    """Instantaneous RF [W m-2 kg-1] at times ``t`` (yr) after the pulse."""
    return _sum_over_components(spec, t, _rf_direct, 0)


def agwp(spec: GasSpec, horizon) -> np.ndarray:
    """Absolute GWP [W m-2 yr kg-1] at ``horizon`` (yr)."""
    return _sum_over_components(spec, horizon, _agwp_direct, 1)


def agtp(spec: GasSpec, horizon) -> np.ndarray:
    """Absolute GTP [K kg-1] at ``horizon`` (yr)."""
    return _sum_over_components(spec, horizon, _agtp_direct, 2)


@lru_cache(maxsize=None)
def annual_cumulative(spec: GasSpec, period: int) -> np.ndarray:
    """
    Cumulative forcing at the end of each year after the pulse: ``AGWP(k + 1)`` for
    ``k = 0, ..., period - 1`` (read-only).

    Element ``k`` belongs to the year that starts ``k`` years after the emission, so
    ``np.diff(..., prepend=0)`` gives the forcing integrated over each year and the annual
    values add up to ``AGWP(period)``.
    """
    arr = np.asarray(agwp(spec, np.arange(1, period + 1, dtype="float64")), dtype="float64")
    arr.setflags(write=False)
    return arr
