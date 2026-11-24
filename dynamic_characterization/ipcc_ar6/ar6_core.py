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
from typing import Optional, Tuple, Union

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
class CO2Release:
    """
    ``kg_co2`` of CO2 released into the atmosphere at the rate ``exp(-t / lifetime) / lifetime``
    (e.g. from the oxidation of a gas with that lifetime), or as a pulse at t = 0 if
    ``lifetime`` is 0. Its forcing follows the CO2 response in use: the AR6 impulse response
    (analytic) by default, or a scenario-dependent one (see :class:`TabulatedCO2Response`).
    """

    kg_co2: float
    lifetime: float = 0.0

    def scaled(self, factor: float) -> "CO2Release":
        return CO2Release(kg_co2=self.kg_co2 * factor, lifetime=self.lifetime)


@dataclass(frozen=True)
class TabulatedCO2Response:
    """
    RF of a 1 kg CO2 pulse [W m-2] at t = 0, 1, ..., n - 1 yr, e.g. from a scenario-dependent
    radiative efficiency and impulse response. Values in between are interpolated linearly;
    horizons beyond n - 1 yr are not available.
    """

    name: str
    rf: Tuple[float, ...]

    @property
    def max_years(self) -> int:
        return len(self.rf) - 1


@dataclass(frozen=True)
class GasSpec:
    """Response of a 1 kg pulse emission, as a sum of :class:`Component` s and :class:`CO2Release` s."""

    name: str
    components: Tuple[Union[Component, CO2Release], ...] = ()

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


@lru_cache(maxsize=None)
def _ar6_co2_component() -> Component:
    """RF of a 1 kg CO2 pulse with the AR6 impulse response (Table 7.SM.7 CO2 row)."""
    p = load_parameters()["co2"]
    A = p["radiative_efficiency_W_m2_ppm"] / (1e3 * kg_per_ppb(load_parameters()["constants"]["M_co2"]))
    a0, *a = p["a"]
    return Component(constant=A * a0, terms=tuple((A * ai, taui) for ai, taui in zip(a, p["tau"])))


def co2_spec() -> GasSpec:
    """1 kg CO2 pulse (no CCF by construction: the impulse response contains the carbon cycle)."""
    return GasSpec(name="CO2", components=(CO2Release(kg_co2=1.0),))


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
    return GasSpec(name=name, components=(CO2Release(kg_co2=kg_co2, lifetime=lifetime),))


@lru_cache(maxsize=None)
def _release_as_component(release: CO2Release) -> Component:
    """A CO2 release with the AR6 CO2 response, as an analytic sum of exponentials."""
    co2 = _ar6_co2_component()
    if release.lifetime == 0:
        return co2.scaled(release.kg_co2)
    lifetime, kg_co2 = release.lifetime, release.kg_co2
    terms = [(-co2.constant * kg_co2, lifetime)]
    constant = co2.constant * kg_co2
    for c, tau in co2.terms:
        factor = c * kg_co2 * tau / (tau - lifetime)
        terms += [(factor, tau), (-factor, lifetime)]
    return Component(terms=tuple(terms), constant=constant)


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


def _grid(years: int) -> np.ndarray:
    return np.linspace(0, years, years * STEPS_PER_YEAR + 1)


def _agtp_from_rf_grid(rf: np.ndarray, t: np.ndarray) -> np.ndarray:
    """
    AGTP on the grid from an RF on the grid: the temperature response is integrated exactly
    over each step and combined with the step-average RF (second-order accurate).
    """
    dt = 1.0 / STEPS_PER_YEAR
    q, d = _thermal()
    edges = np.append(t, t[-1] + dt)
    kernel = sum(qj * -np.diff(np.exp(-edges / dj)) for qj, dj in zip(q, d))[:-1]
    rf_step = (rf[1:] + rf[:-1]) / 2
    agtp = np.zeros_like(t)
    agtp[1:] = _convolve(rf_step, kernel[: rf_step.size])
    return agtp


@lru_cache(maxsize=None)
def _co2_grids(co2_response: Optional[TabulatedCO2Response], years: int) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """(rf, agwp, agtp) of a 1 kg CO2 pulse on the grid t = 0, dt, ..., years."""
    t = _grid(years)
    if co2_response is None:
        co2 = _ar6_co2_component()
        grids = (_rf_direct(co2, t), _agwp_direct(co2, t), _agtp_direct(co2, t))
    else:
        if years > co2_response.max_years:
            raise ValueError(
                f"The CO2 response {co2_response.name!r} covers {co2_response.max_years} years; "
                f"{years} were requested."
            )
        dt = 1.0 / STEPS_PER_YEAR
        rf = np.interp(t, np.arange(len(co2_response.rf), dtype="float64"), np.asarray(co2_response.rf))
        agwp = np.concatenate(([0.0], np.cumsum((rf[1:] + rf[:-1]) / 2 * dt)))
        agtp = _agtp_from_rf_grid(rf, t)
        grids = (rf, agwp, agtp)
    for arr in grids:
        arr.setflags(write=False)
    return grids


@lru_cache(maxsize=None)
def _ccf_on_grid(
    spec: Component, years: int, co2_response: Optional[TabulatedCO2Response] = None
) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """(rf_cc, agwp_cc, agtp_cc) on the grid t = 0, dt, ..., years (AR6 carbon_cycle_adjustment)."""
    t = _grid(years)
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
    co2_rf, co2_agwp, co2_agtp = _co2_grids(co2_response, years)
    rf_cc = _convolve(flux_co2, co2_rf) * dt
    agwp_cc = _convolve(flux_co2, co2_agwp) * dt
    agtp_cc = _convolve(flux_co2, co2_agtp) * dt
    for arr in (rf_cc, agwp_cc, agtp_cc):
        arr.setflags(write=False)
    return rf_cc, agwp_cc, agtp_cc


@lru_cache(maxsize=None)
def _release_on_grid(
    release: CO2Release, years: int, co2_response: TabulatedCO2Response
) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """(rf, agwp, agtp) of a CO2 release with a tabulated CO2 response, on the grid."""
    t = _grid(years)
    grids = _co2_grids(co2_response, years)
    if release.lifetime == 0:
        out = tuple(release.kg_co2 * g for g in grids)
    else:
        # fraction released during each grid step; trapezoid over the release time within the step
        fraction = -np.diff(np.exp(-np.append(t, t[-1] + 1.0 / STEPS_PER_YEAR) / release.lifetime))
        out = tuple(
            release.kg_co2 * _convolve(fraction, (g + np.concatenate(([0.0], g[:-1]))) / 2) for g in grids
        )
    for arr in out:
        arr.setflags(write=False)
    return out


def _grid_index(horizons: np.ndarray) -> np.ndarray:
    idx = np.rint(horizons * STEPS_PER_YEAR).astype(int)
    if not np.allclose(idx, horizons * STEPS_PER_YEAR):
        raise ValueError(
            "With the climate-carbon feedback or a tabulated CO2 response, time horizons must be "
            "multiples of 0.1 yr."
        )
    return idx


def _grid_years(horizons: np.ndarray) -> int:
    return max(int(np.ceil(horizons.max())) if horizons.size else 0, 1)


def _ccf_at(spec: Component, horizons, which: int, co2_response=None) -> np.ndarray:
    horizons = np.atleast_1d(np.asarray(horizons, dtype="float64"))
    if spec.ccf is None:
        return np.zeros_like(horizons)
    idx = _grid_index(horizons)
    return _ccf_on_grid(spec, _grid_years(horizons), co2_response)[which][idx]


# --- public metric functions ---------------------------------------------------------------

_DIRECT = (_rf_direct, _agwp_direct, _agtp_direct)


def _sum_over_components(spec: GasSpec, x, which: int, co2_response=None) -> np.ndarray:
    shape = np.shape(x)
    total = np.zeros(shape, dtype="float64")
    for component in spec.components:
        if isinstance(component, CO2Release):
            if co2_response is None:
                total = total + _DIRECT[which](_release_as_component(component), x)
            else:
                horizons = np.atleast_1d(np.asarray(x, dtype="float64"))
                grid = _release_on_grid(component, _grid_years(horizons), co2_response)[which]
                total = total + grid[_grid_index(horizons)].reshape(shape)
        else:
            total = total + _DIRECT[which](component, x) + _ccf_at(component, x, which, co2_response).reshape(shape)
    return total


def radiative_forcing(spec: GasSpec, t, co2_response: Optional[TabulatedCO2Response] = None) -> np.ndarray:
    """Instantaneous RF [W m-2 kg-1] at times ``t`` (yr) after the pulse."""
    return _sum_over_components(spec, t, 0, co2_response)


def agwp(spec: GasSpec, horizon, co2_response: Optional[TabulatedCO2Response] = None) -> np.ndarray:
    """Absolute GWP [W m-2 yr kg-1] at ``horizon`` (yr)."""
    return _sum_over_components(spec, horizon, 1, co2_response)


def agtp(spec: GasSpec, horizon, co2_response: Optional[TabulatedCO2Response] = None) -> np.ndarray:
    """
    Absolute GTP [K kg-1] at ``horizon`` (yr), with the AR6 temperature response. With a
    tabulated CO2 response, its CO2 parts are convolved numerically with that temperature
    response.
    """
    return _sum_over_components(spec, horizon, 2, co2_response)


@lru_cache(maxsize=None)
def annual_cumulative(spec: GasSpec, period: int, co2_response: Optional[TabulatedCO2Response] = None) -> np.ndarray:
    """
    Cumulative forcing at the end of each year after the pulse: ``AGWP(k + 1)`` for
    ``k = 0, ..., period - 1`` (read-only).

    Element ``k`` belongs to the year that starts ``k`` years after the emission, so
    ``np.diff(..., prepend=0)`` gives the forcing integrated over each year and the annual
    values add up to ``AGWP(period)``.
    """
    arr = np.asarray(
        agwp(spec, np.arange(1, period + 1, dtype="float64"), co2_response=co2_response), dtype="float64"
    )
    arr.setflags(write=False)
    return arr
