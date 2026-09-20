"""Locate FAIR calibration data and run the model (requires `fair`).

The calibrated, constrained AR6 ensemble (``calibration1.4.1``, 841 members)
is not bundled inside the ``fair`` package itself. We first look for the
calibration CSVs inside the installed ``fair`` package (some installs vendor
them); if they are not found, we download the matched 1.4.1 files from the
FaIR example data on GitHub and cache them with ``pooch``.

FAIR is run from 1750 so the climate state spins up correctly; a dynamic-LCA
inventory is applied as a per-species emission perturbation on top of the SSP
marker background (``fill_from_rcmip``). Emissions live on the ``timepoints``
axis (year midpoints), while temperature/forcing live on ``timebounds``
(calendar years); both are indexed here by ``year - 1750``.

A dynamic-LCA characterization needs one model run per perturbed species (plus
a baseline), and a naive implementation redoes the same expensive work every
time. Three things make that cheap here:

* **Template caching** - reading the calibration, applying the 841 parameter
  sets and building the energy-balance matrices takes ~7 s and does not depend
  on the marker; filling the RCMIP background takes ~1 s per marker. Both are
  cached separately and copied into each run.
* **Spin-up reuse** - 1750 to the first inventory year is identical for the
  baseline and every perturbation. It is run once and cached; the actual runs
  restart from that state, so they only simulate the reporting period.
* **Scenario batching** - FAIR vectorizes over its scenario axis, so several
  perturbations are run in one call, in chunks sized to a memory budget.

Results are memoized per (marker, period, perturbation), so asking for
temperature after radiative forcing (or re-running the same inventory) costs
nothing: both outputs come from the same run.
"""

import glob
import hashlib
import os
import threading
import warnings
from collections import OrderedDict
from functools import lru_cache
from typing import Dict, List, NamedTuple, Optional, Sequence, Tuple

import numpy as np

from . import FAIR_IMPORT_ERROR_MSG

_RUN_LOCK = threading.Lock()

# Simulation always starts here so the climate response spins up from
# pre-industrial; reporting years are sliced out afterwards.
_SIM_START_YEAR = 1750

OUTPUTS = ("radiative_forcing", "temperature")

# The calibrated 1.4.1 ensemble is calibrated for - and distributed with - the
# Thornhill (2021) methane-lifetime chemistry, which is what FaIR's own
# calibrated-constrained example runs. FaIR itself defaults to "leach2021",
# under which the per-species `ch4_lifetime_chemical_sensitivity` values in the
# calibration are simply ignored: a NOx perturbation then leaves CH4 untouched
# and the ozone precursors come out with the wrong sign.
_CH4_METHOD = "thornhill2021"

# FaIR keeps roughly a dozen (timebounds, scenario, config, specie) float64
# arrays alive per run; used to size scenario batches against a memory budget.
_ARRAYS_PER_SCENARIO = 12
_BATCH_MEMORY_BUDGET = 1.5e9  # bytes

# Number of (marker, period, perturbation) results kept around.
_RESULT_CACHE_SIZE = 256

# An LCA-scale perturbation (kilograms against a gigatonne background) is so
# small that the perturbed-minus-baseline difference loses most of its
# significant digits to float64 round-off: at real inventory magnitudes the
# CO2 part of the response came out ~20% off. The response is linear in the
# perturbation far beyond this range, so a perturbation is scaled up to this
# fraction of the species' peak background emission before the run, and the
# difference is scaled back down afterwards. Verified linear (<0.1% deviation)
# from this target over four more orders of magnitude.
_LINEARIZATION_TARGET = 1e-3

# FaIR example data (calibration 1.4.1). Pinned by content hash.
_FAIR_DATA_BASE = (
    "https://raw.githubusercontent.com/OMS-NetZero/FAIR/master/"
    "examples/data/calibrated_constrained_ensemble"
)
_PARAMS_FILENAME = "calibrated_constrained_parameters_calibration1.4.1.csv"
_PROPS_FILENAME = "species_configs_properties_calibration1.4.1.csv"
_PARAMS_HASH = (
    "sha256:7b6c5d9fa0b0b0d3eb47189bf5d63cbf77e752ddac682947abee5ff529206780"
)
_PROPS_HASH = (
    "sha256:42d04aa1a8f385cc53eae22beab26f857c535c5aa7dfdb98176d712bfc0c95a0"
)

_PARAM_CANDIDATES = (
    _PARAMS_FILENAME,
    "calibrated_constrained_parameters.csv",
)
_PROP_CANDIDATES = (
    _PROPS_FILENAME,
    "species_configs_properties.csv",
)

# kg -> FAIR emission unit. FaIR wants CO2 in Gt/yr, the other major gases and
# the short-lived forcers in Mt/yr, and every halogenated species in kt/yr -
# getting that last group wrong is a factor of 1000. The per-species units are
# read from FaIR itself; these constants are the fallback.
_KG_TO_FAIR_CO2 = 1e-12
_KG_TO_FAIR_OTHER = 1e-9
_KG_PER_UNIT_PREFIX = {"Gt": 1e-12, "Mt": 1e-9, "kt": 1e-6, "t": 1e-3}


def require_fair():
    """Import and return the `fair` module, or raise a clear ImportError."""
    try:
        import fair  # noqa: WPS433 (lazy, optional)
    except ImportError as exc:
        raise ImportError(FAIR_IMPORT_ERROR_MSG) from exc
    return fair


def _find_in_fair(fair_module, filenames) -> Optional[str]:
    root = os.path.dirname(os.path.abspath(fair_module.__file__))
    for name in filenames:
        hits = glob.glob(os.path.join(root, "**", name), recursive=True)
        if hits:
            return sorted(hits)[0]
    return None


def find_calibration_files() -> Tuple[Optional[str], Optional[str]]:
    """Locate (parameters_csv, properties_csv) inside the installed fair pkg.

    Either entry may be ``None`` if that file is not vendored in the install
    (the stock ``fair`` package ships only the species-properties defaults).
    Use :func:`get_calibration_files` for paths guaranteed to exist.
    """
    fair = require_fair()
    return _find_in_fair(fair, _PARAM_CANDIDATES), _find_in_fair(
        fair, _PROP_CANDIDATES
    )


@lru_cache(maxsize=1)
def get_calibration_files() -> Tuple[str, str]:
    """Return existing (parameters_csv, properties_csv) for calibration 1.4.1.

    Prefers files vendored in the installed ``fair`` package; otherwise
    downloads the matched 1.4.1 pair from the FaIR example data and caches
    them with ``pooch``. Always returns two real file paths.
    """
    require_fair()
    params, props = find_calibration_files()
    if params is not None and props is not None:
        return params, props

    import pooch  # ships as a fair dependency

    params = pooch.retrieve(
        url=f"{_FAIR_DATA_BASE}/{_PARAMS_FILENAME}", known_hash=_PARAMS_HASH
    )
    props = pooch.retrieve(
        url=f"{_FAIR_DATA_BASE}/{_PROPS_FILENAME}", known_hash=_PROPS_HASH
    )
    return params, props


@lru_cache(maxsize=2)
def _config_index(params_csv: str):
    import pandas as pd

    return list(pd.read_csv(params_csv, index_col=0).index)


@lru_cache(maxsize=1)
def _emission_unit_prefixes() -> Dict[str, str]:
    """Species -> unit prefix ('Gt', 'Mt', 'kt') as declared by FaIR."""
    try:
        from fair.structure.units import desired_emissions_units
    except ImportError:
        return {}
    return {
        species: unit.split(" ", 1)[0]
        for species, unit in desired_emissions_units.items()
    }


def _unit_factor(species: str) -> float:
    """kg -> FAIR emission unit for a species."""
    prefix = _emission_unit_prefixes().get(species)
    if prefix in _KG_PER_UNIT_PREFIX:
        return _KG_PER_UNIT_PREFIX[prefix]
    return _KG_TO_FAIR_CO2 if species.startswith("CO2") else _KG_TO_FAIR_OTHER


class _Template(NamedTuple):
    """Everything a FAIR run needs that does not depend on the perturbation."""

    species: list
    properties: dict
    configs: list
    specie_axis: list
    emissions: np.ndarray  # (timepoints, specie), 1750..end
    forcing: np.ndarray  # (timebounds, specie)
    species_configs: object  # xr.Dataset
    climate_configs: object  # xr.Dataset
    ebms: object  # xr.Dataset, energy-balance matrices + seeded noise


def _override_defaults(f, params_csv: str) -> None:
    """Vectorized stand-in for ``FAIR.override_defaults``.

    FaIR's own version fills the calibration parameters with a Python loop
    over configs x columns (841 x 86 xarray writes, ~6 s). The columns are
    ``name`` or ``name[index]``, where the index is a layer for climate
    configs and a specie for species configs; here each column is written in
    one array assignment instead. Verified to give the same values as
    ``override_defaults`` for calibration 1.4.1.
    """
    import pandas as pd
    from fair.io.param_sets import energy_balance_parameters

    df_configs = pd.read_csv(params_csv, index_col=0).reindex(f.configs)
    specie_position = {name: i for i, name in enumerate(f.species)}

    for col in df_configs.columns:
        param_name, _, param_index = col.partition("[")
        param_index = param_index[:-1] if param_index else None
        values = df_configs[col].to_numpy()

        if param_name in energy_balance_parameters:
            target = f.climate_configs[param_name]
        else:
            if param_index is not None and param_index not in specie_position:
                continue
            target = f.species_configs[param_name]

        if param_index is None:
            target.values[...] = values
        elif param_name in energy_balance_parameters:
            target.values[:, int(param_index)] = values  # (config, layer)
        else:
            target.values[:, specie_position[param_index]] = values  # (config, specie)


class _Configs(NamedTuple):
    species: list
    properties: dict
    configs: list
    species_configs: object
    climate_configs: object
    ebms: object


@lru_cache(maxsize=2)
def _calibrated_configs(end_year: int) -> _Configs:
    """Calibrated species/climate configs and energy-balance matrices.

    ``override_defaults`` is a Python loop over 841 configs (~6 s) and
    ``_make_ebms`` builds one energy-balance model per config. Neither depends
    on the SSP marker, so switching background scenarios does not redo them.
    """
    fair = require_fair()
    params_csv, props_csv = get_calibration_files()
    configs = _config_index(params_csv)
    species, properties = fair.io.read_properties(filename=props_csv)

    f = fair.FAIR(ch4_method=_CH4_METHOD)
    f.define_time(_SIM_START_YEAR, end_year, 1)
    f.define_scenarios(["calibration"])
    f.define_configs(configs)
    f.define_species(species, properties)
    f.allocate()
    f.fill_species_configs(props_csv)
    _override_defaults(f, params_csv)  # calibrated climate/species configs
    with warnings.catch_warnings():
        # Same filter FaIR applies around this call inside `run()`.
        warnings.filterwarnings(
            "ignore", category=RuntimeWarning, module="scipy.stats._multivariate"
        )
        f._make_ebms()  # energy-balance matrices and seeded internal variability
    return _Configs(
        species=species,
        properties=properties,
        configs=configs,
        species_configs=f.species_configs,
        climate_configs=f.climate_configs,
        ebms=f.ebms,
    )


@lru_cache(maxsize=8)
def _background(marker: str, end_year: int):
    """RCMIP background emissions and prescribed forcing for one marker.

    The RCMIP background is the same for every config, so only one config's
    slice is kept and broadcast when a run is set up. Concentrations are not
    taken from RCMIP - see :func:`_initialise_preindustrial`.
    """
    fair = require_fair()
    cfg = _calibrated_configs(end_year)
    f = fair.FAIR(ch4_method=_CH4_METHOD)
    f.define_time(_SIM_START_YEAR, end_year, 1)
    f.define_scenarios([marker])
    f.define_configs(cfg.configs)
    f.define_species(cfg.species, cfg.properties)
    f.allocate()
    f.fill_from_rcmip()  # SSP marker background emissions
    return (
        f.emissions.values[:, 0, 0, :].copy(),
        f.forcing.values[:, 0, 0, :].copy(),
        list(f.emissions.specie.values),
    )


def _template(marker: str, end_year: int) -> _Template:
    """Everything a run needs: calibrated configs plus the marker background.

    Both halves are cached, so the expensive setup happens once per process
    (and the marker-independent half survives a scenario switch).
    """
    cfg = _calibrated_configs(end_year)
    emissions, forcing, specie_axis = _background(marker, end_year)
    return _Template(
        species=cfg.species,
        properties=cfg.properties,
        configs=cfg.configs,
        specie_axis=specie_axis,
        emissions=emissions,
        forcing=forcing,
        species_configs=cfg.species_configs,
        climate_configs=cfg.climate_configs,
        ebms=cfg.ebms,
    )


def _noop():
    """Stand-in for ``FAIR._make_ebms`` when the matrices are already set."""


def _new_run(tmpl: _Template, n_scenarios: int, start_year: int, end_year: int):
    """Allocate a FAIR object and fill it from the cached template."""
    fair = require_fair()
    f = fair.FAIR(ch4_method=_CH4_METHOD)
    f.define_time(start_year, end_year, 1)
    f.define_scenarios([f"run{i}" for i in range(n_scenarios)])
    f.define_configs(tmpl.configs)
    f.define_species(tmpl.species, tmpl.properties)
    f.allocate()

    off = start_year - _SIM_START_YEAR
    n_tp = f.emissions.shape[0]
    n_tb = f.forcing.shape[0]
    f.emissions.values[...] = tmpl.emissions[off : off + n_tp, None, None, :]
    f.forcing.values[...] = tmpl.forcing[off : off + n_tb, None, None, :]
    f.species_configs = tmpl.species_configs.copy(deep=True)
    f.climate_configs = tmpl.climate_configs.copy(deep=True)
    # Reuse the cached energy-balance matrices, including the seeded internal
    # variability draws: sliced to this window they are the very same noise
    # realization the full 1750-onwards run would have used, so it still
    # cancels when perturbed and baseline runs are differenced.
    f.ebms = tmpl.ebms.isel(timebounds=slice(off, off + n_tb))
    f._make_ebms = _noop
    return f


def _initialise_preindustrial(f) -> None:
    fair = require_fair()
    # Every specie is emissions-driven or calculated, so the concentration
    # array is an output; only its first timebound is an initial condition,
    # and it comes from the calibration (RCMIP leaves calculated species NaN).
    fair.interface.initialise(f.concentration, f.species_configs["baseline_concentration"])
    fair.interface.initialise(f.forcing, 0)
    fair.interface.initialise(f.temperature, 0)
    fair.interface.initialise(f.cumulative_emissions, 0)
    fair.interface.initialise(f.airborne_emissions, 0)


@lru_cache(maxsize=4)
def _spinup_state(marker: str, end_year: int, start_year: int) -> Optional[Dict]:
    """Run 1750 -> ``start_year`` once and cache the climate state there.

    The spin-up is identical for the baseline and every perturbation, so the
    actual runs restart from this state and only simulate the reporting
    period. FAIR does not expose the internal Cummins state vector, so a
    restarted run re-initialises it from the forcing at ``start_year``; that
    leaves a small transient in the *absolute* trajectory which cancels out in
    the perturbed-minus-baseline differences this module reports (verified at
    <1e-3 relative on the response).
    """
    if start_year <= _SIM_START_YEAR:
        return None
    tmpl = _template(marker, end_year)
    f = _new_run(tmpl, 1, _SIM_START_YEAR, start_year)
    _initialise_preindustrial(f)
    with _RUN_LOCK:
        f.run(progress=False)
    return {
        "concentration": f.concentration.values[-1, 0].copy(),
        "forcing": f.forcing.values[-1, 0].copy(),
        "temperature": f.temperature.values[-1, 0].copy(),
        "airborne_emissions": f.airborne_emissions.values[-1, 0].copy(),
        "cumulative_emissions": f.cumulative_emissions.values[-1, 0].copy(),
        "gas_partitions": f.gas_partitions.values[0].copy(),
    }


def _apply_state(f, state: Optional[Dict]) -> None:
    if state is None:
        _initialise_preindustrial(f)
        return
    for i in range(f.emissions.shape[1]):  # every scenario starts here
        f.concentration.values[0, i] = state["concentration"]
        f.forcing.values[0, i] = state["forcing"]
        f.temperature.values[0, i] = state["temperature"]
        f.airborne_emissions.values[0, i] = state["airborne_emissions"]
        f.cumulative_emissions.values[0, i] = state["cumulative_emissions"]
        f.gas_partitions.values[i] = state["gas_partitions"]


def _perturbation_scale(tmpl, perturbation, years) -> float:
    """Factor that lifts a perturbation out of the float64 noise floor.

    Each species is scaled to at most :data:`_LINEARIZATION_TARGET` of its own
    peak background emission over the run window, and one common factor (the
    smallest) is used so a multi-species perturbation stays self-consistent.
    Perturbations that are already that large are left alone.
    """
    if not perturbation:
        return 1.0
    offset = int(years[0]) - _SIM_START_YEAR
    window = slice(offset, offset + len(years))
    factor = np.inf
    for species_name, delta in perturbation.items():
        if species_name not in tmpl.specie_axis:
            continue
        peak = np.max(np.abs(np.asarray(delta, dtype="float64"))) * _unit_factor(
            species_name
        )
        if peak == 0:
            continue
        si = tmpl.specie_axis.index(species_name)
        background = np.nanmax(np.abs(tmpl.emissions[window, si]))
        if not np.isfinite(background) or background == 0:
            continue
        factor = min(factor, _LINEARIZATION_TARGET * background / peak)
    if not np.isfinite(factor):
        return 1.0
    return max(float(factor), 1.0)


def _apply_perturbation(f, tmpl, scenario_index, perturbation, years, scale) -> None:
    """Add signed per-species emission deltas (kg) to one scenario."""
    if not perturbation:
        return
    start_year = int(years[0])
    n_timepoints = f.emissions.shape[0]
    emissions = f.emissions.values  # (timepoints, scenario, config, specie)
    for species_name, delta in perturbation.items():
        if species_name not in tmpl.specie_axis:
            continue
        si = tmpl.specie_axis.index(species_name)
        scaled = np.asarray(delta, dtype="float64") * (
            _unit_factor(species_name) * scale
        )
        tp = np.asarray(years, dtype="int64") - start_year
        inside = (tp >= 0) & (tp < n_timepoints)
        np.add.at(
            emissions[:, scenario_index, :, si],
            tp[inside],
            scaled[inside][:, None],
        )


def _batch_size(n_timebounds: int, n_configs: int, n_species: int) -> int:
    per_scenario = n_timebounds * n_configs * n_species * 8 * _ARRAYS_PER_SCENARIO
    return max(1, int(_BATCH_MEMORY_BUDGET // max(per_scenario, 1)))


_RESULT_CACHE: "OrderedDict[tuple, Dict[str, np.ndarray]]" = OrderedDict()


def _cache_key(marker, years, perturbation) -> tuple:
    if not perturbation:
        return (marker, int(years[0]), int(years[-1]), None)
    digest = hashlib.blake2b(digest_size=16)
    for species_name in sorted(perturbation):
        digest.update(species_name.encode())
        digest.update(
            np.ascontiguousarray(perturbation[species_name], dtype="float64").tobytes()
        )
    return (marker, int(years[0]), int(years[-1]), digest.hexdigest())


def _cache_put(key, value) -> None:
    _RESULT_CACHE[key] = value
    while len(_RESULT_CACHE) > _RESULT_CACHE_SIZE:
        _RESULT_CACHE.popitem(last=False)


def clear_caches() -> None:
    """Drop all cached FAIR templates, spin-up states and results."""
    _calibrated_configs.cache_clear()
    _background.cache_clear()
    _spinup_state.cache_clear()
    _config_index.cache_clear()
    get_calibration_files.cache_clear()
    _RESULT_CACHE.clear()


def run_perturbations(
    marker: str,
    years: np.ndarray,
    perturbations: Sequence[Optional[Dict[str, np.ndarray]]],
    max_batch: Optional[int] = None,
) -> List[Dict[str, np.ndarray]]:
    """Run FAIR for several emission perturbations of the same background.

    Parameters
    ----------
    marker : str
        FAIR SSP marker scenario (e.g. ``"ssp245"``).
    years : np.ndarray
        Reporting calendar years (ascending, consecutive). The model restarts
        from the cached ``years[0]`` state and runs to ``years[-1]``.
    perturbations : sequence of (dict | None)
        One entry per run. Each maps a FAIR species name to a
        ``(len(years),)`` array of signed emission deltas **in kg**, aligned
        to ``years``. ``None`` or ``{}`` is the unperturbed background.
    max_batch : int, optional
        Cap on scenarios simulated in a single FAIR call. Defaults to what
        fits in the internal memory budget.

    Returns
    -------
    list of dict
        Aligned with ``perturbations``; each dict holds both
        ``"radiative_forcing"`` and ``"temperature"`` as
        ``(n_configs, len(years))`` ensemble arrays, plus ``"scale"``: the
        factor the perturbation was multiplied by for numerical conditioning
        (see :data:`_LINEARIZATION_TARGET`). Divide the
        perturbed-minus-baseline difference by it to get the response to the
        perturbation as given.
    """
    require_fair()
    years = np.asarray(years)
    end_year = int(years[-1])
    start_year = int(years[0])

    keys = [_cache_key(marker, years, p) for p in perturbations]
    results: List[Optional[Dict[str, np.ndarray]]] = [
        _RESULT_CACHE.get(key) for key in keys
    ]
    todo = [i for i, res in enumerate(results) if res is None]
    if not todo:
        return results

    tmpl = _template(marker, end_year)
    state = _spinup_state(marker, end_year, start_year)
    n_timebounds = len(years)
    batch = max_batch or _batch_size(
        n_timebounds, len(tmpl.configs), len(tmpl.specie_axis)
    )
    rows = np.asarray(years, dtype="int64") - start_year

    for chunk_start in range(0, len(todo), batch):
        chunk = todo[chunk_start : chunk_start + batch]
        f = _new_run(tmpl, len(chunk), start_year, end_year)
        _apply_state(f, state)
        scales = [
            _perturbation_scale(tmpl, perturbations[pert_index], years)
            for pert_index in chunk
        ]
        for i, pert_index in enumerate(chunk):
            _apply_perturbation(
                f, tmpl, i, perturbations[pert_index], years, scales[i]
            )
        with _RUN_LOCK:
            f.run(progress=False)
        temperature = f.temperature.values  # (timebounds, scenario, config, layer)
        forcing_sum = f.forcing_sum.values  # (timebounds, scenario, config)
        for i, pert_index in enumerate(chunk):
            result = {
                "temperature": np.asarray(
                    temperature[rows, i, :, 0], dtype="float64"
                ).T,
                "radiative_forcing": np.asarray(
                    forcing_sum[rows, i, :], dtype="float64"
                ).T,
                "scale": scales[i],
            }
            results[pert_index] = result
            _cache_put(keys[pert_index], result)
        del f

    return results


def run_fair(
    marker: str,
    perturbation_by_species: Optional[Dict[str, np.ndarray]],
    years: np.ndarray,
    output: str = "radiative_forcing",
) -> np.ndarray:
    """
    Run FAIR for a marker scenario, optionally adding a per-species emission
    perturbation, and return an ``(n_configs, n_years)`` array of the
    requested output for the reporting ``years``.

    Thin wrapper around :func:`run_perturbations` for a single run; see there
    for how the background, spin-up and results are cached. Note that the
    absolute output belongs to a perturbation scaled by the run's ``"scale"``
    factor; :func:`run_perturbations` is the better entry point when the
    difference against the baseline is what you need.
    """
    if output not in OUTPUTS:
        raise ValueError(
            f"output must be 'radiative_forcing' or 'temperature', not {output!r}"
        )
    return run_perturbations(marker, years, [perturbation_by_species])[0][output]
