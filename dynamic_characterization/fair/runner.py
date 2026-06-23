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
"""

import glob
import os
import threading
from functools import lru_cache
from typing import Dict, Optional, Tuple

import numpy as np

from . import FAIR_IMPORT_ERROR_MSG

_RUN_LOCK = threading.Lock()

# Simulation always starts here so the climate response spins up from
# pre-industrial; reporting years are sliced out afterwards.
_SIM_START_YEAR = 1750

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

# kg -> FAIR emission unit. CO2 species are in Gt CO2/yr (1 Gt = 1e12 kg);
# every other supported species is in Mt/yr (1 Mt = 1e9 kg).
_KG_TO_FAIR_CO2 = 1e-12
_KG_TO_FAIR_OTHER = 1e-9


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


def _unit_factor(species: str) -> float:
    """kg -> FAIR emission unit for a species."""
    return _KG_TO_FAIR_CO2 if species.startswith("CO2") else _KG_TO_FAIR_OTHER


def run_fair(
    marker: str,
    perturbation_by_species: Optional[Dict[str, np.ndarray]],
    years: np.ndarray,
    output: str = "radiative_forcing",
) -> np.ndarray:
    """
    Run FAIR for a marker scenario over 1750..max(years), optionally adding a
    per-species emission perturbation, and return an ``(n_configs, n_years)``
    array of the requested output for the reporting ``years``.

    Parameters
    ----------
    marker : str
        FAIR SSP marker scenario (e.g. ``"ssp245"``).
    perturbation_by_species : dict | None
        Maps a FAIR species name to a ``(len(years),)`` array of signed
        emission deltas **in kg**, aligned to ``years``. ``None`` runs the
        unperturbed background.
    years : np.ndarray
        Reporting calendar years (ascending). The simulation runs from 1750
        to ``years[-1]``; outputs are returned for exactly these years.
    output : str
        ``"radiative_forcing"`` (total forcing, W/m^2) or ``"temperature"``
        (surface layer anomaly, K).

    Returns
    -------
    np.ndarray
        ``(n_configs, len(years))`` ensemble output.
    """
    if output not in ("radiative_forcing", "temperature"):
        raise ValueError(
            f"output must be 'radiative_forcing' or 'temperature', not {output!r}"
        )
    fair = require_fair()
    params_csv, props_csv = get_calibration_files()

    end_year = int(years[-1])

    f = fair.FAIR()
    f.define_time(_SIM_START_YEAR, end_year, 1)
    f.define_scenarios([marker])
    f.define_configs(_config_index(params_csv))

    species, properties = fair.io.read_properties(filename=props_csv)
    f.define_species(species, properties)
    f.allocate()
    f.fill_species_configs(props_csv)
    f.fill_from_rcmip()  # SSP marker background emissions
    f.override_defaults(params_csv)  # calibrated climate/species configs

    if perturbation_by_species:
        specie_axis = list(f.emissions.specie.values)
        n_timepoints = f.emissions.sizes["timepoints"]
        emissions = f.emissions.values  # (timepoints, scenario, config, specie)
        for species_name, delta in perturbation_by_species.items():
            if species_name not in specie_axis:
                continue
            si = specie_axis.index(species_name)
            factor = _unit_factor(species_name)
            for i, year in enumerate(years):
                tp = int(year) - _SIM_START_YEAR
                if 0 <= tp < n_timepoints:
                    emissions[tp, 0, :, si] += float(delta[i]) * factor

    fair.interface.initialise(f.forcing, 0)
    fair.interface.initialise(f.temperature, 0)
    fair.interface.initialise(f.cumulative_emissions, 0)
    fair.interface.initialise(f.airborne_emissions, 0)

    with _RUN_LOCK:
        f.run(progress=False)

    if output == "temperature":
        arr = f.temperature.sel(scenario=marker, layer=0).values  # (tb, config)
    else:
        arr = f.forcing_sum.sel(scenario=marker).values  # (tb, config)

    # Slice the timebounds (calendar years) for the requested reporting years.
    rows = [int(year) - _SIM_START_YEAR for year in years]
    arr = np.asarray(arr, dtype="float64")[rows, :]
    return arr.T  # (n_configs, n_years)
