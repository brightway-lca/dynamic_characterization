"""Locate FAIR calibration data and run the model (requires `fair`)."""

import glob
import os
import threading
from typing import Dict, Optional, Tuple

import numpy as np

from . import FAIR_IMPORT_ERROR_MSG

_RUN_LOCK = threading.Lock()

_PARAM_CANDIDATES = (
    "calibrated_constrained_parameters_calibration1.4.1.csv",
    "calibrated_constrained_parameters.csv",
)
_PROP_CANDIDATES = (
    "species_configs_properties_calibration1.4.1.csv",
    "species_configs_properties.csv",
)


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
    """Locate (parameters_csv, properties_csv) inside the installed fair pkg."""
    fair = require_fair()
    return _find_in_fair(fair, _PARAM_CANDIDATES), _find_in_fair(
        fair, _PROP_CANDIDATES
    )


def run_fair(
    marker: str,
    perturbation_by_species: Optional[Dict[str, np.ndarray]],
    years: np.ndarray,
    output: str = "radiative_forcing",
) -> np.ndarray:
    """
    Run FAIR for a marker scenario, optionally adding a per-species
    emission perturbation, and return an (n_configs, n_years) array of the
    requested output ('radiative_forcing' total forcing or 'temperature').

    The actual f.run() is serialized under a module lock.
    """
    if output not in ("radiative_forcing", "temperature"):
        raise ValueError(
            f"output must be 'radiative_forcing' or 'temperature', not {output!r}"
        )
    fair = require_fair()
    params_csv, props_csv = find_calibration_files()

    f = fair.FAIR()
    f.define_time(int(years[0]), int(years[-1]), 1)
    f.define_scenarios([marker])

    if props_csv is not None:
        species, properties = fair.io.read_properties(filename=props_csv)
    else:  # pragma: no cover - depends on fair install
        species, properties = fair.io.read_properties()
    f.define_species(species, properties)

    if params_csv is not None:
        import pandas as pd

        cfg = pd.read_csv(params_csv, index_col=0)
        f.define_configs(list(cfg.index))
    else:  # pragma: no cover
        f.define_configs(["default"])

    f.allocate()
    if props_csv is not None:
        f.fill_species_configs(filename=props_csv)
    else:  # pragma: no cover
        f.fill_species_configs()
    # Background SSP emissions for the marker.
    f.fill_from_rcmip()
    if params_csv is not None:
        f.override_defaults(params_csv)

    if perturbation_by_species:
        for sp, delta in perturbation_by_species.items():
            # delta is (n_years,) in fair's emission units; broadcast over configs.
            f.emissions.loc[
                dict(specie=sp, scenario=marker)
            ] = f.emissions.loc[dict(specie=sp, scenario=marker)] + delta[
                :, None, None
            ][: f.emissions.loc[dict(specie=sp, scenario=marker)].shape[0]]

    fair.interface.initialise(f.temperature, 0)
    fair.interface.initialise(f.forcing, 0)

    with _RUN_LOCK:
        f.run(progress=False)

    if output == "temperature":
        arr = f.temperature.sel(scenario=marker, layer=0).values  # (time, config)
    else:
        arr = f.forcing_sum.sel(scenario=marker).values  # (time, config)
    return np.asarray(arr, dtype="float64").T  # (n_configs, n_years)
