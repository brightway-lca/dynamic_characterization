# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]
* Add optional FAIR climate-model path with `fair_radiative_forcing` (ΔRF, W/m²) and `fair_temperature` (ΔT, K) metrics
* FAIR metrics return probabilistic ensemble results with per-flow/activity attribution and configurable quantiles (default: 2.5, 25, 50, 75, 97.5 percentiles)
* Select a FAIR-native SSP background with `set_fair_scenario(ssp, rcp)` (8 markers, IAM-agnostic); the 4 IAM scenarios that match a FAIR marker also support FAIR via `set_scenario`
* Add `available_fair_scenarios()` and `available_scenarios(metric)` to programmatically list scenarios supporting a given metric or metric family (`"fair"`, `"prospective"`, `"pGWP"`, etc.)
* FAIR runs use the calibrated, constrained AR6 ensemble (`calibration1.4.1`, 841 members); the calibration data is downloaded and cached (via `pooch`) on first use when not vendored in the `fair` install (verified against `fair` 2.2.4)
* Make the FAIR path fast enough for real dynamic inventories: the calibrated ensemble setup and the SSP background are cached separately, the 1750-onwards spin-up is run once and reused (runs restart from it), species perturbations are batched over FAIR's scenario axis, and results are memoized per (marker, period, perturbation) so both ΔRF and ΔT come from the same runs - a 80k-row `bw_timex` inventory went from ~92 s to ~7 s per metric (and ~0.4 s for the second one)
* Vectorize the FAIR inventory aggregation and per-flow attribution (no Python loop over inventory rows; ensemble quantiles are taken once per species and scaled per flow, mirrored for uptake flows)
* Add `dynamic_characterization.fair.runner.clear_caches()` to drop cached FAIR templates, spin-up states and results
* Fix a numerical-conditioning error in the FAIR metrics: at LCA magnitudes (kilograms against a gigatonne background) the perturbed-minus-baseline difference lost most of its significant digits to float64 round-off, which put the CO2 part of the response up to ~20% off. Perturbations are now scaled to a fixed small fraction of the species' background emission and the response scaled back (verified linear over four further orders of magnitude; results now within ~0.3% of the converged value)
* Rewrite the FAIR species map and match flows exactly instead of by substring. Substring matching routed every `*methane` flow (Tetrafluoromethane, Trifluoromethane, Dichloromethane, Dichlorodifluoromethane, Bromotrifluoromethane, ...) and `NMVOC, non-methane volatile organic compounds` to CH4, so halocarbons were characterized as methane and the VOC species was never used
* Map the F-gases and ozone-depleting substances (SF6, NF3, CF4, C2F6, the HFCs, CFCs, HCFCs and halons) to their own FAIR species, with the flow's CAS number as a fallback when the name is unknown (other ecoinvent versions, SimaPro). FAIR expects halogenated species in kt/yr rather than Mt/yr; the unit per species is now taken from FAIR itself instead of assumed
* Route biogenic CO2 to `CO2 AFOLU` instead of `CO2 FFI` and give uptake flows their proper sign: `Carbon dioxide, non-fossil, resource correction` was counted as an emission, which (with biogenic CO2 mixed into the fossil pool) overstated net CO2 by ~8% for a typical ecoinvent inventory
* Convert flows reported on a different mass basis than FAIR expects (nitric oxide to NO2, sulfur trioxide to SO2)
* Log the flows that no FAIR species covers, so silent under-coverage is visible

## [1.4.2] - (2026-08-14)
* The error raised when a prospective metric is calculated without a scenario now explains what to do: full import path, a copy-pasteable `set_scenario` call, and the available IAM-SSP-RCP combinations.
* Fixed `TypeError: ... got an unexpected keyword argument 'time_varying_re'` for the `pGWP` and `pGTP` metrics: `time_varying_re` was passed to every characterization function, including the IPCC AR6 fallback functions (CO and the GHGs from `decay_multipliers.json`), which don't accept it.

## [1.4.1] - (2026-08-03)
* Fixed a regression introduced in 1.4.0 where ordinary CO2 emission flows (e.g. `Carbon dioxide, fossil` emitted to any air subcategory) got no default characterization function and were silently skipped, understating dynamic climate scores by roughly an order of magnitude. `characterize_uptake=False` now only suppresses the uptake functions instead of all CO2 characterization.

## [1.4.0] - (2026-05-17)
* Add caching
* Vectorize radiative forcing calculations

## [1.3.1] - (2026-01-29)
* Fixed packaging issue

## [1.3.0] - (2026-01-29)
* Add new prospective module for prospective-dynamic LCIA based on Watanabe et al. (2026)

## [1.2.0] - (2025-10-08)
* Instead of not characterizing non-fossil methane, treat it the same as fossil methane. Discussion at https://github.com/brightway-lca/dynamic_characterization/issues/15.

## [1.1.1] - (2025-02-27)
* Use loguru for logging.

## [1.1.0] - (2025-02-07)
* Renamed `characterization_function_dict` to `characterization_functions`.

## [1.0.3] - (2024-09-27)
* Fixed path to default characterization functions

## [1.0.2] - (2024-09-19)
* Fixed paths for data files

## [1.0.1] - (2024-09-19)
* Fixed packaging issue

## [1.0.0] - (2024-09-19)
* Renamed submodules:
    * `dynamic_characterization.timex` -> `dynamic_characterization.ipcc_ar6`
    * `dynamic_characterization.temporalis` -> `dynamic_characterization.original_temporalis_functions`
* Renamed function to apply characterization functions to `dynamic_characterization.characterize()`

## [0.0.4] - (2024-07-17)
* Fixed an issue with the path of data files

## [0.0.3] - (2024-07-17)
* Added characterize_dynamic_inventory function that can directly apply functions to a dynamic inventory dataframe. This also includes a function to add a set of default characterization functions based on the CAS-numbers of bioflows
* Improves computational perfomance by using arrays and namedtuple instead of pd.Series

## [0.0.2] - (2024-07-11)
* Version bump to harmonize with conda package version.

## [0.0.1] - (2024-06-17)
* Initial release
