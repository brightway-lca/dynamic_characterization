# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

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
