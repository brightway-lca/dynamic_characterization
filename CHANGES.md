# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

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
