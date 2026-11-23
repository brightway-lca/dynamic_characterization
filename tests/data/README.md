# Test reference data

- `legacy_golden_v1.4.3.npz`: output of the legacy IPCC AR6 functions of v1.4.3, written by
  `dev/make_legacy_golden.py` (see `tests/test_legacy_unchanged.py`).
- `ar6_metrics_supplement_cleaned.csv`: IPCC AR6 WG1 Table 7.SM.7 (emission metrics of CO2,
  CH4, N2O and the halocarbons), file `data_output/7sm/metrics_supplement_cleaned.csv` of the
  AR6 code repository https://github.com/chrisroadmap/ar6 at commit
  24f279260fc2d1ebd09b2bd4114f2be11859e925. MIT licence, Copyright (c) Chris Smith. Used
  unchanged as reference values in `tests/test_ar6_core.py`.
