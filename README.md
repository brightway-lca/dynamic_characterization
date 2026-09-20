# dynamic_characterization

[![Read the Docs](https://img.shields.io/readthedocs/dynamic-characterization?label=documentation)](https://dynamic-characterization.readthedocs.io/en/latest/)
[![PyPI - Version](https://img.shields.io/pypi/v/dynamic-characterization?color=%2300549f)](https://pypi.org/project/dynamic-characterization/)
[![Conda Version](https://img.shields.io/conda/v/diepers/dynamic_characterization?label=conda)](https://anaconda.org/diepers/dynamic_characterization)
![Conda - License](https://img.shields.io/conda/l/diepers/dynamic_characterization)
[![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/brightway-lca/dynamic_characterization/main?labpath=notebooks%2Fdynamic_characterization_demo.ipynb)

This is a package for the dynamic characterization of Life Cycle Inventories with temporal information. It includes a collection of dynamic characterization functions for various environmental flows. We also provide a simple interface to apply these functions to an existing dynamic LCI (coming from, e.g., [bw_temporalis](https://github.com/brightway-lca/bw_temporalis), [bw_timex](https://github.com/brightway-lca/bw_timex) or [optimex](https://github.com/RWTH-LTT/optimex).

The following dynamic characterization functions are currently included:

| module | impact category | metric | covered emissions | source |
|--------|-----------------|--------|-------------------|--------|
| `ipcc_ar6` | climate change | radiative forcing, GWP | 247 GHGs | radiative efficiencies & lifetimes from [IPCC AR6 Ch.7](https://www.ipcc.ch/report/ar6/wg1/chapter/chapter-7/) |
| `prospective` | climate change | prospective radiative forcing, pGWP, pGTP | CO2, CH4, N2O | [Watanabe et al. (2026)](https://pubs.acs.org/doi/10.1021/acs.est.5c12391) |
| `original_temporalis_functions` | climate change | radiative forcing | CO2, CH4 | [bw_temporalis](https://github.com/brightway-lca/bw_temporalis/tree/main) |

## Installation

```console
$ pip install dynamic-characterization
```

The same thing with [uv](https://docs.astral.sh/uv/) is `uv add dynamic-characterization`, and with conda `conda install -c conda-forge -c diepers dynamic_characterization`. See the [installation guide](https://dynamic-characterization.readthedocs.io/en/latest/content/installation/) for the details.

## Quick start

You bring a dynamic inventory - a DataFrame of `date`, `amount`, `flow`, `activity`, one row per emission at one point in time - and say which metric you want:

```python
from dynamic_characterization import characterize
from dynamic_characterization.ipcc_ar6 import characterize_co2, characterize_ch4

df_characterized = characterize(
    dynamic_inventory_df,
    metric="radiative_forcing",  # or GWP, pGWP, pGTP, prospective_radiative_forcing
    characterization_functions={1: characterize_co2, 3: characterize_ch4},
    time_horizon=100,
)
```

Each function turns one row of the inventory into the characterized time series it causes, so the result carries a value per point in time rather than a single score.

If you work with [Brightway](https://docs.brightway.dev/en/latest/), pass an impact assessment method as `base_lcia_method` instead of the function dictionary, and the flows it characterizes are matched to our functions automatically, by name or CAS number.

## Documentation

The full documentation is at [dynamic-characterization.readthedocs.io](https://dynamic-characterization.readthedocs.io/en/latest/):

- [Usage](https://dynamic-characterization.readthedocs.io/en/latest/content/usage/) - the workflow in full, and every metric you can ask for
- [Prospective characterization](https://dynamic-characterization.readthedocs.io/en/latest/content/prospective/) - scenario-based characterization factors: choosing an IAM-SSP-RCP scenario, the IPCC fallback for unsupported GHGs, and time-varying radiative efficiency
- [Examples](https://dynamic-characterization.readthedocs.io/en/latest/content/examples/) - the demo notebooks, rendered
- [API reference](https://dynamic-characterization.readthedocs.io/en/latest/api/) - including how to write your own characterization function

## Contributing

Contributions are very welcome.
To learn more, see the [Contributor Guide][Contributor Guide].

## License

Distributed under the terms of the [BSD 3-Clause license][License],
_dynamic_characterization_ is free and open source software.

## Issues

If you encounter any problems,
please [file an issue][Issue Tracker] along with a detailed description.

## Support

If you have any questions or need help, do not hesitate to contact Timo Diepers ([timo.diepers@ltt.rwth-aachen.de](mailto:timo.diepers@ltt.rwth-aachen.de))

<!-- github-only -->

[License]: https://github.com/brightway-lca/dynamic_characterization/blob/main/LICENSE
[Contributor Guide]: https://github.com/brightway-lca/dynamic_characterization/blob/main/CONTRIBUTING.md
[Issue Tracker]: https://github.com/brightway-lca/dynamic_characterization/issues
