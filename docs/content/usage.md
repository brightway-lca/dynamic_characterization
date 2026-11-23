---
icon: lucide/footprints
tags:
  - usage
  - dynamic characterization
---

# Usage

The workflow to use this package could look like this:

```python
import pandas as pd
from dynamic_characterization import characterize
from dynamic_characterization.ipcc_ar6 import characterize_co2, characterize_ch4

# defining a dummy dynamic inventory that you somehow got
dynamic_inventory_df = pd.DataFrame(
        data={
            "date": pd.Series(
                data=[
                    "15-12-2020",
                    "20-12-2020",
                    "25-05-2022",
                ],
                dtype="datetime64[s]",
            ),
            "amount": pd.Series(data=[10.0, 20.0, 50.0], dtype="float64"),
            "flow": pd.Series(data=[1, 1, 3], dtype="int"),
            "activity": pd.Series(data=[2, 2, 4], dtype="int"),
        }
    )

df_characterized = characterize(
        dynamic_inventory_df,
        metric="radiative_forcing", # could also be GWP
        characterization_functions={
            1: characterize_co2,
            3: characterize_ch4,
        },
        time_horizon=2,
    )
```

If you use this package with [Brightway](https://docs.brightway.dev/en/latest/), stuff can get even easier: if you have an impact assessment method at hand, you can pass it to the characterize function via the `base_lcia_method` attribute and we'll try to automatically match the flows that are characterized in that method to the flows we have characterization functions for. This matching is based on the names or the CAS numbers, depending on the flow. The function call could look like this then:

```python
method = ('EF v3.1', 'climate change', 'global warming potential (GWP100)')

df_characterized = characterize(
        dynamic_inventory_df,
        metric="radiative_forcing", # could also be GWP
        base_lcia_method=method,
        time_horizon=2,
)
```

## Available Metrics

The `characterize()` function supports the following metrics:

| Metric | Description | Module |
|--------|-------------|--------|
| `radiative_forcing` | Radiative forcing time series (W/m²) | ipcc_ar6 |
| `GWP` | Global Warming Potential (kg CO2eq) | ipcc_ar6 |
| `prospective_radiative_forcing` | Prospective radiative forcing time series (W/m²) | prospective |
| `pGWP` | Prospective GWP with scenario-based RE and CO2 IRF (kg CO2eq) | prospective |
| `pGTP` | Prospective GTP with scenario-based RE and CO2 IRF (kg CO2eq) | prospective |

For the prospective metrics (`pGWP`, `pGTP`), you need to set a background scenario first. See [Prospective Characterization](prospective.md) for details.

## IPCC AR6 parameter set

By default, the characterization functions are those of dynamic_characterization ≤ 1.4
(`parameter_set="legacy"`). They give lower GWPs than the IPCC AR6 metrics, because they leave
out the climate–carbon feedback of non-CO2 gases and the CO2 from fossil-methane oxidation, and
because their annual forcing only adds up to AGWP(H − 1). The `"ar6"` parameter set reproduces
the AR6 metrics (WG1 Ch. 7, Tables 7.15 and 7.SM.7):

```python
characterize(
    dynamic_inventory_df,
    metric="GWP",
    base_lcia_method=method,
    characterization_options={"parameter_set": "ar6"},
)
```

| flow | static (IPCC 2021) | `legacy` | `ar6` |
|---|---|---|---|
| Methane, fossil | 29.8 | 26.7 | 29.8 |
| Methane, non-fossil | 27.0 | 26.7 | 27.0 |
| Dinitrogen monoxide | 273 | 263.6 | 273 |
| Tetrafluoromethane | 7380 | 7058 | 7380 |
| Carbon monoxide, fossil | 4.06 | 1.57 | 3.20 |
| Carbon monoxide, non-fossil | 2.49 | 1.57 | 1.63 |
| Hydrogen | 11.6 | – | 11.4 |
| VOC, volatile organic compounds | 4.5 | – | 4.50 |

GWP100 per kg. All halocarbons of AR6 Table 7.SM.7 are reproduced within the table's rounding,
except SF6, which follows Table 7.15 (lifetime 1000 yr, GWP100 24 300) like the static methods.
Against the ecoinvent 3.12 "IPCC 2021, climate change: total (excl. biogenic CO2)" GWP100 method,
184 of its 189 flows agree within 1 %; the other five have factors rounded to one significant
digit (butane 0.006, tetrafluoroethene 0.004).
In the `"ar6"` set:

- non-CO2 gases include the climate–carbon feedback (Gasser et al. 2017);
- fossil CH4 adds the CO2 from its oxidation; non-fossil CH4 follows the AR6 convention;
- N2O includes its effects on ozone and on methane;
- CFC-11 and CFC-12 get the +12 % tropospheric adjustment;
- the forcing stamped at `date + k` years is the forcing of the year starting then, so the
  annual values add up to AGWP(H). The CO2 reference of the GWP metric uses the same convention.

Flows are matched by CAS number first and then by exact name. "Methane, from soil or biomass
stock" counts as fossil, as in the ecoinvent IPCC 2021 method.

Options (all are keys of `characterization_options`; see `CharacterizationOptions`):

| option | default | meaning |
|---|---|---|
| `parameter_set` | `"legacy"` | `"legacy"` or `"ar6"` |
| `climate_carbon_feedback` | `True` | climate–carbon feedback for non-CO2 gases |
| `ccf_parameters` | `"gasser2017_mean"` | or `"gasser2017_upper"` / `"gasser2017_lower"` (±1 SD) |
| `biogenic_oxidation_co2` | `"ar6"` | non-fossil CH4: `"ar6"` (27.0), `"neutral"` (CH4 pulse only, 27.9) or `"explicit"` (oxidation CO2 counted, uptake modelled in the inventory, 29.8) |
| `ch4_oxidation_yield` | 0.75 | fraction of CH4 oxidised to CO2 |
| `unspecified_origin` | `"fossil"` | origin of CH4 flows without a fossil/non-fossil suffix |
| `flow_overrides` | `{}` | `{flow id or (database, code): {"species": ..., "origin": ...}}` |

`compare_with_static(method, time_horizon, characterization_options)` returns a table of the
static factor and the dynamic GWP of every flow of a method.

Without `characterization_options`, a `FutureWarning` announces that the default will change to
`"ar6"` in version 2.0. With the prospective metrics, the `"ar6"` set follows the chosen
scenario; see [Prospective Characterization](prospective.md#with-the-ar6-parameter-set).

### Indirect gases (CO, VOC, H2, NOx)

AR6 gives no per-kg metrics for ozone precursors, so the `"ar6"` set takes published parameter
sets (`ipcc_ar6/data/precursor_parameters.json`; every value with its table and row, nothing is
fitted). Each forcing component (short- and long-term ozone, CH4, sulfate, stratospheric water
vapour, nitrate) follows the two-time-scale form of Fuglestvedt et al. (2010), with the precursor
lifetime of Warwick et al. (2023) where it matters (H2).

| species | default set | other sets | GWP100 per kg (default) |
|---|---|---|---|
| CO | `fry2012_4regions` (11 models) | `derwent2001`, `berntsen2005_uio`, `berntsen2005_lmdz`, `fuglestvedt2010_mean3` | 1.63 (+1.57 oxidation CO2 if fossil) |
| VOC | `fry2012_4regions` (per kg C; 0.845 kg C per kg VOC) | `collins2002` (per kg VOC) | 4.50 |
| H2 | `warwick2023` | `sand2023` | 11.4 (Sand: 11.7) |
| NOx (opt-in) | `fry2012_4regions` (per kg N) | – | −3.8 per kg NO, −2.5 per kg NO2 |
| nitrate (opt-in) | `bauer2007` (per kg N) | – | −10.4 per kg NO additionally |

- The default CO value is lower than the static 2.49 of the IPCC 2021 methods, which is the
  midpoint of the older single-model range of AR5 (Fuglestvedt et al. 2010). The multi-model mean
  of Fry et al. (2012) is the more robust basis; `fuglestvedt2010_mean3` gives 2.45.
- Non-fossil CO is treated as carbon neutral (no oxidation CO2), as in the static factors, unless
  `biogenic_oxidation_co2="explicit"`.
- NO, NO2 and "Nitrogen oxides" are characterized only with `include_nox=True`: the sign and size
  of the NOx metric are very uncertain (AR5 Table 8.A.3). The static −10.8 of "Nitric oxide"
  corresponds to the per-kg-N literature value; per kg NO it is about −4.
- The sources do not add the climate–carbon feedback; `precursor_ccf=True` adds it (+5–6 %).

| option | default | meaning |
|---|---|---|
| `include_precursors` | `["CO", "VOC", "H2"]` | precursors with indirect effects; CO not listed falls back to oxidation CO2 only |
| `precursor_parameter_sets` | see above | `{species: set}` for "CO", "VOC", "H2", "NOx", "nitrate" |
| `include_nox` | `False` | characterize NO, NO2, "Nitrogen oxides" |
| `nox_nitrate` | `False` | add nitrate aerosol (needs `include_nox`) |
| `precursor_ccf` | `False` | climate–carbon feedback for the precursors |
