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

## Default characterization functions (IPCC AR6)

The default characterization functions reproduce the IPCC AR6 emission metrics (WG1 Ch. 7,
Tables 7.15 and 7.SM.7); no option is needed:

```python
characterize(dynamic_inventory_df, metric="GWP", base_lcia_method=method)
```

| flow | static (IPCC 2021) | ≤ 1.x | 2.0 |
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
digit (butane 0.006, tetrafluoroethene 0.004). The functions include:

- the climate–carbon feedback of non-CO2 gases (Gasser et al. 2017);
- the CO2 from the oxidation of fossil CH4; non-fossil CH4 follows the AR6 convention;
- the effects of N2O on ozone and methane, and the +12 % adjustment of CFC-11 and CFC-12;
- the indirect effects of CO, VOC and H2 (see below);
- annual forcing in the convention: the value stamped at `date + k` years is the forcing of the
  year starting then, so the annual values add up to AGWP(H). The CO2 reference of the GWP metric
  uses the same convention.

Flows are matched by CAS number first and then by exact name. "Methane, from soil or biomass
stock" counts as fossil, as in the ecoinvent IPCC 2021 method.

Options (all are keys of `characterization_options`; see `CharacterizationOptions`):

| option | default | meaning |
|---|---|---|
| `climate_carbon_feedback` | `True` | climate–carbon feedback for non-CO2 gases |
| `ccf_parameters` | `"gasser2017_mean"` | or `"gasser2017_upper"` / `"gasser2017_lower"` (±1 SD) |
| `biogenic_oxidation_co2` | `"ar6"` | non-fossil CH4: `"ar6"` (27.0), `"neutral"` (CH4 pulse only, 27.9) or `"explicit"` (oxidation CO2 counted, uptake modelled in the inventory, 29.8) |
| `ch4_oxidation_yield` | 0.75 | fraction of CH4 oxidised to CO2 |
| `unspecified_origin` | `"fossil"` | origin of CH4 flows without a fossil/non-fossil suffix |
| `flow_overrides` | `{}` | `{flow id or (database, code): {"species": ..., "origin": ...}}` |

`compare_with_static(method, time_horizon, characterization_options)` returns a table of the
static factor and the dynamic GWP of every flow of a method. With the prospective metrics, the
functions follow the chosen scenario; see [Prospective Characterization](prospective.md).

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

## Migrating to 2.0

Version 2.0 keeps the call syntax of 1.x (`characterize()` arguments, metric names, and the public
function names) but changes the results:

- **Default functions**: the IPCC AR6 functions above replace the simplified ones of 1.x (no
  climate–carbon feedback, no oxidation CO2 of fossil CH4, CO as oxidation only, annual forcing
  summed to H − 1). GWP100 rises, e.g. fossil CH4 26.7 → 29.8, N2O 263.6 → 273, CF4 7058 → 7380;
  H2 and VOC are now characterized, NO/NOx on request (`include_nox=True`).
- **Time convention**: the forcing of the first year is stamped at the emission date (1.x: zero at
  the emission date, the first year's forcing one year later).
- **Public names** keep working and point to the new functions: `ipcc_ar6.characterize_co2`,
  `characterize_co2_uptake`, `characterize_n2o`, `characterize_ch4` (fossil CH4) and
  `characterize_co` (fossil CO incl. indirect effects). `create_generic_characterization_function`
  is unchanged.
- **Prospective metrics**: CO2, CH4 and N2O follow Barbosa Watanabe & Cherubini (2026), with
  indirect effects and climate–carbon feedback (1.x: direct forcing only; CH4 pGWP100 rises from
  about 18–20 to about 30). `prospective.agwp`, `prospective.agtp` and
  `prospective.radiative_forcing` keep their signatures. Horizons beyond 100 years are truncated
  with a warning, as before.
- **Removed**: the 1.x functions and their data (`ipcc_ar6/data/decay_multipliers.json`).
- **New**: `characterization_options` (all optional; see the tables above),
  `compare_with_static`, the modules `ipcc_ar6.ar6_core`, `ipcc_ar6.gases`,
  `prospective.watanabe2026`, and fixes of the CAS lookup and of pGTP for gases without scenario
  data.
