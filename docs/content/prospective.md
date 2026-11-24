---
icon: lucide/telescope
tags:
  - prospective
  - scenario
  - gwp
---

# Prospective Characterization

The `prospective` module provides **prospective characterization factors (pCFs)** based on [Watanabe et al. (2026)](https://pubs.acs.org/doi/10.1021/acs.est.5c12391). Unlike conventional characterization factors, which use a fixed radiative efficiency and a fixed impulse response function, prospective CFs account for how atmospheric conditions change over time under different climate scenarios.

## Why Use Prospective Characterization?

Traditional characterization factors (like GWP100) assume constant atmospheric concentrations when calculating radiative forcing. However, greenhouse gas concentrations are projected to change significantly under different climate scenarios. This means:

- **Higher CO2 concentrations** → Lower radiative efficiency per additional kg CO2 (saturation effect)
- **Different carbon cycle feedbacks** → Different CO2 impulse response function → More of an emitted kg stays airborne under high-emission scenarios
- **Different scenarios** → Different atmospheric evolution → Different characterization factors

Both terms of the AGWP integral therefore depend on the scenario, not just the radiative efficiency. For CO2 the IRF is RCP-specific: the airborne fraction 100 years after the emission ranges from 0.46 under RCP2.6 to 0.66 under RCP8.5. CH4 and N2O keep fixed atmospheric lifetimes (11.8 and 109 years), so only their RE is scenario-dependent.

For LCA studies looking at future systems (e.g., infrastructure with 50+ year lifetimes), prospective characterization provides more scenario-consistent results.

## Available Metrics

| Metric | Output | Description |
|--------|--------|-------------|
| `prospective_radiative_forcing` | W/m² time series | Radiative forcing using scenario-based radiative efficiencies (and, for CO2, a scenario-based IRF) |
| `pGWP` | kg CO2eq | Prospective Global Warming Potential - integrated radiative forcing relative to CO2 |
| `pGTP` | kg CO2eq | Prospective Global Temperature Potential - endpoint temperature change relative to CO2 |

## Background Scenarios

The Watanabe module uses Integrated Assessment Model (IAM) scenarios combining SSP (Shared Socioeconomic Pathways) and RCP (Representative Concentration Pathways):

### Available Combinations

| IAM | SSP | Available RCPs |
|-----|-----|----------------|
| IMAGE | SSP1 | 2.6, 4.5, 8.5 |
| MESSAGE | SSP2 | 2.6, 4.5, 6.0, 8.5 |
| AIM | SSP3 | 4.5, 6.0, 8.5 |
| GCAM4 | SSP4 | 2.6, 4.5, 6.0, 8.5 |
| REMIND | SSP5 | 2.6, 4.5, 6.0, 8.5 |

Each IAM comes with exactly one SSP. The full list is also available at runtime as
`dynamic_characterization.prospective.VALID_SCENARIOS`.

### Setting a Scenario

You **must** set a scenario before using pGWP or pGTP metrics:

```python
import dynamic_characterization.prospective as prospective

# Set scenario - this is required before using pGWP or pGTP
prospective.set_scenario(iam="IMAGE", ssp="SSP1", rcp="2.6")

# Check current scenario
scenario = prospective.get_scenario()
print(scenario)  # {'iam': 'IMAGE', 'ssp': 'SSP1', 'rcp': '2.6'}

# Reset scenario (clears current setting)
prospective.reset_scenario()
```

### Choosing a Scenario

- **Low emissions future** (Paris Agreement compatible): IMAGE-SSP1-RCP2.6
- **Middle-of-the-road**: MESSAGE-SSP2-RCP4.5 or GCAM4-SSP4-RCP4.5
- **High emissions**: REMIND-SSP5-RCP8.5

The choice of scenario should align with your study's assumptions about future climate policy and socioeconomic development.

## Using Prospective Metrics

### Basic Usage

```python
import dynamic_characterization.prospective as prospective
from dynamic_characterization import characterize

# 1. Set the scenario first
prospective.set_scenario(iam="IMAGE", ssp="SSP1", rcp="2.6")

# 2. Characterize with prospective radiative forcing (W/m² time series)
df_prf = characterize(
    dynamic_inventory_df,
    metric="prospective_radiative_forcing",
    base_lcia_method=method,
    time_horizon=100,
)

# Or use prospective GWP (kg CO2eq)
df_pgwp = characterize(
    dynamic_inventory_df,
    metric="pGWP",
    base_lcia_method=method,
    time_horizon=100,
)

# Or use prospective GTP (kg CO2eq)
df_pgtp = characterize(
    dynamic_inventory_df,
    metric="pGTP",
    base_lcia_method=method,
    time_horizon=100,
)
```

### With Custom Characterization Functions

```python
from dynamic_characterization.prospective import (
    characterize_co2,
    characterize_ch4,
    characterize_n2o,
)

prospective.set_scenario(iam="IMAGE", ssp="SSP1", rcp="2.6")

df_characterized = characterize(
    dynamic_inventory_df,
    metric="pGWP",
    characterization_functions={
        co2_flow_id: characterize_co2,
        ch4_flow_id: characterize_ch4,
        n2o_flow_id: characterize_n2o,
    },
    time_horizon=100,
)
```

## Fallback to IPCC for Unsupported GHGs

The Watanabe module provides prospective characterization factors only for **CO2, CH4, and N2O**. For other greenhouse gases (like CO and the 244 other GHGs from IPCC AR6), you can choose to either:

1. **Fall back to IPCC** (default): Use standard IPCC AR6 characterization for unsupported GHGs
2. **Skip unsupported GHGs**: Only characterize CO2, CH4, and N2O

```python
# Default behavior: use IPCC for GHGs not in Watanabe
df_prf = characterize(
    dynamic_inventory_df,
    metric="prospective_radiative_forcing",
    base_lcia_method=method,
    fallback_to_ipcc=True,  # default
)

# Strict mode: only characterize CO2, CH4, N2O
df_prf_strict = characterize(
    dynamic_inventory_df,
    metric="prospective_radiative_forcing",
    base_lcia_method=method,
    fallback_to_ipcc=False,
)
```

| `fallback_to_ipcc` | CO2, CH4, N2O | Other GHGs (CO, SF6, etc.) |
|--------------------|---------------|----------------------------|
| `True` (default)   | Watanabe prospective | IPCC AR6 (non-prospective) |
| `False`            | Watanabe prospective | Not characterized |

This option applies to all prospective metrics: `prospective_radiative_forcing`, `pGWP`, and `pGTP`.

## Time-Varying Radiative Efficiency

By default, the radiative efficiency (RE) is fixed at the emission year value (consistent with IPCC methodology). However, you can enable **time-varying RE** for more physical accuracy:

```python
# Works with all prospective metrics
df_prf = characterize(
    dynamic_inventory_df,
    metric="prospective_radiative_forcing",
    base_lcia_method=method,
    time_varying_re=True,  # RE evolves as the gas decays
)

df_pgwp = characterize(
    dynamic_inventory_df,
    metric="pGWP",
    base_lcia_method=method,
    time_varying_re=True,
)
```

### Comparison

| Option | Behavior | Use Case |
|--------|----------|----------|
| `time_varying_re=False` (default) | RE fixed at emission year | IPCC-consistent, comparable results |
| `time_varying_re=True` | RE evolves each year | More physically accurate for long time horizons |

The difference is typically small (<5%) for most scenarios but can be larger for high-emissions scenarios where atmospheric concentrations change rapidly.

## Emission Year Handling

The prospective characterization functions use the emission date to look up the appropriate radiative efficiency:

The bounds come from the RE data itself, which currently spans 2020-2150:

- **Emission years 2020-2150**: Uses exact year data from Watanabe SI tables
- **Emission years < 2020**: Clamped to 2020 with a warning
- **Emission years > 2150**: Clamped to 2150 with a warning

```python
# Emissions at year 2050 - uses 2050 RE data
series_2050 = MockSeries(date="2050-06-15", amount=1.0, flow="CO2")

# Emissions at year 2010 - clamped to 2020 (with warning)
series_2010 = MockSeries(date="2010-01-01", amount=1.0, flow="CO2")
```

Note that the Watanabe paper only reports characterization factors for 2030-2100.
Emission years outside that window use the same underlying RE data but are not
covered by the published tables. With `time_varying_re=True`, an emission year
close to 2150 also has its RE held constant once the decay period runs past the
end of the data.

## Direct AGWP/AGTP Access

For advanced use, you can access the underlying AGWP and AGTP calculation functions directly:

```python
from dynamic_characterization.prospective import agwp, agtp

prospective.set_scenario(iam="IMAGE", ssp="SSP1", rcp="2.6")

# Calculate AGWP for 1 kg CO2 emitted in 2030
agwp_co2 = agwp.agwp_co2(emission_year=2030, time_horizon=100)

# Calculate AGWP for 1 kg CH4
agwp_ch4 = agwp.agwp_ch4(emission_year=2030, time_horizon=100)

# Calculate pGWP
pgwp_ch4 = agwp_ch4 / agwp_co2
print(f"pGWP100 CH4: {pgwp_ch4:.1f}")  # ~27-32 for 2030, depending on scenario

# Same for GTP
agtp_co2 = agtp.agtp_co2(emission_year=2030, time_horizon=100)
agtp_ch4 = agtp.agtp_ch4(emission_year=2030, time_horizon=100)
pgtp_ch4 = agtp_ch4 / agtp_co2
print(f"pGTP100 CH4: {pgtp_ch4:.1f}")  # ~5-8 for 2030, depending on scenario
```

## Method

The prospective metrics cover all gases of the method, evaluated with the scenario for the year of
each emission. Two methods are available (option `prospective_method`):

### Method `"watanabe2026"` (default)

CO2, CH4 and N2O are computed as in Barbosa Watanabe & Cherubini (2026), following the paper's
code (SI `es5c12391_si_028.pdf`, Code #3–#10, ported in `prospective.watanabe2026`), with two
errors of that code corrected:

- CH4: scenario RE plus the AR6 ozone and stratospheric water vapour terms (1.4e-4 and 0.4e-4
  W m⁻² ppb⁻¹); N2O: scenario RE plus ozone (5.5e-4) and the CH4 effect (−1.7 × CH4);
- a climate–carbon feedback (Gasser et al. 2017) driven by an annual two-layer temperature model,
  which also gives the AGTP;
- annual trapezoidal integration; GWP_H and GTP_H are the gas value divided by the CO2 value at
  year H. The paper does not distinguish fossil and non-fossil CH4.

Corrected errors of the published code:

- in the CH4 scripts and the N2O AGWP script, the second and third feedback boxes decayed from
  the previous value of the first box; here every box decays from its own value (as in the N2O
  AGTP script);
- the CO2 impulse response was chosen by the position of the scenario among the IAM's scenarios,
  not by its RCP, which picked the wrong response for AIM and IMAGE (e.g. RCP6.0 for
  IMAGE-SSP1-8.5); here the scenario's own RCP is used.

With both errors restored, the port reproduces every value of the paper's SI tables within their
rounding (this was checked when it was written). With the corrections, the values differ from the
tables where the errors act: CH4 pGWP100 by −2 to +7 % and CH4 pGTP100 by +16 to +34 % in the
scenarios checked (2030); N2O pGTP of GCAM4, MESSAGE and REMIND, which neither error affects, still
equals the tables.

All other gases (halocarbons, CO, VOC, H2, optionally NOx) keep their AR6 or literature forcing;
their CO2 terms follow the same scenario CO2 response, and their AGTP uses the same temperature
scheme, so that the CO2 reference is consistent.

### Method `"ar6"`

`characterization_options={"prospective_method": "ar6"}` applies the AR6 recipe to the scenario
inputs: CH4 and N2O with the AR6 tropospheric adjustments (−14 %, +7 %), the AR6 climate–carbon
feedback and temperature response, fossil/non-fossil CH4 with oxidation CO2, and the CO2 response
of the scenario's own RCP.

### Examples

pGWP100 per kg emitted in 2030 / 2070:

| scenario | CH4, fossil | N2O | CF4 | H2 |
|---|---|---|---|---|
| **`watanabe2026`** | | | | |
| IMAGE-SSP1-2.6 | 32.3 / 38.8 | 266 / 267 | 7319 / 7917 | 11.3 / 12.2 |
| MESSAGE-SSP2-4.5 | 30.4 / 38.8 | 259 / 296 | 7094 / 8746 | 10.9 / 13.6 |
| REMIND-SSP5-8.5 | 27.7 / 37.1 | 242 / 362 | 6619 / 10208 | 10.2 / 16.0 |
| **`ar6`** | | | | |
| IMAGE-SSP1-2.6 | 29.0 / 34.5 | 285 / 291 | 7323 / 7921 | 11.3 / 12.2 |
| MESSAGE-SSP2-4.5 | 27.5 / 35.0 | 278 / 322 | 7098 / 8751 | 10.9 / 13.6 |
| REMIND-SSP5-8.5 | 25.4 / 34.3 | 260 / 392 | 6622 / 10213 | 10.2 / 16.0 |
| AR6 (no scenario) | 29.8 | 273 | 7380 | 11.4 |

### Limitations

- Time horizons up to 100 years (the length of the scenario CO2 impulse response); longer ones are
  truncated to 100 years, with a warning.
- `time_varying_re=True` (see above) is supported by `"watanabe2026"`; the `"ar6"` method uses the
  radiative efficiencies of the emission year and warns.
- REMIND has CH4 radiative efficiencies up to 2100 only; later years use the 2100 value.
- Emission years outside 2020–2150 are clamped, with a warning.

## Reference

Barbosa Watanabe, M. D., & Cherubini, F. (2026). *Environmental Science & Technology*, 60, 3202-3215. https://doi.org/10.1021/acs.est.5c12391 (method and data; the metric code of its Supporting Information is ported in `prospective.watanabe2026`, with two errors of that code corrected).

### Data Sources

The implementation uses data from the paper's Supporting Information:
- **IRF (Impulse Response Functions)**: CO2 (RCP-specific), CH4 (τ=11.8 years), N2O (τ=109 years)
- **RE (Radiative Efficiencies)**: Time series from 2020-2150 for each IAM-SSP-RCP combination
