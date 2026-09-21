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
print(f"pGWP100 CH4: {pgwp_ch4:.1f}")  # ~21-28 depending on scenario

# Same for GTP
agtp_co2 = agtp.agtp_co2(emission_year=2030, time_horizon=100)
agtp_ch4 = agtp.agtp_ch4(emission_year=2030, time_horizon=100)
pgtp_ch4 = agtp_ch4 / agtp_co2
print(f"pGTP100 CH4: {pgtp_ch4:.1f}")  # ~3-6 depending on scenario
```

## Reference

Watanabe, K., Nansai, K., Nakajima, K., & Gibon, T. (2015). **Prospective Characterization Factors for Assessing Climate Change Impacts in Life Cycle Assessments**. *Environmental Science & Technology*, 49(5), 2811-2819. https://doi.org/10.1021/acs.est.5b01118

### Data Sources

The implementation uses data from the paper's Supporting Information:
- **IRF (Impulse Response Functions)**: CO2 (RCP-specific), CH4 (τ=11.8 years), N2O (τ=109 years)
- **RE (Radiative Efficiencies)**: Time series from 2020-2150 for each IAM-SSP-RCP combination
