---
icon: lucide/thermometer
tags:
  - api
  - fair
---

# FAIR

The optional FAIR climate-model path: instead of convolving a fixed impulse response, the inventory is run through the [FaIR](https://github.com/OMS-/FAIR) simple climate model on top of a scenario background. Requires the `fair` extra:

```bash
pip install dynamic_characterization[fair]
```

## Core

::: dynamic_characterization.fair.core

## Runner

::: dynamic_characterization.fair.runner

## Species map

::: dynamic_characterization.fair.species_map

## Allocation

::: dynamic_characterization.fair.allocation
