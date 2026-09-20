---
tags:
  - api
---

# API Reference

This section contains the API documentation generated from the `dynamic_characterization` source code docstrings.

The main user-facing function is [`characterize`](characterize.md). It takes a dynamic inventory and a metric, and applies the matching characterization functions:

- [`characterize`](characterize.md) — the entry point: applies characterization functions to a dynamic inventory.
- [`ipcc_ar6`](ipcc_ar6.md) — radiative forcing and GWP based on IPCC AR6 radiative efficiencies and lifetimes.
- [`prospective`](prospective.md) — scenario-based prospective characterization factors (pGWP, pGTP) after Watanabe et al.
- [`original_temporalis_functions`](original_temporalis_functions.md) — the legacy functions inherited from `bw_temporalis`.
- [`classes`](classes.md) — the data structures a characterization function returns.
