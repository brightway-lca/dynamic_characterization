---
icon: lucide/calculator
tags:
  - api
---

# characterize

The entry point of the package. `characterize` takes a dynamic inventory DataFrame (`date`, `amount`, `flow`, `activity`), works out which characterization function applies to each flow — either from the `characterization_functions` you pass or automatically from a `base_lcia_method` — and returns the characterized inventory for the requested metric.

::: dynamic_characterization.dynamic_characterization
