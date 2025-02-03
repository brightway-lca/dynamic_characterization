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
