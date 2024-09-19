import numpy as np
import pandas as pd
from collections import namedtuple
from dynamic_characterization import characterize


def define_dataframes() -> tuple[pd.DataFrame, pd.DataFrame]:
    """
    df_input:

    | date       | amount | flow | activity |
    |------------|--------|------|----------|
    | 15-12-2020 | 10     | 1    | 2        |
    | 20-12-2020 | 20     | 1    | 2        |
    | 25-05-2022 | 50     | 3    | 4        |

    df_expected_characterize:

    | date       | amount | flow | activity |
    |------------|--------|------|----------|
    | 15-12-2020 | 10     | 1    | 2        |
    | 16-12-2020 | 9      | 1    | 2        |
    | 20-12-2020 | 20     | 1    | 2        |
    | 21-12-2020 | 19     | 1    | 2        |
    | 25-05-2022 | 50     | 3    | 4        |
    | 26-05-2022 | 49     | 3    | 4        |

    """

    df_input = pd.DataFrame(
        data={
            "date": pd.Series(
                data=[
                    "15-12-2020",
                    "20-12-2020",
                    "25-05-2022",
                ],
                dtype="datetime64[s]",
            ),
            "amount": pd.Series(
                data=[10.0, 20.0, 50.0], dtype="float64"
            ),
            "flow": pd.Series(data=[1, 1, 3], dtype="int"),
            "activity": pd.Series(data=[2, 2, 4], dtype="int"),
        }
    )

    df_expected_characterize = pd.DataFrame(
        data={
            "date": pd.Series(
                data=[
                    "15-12-2020",
                    "16-12-2020",
                    "20-12-2020",
                    "21-12-2020",
                    "25-05-2022",
                    "26-05-2022",
                ],
                dtype="datetime64[s]",
            ),
            "amount": pd.Series(
                data=[10.0, 9.0, 20.0, 19.0, 50.0, 49.0], dtype="float64"
            ),
            "flow": pd.Series(data=[1, 1, 1, 1, 3, 3], dtype="int"),
            "activity": pd.Series(data=[2, 2, 2, 2, 4, 4], dtype="int"),
        }
    )

    return (df_input, df_expected_characterize)


def function_characterization_test(series: namedtuple, period: int = 2) -> namedtuple:
    date_beginning: np.datetime64 = series.date.to_numpy()
    dates_characterized: np.ndarray = date_beginning + np.arange(
        start=0, stop=period, dtype="timedelta64[D]"
    ).astype("timedelta64[s]")

    amount_beginning: float = series.amount
    amount_characterized: np.ndarray = amount_beginning - np.arange(
        start=0, stop=period, dtype="int"
    )

    return namedtuple("CharacterizedRow", ["date", "amount", "flow", "activity"])(
        date=np.array(dates_characterized, dtype="datetime64[s]"),
        amount=amount_characterized,
        flow=series.flow,
        activity=series.activity,
    )


def test_characterize_dynamic_inventory():
    df_input, df_expected_characterize = define_dataframes()
    df_characterized = characterize(
        df_input,
        metric="radiative_forcing",
        characterization_function_dict={
            1: function_characterization_test,
            3: function_characterization_test,
        },
        time_horizon=2,
    )

    pd.testing.assert_frame_equal(df_characterized, df_expected_characterize)
    