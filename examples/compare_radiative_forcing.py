"""
Compare radiative forcing between IPCC AR6 and prospective characterization.

This example demonstrates the difference between:
- Standard IPCC AR6 radiative forcing (fixed radiative efficiencies)
- Prospective radiative forcing (scenario-based, time-varying radiative efficiencies)

The prospective approach accounts for how atmospheric GHG concentrations change
over time under different climate scenarios, which affects radiative efficiency.
"""

import sys
import os
import importlib.util
import types

# Set up path for imports
_base_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
_dc_dir = os.path.join(_base_dir, "dynamic_characterization")

# Create fake package structure to enable imports without triggering __init__.py
sys.modules["dynamic_characterization"] = types.ModuleType("dynamic_characterization")
sys.modules["dynamic_characterization"].__path__ = [_dc_dir]

# Load classes module first (needed by radiative_forcing modules)
_classes_path = os.path.join(_dc_dir, "classes.py")
_classes_spec = importlib.util.spec_from_file_location("dynamic_characterization.classes", _classes_path)
classes = importlib.util.module_from_spec(_classes_spec)
sys.modules["dynamic_characterization.classes"] = classes
_classes_spec.loader.exec_module(classes)

# Load ipcc_ar6 radiative_forcing module
_ipcc_rf_path = os.path.join(_dc_dir, "ipcc_ar6", "radiative_forcing.py")
_ipcc_rf_spec = importlib.util.spec_from_file_location("dynamic_characterization.ipcc_ar6.radiative_forcing", _ipcc_rf_path)
ipcc_rf = importlib.util.module_from_spec(_ipcc_rf_spec)
sys.modules["dynamic_characterization.ipcc_ar6.radiative_forcing"] = ipcc_rf
_ipcc_rf_spec.loader.exec_module(ipcc_rf)

# Set up prospective package structure
_prosp_dir = os.path.join(_dc_dir, "prospective")
_prosp_pkg = types.ModuleType("dynamic_characterization.prospective")
_prosp_pkg.__path__ = [_prosp_dir]
_prosp_pkg.__package__ = "dynamic_characterization.prospective"
sys.modules["dynamic_characterization.prospective"] = _prosp_pkg

# Load prospective config module
_config_path = os.path.join(_prosp_dir, "config.py")
_config_spec = importlib.util.spec_from_file_location("dynamic_characterization.prospective.config", _config_path)
config = importlib.util.module_from_spec(_config_spec)
sys.modules["dynamic_characterization.prospective.config"] = config
_config_spec.loader.exec_module(config)

# Load prospective data_loader module
_dl_path = os.path.join(_prosp_dir, "data_loader.py")
_dl_spec = importlib.util.spec_from_file_location("dynamic_characterization.prospective.data_loader", _dl_path)
data_loader = importlib.util.module_from_spec(_dl_spec)
sys.modules["dynamic_characterization.prospective.data_loader"] = data_loader
_dl_spec.loader.exec_module(data_loader)

# Load prospective radiative_forcing module
_prosp_rf_path = os.path.join(_prosp_dir, "radiative_forcing.py")
_prosp_rf_spec = importlib.util.spec_from_file_location("dynamic_characterization.prospective.radiative_forcing", _prosp_rf_path)
prosp_rf = importlib.util.module_from_spec(_prosp_rf_spec)
sys.modules["dynamic_characterization.prospective.radiative_forcing"] = prosp_rf
_prosp_rf_spec.loader.exec_module(prosp_rf)

import matplotlib.pyplot as plt
import numpy as np

# Now use the loaded modules
ipcc_characterize_co2 = ipcc_rf.characterize_co2
ipcc_characterize_ch4 = ipcc_rf.characterize_ch4
ipcc_characterize_n2o = ipcc_rf.characterize_n2o

set_scenario = config.set_scenario
reset_scenario = config.reset_scenario

prospective_characterize_co2 = prosp_rf.characterize_co2
prospective_characterize_ch4 = prosp_rf.characterize_ch4
prospective_characterize_n2o = prosp_rf.characterize_n2o


class MockSeries:
    """Mock series object mimicking a row from dynamic inventory."""

    def __init__(self, date: str, amount: float, flow: str = "test", activity: str = "test"):
        self._date = date
        self.amount = amount
        self.flow = flow
        self.activity = activity

    @property
    def date(self):
        class DateWrapper:
            def __init__(self, date_str):
                self._date = np.datetime64(date_str, "s")

            def to_numpy(self):
                return self._date

        return DateWrapper(self._date)


def compare_co2_radiative_forcing(
    emission_year: int = 2050,
    time_horizon: int = 100,
    scenario: dict = None,
):
    """
    Compare CO2 radiative forcing between IPCC and prospective approaches.

    Parameters
    ----------
    emission_year : int
        Year of emission (default: 2050)
    time_horizon : int
        Time horizon in years (default: 100)
    scenario : dict
        Scenario settings with keys 'iam', 'ssp', 'rcp'.
        Default: IMAGE-SSP1-RCP2.6
    """
    if scenario is None:
        scenario = {"iam": "IMAGE", "ssp": "SSP1", "rcp": "2.6"}

    # Set prospective scenario
    reset_scenario()
    set_scenario(**scenario)

    # Create mock emission (1 kg CO2)
    series = MockSeries(date=f"{emission_year}-01-01", amount=1.0, flow="CO2")

    # Characterize with IPCC AR6
    ipcc_result = ipcc_characterize_co2(series, period=time_horizon)

    # Characterize with prospective (fixed RE at emission year)
    prospective_result = prospective_characterize_co2(
        series, period=time_horizon, time_varying_re=False
    )

    # Characterize with prospective (time-varying RE)
    prospective_varying_result = prospective_characterize_co2(
        series, period=time_horizon, time_varying_re=True
    )

    return {
        "ipcc": ipcc_result,
        "prospective_fixed": prospective_result,
        "prospective_varying": prospective_varying_result,
        "years": np.arange(emission_year, emission_year + time_horizon),
    }


def compare_ch4_radiative_forcing(
    emission_year: int = 2050,
    time_horizon: int = 100,
    scenario: dict = None,
):
    """Compare CH4 radiative forcing between IPCC and prospective approaches."""
    if scenario is None:
        scenario = {"iam": "IMAGE", "ssp": "SSP1", "rcp": "2.6"}

    reset_scenario()
    set_scenario(**scenario)

    series = MockSeries(date=f"{emission_year}-01-01", amount=1.0, flow="CH4")

    ipcc_result = ipcc_characterize_ch4(series, period=time_horizon)
    prospective_result = prospective_characterize_ch4(
        series, period=time_horizon, time_varying_re=False
    )
    prospective_varying_result = prospective_characterize_ch4(
        series, period=time_horizon, time_varying_re=True
    )

    return {
        "ipcc": ipcc_result,
        "prospective_fixed": prospective_result,
        "prospective_varying": prospective_varying_result,
        "years": np.arange(emission_year, emission_year + time_horizon),
    }


def compare_n2o_radiative_forcing(
    emission_year: int = 2050,
    time_horizon: int = 100,
    scenario: dict = None,
):
    """Compare N2O radiative forcing between IPCC and prospective approaches."""
    if scenario is None:
        scenario = {"iam": "IMAGE", "ssp": "SSP1", "rcp": "2.6"}

    reset_scenario()
    set_scenario(**scenario)

    series = MockSeries(date=f"{emission_year}-01-01", amount=1.0, flow="N2O")

    ipcc_result = ipcc_characterize_n2o(series, period=time_horizon)
    prospective_result = prospective_characterize_n2o(
        series, period=time_horizon, time_varying_re=False
    )
    prospective_varying_result = prospective_characterize_n2o(
        series, period=time_horizon, time_varying_re=True
    )

    return {
        "ipcc": ipcc_result,
        "prospective_fixed": prospective_result,
        "prospective_varying": prospective_varying_result,
        "years": np.arange(emission_year, emission_year + time_horizon),
    }


def plot_comparison(results: dict, gas: str, scenario_name: str, ax=None):
    """
    Plot radiative forcing comparison.

    Parameters
    ----------
    results : dict
        Results from compare_*_radiative_forcing functions
    gas : str
        Gas name for title
    scenario_name : str
        Scenario name for title
    ax : matplotlib.axes.Axes, optional
        Axes to plot on. If None, creates new figure.
    """
    if ax is None:
        fig, ax = plt.subplots(figsize=(10, 6))

    years = results["years"]

    # Plot IPCC (reference)
    ax.plot(
        years,
        results["ipcc"].amount,
        label="IPCC AR6 (fixed RE)",
        color="black",
        linewidth=2,
        linestyle="-",
    )

    # Plot prospective with fixed RE
    ax.plot(
        years,
        results["prospective_fixed"].amount,
        label="Prospective (fixed RE at emission year)",
        color="blue",
        linewidth=2,
        linestyle="--",
    )

    # Plot prospective with time-varying RE
    ax.plot(
        years,
        results["prospective_varying"].amount,
        label="Prospective (time-varying RE)",
        color="red",
        linewidth=2,
        linestyle=":",
    )

    ax.set_xlabel("Year", fontsize=12)
    ax.set_ylabel("Radiative Forcing (W/m² per kg)", fontsize=12)
    ax.set_title(f"{gas} Radiative Forcing Comparison\n{scenario_name}", fontsize=14)
    ax.legend(loc="upper right", fontsize=10)
    ax.grid(True, alpha=0.3)
    ax.set_xlim(years[0], years[-1])

    return ax


def plot_cumulative_comparison(results: dict, gas: str, scenario_name: str, ax=None):
    """Plot cumulative radiative forcing (AGWP) comparison."""
    if ax is None:
        fig, ax = plt.subplots(figsize=(10, 6))

    years = results["years"]

    # Calculate cumulative forcing
    ipcc_cumulative = np.cumsum(results["ipcc"].amount)
    prospective_fixed_cumulative = np.cumsum(results["prospective_fixed"].amount)
    prospective_varying_cumulative = np.cumsum(results["prospective_varying"].amount)

    ax.plot(
        years,
        ipcc_cumulative,
        label="IPCC AR6 (fixed RE)",
        color="black",
        linewidth=2,
        linestyle="-",
    )
    ax.plot(
        years,
        prospective_fixed_cumulative,
        label="Prospective (fixed RE at emission year)",
        color="blue",
        linewidth=2,
        linestyle="--",
    )
    ax.plot(
        years,
        prospective_varying_cumulative,
        label="Prospective (time-varying RE)",
        color="red",
        linewidth=2,
        linestyle=":",
    )

    ax.set_xlabel("Year", fontsize=12)
    ax.set_ylabel("Cumulative RF (W·yr/m² per kg)", fontsize=12)
    ax.set_title(f"{gas} Cumulative Radiative Forcing (AGWP)\n{scenario_name}", fontsize=14)
    ax.legend(loc="lower right", fontsize=10)
    ax.grid(True, alpha=0.3)
    ax.set_xlim(years[0], years[-1])

    return ax


def main():
    """Generate comparison plots for CO2, CH4, and N2O."""

    # Settings
    emission_year = 2050
    time_horizon = 100
    scenario = {"iam": "IMAGE", "ssp": "SSP1", "rcp": "2.6"}
    scenario_name = f"IMAGE-SSP1-RCP2.6, emission year {emission_year}"

    # Create figure with subplots for all three gases
    fig, axes = plt.subplots(3, 2, figsize=(14, 15))
    fig.suptitle(
        "Comparison of IPCC AR6 vs Prospective Radiative Forcing",
        fontsize=16,
        fontweight="bold",
    )

    # CO2
    print("Calculating CO2 radiative forcing...")
    co2_results = compare_co2_radiative_forcing(emission_year, time_horizon, scenario)
    plot_comparison(co2_results, "CO2", scenario_name, ax=axes[0, 0])
    plot_cumulative_comparison(co2_results, "CO2", scenario_name, ax=axes[0, 1])

    # CH4
    print("Calculating CH4 radiative forcing...")
    ch4_results = compare_ch4_radiative_forcing(emission_year, time_horizon, scenario)
    plot_comparison(ch4_results, "CH4", scenario_name, ax=axes[1, 0])
    plot_cumulative_comparison(ch4_results, "CH4", scenario_name, ax=axes[1, 1])

    # N2O
    print("Calculating N2O radiative forcing...")
    n2o_results = compare_n2o_radiative_forcing(emission_year, time_horizon, scenario)
    plot_comparison(n2o_results, "N2O", scenario_name, ax=axes[2, 0])
    plot_cumulative_comparison(n2o_results, "N2O", scenario_name, ax=axes[2, 1])

    plt.tight_layout()
    output_path = os.path.join(os.path.dirname(__file__), "radiative_forcing_comparison.png")
    plt.savefig(output_path, dpi=150, bbox_inches="tight")
    print(f"\nPlot saved to '{output_path}'")
    plt.show()

    # Print summary statistics
    print("\n" + "=" * 60)
    print("SUMMARY: Cumulative RF after 100 years (AGWP)")
    print("=" * 60)

    for gas, results in [("CO2", co2_results), ("CH4", ch4_results), ("N2O", n2o_results)]:
        ipcc_agwp = np.sum(results["ipcc"].amount)
        prosp_fixed_agwp = np.sum(results["prospective_fixed"].amount)
        prosp_varying_agwp = np.sum(results["prospective_varying"].amount)

        print(f"\n{gas}:")
        print(f"  IPCC AR6:              {ipcc_agwp:.4e} W·yr/m²/kg")
        print(f"  Prospective (fixed):   {prosp_fixed_agwp:.4e} W·yr/m²/kg")
        print(f"  Prospective (varying): {prosp_varying_agwp:.4e} W·yr/m²/kg")
        print(f"  Ratio (prosp_fixed/IPCC): {prosp_fixed_agwp/ipcc_agwp:.3f}")


def compare_scenarios():
    """Compare radiative forcing across different scenarios for CO2."""

    emission_year = 2050
    time_horizon = 100

    scenarios = [
        {"iam": "IMAGE", "ssp": "SSP1", "rcp": "2.6"},
        {"iam": "MESSAGE", "ssp": "SSP2", "rcp": "4.5"},
        {"iam": "REMIND", "ssp": "SSP5", "rcp": "8.5"},
    ]
    scenario_names = ["IMAGE-SSP1-RCP2.6", "MESSAGE-SSP2-RCP4.5", "REMIND-SSP5-RCP8.5"]

    fig, axes = plt.subplots(1, 2, figsize=(14, 6))
    fig.suptitle(
        f"CO2 Radiative Forcing Across Scenarios (emission year {emission_year})",
        fontsize=14,
        fontweight="bold",
    )

    colors = ["green", "orange", "red"]

    # First, get IPCC reference (same for all scenarios)
    reset_scenario()
    set_scenario(**scenarios[0])
    series = MockSeries(date=f"{emission_year}-01-01", amount=1.0, flow="CO2")
    ipcc_result = ipcc_characterize_co2(series, period=time_horizon)
    years = np.arange(emission_year, emission_year + time_horizon)

    # Plot IPCC reference
    axes[0].plot(
        years, ipcc_result.amount, label="IPCC AR6", color="black", linewidth=2, linestyle="-"
    )
    axes[1].plot(
        years,
        np.cumsum(ipcc_result.amount),
        label="IPCC AR6",
        color="black",
        linewidth=2,
        linestyle="-",
    )

    # Plot prospective for each scenario
    for scenario, name, color in zip(scenarios, scenario_names, colors):
        reset_scenario()
        set_scenario(**scenario)

        prosp_result = prospective_characterize_co2(
            series, period=time_horizon, time_varying_re=False
        )

        axes[0].plot(
            years, prosp_result.amount, label=f"Prospective {name}", color=color, linewidth=1.5
        )
        axes[1].plot(
            years,
            np.cumsum(prosp_result.amount),
            label=f"Prospective {name}",
            color=color,
            linewidth=1.5,
        )

    axes[0].set_xlabel("Year", fontsize=12)
    axes[0].set_ylabel("Radiative Forcing (W/m² per kg)", fontsize=12)
    axes[0].set_title("Marginal Radiative Forcing", fontsize=12)
    axes[0].legend(fontsize=9)
    axes[0].grid(True, alpha=0.3)

    axes[1].set_xlabel("Year", fontsize=12)
    axes[1].set_ylabel("Cumulative RF (W·yr/m² per kg)", fontsize=12)
    axes[1].set_title("Cumulative Radiative Forcing (AGWP)", fontsize=12)
    axes[1].legend(fontsize=9)
    axes[1].grid(True, alpha=0.3)

    plt.tight_layout()
    output_path = os.path.join(os.path.dirname(__file__), "radiative_forcing_scenarios.png")
    plt.savefig(output_path, dpi=150, bbox_inches="tight")
    print(f"\nScenario comparison saved to '{output_path}'")
    plt.show()


if __name__ == "__main__":
    print("Generating radiative forcing comparison plots...")
    print("-" * 60)

    main()

    print("\n" + "-" * 60)
    print("Generating scenario comparison...")
    print("-" * 60)

    compare_scenarios()
