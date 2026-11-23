"""Options of the default characterization functions (``characterization_options``)."""

from dataclasses import dataclass, field, fields
from typing import Mapping, Optional, Tuple

PARAMETER_SETS = ("legacy", "ar6")
_RESERVED_PARAMETER_SETS = ("fair_irf",)
ORIGINS = ("fossil", "non-fossil")
BIOGENIC_OXIDATION_MODES = ("ar6", "neutral", "explicit")
CCF_PARAMETER_SETS = ("gasser2017_mean", "gasser2017_upper", "gasser2017_lower")
_OVERRIDE_KEYS = ("species", "origin")


@dataclass(frozen=True)
class CharacterizationOptions:
    """
    Options for the default characterization functions.

    Pass them as a dict, e.g. ``characterize(..., characterization_options={"parameter_set":
    "ar6"})``; unknown keys raise a ``ValueError``.

    Attributes
    ----------
    parameter_set : {"legacy", "ar6"}
        ``"legacy"`` (default) keeps the functions of dynamic_characterization <= 1.4.
        ``"ar6"`` reproduces the IPCC AR6 metrics (WG1 Ch. 7, Tables 7.15 and 7.SM.7):
        climate-carbon feedback for non-CO2 gases, separate fossil and non-fossil CH4
        (with oxidation CO2), N2O with its ozone and CH4 effects, the CFC-11/12 adjustment, and
        annual forcing that adds up to AGWP(H).
    climate_carbon_feedback : bool
        ``"ar6"`` only. Add the climate-carbon feedback to non-CO2 gases. Default True (AR6).
    ccf_parameters : {"gasser2017_mean", "gasser2017_upper", "gasser2017_lower"}
        ``"ar6"`` only. Gasser et al. (2017) parameter set; mean is AR6, upper/lower are +-1 SD.
    biogenic_oxidation_co2 : {"ar6", "neutral", "explicit"}
        ``"ar6"`` only. Treatment of the carbon in non-fossil CH4 (see
        ``ipcc_ar6.gases.ch4``). Default ``"ar6"`` (GWP100 27.0).
    ch4_oxidation_yield : float, optional
        ``"ar6"`` only. Fraction of CH4 oxidised to CO2. Default None = AR6 value 0.75.
    unspecified_origin : {"fossil", "non-fossil"}
        ``"ar6"`` only. Origin assumed for CH4 flows without a fossil/non-fossil suffix.
    flow_overrides : mapping
        ``"ar6"`` only. ``{flow id or (database, code): {"species": ..., "origin": ...}}`` to
        override the classification of individual flows. ``species`` is one of "CO2",
        "CO2 uptake", "CH4", "N2O", "CO", a CAS number, or None to leave the flow
        uncharacterized.
    """

    parameter_set: str = "legacy"
    climate_carbon_feedback: bool = True
    ccf_parameters: str = "gasser2017_mean"
    biogenic_oxidation_co2: str = "ar6"
    ch4_oxidation_yield: Optional[float] = None
    unspecified_origin: str = "fossil"
    flow_overrides: Tuple[Tuple[object, Tuple[Tuple[str, object], ...]], ...] = field(default=())

    def __post_init__(self):
        if self.parameter_set in _RESERVED_PARAMETER_SETS:
            raise NotImplementedError(
                f"parameter_set={self.parameter_set!r} is planned but not available yet."
            )
        _check("parameter_set", self.parameter_set, PARAMETER_SETS)
        _check("ccf_parameters", self.ccf_parameters, CCF_PARAMETER_SETS)
        _check("biogenic_oxidation_co2", self.biogenic_oxidation_co2, BIOGENIC_OXIDATION_MODES)
        _check("unspecified_origin", self.unspecified_origin, ORIGINS)
        if not isinstance(self.climate_carbon_feedback, bool):
            raise ValueError("climate_carbon_feedback must be True or False")
        if self.ch4_oxidation_yield is not None and not 0 <= self.ch4_oxidation_yield <= 1:
            raise ValueError("ch4_oxidation_yield must be between 0 and 1")
        for _, override in self.flow_overrides:
            for key, value in override:
                if key not in _OVERRIDE_KEYS:
                    raise ValueError(f"Unknown flow override key {key!r}; use {_OVERRIDE_KEYS}")
                if key == "origin" and value is not None:
                    _check("flow_overrides origin", value, ORIGINS)

    @classmethod
    def from_value(cls, value) -> "CharacterizationOptions":
        """Build from None (defaults), a dict, or an existing instance."""
        if value is None:
            return cls()
        if isinstance(value, cls):
            return value
        if not isinstance(value, Mapping):
            raise TypeError("characterization_options must be a dict or CharacterizationOptions")
        known = {f.name for f in fields(cls)}
        unknown = set(value) - known
        if unknown:
            raise ValueError(
                f"Unknown characterization option(s) {sorted(unknown)}; available: {sorted(known)}"
            )
        kwargs = dict(value)
        if "flow_overrides" in kwargs:
            kwargs["flow_overrides"] = tuple(
                sorted(
                    ((_hashable(flow), tuple(sorted(dict(override).items()))) for flow, override in dict(kwargs["flow_overrides"]).items()),
                    key=repr,
                )
            )
        return cls(**kwargs)

    def override_for(self, node) -> Optional[dict]:
        """The override for a biosphere node (matched by id or (database, code)), if any."""
        keys = {getattr(node, "id", None), (node.get("database"), node.get("code"))}
        for flow, override in self.flow_overrides:
            if flow in keys:
                return dict(override)
        return None


def _hashable(value):
    return tuple(value) if isinstance(value, list) else value


def _check(name, value, allowed):
    if value not in allowed:
        raise ValueError(f"{name} must be one of {allowed}, not {value!r}")
