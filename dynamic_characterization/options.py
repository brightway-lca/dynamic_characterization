"""Options of the default characterization functions (``characterization_options``)."""

from dataclasses import dataclass, field, fields
from typing import Mapping, Optional, Tuple

from dynamic_characterization.ipcc_ar6.gases import precursor_parameter_set

PARAMETER_SETS = ("legacy", "ar6")
_RESERVED_PARAMETER_SETS = ("fair_irf",)
ORIGINS = ("fossil", "non-fossil")
BIOGENIC_OXIDATION_MODES = ("ar6", "neutral", "explicit")
CCF_PARAMETER_SETS = ("gasser2017_mean", "gasser2017_upper", "gasser2017_lower")
_OVERRIDE_KEYS = ("species", "origin")
PRECURSORS = ("CO", "H2", "VOC")
PROSPECTIVE_METHODS = ("watanabe2026", "ar6")
_PARAMETER_SET_SPECIES = ("CO", "VOC", "H2", "NOx", "nitrate")


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
    include_precursors : collection of {"CO", "VOC", "H2"}
        ``"ar6"`` only. Precursors characterized with their indirect forcing (ozone, CH4,
        sulfate, stratospheric water vapour) from published parameter sets. Default: all three.
        CO not listed falls back to its oxidation CO2 only; VOC and H2 not listed are not
        characterized.
    precursor_parameter_sets : mapping
        ``"ar6"`` only. ``{species: set name}`` for "CO", "VOC", "H2", "NOx", "nitrate"; see
        ``ipcc_ar6/data/precursor_parameters.json``. Defaults: CO and VOC "fry2012_4regions",
        H2 "warwick2023", NOx "fry2012_4regions", nitrate "bauer2007".
    include_nox : bool
        ``"ar6"`` only. Characterize NO, NO2 and "Nitrogen oxides" (default False: sign and
        size of the NOx metric are very uncertain).
    nox_nitrate : bool
        ``"ar6"`` only, with ``include_nox``. Add the direct forcing of nitrate aerosol.
    precursor_ccf : bool
        ``"ar6"`` only. Add the climate-carbon feedback to the precursors (the published sets
        do not include it; AR6 adds it to all its metrics). Default False.
    prospective_method : {"watanabe2026", "ar6"}
        ``"ar6"`` with the prospective metrics only. ``"watanabe2026"`` (default) computes CO2,
        CH4 and N2O as Barbosa Watanabe & Cherubini (2026), with two errors of the published code
        corrected (see ``prospective.watanabe2026``);
        ``"ar6"`` applies the AR6 recipe (tropospheric adjustments, AR6 climate-carbon feedback and
        temperature response) to the scenario inputs. See ``ar6_prospective``.
    """

    parameter_set: str = "legacy"
    climate_carbon_feedback: bool = True
    ccf_parameters: str = "gasser2017_mean"
    biogenic_oxidation_co2: str = "ar6"
    ch4_oxidation_yield: Optional[float] = None
    unspecified_origin: str = "fossil"
    flow_overrides: Tuple[Tuple[object, Tuple[Tuple[str, object], ...]], ...] = field(default=())
    include_precursors: Tuple[str, ...] = PRECURSORS
    precursor_parameter_sets: Tuple[Tuple[str, str], ...] = ()
    include_nox: bool = False
    nox_nitrate: bool = False
    precursor_ccf: bool = False
    prospective_method: str = "watanabe2026"

    def __post_init__(self):
        if self.parameter_set in _RESERVED_PARAMETER_SETS:
            raise NotImplementedError(
                f"parameter_set={self.parameter_set!r} is planned but not available yet."
            )
        _check("parameter_set", self.parameter_set, PARAMETER_SETS)
        _check("ccf_parameters", self.ccf_parameters, CCF_PARAMETER_SETS)
        _check("biogenic_oxidation_co2", self.biogenic_oxidation_co2, BIOGENIC_OXIDATION_MODES)
        _check("unspecified_origin", self.unspecified_origin, ORIGINS)
        _check("prospective_method", self.prospective_method, PROSPECTIVE_METHODS)
        for name in ("climate_carbon_feedback", "include_nox", "nox_nitrate", "precursor_ccf"):
            if not isinstance(getattr(self, name), bool):
                raise ValueError(f"{name} must be True or False")
        unknown = set(self.include_precursors) - set(PRECURSORS)
        if unknown:
            raise ValueError(f"include_precursors must be a subset of {PRECURSORS}, not {sorted(unknown)}")
        if self.nox_nitrate and not self.include_nox:
            raise ValueError("nox_nitrate requires include_nox=True")
        for species, name in self.precursor_parameter_sets:
            if species not in _PARAMETER_SET_SPECIES:
                raise ValueError(f"precursor_parameter_sets: unknown species {species!r}; use {_PARAMETER_SET_SPECIES}")
            precursor_parameter_set(species, name)  # raises for unknown set names
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
        if "include_precursors" in kwargs:
            kwargs["include_precursors"] = tuple(sorted(set(kwargs["include_precursors"])))
        if "precursor_parameter_sets" in kwargs:
            kwargs["precursor_parameter_sets"] = tuple(sorted(dict(kwargs["precursor_parameter_sets"]).items()))
        if "flow_overrides" in kwargs:
            kwargs["flow_overrides"] = tuple(
                sorted(
                    ((_hashable(flow), tuple(sorted(dict(override).items()))) for flow, override in dict(kwargs["flow_overrides"]).items()),
                    key=repr,
                )
            )
        return cls(**kwargs)

    def precursor_parameter_set(self, species: str) -> Optional[str]:
        """Chosen parameter set of a precursor species, or None for the default."""
        return dict(self.precursor_parameter_sets).get(species)

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
