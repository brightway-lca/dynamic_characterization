"""Helpers for matching biosphere flows to characterization data."""

import re
from typing import Dict, Iterable, Optional

from loguru import logger

_CAS_PATTERN = re.compile(r"^(\d+)-(\d{2})-(\d)$")

# Known CAS typos in the bundled characterization data, as {correct CAS: CAS in the data}.
# Both forms resolve to the same data entry.
KNOWN_CAS_TYPOS = {
    # HFO-1234ze(E): the Hodnebrog et al. (2020) table lists 29188-24-9, which fails the
    # CAS check digit; the registry number is 29118-24-9.
    "29118-24-9": "29188-24-9",
}


def _cas_check_digit(digits: str) -> int:
    """Check digit of a CAS number, given all digits before it (CAS registry rule)."""
    return sum(i * int(d) for i, d in enumerate(reversed(digits), start=1)) % 10


def normalize_cas(cas: Optional[str]) -> Optional[str]:
    """
    Return a CAS registry number in canonical form, or None if it is not CAS-formatted.

    The canonical form has no whitespace and no leading zeros in the first group, so that the
    zero-padded form used by ecoinvent (``000075-38-7``) and the plain form (``75-38-7``) match.
    Whitespace anywhere in the string is removed (e.g. ``"0 5595-10-8"``, ``"000375-17-7 "``).
    The check digit is not enforced here, so that a number is matched as written; use
    :func:`is_valid_cas` to validate it.

    Examples
    --------
    >>> normalize_cas("000124-38-9")
    '124-38-9'
    >>> normalize_cas(" 0 5595-10-8")
    '5595-10-8'
    >>> normalize_cas("not a CAS") is None
    True
    """
    if not isinstance(cas, str):
        return None
    match = _CAS_PATTERN.match(re.sub(r"\s", "", cas))
    if match is None:
        return None
    first, second, check = match.groups()
    first = first.lstrip("0")
    if not 2 <= len(first) <= 7:
        return None
    return f"{first}-{second}-{check}"


def is_valid_cas(cas: Optional[str]) -> bool:
    """True if ``cas`` is CAS-formatted (see :func:`normalize_cas`) and its check digit is correct."""
    canonical = normalize_cas(cas)
    if canonical is None:
        return False
    first, second, check = canonical.split("-")
    return _cas_check_digit(first + second) == int(check)


def build_cas_index(keys: Iterable[str]) -> Dict[str, str]:
    """
    Map the canonical CAS number of each key to the key itself.

    Keys that are not CAS-formatted are skipped with a warning. If two keys normalise to the same
    CAS number, the first one wins and a warning is logged. Keys listed in
    :data:`KNOWN_CAS_TYPOS` are also indexed under their correct CAS number.
    """
    index: Dict[str, str] = {}
    for key in keys:
        canonical = normalize_cas(key)
        if canonical is None:
            logger.warning("Skipping {!r} in characterization data: not a CAS number.", key)
            continue
        if canonical in index:
            logger.warning(
                "CAS numbers {!r} and {!r} are the same after normalisation; using {!r}.",
                index[canonical],
                key,
                index[canonical],
            )
            continue
        index[canonical] = key
    for correct, typo in KNOWN_CAS_TYPOS.items():
        if typo in index and correct not in index:
            index[correct] = index[typo]
    return index
