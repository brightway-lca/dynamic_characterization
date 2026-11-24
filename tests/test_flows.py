"""CAS normalisation and the CAS-based lookup of the default characterization functions."""

import csv
import os

import bw2data as bd
import pytest
from bw2data.tests import bw2test

from dynamic_characterization.dynamic_characterization import (
    clear_characterization_function_cache,
    create_characterization_functions_from_method,
)
from dynamic_characterization.flows import (
    KNOWN_CAS_TYPOS,
    build_cas_index,
    is_valid_cas,
    normalize_cas,
)

HALOCARBON_TABLE = os.path.join(
    os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
    "dynamic_characterization",
    "ipcc_ar6",
    "data",
    "hodnebrog20.csv",
)


def _table_cas_numbers():
    """CAS numbers of the gases with a lifetime in the halocarbon table, as written there."""
    with open(HALOCARBON_TABLE, newline="", encoding="utf-8-sig") as f:
        return [row["CASRN"] for row in csv.DictReader(f) if row["CASRN"] and row["Lifetime (yr)"]]


@pytest.mark.parametrize(
    "raw, expected",
    [
        ("000124-38-9", "124-38-9"),  # ecoinvent zero padding
        ("124-38-9", "124-38-9"),
        ("75-38-7", "75-38-7"),  # vinylidene fluoride, unpadded in ecoinvent 3.12
        ("000375-17-7 ", "375-17-7"),  # trailing space
        ("0 5595-10-8", "5595-10-8"),  # inner space
        ("1309353-34-1", "1309353-34-1"),  # seven-digit first group
        ("0001309353-34-1", "1309353-34-1"),
    ],
)
def test_normalize_cas(raw, expected):
    assert normalize_cas(raw) == expected


@pytest.mark.parametrize("raw", [None, "", "abc", "124-38", "1-38-9", "12345678-38-9", 12438])
def test_normalize_cas_rejects_non_cas(raw):
    assert normalize_cas(raw) is None


def test_is_valid_cas():
    assert is_valid_cas("000124-38-9")
    assert not is_valid_cas("124-38-8")
    assert not is_valid_cas("29188-24-9")  # typo in the bundled data
    assert is_valid_cas("29118-24-9")


def test_every_table_cas_number_is_indexed():
    keys = _table_cas_numbers()
    index = build_cas_index(keys)
    assert set(index.values()) == set(keys)
    assert len(index) == len(keys) + len(KNOWN_CAS_TYPOS)


def test_table_cas_numbers_have_valid_check_digits_except_known_typos():
    invalid = {normalize_cas(k) for k in _table_cas_numbers() if not is_valid_cas(k)}
    assert invalid == set(KNOWN_CAS_TYPOS.values())


def test_known_typo_resolves_under_correct_cas():
    index = build_cas_index(["029188-24-9"])
    assert index["29118-24-9"] == index["29188-24-9"] == "029188-24-9"


METHOD = ("test", "climate change", "CAS lookup")
BIOFLOWS = {
    # (name, CAS as stored on the flow)
    "vinylidene_fluoride": ("Vinylidene fluoride", "75-38-7"),
    "cf4_padded": ("Tetrafluoromethane", "000075-73-0"),
    "hfc329p_unpadded": ("HFC-329p", "375-17-7"),
    "hfo1225ye_padded": ("HFO-1225ye(E)", "005595-10-8"),
    "hfo1234ze_correct_cas": ("HFO-1234ze(E)", "29118-24-9"),
    "no_cas": ("Some flow", None),
    "unknown_cas": ("Water", "007732-18-5"),
}


def setup_project():
    clear_characterization_function_cache()
    biosphere = bd.Database(bd.config.biosphere)
    biosphere.write(
        {
            (bd.config.biosphere, code): {
                "code": code,
                "name": name,
                "categories": ("air",),
                "type": "emission",
            }
            | ({"CAS number": cas} if cas else {})
            for code, (name, cas) in BIOFLOWS.items()
        }
    )
    bd.Method(METHOD).write([((bd.config.biosphere, code), 1.0) for code in BIOFLOWS])
    return {code: biosphere.get(code).id for code in BIOFLOWS}


@bw2test
def test_lookup_ignores_cas_formatting():
    ids = setup_project()
    functions = create_characterization_functions_from_method(METHOD)
    for code in (
        "vinylidene_fluoride",
        "cf4_padded",
        "hfc329p_unpadded",
        "hfo1225ye_padded",
        "hfo1234ze_correct_cas",
    ):
        assert ids[code] in functions, code
    assert ids["no_cas"] not in functions
    assert ids["unknown_cas"] not in functions
