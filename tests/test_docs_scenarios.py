"""Guard that README documents the FAIR optional install extra and both metrics."""

import os


def test_readme_documents_fair():
    root = os.path.join(os.path.dirname(__file__), "..")
    with open(os.path.join(root, "README.md")) as fh:
        text = fh.read()
    assert "dynamic_characterization[fair]" in text
    assert "fair_radiative_forcing" in text
    assert "fair_temperature" in text
    assert "available_scenarios" in text
