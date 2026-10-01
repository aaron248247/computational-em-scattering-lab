import pytest

from em_rf_lab.simulation import meep_available


def pytest_collection_modifyitems(config, items):
    if meep_available():
        return
    skip = pytest.mark.skip(reason="MEEP (pymeep) not installed")
    for item in items:
        if "meep" in item.keywords:
            item.add_marker(skip)
