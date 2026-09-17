"""Pytest options for the optional large-basis reproduction tests."""
import pytest


def pytest_addoption(parser):
    parser.addoption(
        "--run-slow",
        action="store_true",
        default=False,
        help="Recompute every published spectrum, including the 771,425-state case.",
    )


def pytest_configure(config):
    config.addinivalue_line(
        "markers",
        "slow: full spectrum reproduction; may use several GB of memory",
    )


def pytest_collection_modifyitems(config, items):
    if config.getoption("--run-slow"):
        return
    skip_slow = pytest.mark.skip(reason="Full reproduction requires --run-slow.")
    for item in items:
        if "slow" in item.keywords:
            item.add_marker(skip_slow)
