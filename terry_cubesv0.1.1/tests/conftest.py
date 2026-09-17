"""Pytest options for the optional large-basis reproduction tests."""
import pytest


def pytest_addoption(parser):
    parser.addoption("--run-research", action="store_true", default=False,
                     help="Run independent tensors, extended references and exact certificates.")
    parser.addoption(
        "--run-slow",
        action="store_true",
        default=False,
        help="Recompute every published spectrum, including the 771,425-state case.",
    )


def pytest_configure(config):
    config.addinivalue_line(
        "markers",
        "research: independent tensors, extended references and exact certificates",
    )
    config.addinivalue_line(
        "markers",
        "slow: full spectrum reproduction; may use several GB of memory",
    )


def pytest_collection_modifyitems(config, items):
    skip_slow = pytest.mark.skip(reason="Full reproduction requires --run-slow.")
    skip_research = pytest.mark.skip(reason="Independent and certificate checks require --run-research.")
    for item in items:
        if "slow" in item.keywords and not config.getoption("--run-slow"):
            item.add_marker(skip_slow)
        if "research" in item.keywords and not config.getoption("--run-research"):
            item.add_marker(skip_research)
