"""Numerical regressions against the unchanged joined-cube reference spectra.

The single-cube interval check tests agreement with a separately established
certificate. It does not rerun that certificate's exact-arithmetic verifier.
"""
from functools import lru_cache
import json
from pathlib import Path

import numpy as np
import pytest

from terry_cubes import JoinedCubesModel, SpectrumSolver


DATA_DIR = Path(__file__).resolve().parents[1] / "data"
ENERGY_ATOL = 1e-9
CERTIFIED_SINGLE_CUBE_INTERVAL = (2.941733, 2.958693)


@lru_cache(maxsize=3)
def _reference_records(num_cubes):
    path = DATA_DIR / f"SU2_Joined_Cubes_N{num_cubes}_Spectra.json"
    with path.open(encoding="utf-8") as stream:
        document = json.load(stream)
    assert document["model"]["num_cubes"] == num_cubes
    assert document["units"] == "kappa=1"
    return tuple(document["spectra"])


def _reference(num_cubes, cutoff, ratio):
    matches = [
        row
        for row in _reference_records(num_cubes)
        if (
            row["num_cubes"],
            row["energy_cap"],
            row["nu_over_kappa"],
        ) == (num_cubes, cutoff, ratio)
    ]
    assert len(matches) == 1, "A baseline must be identified by its model and cutoff."
    return matches[0]


@lru_cache(maxsize=None)
def _compute(num_cubes, cutoff, ratio):
    model = JoinedCubesModel(num_cubes=num_cubes, kappa=1.0, nu=ratio)
    return SpectrumSolver(model, electric_cutoff=cutoff).compute_mass_gap()


def _assert_matches_reference(result, reference):
    assert result.dimension == reference["dimension"]
    assert result.converged
    np.testing.assert_allclose(
        [result.e0, result.e1, result.gap],
        [reference["E0"], reference["E1"], reference["gap"]],
        atol=ENERGY_ATOL,
        rtol=0,
    )
    assert result.e1 >= result.e0
    assert result.gap == pytest.approx(result.e1 - result.e0, abs=1e-12)
    np.testing.assert_allclose(
        result.eigenvalues[:2], [result.e0, result.e1], atol=1e-12, rtol=0
    )
    residuals = np.asarray(result.residual_norms)
    assert residuals.size >= 2
    assert np.all(np.isfinite(residuals))
    assert np.all(residuals <= ENERGY_ATOL)


FAST_CASES = [
    pytest.param(cubes, cutoff, ratio, id=f"N{cubes}-cap{cutoff}-nu{ratio}")
    for cubes, cutoff in ((1, 12), (2, 8), (3, 6))
    for ratio in (0.25, 1.0)
]


@pytest.mark.parametrize("num_cubes,cutoff,ratio", FAST_CASES)
def test_reference_spectrum(num_cubes, cutoff, ratio):
    """Exercise higher-spin links and the two shared planes of three cubes."""
    result = _compute(num_cubes, cutoff, ratio)
    _assert_matches_reference(result, _reference(num_cubes, cutoff, ratio))


def test_single_cube_agrees_with_previously_certified_interval():
    """A finite-basis estimate agrees with the prior full-cube enclosure.

    At kappa=nu=1 the separate single-cube certificate bounds the untruncated
    gap. The retained cap-12 estimate lies inside that interval. This numerical
    containment assertion is a regression check, not a new certificate.
    """
    result = _compute(1, 12, 1.0)
    lower, upper = CERTIFIED_SINGLE_CUBE_INTERVAL
    assert lower <= result.gap <= upper


@pytest.mark.parametrize("num_cubes", (1, 2, 3))
def test_common_coupling_scale_changes_energies_not_basis(num_cubes):
    cutoff = 6
    ratio = 0.25
    scale = 2.5
    base = _compute(num_cubes, cutoff, ratio)
    model = JoinedCubesModel(
        num_cubes=num_cubes, kappa=scale, nu=scale * ratio
    )
    scaled = SpectrumSolver(model, electric_cutoff=cutoff).compute_mass_gap()
    assert scaled.dimension == base.dimension
    assert scaled.converged
    np.testing.assert_allclose(
        [scaled.e0, scaled.e1, scaled.gap],
        scale * np.array([base.e0, base.e1, base.gap]),
        atol=ENERGY_ATOL,
        rtol=0,
    )


@pytest.mark.parametrize("num_cubes", (1, 2, 3))
def test_pure_electric_gap_includes_degenerate_first_excitation(num_cubes):
    """The shortest fundamental loop costs four times j(j+1) = 3."""
    kappa = 2.0
    model = JoinedCubesModel(num_cubes=num_cubes, kappa=kappa, nu=0.0)
    result = SpectrumSolver(model, electric_cutoff=6).compute_mass_gap()
    assert result.converged
    assert result.e0 == pytest.approx(0.0, abs=1e-12)
    assert result.e1 == pytest.approx(3.0 * kappa, abs=1e-12)
    assert result.gap == pytest.approx(3.0 * kappa, abs=1e-12)
    np.testing.assert_allclose(
        result.eigenvalues[:3], [0.0, 3.0 * kappa, 3.0 * kappa],
        atol=1e-12, rtol=0,
    )


@pytest.mark.parametrize(
    "kwargs",
    [
        {"num_cubes": 0},
        {"num_cubes": -1},
        {"num_cubes": 1.5},
        {"num_cubes": 1, "kappa": 0},
        {"num_cubes": 1, "kappa": -1},
        {"num_cubes": 1, "kappa": float("inf")},
        {"num_cubes": 1, "nu": -0.25},
        {"num_cubes": 1, "nu": float("nan")},
        {"num_cubes": 1, "pairing": "unknown"},
    ],
)
def test_invalid_model_parameters_raise_value_error(kwargs):
    with pytest.raises(ValueError):
        JoinedCubesModel(**kwargs)


@pytest.mark.parametrize(
    "kwargs",
    [
        {"electric_cutoff": -1},
        {"electric_cutoff": float("inf")},
        {"electric_cutoff": float("nan")},
        {"electric_cutoff": 6, "tolerance": 0},
        {"electric_cutoff": 6, "max_iterations": 0},
        {"electric_cutoff": 6, "levels": 1},
        {"electric_cutoff": 6, "backend": "unknown"},
    ],
)
def test_invalid_solver_parameters_raise_value_error(kwargs):
    with pytest.raises(ValueError):
        SpectrumSolver(JoinedCubesModel(num_cubes=1), **kwargs).compute_mass_gap()


def test_cutoff_with_only_vacuum_cannot_define_gap():
    with pytest.raises(ValueError):
        SpectrumSolver(
            JoinedCubesModel(num_cubes=1), electric_cutoff=0
        ).compute_mass_gap()


FULL_CASES = [
    pytest.param(
        num_cubes, row["energy_cap"], row["nu_over_kappa"],
        id=f"N{num_cubes}-cap{row['energy_cap']}-nu{row['nu_over_kappa']}",
    )
    for num_cubes in (1, 2, 3)
    for row in _reference_records(num_cubes)
]


@pytest.mark.slow
@pytest.mark.parametrize("num_cubes,cutoff,ratio", FULL_CASES)
def test_all_published_spectra(num_cubes, cutoff, ratio):
    """Opt in with pytest --run-slow; largest case has 771,425 states."""
    result = _compute(num_cubes, cutoff, ratio)
    _assert_matches_reference(result, _reference(num_cubes, cutoff, ratio))
