"""CLI and solver failure paths, in addition to physics regression tests."""
import json
import subprocess
import sys

import numpy as np
import pytest

from terry_cubes import JoinedCubesModel, SpectrumSolver, SpectrumConvergenceError
from terry_cubes.cli import main
from terry_cubes.data import load_spectra
from terry_cubes.solver import SymmetricOperator, smallest_ritz


def test_cli_json_file_matches_api(tmp_path, capsys):
    path = tmp_path / "nested" / "results.json"
    assert main(["run", "--cubes", "2", "--nu", ".25", "--cutoff", "6",
                 "--output", str(path)]) == 0
    saved = json.loads(path.read_text(encoding="utf-8"))
    baseline = next(row for row in load_spectra(2)["spectra"]
                    if row["energy_cap"] == 6 and row["nu_over_kappa"] == .25)
    assert saved["e0"] == pytest.approx(baseline["E0"], abs=1e-9, rel=0)
    assert saved["e1"] == pytest.approx(baseline["E1"], abs=1e-9, rel=0)
    assert saved["certified"] is False
    assert saved["converged"]
    assert saved["num_cubes"] == 2


def test_cli_stdout_and_module_entrypoint():
    process = subprocess.run(
        [sys.executable, "-m", "terry_cubes", "run", "--cubes", "1",
         "--nu", "0", "--cutoff", "3"],
        text=True, capture_output=True, check=True,
    )
    result = json.loads(process.stdout)
    assert result["gap"] == 3.0


def test_invalid_cli_keeps_existing_output(tmp_path, capsys):
    path = tmp_path / "results.json"
    path.write_text("previous result", encoding="utf-8")
    assert main(["run", "--cubes", "0", "--output", str(path)]) == 2
    assert path.read_text(encoding="utf-8") == "previous result"
    assert "num_cubes" in capsys.readouterr().err


def test_unconverged_solver_raises_and_cli_writes_no_result(tmp_path, capsys):
    model = JoinedCubesModel(2)
    with pytest.raises(SpectrumConvergenceError) as caught:
        SpectrumSolver(model, 6, max_iterations=3).compute_mass_gap()
    assert caught.value.result is not None
    assert not caught.value.result.converged
    path = tmp_path / "failed.json"
    assert main(["run", "--cubes", "2", "--cutoff", "6",
                 "--max-iterations", "3", "--output", str(path)]) == 2
    assert not path.exists()


def test_arpack_reproduces_the_same_retained_hamiltonian():
    model = JoinedCubesModel(2, nu=.25)
    original = SpectrumSolver(model, 8).compute_mass_gap()
    arpack = SpectrumSolver(model, 8, backend="arpack").compute_mass_gap()
    np.testing.assert_allclose(
        [arpack.e0, arpack.e1, arpack.gap],
        [original.e0, original.e1, original.gap], atol=1e-9, rtol=0,
    )
    assert max(arpack.residual_norms) <= 1e-10


def test_original_sparse_lanczos_against_analytic_spectrum():
    n = 180
    rows = np.concatenate([np.arange(n - 1), np.arange(1, n)])
    columns = np.concatenate([np.arange(1, n), np.arange(n - 1)])
    operator = SymmetricOperator(np.full(n, 2.), rows, columns, np.full(len(rows), -1.))
    values, _, info = smallest_ritz(operator, max_iterations=n)
    expected = 2 - 2 * np.cos(np.arange(1, 4) * np.pi / (n + 1))
    np.testing.assert_allclose(values, expected, atol=1e-11, rtol=0)
    assert info["converged"]


@pytest.mark.parametrize("cubes", [1, 2, 3])
def test_installed_dataset_loader(cubes):
    data = load_spectra(cubes)
    assert data["spectra"][0]["num_cubes"] == cubes
    assert "Numerical finite-basis" in data["status"]



def test_changing_notebook_parameters_rebuilds_the_hamiltonian():
    solver = SpectrumSolver(JoinedCubesModel(1), 6)
    small = solver.compute_mass_gap()
    solver.electric_cutoff = 8
    larger = solver.compute_mass_gap()
    assert small.dimension == 32
    assert larger.dimension == 86
    assert larger.gap != pytest.approx(small.gap, abs=1e-9, rel=0)
    assert larger.electric_cutoff == "8"
    solver.model = JoinedCubesModel(1, kappa=2, nu=2)
    scaled = solver.compute_mass_gap()
    assert scaled.dimension == larger.dimension
    assert scaled.gap == pytest.approx(2 * larger.gap, abs=1e-9, rel=0)


def test_mutated_solver_settings_are_revalidated():
    solver = SpectrumSolver(JoinedCubesModel(1), 6)
    solver.backend = "unknown"
    with pytest.raises(ValueError, match="backend"):
        solver.compute_mass_gap()


def test_result_json_handles_unrepresentable_normalized_ratios():
    result = SpectrumSolver(JoinedCubesModel(1, kappa=5e-324, nu=1), 3).compute_mass_gap()
    saved = result.to_dict()
    json.dumps(saved, allow_nan=False)
    assert saved["nu_over_kappa"] is None
    assert np.isfinite(saved["e0"])
    assert np.isfinite(saved["e1"])
