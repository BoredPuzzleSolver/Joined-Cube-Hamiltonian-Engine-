"""Regression and independent checks for the added research interfaces."""
import ast
import json
from hashlib import sha256
from pathlib import Path

import numpy as np
import pytest

from terry_cubes import (JoinedCubesModel, SpectrumSolver, TwoPlaquetteModel,
                         connected_correlator, cutoff_study, plaquette_operator, __version__)
from terry_cubes.cli import main
from terry_cubes.data import load_spectra, load_research_dataset, list_research_datasets
from terry_cubes.certification import verify_certificate, CertificateVerificationError


def test_original_math_files_and_extra_datasets_are_preserved():
    root = Path(__file__).resolve().parents[1]
    manifest = json.loads((root / "RELEASE_PROVENANCE.json").read_text())
    import terry_cubes
    package = Path(terry_cubes.__file__).parent
    for name, expected in manifest["unchanged_runtime_lf_sha256"].items():
        assert sha256((package/name).read_bytes().replace(b"\r\n", b"\n")).hexdigest() == expected
    for name, expected in manifest["added_reference_data_lf_sha256"].items():
        assert sha256((root/"data"/name).read_bytes().replace(b"\r\n", b"\n")).hexdigest() == expected
    # Python 3.12 added AST fields that do not exist in Python 3.10/3.11.
    # Parse both sources with this interpreter instead of comparing ast.dump()
    # against a fingerprint produced by a different Python version.
    reference_name = "su2_sparse_spectrum.py"
    reference_bytes = (root/"tests"/"reference"/reference_name).read_bytes()
    source_manifest = json.loads((root/"SOURCE_PROVENANCE.json").read_text())
    expected_reference = source_manifest["source_lf_normalized_sha256"][reference_name]
    assert sha256(reference_bytes.replace(b"\r\n", b"\n")).hexdigest() == expected_reference
    tree = ast.parse((package/"solver.py").read_text(encoding="utf-8"))
    reference_tree = ast.parse(reference_bytes.decode("utf-8"))
    names = set(manifest["unchanged_solver_ast_sha256"])
    actual_nodes = {node.name: node for node in tree.body
                    if getattr(node, "name", None) in names}
    reference_nodes = {node.name: node for node in reference_tree.body
                       if getattr(node, "name", None) in names}
    assert set(actual_nodes) == names, "A preserved solver definition is missing"
    assert set(reference_nodes) == names, "An original solver definition is missing"
    for name in sorted(names):
        assert ast.dump(actual_nodes[name], include_attributes=False) == ast.dump(
            reference_nodes[name], include_attributes=False
        ), f"Original solver definition changed: {name}"
    assert __version__ == "0.1.1"


@pytest.mark.parametrize("nu", [0.0, 1.0])
def test_eigenvectors_preserve_basis_order_and_direct_equation(nu):
    solver = SpectrumSolver(JoinedCubesModel(1, nu=nu), 6, levels=4, backend="arpack")
    eigen = solver.compute_eigensystem()
    matrix = solver.build_hamiltonian().to_scipy()
    assert eigen.states == solver.build_hamiltonian().states
    np.testing.assert_allclose(matrix @ eigen.eigenvectors,
                               eigen.eigenvectors * eigen.spectrum.eigenvalues, atol=1e-10, rtol=0)
    np.testing.assert_allclose(eigen.eigenvectors.T @ eigen.eigenvectors, np.eye(4), atol=1e-12)
    assert not eigen.complete_retained_spectrum


def test_plaquette_observables_reconstruct_the_existing_hamiltonian():
    model = JoinedCubesModel(2, kappa=2, nu=.5, pairing="mixed")
    solver = SpectrumSolver(model, 6)
    data = solver.build_hamiltonian()
    wilson = sum(plaquette_operator(model, data.states, f) for f in range(model.graph.face_count))
    rebuilt = np.diag(data.diagonal) - model.nu * wilson.toarray()
    np.testing.assert_allclose(rebuilt, data.to_scipy().toarray(), atol=0, rtol=0)


def test_correlator_matches_an_independent_three_level_spectral_sum():
    e = np.array([0., 2., 5.])
    operator = np.array([[0., 1., 2.], [1., 0., 0.], [2., 0., 0.]])
    times = [0., .3, 2.]
    result = connected_correlator(e, np.eye(3), operator, times, complete_spectrum=True)
    for row in result["correlators"]:
        t = row["time"]
        expected = np.exp(-2*t)+4*np.exp(-5*t)
        later = np.exp(-2*(t+.25))+4*np.exp(-5*(t+.25))
        assert row["C_connected"] == pytest.approx(expected, rel=1e-14)
        assert row["effective_gap"] == pytest.approx(np.log(expected/later)/.25)
    partial = connected_correlator(e[:2], np.eye(3)[:, :2], operator, [1.])
    assert not partial["complete_retained_spectrum"]
    assert partial["unrepresented_weight"] == pytest.approx(4.)
    far = connected_correlator(e, np.eye(3), operator, [1000000.])
    assert far["correlators"][0]["C_connected"] == 0
    assert far["correlators"][0]["effective_gap"] == pytest.approx(2.)
    empty = connected_correlator(e, np.eye(3), np.zeros((3, 3)), [0., 1.])
    assert empty["lowest_detected_gap"] is None
    assert empty["correlators"][0]["effective_gap"] is None


@pytest.mark.parametrize("kwargs", [{"dt": 0}, {"weight_threshold": -1}, {"dt": float("nan")}])
def test_correlator_rejects_invalid_parameters(kwargs):
    with pytest.raises(ValueError):
        connected_correlator([0., 1.], np.eye(2), np.eye(2), [0.], **kwargs)


def test_cutoff_study_matches_original_larger_basis_diagnostics():
    result = cutoff_study(JoinedCubesModel(1), [8, 6, 8], embedding_cutoff=18)
    assert [r["electric_cutoff"] for r in result["results"]] == ["6", "8"]
    reference = load_spectra(1)["spectra"]
    for row in result["results"]:
        old = next(x for x in reference if x["energy_cap"] == int(row["electric_cutoff"]) and x["nu_over_kappa"] == 1.)
        np.testing.assert_allclose([row["e0"], row["e1"], row["gap"]],
                                   [old["E0"], old["E1"], old["gap"]], atol=1e-9, rtol=0)
        np.testing.assert_allclose(row["residual_in_largest_test_basis"],
                                   old["residual_in_largest_test_basis"], atol=1e-9, rtol=0)
        assert row["nested_energy_ordering_consistent"] and not row["certified"]
    assert result["results"][-1]["second_seed_max_energy_difference"] < 1e-9


@pytest.mark.parametrize("cutoffs,kwargs", [([], {}), ([0, 3], {}), ([6, 8], {"embedding_cutoff": 6}),
                                          ([6], {"second_seed": 1097})])
def test_invalid_cutoff_study(cutoffs, kwargs):
    with pytest.raises(ValueError):
        cutoff_study(JoinedCubesModel(1), cutoffs, **kwargs)


@pytest.mark.parametrize("max_twice,ratio", [(1, 1), (4, 1), (6, 4), (1, 0)])
def test_two_square_public_model_matches_historical_reference(max_twice, ratio):
    result = TwoPlaquetteModel(nu=ratio).compute_spectrum(max_twice)
    rows = load_research_dataset("two_plaquette_study")["record"]["cutoff_sweep"]
    old = next(r for r in rows if r["jmax"] == max_twice/2 and r["nu_over_kappa"] == ratio)
    np.testing.assert_allclose([result.eigenvalues[0], result.eigenvalues[1], result.gap],
                               [old["E0"], old["E1"], old["gap"]], atol=1e-9, rtol=0)
    assert result.to_dict()["model_id"] == "two_plaquette"
    assert not result.to_dict()["certified"]


def test_two_square_correlators_preserve_the_recorded_symmetry_selection():
    result = TwoPlaquetteModel().compute_spectrum(12)
    actual = result.correlators()
    old = load_research_dataset("two_plaquette_study")["record"]["correlation_study"]
    for channel in old["channels"]:
        new = actual[channel["name"]]
        assert new["lowest_detected_gap"] == pytest.approx(channel["lowest_detected_gap"], abs=1e-9)
        for row, reference in zip(new["correlators"], channel["correlators"]):
            assert row["effective_gap"] == pytest.approx(reference["effective_gap"], abs=1e-9)
    assert actual["symmetric"]["lowest_detected_gap"] == pytest.approx(result.gap, abs=1e-9)
    assert actual["antisymmetric"]["lowest_detected_gap"] > result.gap + .01


def test_research_loaders_distinguish_reading_from_certification():
    assert len(list_research_datasets()) == 7
    for item in list_research_datasets():
        record = load_research_dataset(item["name"])
        assert record["model_id"] in ("single_cube", "two_plaquette")
        assert record["verified_in_this_call"] is False
    with pytest.raises(ValueError):
        load_research_dataset("../../other")


@pytest.mark.parametrize("args", [
    ["study", "--cubes", "1", "--cutoffs", "3", "6"],
    ["two-plaquette", "--max-twice", "4", "--correlators"],
    ["data"], ["data", "--name", "single_cube_energy"],
])
def test_new_cli_commands_produce_strict_json(args, tmp_path):
    path = tmp_path/"result.json"
    assert main(args + ["--output", str(path)]) == 0
    assert isinstance(json.loads(path.read_text()), dict)


def test_unsupported_certificate_cannot_overwrite_a_result(tmp_path):
    path = tmp_path/"result.json"
    path.write_text("previous")
    assert main(["certify", "single-cube", "--cutoff", "14", "--output", str(path)]) == 2
    assert path.read_text() == "previous"
    with pytest.raises(ValueError):
        verify_certificate("joined_cubes")
    with pytest.raises(ValueError):
        verify_certificate("two_plaquette", cutoff=12)


def test_small_certificate_is_fresh_and_checks_omitted_couplings():
    result = verify_certificate(cutoff=8)
    assert result["certified"] and result["package_bridge_verified"]
    assert result["comparison_dimension"] == 508
    assert result["source_manifest_files_verified"] == 42
    assert result["gap_interval_exact"] == ["2236257/1000000", "3068193/1000000"]
    assert result["package_transitions_checked_including_omitted_targets"] > 0
    assert result["package_electric_energies_checked"] == result["comparison_dimension"]


def test_certificate_timeout_fails_without_claiming_success():
    with pytest.raises(CertificateVerificationError, match="timeout"):
        verify_certificate(cutoff=8, timeout=.001)
