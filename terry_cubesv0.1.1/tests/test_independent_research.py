"""Opt-in checks independent of the package's recoupling implementation."""
from fractions import Fraction
from pathlib import Path
import importlib.util
import json
import math
import os

import numpy as np
import pytest

from terry_cubes import JoinedCubesModel, SpectrumSolver, verify_certificate
from terry_cubes import geometry, hamiltonian, recoupling
from terry_cubes.data import load_research_dataset

pytestmark = pytest.mark.research
ROOT = Path(__file__).resolve().parents[1]


def module(name):
    path = ROOT/"validation/reference_joined"/(name+".py")
    spec = importlib.util.spec_from_file_location("independent_"+name, path)
    result = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(result)
    return result


def save_record(name, result):
    if os.environ.get("TERRY_VALIDATION_DIR"):
        path = Path(os.environ["TERRY_VALIDATION_DIR"])/name
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(json.dumps(result, indent=2, allow_nan=False)+"\n", encoding="utf-8")


def test_full_exact_orientation_sweep_and_basis_phase_identity():
    result = recoupling.oriented_self_checks(max_twice=14)
    assert result["exact_vertex_contractions"] == 3368
    for n in (1, 2, 3, 4):
        graph = geometry.joined_graph(n)
        selected = {s.auxiliary_edge for s in graph.splits
                    if graph.real_coordinates[s.original_vertex][1] != graph.real_coordinates[s.original_vertex][2]}
        for (cycle, _), phase in zip(graph.face_data, graph.face_orientation_signs):
            assert (-1)**sum(e in selected for e in cycle) == phase
    save_record("orientation_fresh.json", result)


@pytest.mark.parametrize("n,expected", [(1, 32), (2, 868), (3, 25676)])
def test_independent_cycle_space_counts_match_package_states(n, expected):
    reference = module("su2_joined_cubes_independent_check")
    states, _ = reference.full_spin_half_states(n)
    assert len(states) == expected
    assert states == hamiltonian.spin_basis(geometry.joined_graph(n), 1)


@pytest.mark.parametrize("n,cap,pairing", [(1, 12, "longitudinal"), (2, 8, "longitudinal"),
    (2, 8, "crossed"), (2, 8, "mixed"), (3, 6, "longitudinal")])
def test_independent_casimir_tensor_face_matrices(n, cap, pairing):
    probe = module("su2_joined_cubes_tensor_probe")
    model = geometry.joined_graph(n, pairing)
    states = hamiltonian.energy_basis(model, Fraction(cap))
    transitions = hamiltonian.wilson_transitions(model, states)
    constraints = [[] for _ in states]
    for face, entries in zip(model.faces, transitions):
        independent = probe.tensor_face_matrix(states, model.edges, model.vertices, face)
        actual = np.zeros_like(independent)
        for row, col, phase, square in entries:
            actual[row, col] = phase*math.sqrt(float(square))
            ratio = phase*(1 if independent[row, col] > 0 else -1)
            constraints[row].append((col, ratio))
        np.testing.assert_allclose(np.abs(independent), np.abs(actual), atol=1e-12, rtol=0)
    phases = {}
    for start in range(len(states)):
        if start in phases:
            continue
        phases[start] = 1
        queue = [start]
        for state in queue:
            for target, ratio in constraints[state]:
                wanted = phases[state]*ratio
                if target in phases:
                    assert phases[target] == wanted
                else:
                    phases[target] = wanted
                    queue.append(target)


@pytest.mark.parametrize("n", [2, 3])
def test_iterative_low_levels_against_full_dense_diagonalization(n):
    model = JoinedCubesModel(n)
    solver = SpectrumSolver(model, 8, levels=6, backend="arpack", max_iterations=1000)
    dense = np.linalg.eigvalsh(solver.build_hamiltonian().to_scipy().toarray())
    arpack = solver.compute_mass_gap()
    np.testing.assert_allclose(arpack.eigenvalues, dense[:6], atol=1e-10, rtol=0)
    solver.backend = "lanczos"
    lanczos = solver.compute_mass_gap()
    np.testing.assert_allclose(lanczos.eigenvalues[:2], dense[:2], atol=1e-10, rtol=0)
    # Single-vector Lanczos is deliberately not asserted to recover multiplicities.


@pytest.mark.parametrize("cap,ratio", [(24, 1), (30, 1), (24, 4), (30, 4)])
def test_extended_single_cube_reference(cap, ratio):
    rows = load_research_dataset("single_cube_energy")["record"]["spectra"]
    old = next(r for r in rows if r["electric_energy_cap_over_kappa"] == cap and r["nu_over_kappa"] == ratio)
    result = SpectrumSolver(JoinedCubesModel(1, nu=ratio), cap, backend="arpack", max_iterations=1000).compute_mass_gap()
    assert result.dimension == old["dimension"]
    np.testing.assert_allclose([result.e0, result.e1, result.gap], [old["E0"], old["E1"], old["gap"]], atol=1e-9, rtol=0)
    save_record(f"single_cube_cap{cap}_ratio{ratio}_fresh.json", result.to_dict())


def test_strongest_single_cube_certificate():
    result = verify_certificate(cutoff=12)
    assert result["certified"] and result["package_bridge_verified"]
    assert result["gap_interval_exact"] == ["2941733/1000000", "2958693/1000000"]
    assert result["comparison_dimension"] == 2042
    save_record("single_cube_certificate_fresh.json", result)


def test_all_original_two_plaquette_certificates():
    result = verify_certificate("two_plaquette")
    assert result["certified"] and result["package_bridge_verified"]
    assert len(result["benchmarks"]) == 3
    assert result["benchmarks"][0]["certificate"]["gap_interval_exact"] == ["154597949/50000000", "309195907/100000000"]
    save_record("two_plaquette_certificate_fresh.json", result)
