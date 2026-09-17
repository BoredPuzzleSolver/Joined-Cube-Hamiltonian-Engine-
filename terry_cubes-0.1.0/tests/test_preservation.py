"""Protect original SU(2) indexing, exact amplitudes, and intertwiner channels."""
from __future__ import annotations

from dataclasses import asdict
from fractions import Fraction
import importlib
import math
from pathlib import Path
import sys

import numpy as np
import pytest
from sympy import Rational, sign
from sympy.physics.wigner import wigner_6j

from terry_cubes import geometry, hamiltonian, recoupling


@pytest.fixture(scope="module")
def original_modules():
    """Load the unmodified research scripts with their original import names."""
    reference = Path(__file__).parent / "reference"
    names = (
        "su2_recoupling",
        "su2_cube_model",
        "su2_oriented_vertex",
        "su2_joined_cubes_model",
        "su2_sparse_spectrum",
    )
    missing = [name for name in names if not (reference / (name + ".py")).is_file()]
    assert not missing, f"Missing original reference scripts: {missing}"
    previous = {name: sys.modules.get(name) for name in names}
    with pytest.MonkeyPatch.context() as patch:
        patch.syspath_prepend(str(reference))
        for name in names:
            sys.modules.pop(name, None)
        try:
            yield {name: importlib.import_module(name) for name in names}
        finally:
            for name, module in previous.items():
                if module is None:
                    sys.modules.pop(name, None)
                else:
                    sys.modules[name] = module


@pytest.mark.parametrize("num_cubes", [1, 2, 3])
@pytest.mark.parametrize("pairing", ["longitudinal", "crossed", "mixed"])
def test_joined_geometry_preserves_exact_order_and_phases(original_modules, num_cubes, pairing):
    original = original_modules["su2_joined_cubes_model"]
    actual = geometry.joined_graph(num_cubes, pairing)
    expected = original.joined_graph(num_cubes, pairing)
    assert asdict(actual) == asdict(expected)
    assert geometry.graph_metadata(actual) == original.graph_metadata(expected)
    assert actual.physical_edge_count == 8 * num_cubes + 4
    assert actual.face_count == 5 * num_cubes + 1
    assert len(actual.splits) == 4 * (num_cubes - 1)


def test_cube_constants_preserve_indexing(original_modules):
    original = original_modules["su2_cube_model"]
    for name in ("VERTICES", "EDGES", "FACES", "EDGE_INDEX", "INCIDENT", "FACE_EDGES", "FACE_DATA"):
        assert getattr(geometry, name) == getattr(original, name), name


def _exact_entry_dictionary(entries):
    result = {(row, col): (phase, square) for row, col, phase, square in entries}
    assert len(result) == len(entries), "A single face should not repeat a matrix entry"
    assert all(isinstance(square, Fraction) for phase, square in result.values())
    return result


@pytest.mark.parametrize("pairing", ["longitudinal", "crossed", "mixed"])
def test_joined_basis_and_exact_transitions_match_original(original_modules, pairing):
    original = original_modules["su2_joined_cubes_model"]
    graph = geometry.joined_graph(2, pairing)
    original_graph = original.joined_graph(2, pairing)
    states = hamiltonian.energy_basis(graph, Fraction(6))
    expected_states = original.energy_basis(original_graph, Fraction(6))
    assert states == expected_states  # Tuple order fixes every matrix index.
    actual = hamiltonian.wilson_transitions(graph, states)
    expected = original.wilson_transitions(original_graph, expected_states)
    assert len(actual) == len(expected) == graph.face_count
    for entries, reference_entries in zip(actual, expected):
        lookup = _exact_entry_dictionary(entries)
        assert lookup == _exact_entry_dictionary(reference_entries)
        assert entries == reference_entries  # Preserve traversal order as well.
        for (row, col), amplitude in lookup.items():
            assert lookup[(col, row)] == amplitude
    if pairing == "mixed":
        assert any(len(face) % 2 for face in graph.faces)


@pytest.mark.parametrize("pairing", ["longitudinal", "crossed", "mixed"])
def test_virtual_intertwiner_channels_exceed_physical_spin_cutoff(original_modules, pairing):
    original = original_modules["su2_joined_cubes_model"]
    graph = geometry.joined_graph(2, pairing)
    expected_graph = original.joined_graph(2, pairing)
    states = hamiltonian.spin_basis(graph, 1)
    assert states == original.spin_basis(expected_graph, 1)
    physical_count = graph.physical_edge_count
    assert all(max(state[:physical_count]) <= 1 for state in states)
    higher_channels = [state for state in states if max(state[physical_count:]) == 2]
    assert higher_channels, "Four spin-half links must retain the spin-one coupling channel"
    for state in higher_channels:
        physical_energy = Fraction(sum(q * (q + 2) for q in state[:physical_count]), 4)
        assert hamiltonian.electric_exact(graph, state) == physical_energy
        assert hamiltonian.electric_exact(graph, state) == original.electric_exact(expected_graph, state)


@pytest.mark.parametrize("max_twice", [1, 2])
def test_single_cube_basis_and_exact_wilson_entries_match_original(original_modules, max_twice):
    original = original_modules["su2_cube_model"]
    states, entries = hamiltonian.cube_wilson_transitions(max_twice)
    expected_states, expected_entries = original.cube_wilson_transitions(max_twice)
    assert states == expected_states
    assert entries == expected_entries
    for face_entries in entries:
        lookup = _exact_entry_dictionary(face_entries)
        for (row, col), amplitude in lookup.items():
            assert lookup[(col, row)] == amplitude
    np.testing.assert_array_equal(
        hamiltonian.cube_electric_diagonal(states), original.electric_diagonal(expected_states)
    )


def test_single_cube_complete_electric_cutoff_matches_original(original_modules):
    original = original_modules["su2_cube_model"]
    assert hamiltonian.cube_energy_basis(Fraction(12)) == original.energy_basis(Fraction(12))


def test_sparse_entries_preserve_face_order_and_repeated_entries(original_modules):
    original = original_modules["su2_joined_cubes_model"]
    graph = geometry.joined_graph(2, "mixed")
    states = hamiltonian.energy_basis(graph, Fraction(6))
    expected_graph = original.joined_graph(2, "mixed")
    per_face = original.wilson_transitions(expected_graph, states)
    flattened = [entry for face in per_face for entry in face]
    rows, cols, values = hamiltonian.build_sparse(graph, states)
    np.testing.assert_array_equal(rows, [row for row, col, phase, square in flattened])
    np.testing.assert_array_equal(cols, [col for row, col, phase, square in flattened])
    np.testing.assert_array_equal(
        values, [phase * math.sqrt(float(square)) for row, col, phase, square in flattened]
    )


@pytest.mark.parametrize("arguments", [
    (0, 0, 0, 1, 1, 1),
    (0, 1, 1, 0, 1, 1),
    (2, 2, 2, 2, 2, 2),
    (1, 2, 3, 2, 1, 2),
    (4, 4, 4, 4, 4, 4),
])
def test_exact_six_j_symbols_against_independent_sympy(arguments):
    reference = wigner_6j(*(Rational(q, 2) for q in arguments))
    numerator, denominator = (reference ** 2).as_numer_denom()
    assert recoupling.wigner_6j_squared(*arguments) == Fraction(int(numerator), int(denominator))
    assert recoupling.sixj_sign(*arguments) == int(sign(reference))


def test_oriented_invariant_norms_signs_and_exact_radicals():
    result = recoupling.oriented_self_checks(max_twice=3)
    assert result["exact_vertex_contractions"] == 62
    assert result["unit_3j_norms"]
    assert result["squared_6j_identity"]
    assert result["all_radical_ratios_exact_rational_squares"]
    assert result["exact_local_sign_relation_to_reference_6j_formula"]
    assert result["exact_local_dimension_ratio_to_reference_6j_formula"]
