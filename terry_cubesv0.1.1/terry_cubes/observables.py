"""Normalized plaquette operators and finite-basis connected correlations."""
import math
import numpy as np
from scipy.sparse import coo_matrix
from scipy.special import logsumexp

from .hamiltonian import JoinedCubesModel, state_face_transitions, _integer


def plaquette_operator(model, states, face_index):
    """Return P(Tr_fund U_face/2)P in the supplied, unchanged basis order."""
    if not isinstance(model, JoinedCubesModel):
        raise TypeError("model must be a JoinedCubesModel")
    face_index = _integer(face_index, "face_index", 0)
    if face_index >= model.graph.face_count:
        raise ValueError("face_index is outside the graph")
    states = tuple(states)
    index = {state: i for i, state in enumerate(states)}
    if len(index) != len(states):
        raise ValueError("states must not contain duplicate labels")
    rows, cols, entries = [], [], []
    for col, state in enumerate(states):
        for target, sign, square in state_face_transitions(model.graph, state, face_index):
            row = index.get(target)
            if row is not None:
                rows.append(row)
                cols.append(col)
                entries.append(sign * math.sqrt(float(square)))
    return coo_matrix((entries, (rows, cols)), shape=(len(states), len(states))).tocsr()


def connected_correlator(energies, eigenvectors, operator, times, *, dt=0.25,
                         weight_threshold=1e-18, complete_spectrum=False):
    """Compute sum_{n>0}|<n|O|0>|^2 exp[-(En-E0)t] for supplied modes.

    Imaginary time uses inverse Hamiltonian energy units (hbar=1). A partial
    eigensystem gives only a partial spectral sum. A symmetry-forbidden
    channel can miss the lowest gap even with a complete finite spectrum.
    Small weights are explicitly excluded and reported, as in the original
    two-square study. This is not an independent mass measurement.
    """
    e = np.asarray(energies, dtype=float)
    vectors = np.asarray(eigenvectors)
    times = np.asarray(list(times), dtype=float)
    if e.ndim != 1 or len(e) < 2 or not np.all(np.isfinite(e)) or np.any(np.diff(e) < 0):
        raise ValueError("energies must be a finite ascending sequence with at least two levels")
    if vectors.ndim != 2 or vectors.shape[1] != len(e) or not np.all(np.isfinite(vectors)):
        raise ValueError("eigenvectors must have one finite column per energy")
    if operator.shape != (vectors.shape[0], vectors.shape[0]):
        raise ValueError("operator shape must match the eigenvector basis")
    if not np.allclose(vectors.conj().T @ vectors, np.eye(len(e)), rtol=0, atol=1e-8):
        raise ValueError("eigenvectors must be orthonormal")
    if not np.isfinite(dt) or dt <= 0 or not np.isfinite(weight_threshold) or weight_threshold < 0:
        raise ValueError("dt must be positive and weight_threshold nonnegative, both finite")
    if times.ndim != 1 or not len(times) or not np.all(np.isfinite(times)) or np.any(times < 0):
        raise ValueError("times must be a nonempty finite sequence of nonnegative values")
    if complete_spectrum and vectors.shape[0] != len(e):
        raise ValueError("complete_spectrum requires a full orthonormal basis")
    if e[1] <= e[0]:
        raise ValueError("this correlator requires a unique lowest input energy")
    applied = operator @ vectors[:, 0]
    if not np.all(np.isfinite(applied)):
        raise ValueError("operator action must be finite")
    weights = np.abs(vectors[:, 1:].conj().T @ applied) ** 2
    gaps = e[1:] - e[0]
    selected = weights > weight_threshold
    expectation = np.vdot(vectors[:, 0], applied)
    if abs(expectation.imag) > 1e-9:
        raise ValueError("operator must have a real vacuum expectation")
    rows = []
    for t in times:
        if np.any(selected):
            log_weights = np.log(weights[selected])
            log0 = float(logsumexp(log_weights - gaps[selected] * t))
            log1 = float(logsumexp(log_weights - gaps[selected] * (t + dt)))
            if not np.isfinite(log0) or not np.isfinite(log1):
                raise ValueError("time and energy scales exceed floating-point range")
            c, effective = math.exp(log0), (log0 - log1) / dt
        else:
            c, effective = 0.0, None
        rows.append({"time": float(t), "C_connected": c, "effective_gap": effective})
    return {"vacuum_expectation": float(expectation.real),
            "first_excitation_overlap_squared": float(weights[0]),
            "lowest_detected_gap": float(gaps[selected][0]) if np.any(selected) else None,
            "spectral_weight_threshold": float(weight_threshold),
            "discarded_weight": float(weights[~selected].sum()), "dt": float(dt),
            "complete_retained_spectrum": bool(complete_spectrum),
            "unrepresented_weight": max(0.0, float(np.vdot(applied, applied).real - abs(expectation)**2 - weights.sum())),
            "correlators": rows, "certified": False,
            "scope": "Connected imaginary-time spectral sum in a finite basis; partial modes or symmetry selection may miss the lowest excitation. No physical mass units are inferred."}
