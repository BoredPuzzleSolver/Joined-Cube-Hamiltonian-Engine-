"""Reproducible cutoff and coupling study of two shared SU(2) plaquettes.

All energies use kappa=1. Numerical values are distinguished from the
separate exact rational certificates in verify_su2_cutoff_certificates.py.
Requires Python and NumPy; run from any working directory.
"""
from pathlib import Path
import json
import math
import numpy as np

from su2_recoupling import wilson_matrices, electric_diagonal, self_checks
from su2_cutoff_bounds import numerical_cutoff_bounds


def main():
    cutoffs = [1, 2, 3, 4, 6, 8, 10, 12]  # doubled spin
    couplings = [0, .1, .25, .5, 1, 2, 4, 10, 40]
    states, wa, wb = wilson_matrices(max(cutoffs)+1)
    electric = electric_diagonal(states)
    sizes = np.array([max(s) for s in states])
    rows = []
    for coupling in couplings:
        previous = None
        for cutoff in cutoffs:
            ids = np.flatnonzero(sizes <= cutoff)
            w = (wa+wb)[np.ix_(ids, ids)]
            h = np.diag(electric[ids]+2*coupling)-coupling*w
            values = np.linalg.eigvalsh(h)
            if previous is not None:
                assert np.all(values[:3] <= previous+1e-10)
            row = {
                'jmax': cutoff/2, 'dimension': len(ids), 'nu_over_kappa': coupling,
                'E0': float(values[0]), 'E1': float(values[1]),
                'E2': float(values[2]), 'gap': float(values[1]-values[0]),
                'change_in_E0': None if previous is None else float(values[0]-previous[0]),
                'change_in_E1': None if previous is None else float(values[1]-previous[1]),
            }
            # Omitted-state bounds at representative couplings, all cutoffs.
            if coupling in [1, 4, 40]:
                ext_ids = np.flatnonzero(sizes <= cutoff+1)
                ext_h = (np.diag(electric[ext_ids]+2*coupling)
                         -coupling*(wa+wb)[np.ix_(ext_ids, ext_ids)])
                p = sizes[ext_ids] <= cutoff
                row['analytic_bounds_evaluated_numerically'] = numerical_cutoff_bounds(
                    ext_h, electric[ext_ids], p, cutoff
                )
            rows.append(row)
            previous = values[:3]
        print(f"nu/kappa={coupling:g}: J=6 gap={rows[-1]['gap']:.12f}", flush=True)

    # Connected Euclidean-time correlators on the same J=6 Hamiltonian.
    cutoff = 12
    ids = np.flatnonzero(sizes <= cutoff)
    basis = [states[i] for i in ids]
    a, b = wa[np.ix_(ids, ids)], wb[np.ix_(ids, ids)]
    h = np.diag(electric[ids]+2)-(a+b)
    energies, vectors = np.linalg.eigh(h)
    vacuum = vectors[:, 0]
    lookup = {s:i for i,s in enumerate(basis)}
    exchange = [lookup[(s[1], s[0], s[2])] for s in basis]
    parity = [float(vectors[:, n] @ vectors[exchange, n]) for n in range(6)]
    gaps = energies[1:]-energies[0]
    channels = []
    for name, operator in [('symmetric', (a+b)/math.sqrt(2)),
                           ('antisymmetric', (a-b)/math.sqrt(2))]:
        weights = abs(vectors[:, 1:].T @ (operator @ vacuum))**2
        threshold = 1e-18
        selected = np.flatnonzero(weights > threshold)
        assert len(selected)
        correlator_rows = []
        dt = .25
        for t in [.25, .5, 1, 2, 4, 8]:
            # Tiny forbidden-channel roundoff is excluded explicitly.
            c0 = float(weights[selected] @ np.exp(-gaps[selected]*t))
            c1 = float(weights[selected] @ np.exp(-gaps[selected]*(t+dt)))
            correlator_rows.append({'time': t, 'C_connected': c0,
                                    'effective_gap': math.log(c0/c1)/dt})
        channels.append({
            'name': name, 'vacuum_expectation': float(vacuum @ operator @ vacuum),
            'first_excitation_overlap_squared': float(weights[0]),
            'lowest_detected_gap': float(gaps[selected[0]]),
            'spectral_weight_threshold': threshold,
            'discarded_weight': float(weights[weights <= threshold].sum()),
            'dt': dt, 'correlators': correlator_rows,
        })
    checks = self_checks()
    assert abs(rows[4*len(cutoffs)]['gap']-3.1596028310212443) < 1e-12
    result = {
        'model': 'Open two-square SU(2) lattice: 7 links, 6 vertices, Gauss law at every vertex',
        'hamiltonian': 'H/kappa=sum_link j(j+1)+(nu/kappa)*(2-W_A-W_B); W=Tr_fund(U)/2',
        'boundary': 'Open; no external charges',
        'units': 'kappa=1; Euclidean time measured in inverse kappa; hbar=1',
        'cutoff': 'Every link spin <= jmax; exact orthogonal projection, no cyclic wrap',
        'status': 'Floating-point study; exact fixed-graph certificates are supplied separately',
        'cutoff_sweep': rows,
        'correlation_study': {
            'jmax': 6, 'dimension': len(ids), 'nu_over_kappa': 1,
            'first_six_energies': energies[:6].tolist(),
            'exchange_parities': parity, 'channels': channels,
            'interpretation': 'Correlators from the same Hamiltonian test operator overlap, not independent data',
        },
        'matrix_validation': checks,
    }
    output = Path(__file__).resolve().parent.parent/'outputs'/'Shared_SU2_Cutoff_Study_Results.json'
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(result, indent=2)+'\n', encoding='utf-8')
    print(output)


if __name__ == '__main__':
    main()
