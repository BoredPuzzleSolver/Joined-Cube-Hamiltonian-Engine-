"""Three spatial face orientations: numerical open-cube SU(2) study."""
from pathlib import Path
import json
import math
import time
import numpy as np
from su2_cube_model import cube_wilson_transitions, electric_diagonal, EDGES, FACES
from su2_sparse_spectrum import SymmetricOperator, smallest_ritz

ROOT = Path(__file__).resolve().parent.parent


def main():
    started = time.perf_counter()
    states, face_entries = cube_wilson_transitions(3)
    entries = [entry for face in face_entries for entry in face]
    rows = np.array([e[0] for e in entries],dtype=np.int64)
    cols = np.array([e[1] for e in entries],dtype=np.int64)
    data = np.array([e[2]*math.sqrt(float(e[3])) for e in entries])
    assert not np.any(rows==cols)
    by_pair = {(r,c):(sgn,sq) for r,c,sgn,sq in entries}
    assert len(by_pair)==len(entries), 'Unexpected overlap of distinct face transitions'
    assert all(by_pair[(c,r)]==value for (r,c),value in by_pair.items())
    print(f'J=1.5 basis {len(states)}; directed transitions {len(entries)}; construction {time.perf_counter()-started:.2f}s',flush=True)
    t = electric_diagonal(states)
    cutoffs = np.max(np.array(states),axis=1)
    results = []
    previous = {}
    for cutoff in [1,2,3]:
        mask = cutoffs <= cutoff
        keep = mask[rows]&mask[cols]
        mapping = np.full(len(states),-1,dtype=np.int64)
        mapping[mask] = np.arange(np.sum(mask))
        rr,cc,dd = mapping[rows[keep]],mapping[cols[keep]],data[keep]
        for nu in [.25,1,4]:
            op = SymmetricOperator(t[mask]+6*nu,rr,cc,-nu*dd)
            if op.size <= 1100:
                dense = op.dense()
                values,vectors = np.linalg.eigh(dense)
                eigs = values[:2]
                residuals = [float(np.linalg.norm(dense@vectors[:,i]-eigs[i]*vectors[:,i])) for i in range(2)]
                solver = {'method':'dense NumPy eigh','residual_norms':residuals,'converged':max(residuals)<1e-10}
                low_levels = values[:8].tolist()
            else:
                eigs,vectors,solver = smallest_ritz(op,levels=2,max_iterations=240,tolerance=1e-10)
                assert solver['converged'], solver
                check_values,_,second = smallest_ritz(op,levels=2,max_iterations=240,tolerance=1e-10,seed=4811)
                assert second['converged'], second
                assert np.max(abs(eigs-check_values))<1e-9
                solver['second_seed_max_energy_difference'] = float(np.max(abs(eigs-check_values)))
                solver['second_seed_solver'] = second
                low_levels = eigs.tolist()
            assert solver['converged']
            if nu in previous:
                assert np.all(eigs<=previous[nu]+1e-9)
            previous[nu] = eigs
            row = {'jmax':cutoff/2,'dimension':op.size,'nu_over_kappa':nu,
                   'E0':float(eigs[0]),'E1':float(eigs[1]),'gap':float(eigs[1]-eigs[0]),
                   'first_eigenvalues':low_levels,'directed_magnetic_nonzeros':len(dd),
                   'solver':solver}
            results.append(row)
            print(f'J={cutoff/2:g}, nu={nu:g}: E0={eigs[0]:.12f}, E1={eigs[1]:.12f}, gap={row["gap"]:.12f}',flush=True)
    out = {
        'title':'SU(2) on the open single cube',
        'graph':{'vertices':8,'edges':EDGES,'faces':FACES,'independent_cycles':5,
                 'description':'Twelve links, six square faces in xy/xz/yz planes, all vertices on boundary'},
        'hamiltonian':'H/kappa=sum_links j(j+1)+(nu/kappa)*(6-sum_faces Tr_fund(U_face)/2)',
        'gauge_constraint':'Gauss law at all eight trivalent vertices; no external charges',
        'basis':'All triangle-admissible doubled-link-spin configurations with every link spin <= jmax; unique trivalent intertwiners',
        'units':'kappa=1; no physical energy calibration is inferred',
        'status':'Numerical projected spectra; omitted-state certificates are separate',
        'exact_electric_limit':{'nu_over_kappa':0,'E0':0,'E1':3,'gap':3,'first_excitation_degeneracy':6},
        'matrix_checks':{'exact_signed_square_hermiticity':True,'distinct_face_transition_supports':True,
                         'largest_basis_dimension':len(states),'largest_directed_nonzeros':len(entries)},
        'spectra':results,
    }
    output = ROOT/'outputs'/'SU2_Cube_Study_Results.json'
    output.write_text(json.dumps(out,indent=2)+'\n',encoding='utf-8')
    print(output,flush=True)


if __name__ == '__main__':
    main()
