"""Complete-electric-cutoff spectra for open chains of SU(2) cubes.

Run only with the orientation-aware model and its independent checks.
Finite Ritz gaps are estimates, not certified infinite-Hilbert-space gaps.
"""
from fractions import Fraction
from pathlib import Path
from array import array
import argparse
import json
import math
import time
import numpy as np
from su2_joined_cubes_model import (joined_graph, energy_basis,
    electric_diagonal, state_face_transitions, graph_metadata)
from su2_sparse_spectrum import SymmetricOperator, smallest_ritz

HERE=Path(__file__).resolve().parent
OUT=HERE.parent/'outputs' if HERE.name=='work' else HERE


def build_sparse(model, states):
    index={state:i for i,state in enumerate(states)}
    rr,cc,dd=array('q'),array('q'),array('d')
    for face in range(model.face_count):
        for col,state in enumerate(states):
            for target,sign,square in state_face_transitions(model,state,face):
                row=index.get(target)
                if row is not None:
                    rr.append(row);cc.append(col);dd.append(sign*math.sqrt(float(square)))
    return np.array(rr),np.array(cc),np.array(dd)


def numerical_hermiticity(rows,cols,data,n):
    # Entries on different faces are retained separately; compare their sum.
    keys=rows*n+cols
    unique,inverse=np.unique(keys,return_inverse=True)
    summed=np.bincount(inverse,weights=data)
    transposed=(unique%n)*n+unique//n
    order=np.argsort(transposed)
    assert np.array_equal(unique,transposed[order])
    error=float(np.max(np.abs(summed-summed[order]))) if len(data) else 0.
    assert error<1e-12,error
    return error


def study(cubes,max_cap,caps):
    start=time.perf_counter()
    model=joined_graph(cubes)
    states=energy_basis(model,Fraction(max_cap))
    electric=electric_diagonal(model,states)
    print(f'N={cubes} cap={max_cap}: {len(states)} states; enumerated in {time.perf_counter()-start:.2f}s',flush=True)
    rows,cols,data=build_sparse(model,states)
    hermiticity=numerical_hermiticity(rows,cols,data,len(states))
    print(f'N={cubes}: {len(data)} directed face entries; constructed in {time.perf_counter()-start:.2f}s',flush=True)
    spectra=[];previous={}
    for cap in caps:
        mask=electric<=cap
        remap=np.full(len(states),-1,dtype=np.int64)
        remap[mask]=np.arange(np.sum(mask))
        keep=mask[rows]&mask[cols]
        local_rows,local_cols,local_data=remap[rows[keep]],remap[cols[keep]],data[keep]
        for nu in (0.25,1.0):
            op=SymmetricOperator(electric[mask]+model.face_count*nu,
                local_rows,local_cols,-nu*local_data)
            vals,vecs,info=smallest_ritz(op,levels=2,max_iterations=240,tolerance=1e-10)
            assert info['converged'],info
            if nu in previous:
                assert np.all(vals<=previous[nu]+1e-9),(previous[nu],vals)
            previous[nu]=vals
            result={'num_cubes':cubes,'energy_cap':cap,'dimension':op.size,
                'nu_over_kappa':nu,'E0':float(vals[0]),'E1':float(vals[1]),
                'gap':float(vals[1]-vals[0]),'solver':info}
            # This residual includes transitions into the larger retained basis,
            # when one is available. It is not a full omitted-space bound.
            if cap<max_cap:
                embedded=np.zeros((len(states),2));embedded[mask]=vecs
                full=SymmetricOperator(electric+model.face_count*nu,rows,cols,-nu*data)
                result['residual_in_largest_test_basis']=[float(np.linalg.norm(
                    full.matvec(embedded[:,i])-vals[i]*embedded[:,i])) for i in range(2)]
            else:
                vals2,_,info2=smallest_ritz(op,levels=2,max_iterations=240,
                    tolerance=1e-10,seed=8917)
                assert info2['converged'] and np.max(np.abs(vals-vals2))<1e-9
                result['second_seed_max_energy_difference']=float(np.max(np.abs(vals-vals2)))
                result['second_seed_solver']=info2
            spectra.append(result)
            print(f'N={cubes} cap={cap} nu={nu}: n={op.size}; E0={vals[0]:.12f}; E1={vals[1]:.12f}; gap={result["gap"]:.12f}',flush=True)
    output={'model':graph_metadata(model),'hamiltonian':'H = kappa sum_physical_links j(j+1) + nu sum_faces (1 - Tr(U_face)/2)',
        'units':'kappa=1','cutoff':'Complete total physical electric energy <= cap; all compatible virtual coupling labels included',
        'status':'Numerical finite-basis spectra. No untruncated or infinite-volume gap certificate is asserted.',
        'max_cap':max_cap,'hermiticity_error':hermiticity,'spectra':spectra,
        'elapsed_seconds':time.perf_counter()-start}
    path=OUT/f'SU2_Joined_Cubes_N{cubes}_Spectra.json'
    path.write_text(json.dumps(output,indent=2)+'\n',encoding='utf-8')
    print(path,flush=True)
    return output


if __name__=='__main__':
    parser=argparse.ArgumentParser()
    parser.add_argument('--cubes',type=int,required=True)
    parser.add_argument('--max-cap',type=int,required=True)
    parser.add_argument('--caps',type=int,nargs='+')
    args=parser.parse_args()
    caps=args.caps or [cap for cap in (6,8,10,12,14,16,18) if cap<=args.max_cap]
    if args.max_cap not in caps:caps.append(args.max_cap)
    study(args.cubes,args.max_cap,sorted(set(caps)))
