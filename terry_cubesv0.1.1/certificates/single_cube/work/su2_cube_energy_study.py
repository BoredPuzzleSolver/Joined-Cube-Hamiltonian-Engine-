"""A second, complete cutoff route for the SU(2) cube: total electric energy."""
from fractions import Fraction
from pathlib import Path
import json
import math
import time
import numpy as np
from su2_cube_model import energy_basis,graph_wilson_transitions,electric_diagonal
from su2_sparse_spectrum import SymmetricOperator,smallest_ritz

ROOT=Path(__file__).resolve().parent.parent


def main():
    start=time.perf_counter()
    states=energy_basis(Fraction(30))
    faces=graph_wilson_transitions(states)
    rr=[];cc=[];dd=[]
    for face in faces:
        for row,col,sign,square in face:
            rr.append(row);cc.append(col);dd.append(sign*math.sqrt(float(square)))
    rows=np.array(rr,dtype=np.int64);cols=np.array(cc,dtype=np.int64);data=np.array(dd)
    del faces,rr,cc,dd
    electric=electric_diagonal(states)
    print(f'Electric cap30: {len(states)}states, {len(data)}entries, construction {time.perf_counter()-start:.2f}s',flush=True)
    results=[]
    previous={}
    for cap in [12,18,24,30]:
        mask=electric<=cap
        keep=mask[rows]&mask[cols]
        remap=np.full(len(states),-1,dtype=np.int64)
        remap[mask]=np.arange(np.sum(mask))
        for nu in [1,4]:
            op=SymmetricOperator(electric[mask]+6*nu,remap[rows[keep]],remap[cols[keep]],-nu*data[keep])
            values,vectors,check=smallest_ritz(op,levels=2,max_iterations=240,tolerance=1e-10)
            assert check['converged'],check
            if nu in previous:
                assert np.all(values<=previous[nu]+1e-9)
            previous[nu]=values
            row={'electric_energy_cap_over_kappa':cap,'dimension':op.size,
                 'largest_link_spin_present':max(max(s) for s,t in zip(states,electric) if t<=cap)/2,
                 'nu_over_kappa':nu,'E0':float(values[0]),'E1':float(values[1]),
                 'gap':float(values[1]-values[0]),'solver':check}
            if cap==30:
                check2,_,info2=smallest_ritz(op,levels=2,max_iterations=240,tolerance=1e-10,seed=8917)
                assert info2['converged'] and np.max(abs(values-check2))<1e-9
                row['second_seed_max_energy_difference']=float(np.max(abs(values-check2)))
                row['second_seed_solver']=info2
            results.append(row)
            print(f'cap{cap} nu{nu}: {op.size} states, gap={row["gap"]:.12f}',flush=True)
    out={'model':'Same open SU(2) cube and Hamiltonian as su2_cube_study.py',
         'cutoff':'Complete total electric energy T/kappa <= cap, with no additional spin restriction',
         'units':'kappa=1','status':'Numerical variational spectra; differences alone are not gap enclosures',
         'basis_completeness':'Exact energy-pruned enumeration; each doubled label q obeys q(q+2)<=4*cap, yielding a sufficient finite link bound',
         'spectra':results}
    path=ROOT/'outputs'/'SU2_Cube_Energy_Study_Results.json'
    path.write_text(json.dumps(out,indent=2)+'\n',encoding='utf-8')
    print(path)


if __name__=='__main__':
    main()
