"""Stronger complete-cube bounds from a complete total-electric-energy basis."""
from fractions import Fraction as F
from pathlib import Path
import json
import math
import time
import numpy as np
from su2_cube_model import energy_basis
from su2_cube_rational_certificate import comparison_matrices,inertia_a,inertia_lower
from su2_sparse_spectrum import SymmetricOperator,smallest_ritz
from su2_interval_inertia import interval_inertia

ROOT=Path(__file__).resolve().parent.parent


def main():
    results=[]
    for cap in [8,10,12]:
        start=time.perf_counter()
        retained=energy_basis(F(cap))
        states,diagonal,matrix,eta_a,eta_l=comparison_matrices(retained)
        p=len(retained)
        rr=[];cc=[];dd=[]
        for i,row in enumerate(matrix):
            for j,value in row.items():
                rr.append(i);cc.append(j);dd.append(float(value))
        op=SymmetricOperator([float(x) for x in diagonal],rr,cc,dd)
        values,_,solver=smallest_ritz(op,levels=2,max_iterations=240,tolerance=1e-10)
        assert solver['converged']
        a=np.diag([float(x) for x in diagonal[:p]])
        for i in range(p):
            for j,value in matrix[i].items():
                if j<p:a[i,j]=float(value)
        upper_values=np.linalg.eigvalsh(a)[:2]
        threshold=F(4*cap+1,4) # Every electric energy is a quarter-integer.
        # Rounding and solver error only suggest probes. Exact inertia below
        # decides validity; no floating-point error bound is presumed.
        lower=[F(math.floor(min(float(x)-float(eta_l),float(threshold))*10**6)-1,10**6) for x in values]
        upper=[F(math.ceil((float(x)+float(eta_a))*10**6)+1,10**6) for x in upper_values]
        print(f'cap{cap}: P{p}, L{len(states)}, proposed gap[{float(lower[1]-upper[0]):.8f},{float(upper[1]-lower[0]):.8f}]',flush=True)
        rows=[]
        for j,(lo,up) in enumerate(zip(lower,upper)):
            low=inertia_lower(diagonal,matrix,p,lo+eta_l,counter=interval_inertia)
            high=inertia_a(diagonal,matrix,p,up-eta_a,counter=interval_inertia)
            assert low['negative']<=j and lo<=threshold
            assert high['negative']>=j+1
            rows.append({'index':j,'lower_exact':str(lo),'upper_exact':str(up),
                         'lower_decimal':float(lo),'upper_decimal':float(up),
                         'lower_shift_exact':str(lo+eta_l),'upper_shift_exact':str(up-eta_a),
                         'lower_inertia':low,'upper_inertia':high,'certified':True})
            print(f'cap{cap} level{j}: both exact probes passed after {time.perf_counter()-start:.2f}s',flush=True)
        gap=[lower[1]-upper[0],upper[1]-lower[0]]
        results.append({'electric_energy_cap_over_kappa':cap,'retained_dimension':p,
                        'auxiliary_dimension':len(states),'remaining_electric_threshold_exact':str(threshold),
                        'A_error_bound_exact':str(eta_a),'L_error_bound_exact':str(eta_l),
                        'energy_intervals':rows,'gap_interval_exact':[str(x) for x in gap],
                        'gap_interval_decimal':[float(x) for x in gap],'all_certified':True,
                        'numeric_L_eigenvalues':values.tolist(),'numeric_A_eigenvalues':upper_values.tolist(),
                        'elapsed_seconds':time.perf_counter()-start})
    output={'model':'Complete fixed open SU(2) cube, kappa=nu=1',
            'method':'Exact radical-rounding/operator-norm/Schur comparison as su2_cube_rational_certificate.py; pivot signs certified by outward-rounded integer fixed-point interval elimination at 128 bits',
            'basis_completeness':'energy_basis includes all gauge-admissible states with T<=cap, using a sufficient link bound and exact energy pruning',
            'remaining_threshold_proof':'T is a multiple of 1/4; all states outside the complete integer energy cap have T>=cap+1/4',
            'certificates':results}
    path=ROOT/'outputs'/'SU2_Cube_Energy_Certificates.json'
    path.write_text(json.dumps(output,indent=2)+'\n',encoding='utf-8')
    print(path,flush=True)


if __name__=='__main__':
    main()
