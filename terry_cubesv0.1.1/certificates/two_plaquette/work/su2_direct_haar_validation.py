"""Independent Haar-quadrature validation for two SU(2) loop holonomies.

No Wigner 6j symbols or recoupling formulas are used here. Spin projectors come
from independently diagonalized total-angular-momentum Casimir matrices.
Run with Python + NumPy. Doubled integer spin labels avoid half-integer tests.
"""
from pathlib import Path
import json
import math
import importlib.util
import numpy as np

BASE=Path(__file__).resolve().parent
OUT=BASE.parent/'outputs' if BASE.name=='work' else BASE


def angular_momentum(twice_spin):
    j=twice_spin/2
    d=twice_spin+1
    m=np.arange(j,-j-1,-1,dtype=float)
    plus=np.zeros((d,d),dtype=complex)
    for column in range(1,d):
        mm=m[column]
        plus[column-1,column]=math.sqrt(j*(j+1)-mm*(mm+1))
    minus=plus.conj().T
    matrices=((plus+minus)/2,(plus-minus)/(2j),np.diag(m).astype(complex))
    assert np.max(abs(matrices[0]@matrices[1]-matrices[1]@matrices[0]-1j*matrices[2]))<1e-12
    return matrices


def basis_labels(cutoff):
    return [(left,right,middle) for left in range(cutoff+1)
            for right in range(cutoff+1)
            for middle in range(abs(left-right),min(left+right,cutoff)+1,2)]


def casimir_projectors(left,right):
    jl=angular_momentum(left)
    jr=angular_momentum(right)
    dl,dr=left+1,right+1
    total=[np.kron(jl[k],np.eye(dr))+np.kron(np.eye(dl),jr[k]) for k in range(3)]
    casimir=sum(j@j for j in total)
    values,vectors=np.linalg.eigh(casimir)
    projectors={}
    max_error=0.
    for middle in range(abs(left-right),left+right+1,2):
        j=middle/2
        chosen=abs(values-j*(j+1))<1e-9
        assert np.count_nonzero(chosen)==middle+1
        vv=vectors[:,chosen]
        projector=vv@vv.conj().T
        error=max(np.max(abs(projector@projector-projector)),
                  np.max(abs(casimir@projector-j*(j+1)*projector)))
        max_error=max(max_error,float(error))
        projectors[middle]=projector
    assert np.max(abs(sum(projectors.values())-np.eye(dl*dr)))<1e-12
    return projectors,max_error


def class_angle_rule(order):
    k=np.arange(1,order+1)
    alpha=k*np.pi/(order+1)
    # Normalized SU(2) class measure: (2/pi) sin(alpha)^2 d alpha.
    return np.cos(alpha),2*np.sin(alpha)**2/(order+1)


def build_spin_functions(cutoff,order):
    labels=basis_labels(cutoff)
    x,wxy=class_angle_rule(order)
    u,wu=np.polynomial.legendre.leggauss(order)
    wu=wu/2
    jmat={spin:angular_momentum(spin) for spin in range(cutoff+1)}
    alpha=np.arccos(x)
    a_diagonals={spin:np.exp(2j*alpha[:,None]*np.diag(jmat[spin][2]).real[None,:])
                 for spin in range(cutoff+1)}
    b_matrices={}
    for spin in range(cutoff+1):
        d=spin+1
        matrices=np.zeros((order,order,d,d),dtype=complex)
        for index,dot in enumerate(u):
            generator=math.sqrt(max(0.,1-dot*dot))*jmat[spin][0]+dot*jmat[spin][2]
            vals,vecs=np.linalg.eigh(generator)
            phases=np.exp(2j*alpha[:,None]*vals[None,:])
            matrices[:,index,:,:]=np.einsum('ak,xk,bk->xab',vecs,phases,vecs.conj())
        b_matrices[spin]=matrices
    functions=[]
    projector_errors=[]
    cache={}
    for left,right,middle in labels:
        if (left,right) not in cache:
            projectors,error=casimir_projectors(left,right)
            cache[left,right]=projectors
            projector_errors.append(error)
        projector=cache[left,right][middle]
        dl,dr=left+1,right+1
        p4=projector.reshape(dl,dr,dl,dr)
        # Tr[P (D(A) tensor D(B))], with D(A) diagonal.
        diagonal_blocks=np.asarray([p4[a,:,a,:] for a in range(dl)])
        values=math.sqrt(dl*dr/(middle+1))*np.einsum(
            'xa,abc,yucb->xyu',a_diagonals[left],diagonal_blocks,b_matrices[right])
        functions.append(values.ravel())
    functions=np.asarray(functions).T
    weights=np.einsum('x,y,u->xyu',wxy,wxy,wu).ravel()
    multiplier_a=np.broadcast_to(x[:,None,None],(order,order,order)).ravel()
    multiplier_b=np.broadcast_to(x[None,:,None],(order,order,order)).ravel()
    return labels,functions,weights,multiplier_a,multiplier_b,max(projector_errors)


def quadrature_result(cutoff,order):
    labels,f,w,a,b,projector_error=build_spin_functions(cutoff,order)
    gram=f.conj().T@(w[:,None]*f)
    wa=f.conj().T@((w*a)[:,None]*f)
    wb=f.conj().T@((w*b)[:,None]*f)
    errors={'normalization_error':float(abs(sum(w)-1)),
            'gram_max_error':float(np.max(abs(gram-np.eye(len(labels))))),
            'imaginary_spin_function_max':float(np.max(abs(f.imag))),
            'matrix_imaginary_max':float(max(np.max(abs(wa.imag)),np.max(abs(wb.imag)))),
            'casimir_projector_max_error':projector_error,
            'hermiticity_max_error':float(max(np.max(abs(wa-wa.conj().T)),np.max(abs(wb-wb.conj().T))))}
    assert errors['gram_max_error']<1e-11
    assert errors['matrix_imaginary_max']<1e-11
    lookup={label:i for i,label in enumerate(labels)}
    permutation=[lookup[(right,left,middle)] for left,right,middle in labels]
    errors['loop_exchange_matrix_error']=float(np.max(abs(wb-wa[np.ix_(permutation,permutation)])))
    assert errors['loop_exchange_matrix_error']<1e-11
    selection_error=0.
    for row,(left,right,middle) in enumerate(labels):
        for column,(lp,rp,mp) in enumerate(labels):
            if not(right==rp and abs(left-lp)==1 and abs(middle-mp)==1):
                selection_error=max(selection_error,float(abs(wa[row,column])))
    errors['forbidden_transition_max']=selection_error
    assert selection_error<1e-11
    vacuum=lookup[(0,0,0)]
    assert abs(wa[vacuum,lookup[(1,0,1)]]-.5)<1e-12
    assert abs(wb[vacuum,lookup[(0,1,1)]]-.5)<1e-12
    return {'twice_spin_cutoff':cutoff,'quadrature_order_per_variable':order,
            'quadrature_points':order**3,'basis_dimension':len(labels),
            'basis_doubled_spins':[list(label) for label in labels],
            'W_A_matrix':wa.real.tolist(),'W_B_matrix':wb.real.tolist(),
            'checks':errors}


def main():
    results=[]
    for cutoff in [2,3]:
        coarse=quadrature_result(cutoff,8)
        fine=quadrature_result(cutoff,12)
        error=max(float(np.max(abs(np.asarray(fine[key])-np.asarray(coarse[key]))))
                  for key in ['W_A_matrix','W_B_matrix'])
        fine['checks']['quadrature_8_vs_12_max_error']=error
        assert error<1e-11
        results.append(fine)
    comparisons=[]
    candidate_path=BASE/'su2_recoupling.py'
    if candidate_path.exists():
        # Only the comparison stage calls the separate proposed formula.
        # All quadrature matrices above have already been independently built.
        specification=importlib.util.spec_from_file_location('recoupling_candidate',candidate_path)
        candidate=importlib.util.module_from_spec(specification)
        specification.loader.exec_module(candidate)
        for row in results:
            labels,wa,wb=candidate.wilson_matrices(row['twice_spin_cutoff'])
            order=[labels.index(tuple(label)) for label in row['basis_doubled_spins']]
            errors={key:float(np.max(abs(np.asarray(row[key])-matrix[np.ix_(order,order)])))
                    for key,matrix in [('W_A_matrix',wa),('W_B_matrix',wb)]}
            assert max(errors.values())<1e-11
            comparisons.append({'twice_spin_cutoff':row['twice_spin_cutoff'],
                                'matrix_dimension':len(labels),'max_absolute_entry_errors':errors})
    output={'title':'Independent direct Haar validation of shared-edge SU(2) loop multiplication matrices',
            'normalized_basis':'f_(jL,jR,jM)(A,B)=sqrt((2jL+1)(2jR+1)/(2jM+1)) Tr[P_jM (D^jL(A) tensor D^jR(B))].',
            'projectors':'P_jM constructed by eigendecomposition of (J_L tensor I + I tensor J_R)^2; no Clebsch-Gordan or 6j implementation is called.',
            'measure':'After simultaneous conjugation, A=cos(alpha)I+i sin(alpha)sigma_z and B=cos(beta)I+i sin(beta)(sqrt(1-u^2)sigma_x+u sigma_z). Independent Haar measure reduces to (2/pi)^2 sqrt(1-x^2)sqrt(1-y^2) dx dy du/2, x=cos(alpha), y=cos(beta), u in [-1,1].',
            'observables':'W_A=Tr_fund(A)/2=x; W_B=Tr_fund(B)/2=y.',
            'method':'Gauss-Chebyshev second-kind rules for x,y and Gauss-Legendre for u, each at orders8 and12; spectra of angular-momentum matrices determine representation exponentials.',
            'independence':'This file deliberately does not calculate Wigner6j symbols. The stored matrix entries can be compared with a separate recoupling construction.',
            'candidate_recoupling_comparison':comparisons,
            'candidate_formula':'W_A(initial,final)=(1/2)sqrt(dL*dM*dLprime*dMprime)*{jL,jR,jM;jMprime,1/2,jLprime}^2, jRprime=jR. W_B exchanges L and R.',
            'comparison_note':'When su2_recoupling.py is present alongside this script, its proposed matrices are compared only after the independent direct Haar matrices are complete.',
            'precision_note':'Floating-point quadrature verified at two orders and against orthonormality/selection identities; not interval-certified numerical arithmetic.',
            'results':results,
            'all_assertions_passed':True}
    OUT.mkdir(parents=True,exist_ok=True)
    path=OUT/'SU2_Direct_Haar_Validation.json'
    path.write_text(json.dumps(output,indent=2)+'\n',encoding='utf-8')
    for result in results:
        print(result['twice_spin_cutoff'],result['basis_dimension'],result['checks'])
    print(path)


if __name__=='__main__':
    main()
