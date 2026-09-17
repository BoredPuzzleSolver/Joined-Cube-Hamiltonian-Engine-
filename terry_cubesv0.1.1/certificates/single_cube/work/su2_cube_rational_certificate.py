"""Exact radical rounding and Schur inertia for the complete SU(2) cube.

The retained 32-state basis and its 180 directly coupled omitted states give
a finite lower comparison. Every other omitted state has electric energy >=
6.5. Integer square-root brackets control all radical rounding rigorously.
"""
from fractions import Fraction as F
from pathlib import Path
import json
import math
import numpy as np
from su2_cube_model import cube_basis, FACE_DATA, state_face_transitions
from su2_cutoff_bounds import rational_inertia

ROOT = Path(__file__).resolve().parent.parent


def electric(state):
    return F(sum(x*(x+2) for x in state),4)


def rounded_radical(sign, square, scale):
    """Return a rational signed approximation and an exact absolute error cap."""
    square = F(square)
    assert square >= 0 and sign in (-1,0,1) and scale>0
    numerator = math.isqrt(square.numerator*scale*scale//square.denominator)
    lower = F(numerator,scale)
    assert lower*lower <= square < F(numerator+1,scale)**2
    error = F(0) if lower*lower==square else F(1,scale)
    return sign*lower,error


def comparison_matrices(retained, nu=F(1), scale=10**12):
    """Build A_hat and L_hat exactly and rigorous rounding norm caps.

    The physical A and L each differ from the rounded matrix by a symmetric
    error bounded in operator norm by the largest absolute row-error sum.
    """
    p = len(retained)
    lookup = {s:i for i,s in enumerate(retained)}
    targets = set()
    source_entries = []
    for col,state in enumerate(retained):
        for face in FACE_DATA:
            for target,sign,square in state_face_transitions(state,face):
                targets.add(target)
                source_entries.append((target,col,-sign,nu*nu*square))
    states = list(retained)+sorted(targets-set(retained))
    lookup = {s:i for i,s in enumerate(states)}
    values = {}
    errors = {}
    for target,col,sign,square in source_entries:
        row = lookup[target]
        key = tuple(sorted((row,col)))
        value,error = rounded_radical(sign,square,scale)
        if key in values:
            assert values[key]==value and errors[key]==error
        else:
            values[key],errors[key] = value,error
    diagonal = [electric(s)+(6*nu if i<p else 0) for i,s in enumerate(states)]
    row_errors = [F(0) for _ in states]
    a_errors = [F(0) for _ in retained]
    matrix = [dict() for _ in states]
    for (i,j),value in values.items():
        assert i != j
        matrix[i][j]=matrix[j][i]=value
        row_errors[i]+=errors[(i,j)]
        row_errors[j]+=errors[(i,j)]
        if i<p and j<p:
            a_errors[i]+=errors[(i,j)]
            a_errors[j]+=errors[(i,j)]
    return states,diagonal,matrix,max(a_errors),max(row_errors)


def inertia_a(diagonal,matrix,p,probe,counter=rational_inertia):
    out = [[matrix[i].get(j,F(0)) if i!=j else diagonal[i]-probe
            for j in range(p)] for i in range(p)]
    return counter(out)


def inertia_lower(diagonal,matrix,p,probe,counter=rational_inertia):
    """Eliminate the exactly diagonal omitted block; count its positive signs."""
    out = [[matrix[i].get(j,F(0)) if i!=j else diagonal[i]-probe
            for j in range(p)] for i in range(p)]
    for q in range(p,len(diagonal)):
        denominator = diagonal[q]-probe
        assert denominator>0, 'This certificate needs probes below the omitted electric block'
        neighbors = list(matrix[q])
        assert all(i<p for i in neighbors)
        for offset,i in enumerate(neighbors):
            for j in neighbors[offset:]:
                correction = matrix[q][i]*matrix[q][j]/denominator
                out[i][j]-=correction
                if i!=j:
                    out[j][i]-=correction
    result = counter(out)
    result['positive']+=len(diagonal)-p
    if 'dimension' in result:
        result['schur_dimension']=result['dimension']
        result['dimension']=len(diagonal)
    return result


def main():
    retained = cube_basis(1)
    states,diagonal,matrix,eta_a,eta_l = comparison_matrices(retained)
    p = len(retained)
    assert p==32 and len(states)==212
    # These fixed, exact decimal rationals are checked, never accepted on the
    # basis of a floating-point eigensolver tolerance.
    lower = [F('5.30400030'),F('6.20980338')]
    upper = [F('5.54178379'),F('8.77237429')]
    t_rest = F(13,2)
    rows = []
    for j,(lo,up) in enumerate(zip(lower,upper)):
        low_inertia = inertia_lower(diagonal,matrix,p,lo+eta_l)
        high_inertia = inertia_a(diagonal,matrix,p,up-eta_a)
        assert low_inertia['negative']<=j and lo<=t_rest
        assert high_inertia['negative']>=j+1
        rows.append({'index':j,'lower_exact':str(lo),'upper_exact':str(up),
                     'lower_decimal':float(lo),'upper_decimal':float(up),
                     'rounded_lower_probe_exact':str(lo+eta_l),
                     'rounded_upper_probe_exact':str(up-eta_a),
                     'lower_auxiliary_inertia':low_inertia,
                     'upper_projection_inertia':high_inertia,
                     'certified':True})
    gap = [lower[1]-upper[0],upper[1]-lower[0]]
    dense = np.diag([float(x) for x in diagonal])
    for i,row in enumerate(matrix):
        for j,value in row.items():
            dense[i,j]=float(value)
    result = {
        'model':'Open single SU(2) cube, Gauss law all vertices, kappa=nu=1',
        'scope':'Untruncated physical Hilbert space on the fixed cube, conditional on stated basis and recoupling identities',
        'retained_dimension':p,'reachable_auxiliary_dimension':len(states),
        'remaining_electric_threshold_exact':str(t_rest),
        'threshold_proof':'See su2_cube_bounds.py and accompanying bound note; every state outside link-spin cutoff 1/2 has electric energy at least 13/2',
        'matrix_comparison':'H >= L direct-sum T_rest; H minus comparison is QVQ >= 0',
        'rounding_scale':10**12,'A_operator_error_bound_exact':str(eta_a),
        'L_operator_error_bound_exact':str(eta_l),
        'rounding_method':'Integer-isqrt brackets of exact rational squared amplitudes; symmetric maximum row-error sum bounds spectral norm',
        'inertia_method':'Exact rational Schur complement eliminates positive diagonal omitted block; 32-by-32 integer Bareiss inertia',
        'energy_intervals':rows,'gap_interval_exact':[str(x) for x in gap],
        'gap_interval_decimal':[float(x) for x in gap],'all_certified':True,
        'numerical_comparison_only':{'first_L_levels':np.linalg.eigvalsh(dense)[:4].tolist(),
                                     'first_A_levels':np.linalg.eigvalsh(dense[:p,:p])[:4].tolist()},
    }
    output=ROOT/'outputs'/'SU2_Cube_Gap_Certificate.json'
    output.write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8')
    print(json.dumps({'gap_interval':result['gap_interval_decimal'],'A_error':str(eta_a),'L_error':str(eta_l),'all_certified':True}),flush=True)


if __name__ == '__main__':
    main()
