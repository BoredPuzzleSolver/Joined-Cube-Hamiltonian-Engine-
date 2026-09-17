"""Independent rounding and exact-inertia audit of the cube certificate.

Rebuilds rounded blocks from explicit face-label target enumeration. Radical
rounding uses integer binary search, not isqrt. Inertia uses ordinary exact
Fraction Schur elimination of the whole matrix, with omitted indices first.
"""
from fractions import Fraction as F
from functools import lru_cache
from pathlib import Path
import itertools
import json
from su2_cube_model import EDGES, INCIDENT, FACE_DATA, cube_basis, energy_basis, face_amplitude_exact
from su2_recoupling import admissible
import su2_cube_rational_certificate as candidate

ROOT=Path(__file__).resolve().parent.parent


@lru_cache(maxsize=None)
def binary_radical(square,scale):
    numerator=square.numerator*scale*scale
    denominator=square.denominator
    low=0;high=1
    while high*high*denominator<=numerator:
        high*=2
    while high-low>1:
        middle=(low+high)//2
        if middle*middle*denominator<=numerator:
            low=middle
        else:
            high=middle
    value=F(low,scale)
    error=F(0) if low*low*denominator==numerator else F(1,scale)
    assert value*value<=square< F(low+1,scale)**2
    return value,error


def independent_inertia(matrix,index_order):
    a=[{j:F(value) for j,value in enumerate(row) if value} for row in matrix]
    active=list(index_order)
    positive=negative=zero=0
    while active:
        k=next((i for i in active if a[i].get(i,0)),None)
        if k is None:
            assert not any(a[i].get(j,0) for i in active for j in active)
            zero+=len(active)
            break
        pivot=a[k][k]
        positive+=int(pivot>0);negative+=int(pivot<0)
        neighbors=[i for i in active if i!=k and a[i].get(k,0)]
        column={i:a[i][k] for i in neighbors}
        for offset,i in enumerate(neighbors):
            ratio=column[i]/pivot
            for j in neighbors[offset:]:
                value=a[i].get(j,F(0))-ratio*column[j]
                if value:
                    a[i][j]=a[j][i]=value
                else:
                    a[i].pop(j,None);a[j].pop(i,None)
        for i in neighbors:
            a[i].pop(k,None)
        active.remove(k);a[k]={}
    return {'negative':negative,'zero':zero,'positive':positive}


def independent_energy_basis(cap):
    """Fixed-edge-order enumeration without the candidate's vertex-cost pruning."""
    cap=F(cap)
    budget=(4*cap).numerator//(4*cap).denominator
    state=[-1]*len(EDGES)
    answers=[]
    def visit(edge,remaining):
        if edge==len(EDGES):
            answers.append(tuple(state))
            return
        value=0
        while value*(value+2)<=remaining:
            state[edge]=value
            valid=True
            for vertex in EDGES[edge]:
                labels=[state[e] for e in INCIDENT[vertex]]
                if min(labels)>=0:
                    total=sum(labels)
                    if total%2 or 2*max(labels)>total:
                        valid=False
                        break
            if valid:
                visit(edge+1,remaining-value*(value+2))
            value+=1
        state[edge]=-1
    visit(0,budget)
    return sorted(answers)


def main():
    saved=json.loads((ROOT/'outputs/SU2_Cube_Gap_Certificate.json').read_text())
    scale=saved['rounding_scale']
    retained=cube_basis(1)
    assert len(retained)==32
    p=len(retained)
    all_targets=set()
    entries=[]
    for initial in retained:
        for face in FACE_DATA:
            edges,_=face
            for choices in itertools.product((-1,1),repeat=4):
                target=list(initial)
                for edge,change in zip(edges,choices):
                    target[edge]+=change
                target=tuple(target)
                if min(target)<0 or not all(admissible(*(target[edge] for edge in incident)) for incident in INCIDENT):
                    continue
                sign,square=face_amplitude_exact(initial,target,face)
                if square:
                    assert sign in (-1,1)
                    all_targets.add(target)
                    value,error=binary_radical(square,scale)
                    entries.append((initial,target,-sign*value,error))
    states=retained+sorted(all_targets-set(retained))
    assert len(states)==212
    lookup={state:i for i,state in enumerate(states)}
    electric=[F(sum(value*(value+2) for value in state),4) for state in states]
    diagonal=[energy+(6 if i<p else 0) for i,energy in enumerate(electric)]
    lower=[[F(0) for _ in states] for _ in states]
    errors=[[F(0) for _ in states] for _ in states]
    observed={}
    for initial,target,value,error in entries:
        i,j=lookup[initial],lookup[target]
        key=tuple(sorted((i,j)))
        if key in observed:
            assert observed[key]==(value,error)
        observed[key]=(value,error)
        lower[i][j]=lower[j][i]=value
        errors[i][j]=errors[j][i]=error
    for i in range(len(states)):
        lower[i][i]=diagonal[i]
    assert all(lower[i][j]==0 for i in range(p,len(states)) for j in range(p,len(states)) if i!=j)
    eta_l=max(sum(row) for row in errors)
    eta_a=max(sum(row[:p]) for row in errors[:p])
    assert str(eta_l)==saved['L_operator_error_bound_exact']
    assert str(eta_a)==saved['A_operator_error_bound_exact']
    proposed_states,proposed_diagonal,proposed_matrix,proposed_eta_a,proposed_eta_l=candidate.comparison_matrices(retained,scale=scale)
    assert states==proposed_states and diagonal==proposed_diagonal
    assert eta_a==proposed_eta_a and eta_l==proposed_eta_l
    for i in range(len(states)):
        for j in range(len(states)):
            if i!=j:
                assert lower[i][j]==proposed_matrix[i].get(j,F(0))
    probes=[]
    for row in saved['energy_intervals']:
        lo,up=F(row['lower_exact']),F(row['upper_exact'])
        lower_probe=lo+eta_l;upper_probe=up-eta_a
        assert lower_probe==F(row['rounded_lower_probe_exact'])
        assert upper_probe==F(row['rounded_upper_probe_exact'])
        assert min(electric[p:])>lower_probe
        full_shifted=[line.copy() for line in lower]
        for i in range(len(states)):
            full_shifted[i][i]-=lower_probe
        low_inertia=independent_inertia(full_shifted,list(range(p,len(states)))+list(range(p)))
        upper_shifted=[[lower[i][j]-(upper_probe if i==j else 0) for j in range(p)] for i in range(p)]
        up_inertia=independent_inertia(upper_shifted,range(p))
        assert low_inertia==row['lower_auxiliary_inertia']
        assert up_inertia==row['upper_projection_inertia']
        j=row['index']
        assert low_inertia['negative']<=j and up_inertia['negative']>=j+1
        assert lo<=F(saved['remaining_electric_threshold_exact'])
        probes.append({'level':j,'lower_inertia':low_inertia,'upper_inertia':up_inertia})
    intervals=saved['energy_intervals']
    gap=[F(intervals[1]['lower_exact'])-F(intervals[0]['upper_exact']),
         F(intervals[1]['upper_exact'])-F(intervals[0]['lower_exact'])]
    assert list(map(str,gap))==saved['gap_interval_exact']
    result={'retained_dimension':p,'direct_omitted_dimension':len(states)-p,
            'exact_radical_intervals_match_independent_binary_search':True,
            'complete_one_face_target_enumeration_matches':True,
            'A_norm_error_exact':str(eta_a),'L_norm_error_exact':str(eta_l),
            'minimum_direct_omitted_electric_energy':str(min(electric[p:])),
            'four_exact_probe_inertias_match':probes,
            'gap_interval_exact':list(map(str,gap)),
            'all_checks_passed':True}
    energy_basis_checks=[]
    for cap in [8,10,12]:
        independently_enumerated=independent_energy_basis(cap)
        proposed_basis=energy_basis(F(cap))
        assert independently_enumerated==proposed_basis
        energy_basis_checks.append({'cap':cap,'dimension':len(proposed_basis),
                                    'full_sorted_basis_matches_independent_fixed_edge_enumeration':True})
    result['energy_basis_generalization_checks']=energy_basis_checks
    (ROOT/'outputs/SU2_Cube_Certificate_Independent_Audit.json').write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8')
    print(json.dumps(result,indent=2))


if __name__=='__main__':
    main()
