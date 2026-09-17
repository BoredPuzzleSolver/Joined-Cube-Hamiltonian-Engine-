"""Exact Haar contractions for the J=1/2 physical SU(2) cube.

Independent of Wigner symbols: each spin configuration is a disjoint union
of loops and its normalized wavefunction is the product of fundamental
holonomy traces. Haar integration uses only delta and epsilon contractions.
"""
from pathlib import Path
from fractions import Fraction
import importlib.util
import itertools
import json
import math
import numpy as np

BASE=Path(__file__).resolve().parent
OUT=BASE.parent/'outputs' if BASE.name=='work' else BASE
EDGES=sorted((u,v) for u in range(8) for v in range(u+1,8) if u^v in (1,2,4))
FACES=[(0,1,3,2),(4,5,7,6),(0,1,5,4),(2,3,7,6),(0,2,6,4),(1,3,7,5)]
EDGE_IDS={edge:i for i,edge in enumerate(EDGES)}
INCIDENT=[[i for i,edge in enumerate(EDGES) if vertex in edge] for vertex in range(8)]
FACE_MASKS=[sum(1<<EDGE_IDS[tuple(sorted((face[k],face[(k+1)%4])))] for k in range(4)) for face in FACES]


def cycles(mask):
    adjacency={vertex:[] for vertex in range(8)}
    for i,(u,v) in enumerate(EDGES):
        if mask>>i&1:
            adjacency[u].append(v);adjacency[v].append(u)
    assert all(len(neighbors) in (0,2) for neighbors in adjacency.values())
    unused={vertex for vertex,neighbors in adjacency.items() if neighbors}
    loops=[]
    while unused:
        start=min(unused);previous=None;current=start;loop=[]
        while True:
            loop.append(current);unused.remove(current)
            neighbors=sorted(adjacency[current])
            following=neighbors[0] if neighbors[0]!=previous else neighbors[1]
            previous,current=current,following
            if current==start:
                break
        loops.append(tuple(loop))
    return loops


def haar_loop_integral(loop_factors,coefficient=Fraction(1)):
    """Exactly integrate a product of loops when each used edge occurs twice.

    For a canonical edge U_ab, a reversed path contributes conjugate(U_ab)
    after indices are attached to the canonical endpoints. Mixed conjugation
    gives delta_ac delta_bd/2. Equal conjugation gives epsilon_ac epsilon_bd/2,
    where epsilon_01=1 and epsilon_10=-1.
    """
    occurrences={}
    variable_count=0
    for loop in loop_factors:
        variables={vertex:variable_count+i for i,vertex in enumerate(loop)}
        variable_count+=len(loop)
        for k,vertex in enumerate(loop):
            following=loop[(k+1)%len(loop)]
            u,v=sorted((vertex,following))
            occurrences.setdefault((u,v),[]).append((variables[u],variables[v],vertex!=u))
    parent=list(range(variable_count));offset=[0]*variable_count
    parity_coefficients=[0]*variable_count
    def find(i):
        if parent[i]!=i:
            old=parent[i]
            parent[i],delta=find(old)
            offset[i]^=delta
        return parent[i],offset[i]
    def impose(i,j,parity):
        ri,pi=find(i);rj,pj=find(j)
        if ri==rj:
            return pi^pj==parity
        parent[ri]=rj;offset[ri]=pi^pj^parity
        return True
    for factors in occurrences.values():
        assert len(factors)==2, 'This exact contraction rule expects zero or two factors per edge.'
        a,b,conjugated=factors[0];c,d,other_conjugated=factors[1]
        parity=int(conjugated==other_conjugated)
        if not impose(a,c,parity) or not impose(b,d,parity):
            return Fraction(0)
        if parity:
            parity_coefficients[a]^=1
            parity_coefficients[b]^=1
    free_coefficients={};constant=0
    for variable in range(variable_count):
        root,parity=find(variable)
        free_coefficients.setdefault(root,0)
        if parity_coefficients[variable]:
            free_coefficients[root]^=1;constant^=parity
    if any(free_coefficients.values()):
        return Fraction(0)
    return coefficient*Fraction((-1)**constant*2**len(free_coefficients),2**len(occurrences))


def count_spin_one_basis():
    # Independent direct enumeration of all 3^12 link assignments, in chunks.
    count=0
    for first_six in itertools.product(range(3),repeat=6):
        tail=np.asarray(list(itertools.product(range(3),repeat=6)),dtype=np.int8)
        states=np.concatenate((np.tile(first_six,(len(tail),1)),tail),axis=1)
        valid=np.ones(len(states),dtype=bool)
        for edges in INCIDENT:
            local=states[:,edges]
            total=local.sum(axis=1)
            valid&=(total%2==0)&(2*local.max(axis=1)<=total)
        count+=int(np.count_nonzero(valid))
    return count


def main():
    masks=[mask for mask in range(1<<len(EDGES)) if all(sum(mask>>edge&1 for edge in incident)%2==0 for incident in INCIDENT)]
    assert len(masks)==32
    lookup={mask:i for i,mask in enumerate(masks)}
    configurations={mask:cycles(mask) for mask in masks}
    for loops in configurations.values():
        assert haar_loop_integral(loops+loops)==1
    wilson=[];transitions=[]
    for face_index,(face,face_mask) in enumerate(zip(FACES,FACE_MASKS)):
        matrix=[[Fraction(0) for _ in masks] for _ in masks]
        for initial in masks:
            final=initial^face_mask
            assert final in lookup
            value=haar_loop_integral(configurations[initial]+configurations[final]+[face],Fraction(1,2))
            matrix[lookup[final]][lookup[initial]]=value
            spectator_edges=set(edge for vertex in face for edge in INCIDENT[vertex] if not(face_mask>>edge&1))
            spectator_count=sum(initial>>edge&1 for edge in spectator_edges)
            assert spectator_count%2==0
            assert abs(value)==Fraction(1,2**(1+spectator_count//2))
            transitions.append({'face_index':face_index,'initial_mask':initial,'final_mask':final,
                                'occupied_spectator_edges':spectator_count,'coefficient_exact':str(value)})
        assert all(matrix[i][j]==matrix[j][i] for i in range(32) for j in range(32))
        assert matrix[lookup[face_mask]][lookup[0]]==Fraction(1,2)
        wilson.append(matrix)
    output={'title':'Independent exact Haar validation of physical SU(2) cube at J=1/2',
            'vertices':'0,...,7 with vertex=x+2y+4z','edges':EDGES,'faces':FACES,'face_masks':FACE_MASKS,
            'basis':'Eulerian occupied-edge masks, sorted ascending. Every occupied vertex has degree2, hence a unique invariant and disjoint closed loops.',
            'basis_masks':masks,'basis_loop_cycles':{str(mask):configurations[mask] for mask in masks},
            'basis_normalization':'f_C is the product of unnormalized fundamental traces around its cycles. Each f_C has exact Haar norm1. Distinct masks are orthogonal because some edge carries a single fundamental factor.',
            'haar_method':'Every nonzero Wilson matrix element has exactly two fundamental factors per used edge. Integrate with delta/epsilon identities and count consistent binary-index assignments, retaining all epsilon signs. No6j or Clebsch-Gordan tables are used.',
            'W_face_normalization':'W_face=(1/2)Trfund holonomy; each nonzero transition toggles the four face edges.',
            'wilson_matrices_exact':[[[str(value) for value in row] for row in matrix] for matrix in wilson],
            'transitions':transitions,
            'checks':{'spin_half_basis_dimension':len(masks),'spin_one_basis_dimension':count_spin_one_basis(),
                      'exact_basis_norms_all_one':True,'all_six_matrices_exactly_symmetric':True,
                      'absolute_coefficient_spectator_rule':'2^(-1-number_of_occupied_spectator_edges/2)',
                      'negative_directed_matrix_entries':sum(row['coefficient_exact'].startswith('-') for row in transitions),
                      'nonzero_directed_matrix_entries':len(transitions),
                      'distinct_nonzero_coefficients':sorted(set(row['coefficient_exact'] for row in transitions))},
            'all_assertions_passed':True}
    candidate_path=BASE/'su2_cube_model.py'
    if candidate_path.exists():
        spec=importlib.util.spec_from_file_location('cube_candidate',candidate_path)
        candidate=importlib.util.module_from_spec(spec)
        spec.loader.exec_module(candidate)
        candidate_states,candidate_sparse=candidate.cube_wilson_transitions(1)
        candidate_masks=[sum(value<<edge for edge,value in enumerate(state)) for state in candidate_states]
        assert set(candidate_masks)==set(masks)
        candidates=[]
        sign_relations=[]
        for face_index,entries in enumerate(candidate_sparse):
            entries_by_masks={}
            for row,column,sign,square in entries:
                numerator=math.isqrt(square.numerator);denominator=math.isqrt(square.denominator)
                assert numerator*numerator==square.numerator and denominator*denominator==square.denominator
                value=sign*Fraction(numerator,denominator)
                final,initial=candidate_masks[row],candidate_masks[column]
                entries_by_masks[final,initial]=value
            assert len(entries_by_masks)==32
            matrix=[[Fraction(0) for _ in masks] for _ in masks]
            for (final,initial),value in entries_by_masks.items():
                independent=wilson[face_index][lookup[final]][lookup[initial]]
                assert independent and abs(value)==abs(independent)
                ratio=value/independent
                assert ratio in (-1,1)
                sign_relations.append((initial,final,int(ratio)))
                matrix[lookup[final]][lookup[initial]]=value
            candidates.append(matrix)
        phases={0:1}
        changed=True
        while changed:
            changed=False
            for initial,final,ratio in sign_relations:
                if initial in phases:
                    desired=ratio*phases[initial]
                    if final in phases:
                        assert phases[final]==desired, 'Candidate phases cannot be reconciled by one state-basis rephasing.'
                    else:
                        phases[final]=desired;changed=True
        assert len(phases)==32
        for face_index in range(6):
            for i,final in enumerate(masks):
                for j,initial in enumerate(masks):
                    assert candidates[face_index][i][j]==phases[final]*phases[initial]*wilson[face_index][i][j]
        assert len(candidate.cube_basis(2))==output['checks']['spin_one_basis_dimension']
        output['recoupling_comparison']={
            'method':'Compare candidate exact signed squared magnitudes with independently completed Haar contractions; allow one common diagonal state-basis sign change for all six faces.',
            'all_192_directed_coefficients_match_exactly_after_common_rephasing':True,
            'negative_state_rephasing_count':sum(value<0 for value in phases.values()),
            'state_rephasing_by_mask':{str(mask):phases[mask] for mask in masks},
            'spin_one_basis_count_matches_independent_enumeration':True}
    OUT.mkdir(parents=True,exist_ok=True)
    path=OUT/'SU2_Cube_Independent_Validation.json'
    path.write_text(json.dumps(output,indent=2)+'\n',encoding='utf-8')
    print(json.dumps(output['checks'],indent=2))
    print(path)


if __name__=='__main__':
    main()
