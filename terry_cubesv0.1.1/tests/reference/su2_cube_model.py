"""Gauge-invariant SU(2) on an open cube; all spin labels are doubled integers.

Gauss law is imposed at every vertex, with no external charges. Link irreps
and normalized trivalent intertwiners form an orthonormal spin-network basis.
W_f=Tr_fund(U_f)/2; H=kappa sum_l j_l(j_l+1)+nu sum_f(1-W_f).

The trivalent recoupling formula is Eq.(22) of
https://www.nature.com/articles/s42005-024-01697-4 . For spectator e, forward
edge f and backward edge b at a visited vertex, its unnormalized trace factor
is (-1)^(j_e+j_f+J_b+1/2) sqrt(d_Jf*d_jb)
times {j_e j_f j_b; 1/2 J_b J_f}. Multiply all vertices and divide by two.
This module retains the signs and exact rational squared magnitudes.
The phase convention differs from product-of-loop-characters at J=1/2 by a
single diagonal sign matrix, common to all six cube faces. An independent
exact Haar contraction verified that equivalence for all 192 directed
entries. It differs from the coupled-projector ladder basis by
(-1)^(jL+jR-jM). Negative matrix entries must not be discarded.
"""
from __future__ import annotations

from functools import lru_cache
from fractions import Fraction
import itertools
import json
import math
import time

import numpy as np

from su2_recoupling import admissible, triangle_delta_squared, wigner_6j_squared


VERTICES = tuple(range(8))
EDGES = tuple((u,v) for u in VERTICES for v in VERTICES if u<v and (u^v) in (1,2,4))
FACES = ((0,1,3,2),(4,5,7,6),(0,1,5,4),(2,3,7,6),(0,2,6,4),(1,3,7,5))
EDGE_INDEX = {e:i for i,e in enumerate(EDGES)}
INCIDENT = tuple(tuple(i for i,e in enumerate(EDGES) if v in e) for v in VERTICES)
FACE_EDGES = tuple(tuple(EDGE_INDEX[tuple(sorted((cyc[k],cyc[(k+1)%len(cyc)])))]
                              for k in range(len(cyc))) for cyc in FACES)


def graph_basis(max_twice: int, edges=EDGES, vertices=VERTICES, *, max_four_energy=None):
    """Backtracking with exact triangle/parity constraints; supports degree 2/3."""
    if max_twice < 0 or int(max_twice) != max_twice:
        raise ValueError('cutoff must be a nonnegative doubled-spin integer')
    incident = {v:tuple(i for i,e in enumerate(edges) if v in e) for v in vertices}
    if any(len(es) not in (2,3) for es in incident.values()):
        raise ValueError('only degree-two and degree-three vertices are supported')
    edge_vertices = [tuple(v for v in vertices if i in incident[v]) for i in range(len(edges))]
    assigned = [-1]*len(edges)
    answers = []
    all_values = set(range(max_twice+1))
    costs = [q*(q+2) for q in range(max_twice+1)]
    if max_four_energy is not None:
        completions = {
            2:[(q,q) for q in range(max_twice+1)],
            3:[t for t in itertools.product(range(max_twice+1),repeat=3) if admissible(*t)],
        }
        @lru_cache(maxsize=None)
        def minimum_vertex_cost(values):
            best = math.inf
            for t in completions[len(values)]:
                if all(v<0 or v==q for v,q in zip(values,t)):
                    best = min(best,sum(costs[q] for v,q in zip(values,t) if v<0))
            return best
    def candidates(e):
        allowed = all_values.copy()
        for v in edge_vertices[e]:
            others = [assigned[k] for k in incident[v] if k != e]
            if len(others) == 1:
                if others[0] >= 0:
                    allowed.intersection_update([others[0]])
            elif min(others) >= 0:
                a,b = others
                allowed.intersection_update(range(abs(a-b),min(a+b,max_twice)+1,2))
        return sorted(allowed)
    def recurse(left,energy=0):
        if max_four_energy is not None:
            if energy>max_four_energy:
                return
            # Each unassigned edge is counted at its two endpoints. A sum of
            # independent vertex completion minima divided by two is therefore
            # a rigorous lower bound for the remaining electric cost.
            local_min_sum = sum(minimum_vertex_cost(tuple(sorted(assigned[e] for e in es)))
                                for es in incident.values())
            if 2*energy+local_min_sum>2*max_four_energy:
                return
        if not left:
            answers.append(tuple(assigned))
            return
        best = None
        for e in left:
            vals = candidates(e)
            if not vals:
                return
            known = sum(assigned[k] >= 0 for v in edge_vertices[e] for k in incident[v] if k != e)
            score = (len(vals),-known,e)
            if best is None or score < best[0]:
                best = score,e,vals
        _,e,vals = best
        remaining = [k for k in left if k != e]
        for value in vals:
            assigned[e] = value
            recurse(remaining,energy+costs[value])
        assigned[e] = -1
    recurse(list(range(len(edges))))
    return sorted(answers)


def cube_basis(max_twice: int):
    return graph_basis(max_twice)


def energy_basis(energy_cutoff: Fraction):
    """All physical states with sum j(j+1)<=energy_cutoff, without a spin cutoff.

    Arithmetic pruning uses integer four-times-electric-energy. Positivity
    implies q(q+2)<=4*Ecut for each doubled spin q, hence the sufficient
    per-link bound q<=isqrt(floor(4*Ecut)+1)-1. The returned list is complete
    for this energy cutoff; it is not a subset of a separately fixed J basis.
    """
    cap = Fraction(energy_cutoff)
    if cap<0:
        return []
    max_four_energy = (4*cap).numerator//(4*cap).denominator
    max_twice = math.isqrt(max_four_energy+1)-1
    return graph_basis(max_twice,max_four_energy=max_four_energy)


def electric_diagonal(states,kappa=1.0):
    return np.array([kappa*sum(j*(j+2) for j in s)/4 for s in states],dtype=float)


@lru_cache(maxsize=None)
def sixj_sign(a,b,c,d,e,f):
    if not all(triangle_delta_squared(*t) for t in [(a,b,c),(a,e,f),(d,b,f),(d,e,c)]):
        return 0
    lower = [(a+b+c)//2,(a+e+f)//2,(d+b+f)//2,(d+e+c)//2]
    upper = [(a+b+d+e)//2,(a+c+d+f)//2,(b+c+e+f)//2]
    total = Fraction(0)
    for z in range(max(lower),min(upper)+1):
        den = math.prod(math.factorial(z-x) for x in lower)*math.prod(math.factorial(x-z) for x in upper)
        total += Fraction((-1)**z*math.factorial(z+1),den)
    return (total>0)-(total<0)


@lru_cache(maxsize=None)
def local_vertex_factor(spectator,forward,backward,new_forward,new_backward):
    """Return sign and exact square of one normalized-vertex recoupling factor."""
    args = (spectator,forward,backward,1,new_backward,new_forward)
    square = wigner_6j_squared(*args)*(new_forward+1)*(backward+1)
    if not square:
        return 0,Fraction(0)
    exponent_twice = spectator+forward+new_backward+1
    if exponent_twice % 2:
        raise AssertionError('nonintegral recoupling phase')
    sign = (-1)**(exponent_twice//2)*sixj_sign(*args)
    return sign,square


def graph_face_data(edges,vertices,face_vertices):
    index = {tuple(sorted(e)):i for i,e in enumerate(edges)}
    incident = {v:tuple(i for i,e in enumerate(edges) if v in e) for v in vertices}
    cycle_edges = tuple(index[tuple(sorted((face_vertices[k],face_vertices[(k+1)%len(face_vertices)])))]
                        for k in range(len(face_vertices)))
    data = []
    for k,v in enumerate(face_vertices):
        f,b = cycle_edges[k],cycle_edges[k-1]
        spectator = [i for i in incident[v] if i not in (f,b)]
        if len(spectator)>1:
            raise ValueError('face has nontrivalent vertex')
        data.append((spectator[0] if spectator else None,f,b))
    return cycle_edges,tuple(data)


FACE_DATA = tuple(graph_face_data(EDGES,VERTICES,f) for f in FACES)


def face_amplitude_exact(initial,final,face_data):
    """Return (sign, Fraction squared magnitude) of normalized Wilson trace."""
    cycle_edges,vertices = face_data
    active = set(cycle_edges)
    if any(initial[k] != final[k] for k in range(len(initial)) if k not in active):
        return 0,Fraction(0)
    if any(abs(initial[k]-final[k]) != 1 for k in active):
        return 0,Fraction(0)
    sign,square = 1,Fraction(1,4)
    for spectator,f,b in vertices:
        x = initial[spectator] if spectator is not None else 0
        phase,factor = local_vertex_factor(x,initial[f],initial[b],final[f],final[b])
        if not factor:
            return 0,Fraction(0)
        sign *= phase
        square *= factor
    return sign,square


def amplitude_float(exact):
    sign,square = exact
    return sign*math.sqrt(float(square))


def state_face_transitions(initial,face_data,max_twice=None):
    """Yield final spin tuple, sign, exact squared magnitude, with no basis lookup."""
    cycle_edges,_ = face_data
    for shifts in itertools.product((-1,1),repeat=len(cycle_edges)):
        target = list(initial)
        for e,delta in zip(cycle_edges,shifts):
            target[e] += delta
        if min(target)<0 or (max_twice is not None and max(target)>max_twice):
            continue
        target = tuple(target)
        sign,square = face_amplitude_exact(initial,target,face_data)
        if square:
            yield target,sign,square


def graph_wilson_transitions(states,edges=EDGES,vertices=VERTICES,faces=FACES):
    """Return per-face sparse lists (row,col,sign,Fraction squared magnitude)."""
    lookup = {s:i for i,s in enumerate(states)}
    out = []
    for face in faces:
        data = graph_face_data(edges,vertices,face)
        entries = []
        for col,state in enumerate(states):
            for target,sign,square in state_face_transitions(state,data):
                row = lookup.get(target)
                if row is not None:
                    entries.append((row,col,sign,square))
        out.append(entries)
    return out


def cube_wilson_transitions(max_twice: int):
    states = cube_basis(max_twice)
    return states,graph_wilson_transitions(states)


def dense_wilson_matrices(states,transition_lists):
    matrices = []
    for entries in transition_lists:
        w = np.zeros((len(states),len(states)))
        for row,col,sign,square in entries:
            w[row,col] = sign*math.sqrt(float(square))
        matrices.append(w)
    return matrices


def cube_wilson_matrices(max_twice: int):
    states,entries = cube_wilson_transitions(max_twice)
    return states,dense_wilson_matrices(states,entries)


def self_checks():
    """Finite exact structural checks; no claim of a continuum limit."""
    from su2_recoupling import wilson_matrices as theta_wilson
    result = {'primary_recoupling_source':'https://www.nature.com/articles/s42005-024-01697-4',
              'source_equation':22,'normalization':'W=Tr_fund/2',
              'vertices':VERTICES,'edges':EDGES,'faces':FACES,
              'ladder_comparison':[],'cube_cutoffs':[]}
    edges = ((0,1),(1,2),(3,4),(4,5),(0,3),(1,4),(2,5))
    faces = ((0,1,4,3),(1,2,5,4))
    for cutoff in [1,2,3]:
        states = graph_basis(cutoff,edges,range(6))
        transitions = graph_wilson_transitions(states,edges,range(6),faces)
        mats = dense_wilson_matrices(states,transitions)
        reference_states,wa,wb = theta_wilson(cutoff)
        triples = [(s[0],s[1],s[5]) for s in states]
        order = [reference_states.index(t) for t in triples]
        phase = np.array([(-1)**((a+b-c)//2) for a,b,c in triples])
        errors = [float(np.max(abs(w*np.outer(phase,phase)-ref[np.ix_(order,order)])))
                  for w,ref in zip(mats,[wa,wb])]
        assert max(errors)<1e-14
        result['ladder_comparison'].append({'jmax':cutoff/2,'dimension':len(states),
            'phase':'(-1)^(jL+jR-jM)','WA_error':errors[0],'WB_error':errors[1]})
    for cutoff in [1,2,3]:
        start = time.perf_counter()
        states,transitions = cube_wilson_transitions(cutoff)
        generation_seconds = time.perf_counter()-start
        assert len(states)=={1:32,2:1013,3:14879}[cutoff]
        assert all(admissible(*(s[e] for e in es)) for s in states for es in INCIDENT)
        reverse_checked = 0
        for entries in transitions:
            lookup = {(r,c):(sign,sq) for r,c,sign,sq in entries}
            assert len(lookup)==len(entries)
            for r,c,sign,sq in entries:
                assert lookup[(c,r)]==(sign,sq)
                reverse_checked += 1
        vacuum = (0,)*len(EDGES)
        for data in FACE_DATA:
            candidate = list(vacuum)
            for e in data[0]: candidate[e]=1
            assert face_amplitude_exact(vacuum,tuple(candidate),data)==(1,Fraction(1,4))
        result['cube_cutoffs'].append({'jmax':cutoff/2,'dimension':len(states),
            'directed_entries_per_face':[len(t) for t in transitions],
            'negative_directed_entries':sum(sign<0 for es in transitions for _,_,sign,_ in es),
            'exact_Hermiticity_entries_checked':reverse_checked,
            'vacuum_to_face_amplitude_exact':'1/2',
            'generation_seconds':generation_seconds})
    return result


def main():
    results = self_checks()
    full_spin_two=cube_basis(4)
    full_energies=electric_diagonal(full_spin_two)
    energy_checks=[]
    for cap in [8,10,12,18,24,30]:
        states=energy_basis(Fraction(cap))
        overlap=[s for s in states if max(s)<=4]
        independent_overlap=[s for s,t in zip(full_spin_two,full_energies) if t<=cap]
        assert overlap==independent_overlap
        energy_checks.append({'electric_energy_cap':cap,'dimension':len(states),
                              'largest_link_spin_present':max(max(s) for s in states)/2,
                              'states_beyond_link_spin_two':len(states)-len(overlap),
                              'all_overlapping_states_match_full_spin_two_basis':True})
    results['spin_two_basis_dimension']=len(full_spin_two)
    results['energy_cutoff_checks']=energy_checks
    from pathlib import Path
    output=Path(__file__).resolve().parent.parent/'outputs'/'SU2_Cube_Recoupling_Validation.json'
    output.parent.mkdir(parents=True,exist_ok=True)
    output.write_text(json.dumps(results,indent=2)+'\n',encoding='utf-8')
    print(json.dumps(results,indent=2))
    return results


if __name__ == '__main__':
    main()
