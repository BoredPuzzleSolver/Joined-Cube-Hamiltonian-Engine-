"""Omitted-state bounds on the untruncated, fixed, open SU(2) cube.

Threshold minimization is exact integer arithmetic. Refined spectra evaluated
with NumPy/Ritz methods are numerical; rational endpoint certificates belong
in the separate certificate script.
"""
from fractions import Fraction
from functools import lru_cache
from pathlib import Path
import json
import math
import time
import numpy as np

VERTICES = tuple(range(8))
EDGES = tuple((u,v) for u in VERTICES for v in VERTICES
              if u < v and (u^v) in (1,2,4))
INCIDENT = tuple(tuple(i for i,e in enumerate(EDGES) if v in e)
                 for v in VERTICES)
ROOT = Path(__file__).resolve().parents[1]


def casimir4(a):
    return a*(a+2)


def admissible(a, b, c):
    return abs(a-b) <= c <= a+b and (a+b+c)%2 == 0


def electric4(state):
    return sum(casimir4(a) for a in state)


def endpoint_lower_bound4(marked):
    """Marked edge plus its four distinct neighboring edges, in units kappa/4."""
    a, b = marked//2, (marked+1)//2
    return casimir4(marked)+2*(casimir4(a)+casimir4(b))


def theta_witness(marked):
    """Two adjacent face return paths, split as evenly as possible."""
    a, b = marked//2, (marked+1)//2
    index = {edge: i for i, edge in enumerate(EDGES)}
    state = [0]*len(EDGES)
    state[index[(0,1)]] = marked
    for edge in [(0,2),(2,3),(1,3)]:
        state[index[edge]] = a
    for edge in [(0,4),(4,5),(1,5)]:
        state[index[edge]] = b
    assert all(admissible(*(state[e] for e in inc)) for inc in INCIDENT)
    return tuple(state)


@lru_cache(maxsize=None)
def outside_threshold_data(max_twice):
    """Exact branch-and-bound proof of the minimum outside a per-link cutoff.

    Cube edge transitivity maps any maximal-spin edge onto edge (0,1).
    Search all possible maximal doubled spins q. A monotonically increasing
    bound on the marked edge and its four neighbors ends that search.
    """
    if int(max_twice) != max_twice or max_twice < 0:
        raise ValueError("Use a nonnegative integer doubled-spin cutoff")
    first = int(max_twice)+1
    witness = theta_witness(first)
    best = electric4(witness)
    explored = []
    q = first
    while endpoint_lower_bound4(q) < best:
        assigned = [-1]*len(EDGES)
        assigned[0] = q
        costs = [casimir4(v) for v in range(q+1)]
        nodes = 0
        found = 0

        @lru_cache(maxsize=None)
        def vertex_minimum(values):
            known = [v for v in values if v >= 0]
            if len(known) == 0:
                return 0
            if len(known) == 1:
                x = known[0]
                return costs[x]+costs[x//2]+costs[(x+1)//2]
            if len(known) == 2:
                a, b = known
                return costs[a]+costs[b]+costs[abs(a-b)]
            a, b, c = known
            return costs[a]+costs[b]+costs[c] if admissible(a,b,c) else math.inf

        def values_for(edge, partial_cost):
            allowed = set(range(q+1))
            for vertex in EDGES[edge]:
                others = [assigned[k] for k in INCIDENT[vertex] if k != edge]
                if min(others) >= 0:
                    a, b = others
                    allowed.intersection_update(range(abs(a-b),min(a+b,q)+1,2))
            return [v for v in sorted(allowed) if partial_cost+costs[v] < best]

        def recurse(left, partial_cost):
            nonlocal best, witness, nodes, found
            nodes += 1
            if partial_cost >= best:
                return
            # Every edge occurs at two vertices. Independent local minima,
            # summed and divided by two, give a safe global energy bound.
            local_sum = sum(vertex_minimum(tuple(sorted(assigned[e] for e in inc)))
                            for inc in INCIDENT)
            if local_sum >= 2*best:
                return
            if not left:
                best, witness = partial_cost, tuple(assigned)
                found += 1
                return
            chosen = None
            for edge in left:
                vals = values_for(edge, partial_cost)
                if not vals:
                    return
                known_neighbors = sum(assigned[k] >= 0
                                      for v in EDGES[edge] for k in INCIDENT[v]
                                      if k != edge)
                score = (len(vals),-known_neighbors,edge)
                if chosen is None or score < chosen[0]:
                    chosen = score, edge, vals
            _, edge, vals = chosen
            remaining = [e for e in left if e != edge]
            for value in vals:
                assigned[edge] = value
                recurse(remaining,partial_cost+costs[value])
            assigned[edge] = -1

        recurse(list(range(1,len(EDGES))),costs[q])
        explored.append({"maximal_twice_spin": q, "visited_nodes": nodes,
                         "strict_improvements": found})
        q += 1
    assert max(witness) > max_twice
    assert electric4(witness) == best
    assert all(admissible(*(witness[e] for e in inc)) for inc in INCIDENT)
    return {
        "cutoff_twice": max_twice,
        "minimum_electric_energy_over_kappa_exact": str(Fraction(best,4)),
        "minimum_electric_energy_over_kappa": best/4,
        "minimum_energy_integer_units_kappa_over_4": best,
        "minimizing_doubled_spin_state": list(witness),
        "searched_maximal_spins": explored,
        "all_larger_maximal_spins_excluded_from": q,
        "endpoint_lower_bound_at_exclusion_integer_units": endpoint_lower_bound4(q),
        "proof_method": (
            "Exact integer constraint enumeration with admissibility and "
            "energy pruning; edge transitivity reduces the maximal edge "
            "to (0,1); the increasing marked-edge endpoint bound excludes "
            "every larger maximal spin."
        ),
    }


def outside_electric_threshold(max_twice, kappa=Fraction(1)):
    kappa = kappa if isinstance(kappa,Fraction) else Fraction(str(kappa))
    if kappa <= 0:
        raise ValueError("The electric coefficient must be positive")
    return kappa*Fraction(outside_threshold_data(max_twice)[
        "minimum_energy_integer_units_kappa_over_4"],4)


def numerical_boundary_comparison(max_twice, kappa=1.0, nu=1.0):
    """Keep only actual P-to-Q neighbors in the finite auxiliary matrix.

    All other states retain their electric operator in the lower comparison.
    The safe remainder threshold is the minimum outside P, even though the
    reachable portion of that complement has been explicitly retained.
    """
    from su2_cube_model import (
        cube_basis, FACE_DATA, state_face_transitions, electric_diagonal
    )
    from su2_sparse_spectrum import SymmetricOperator, smallest_ritz
    started = time.perf_counter()
    p_states = cube_basis(max_twice)
    p_set = set(p_states)
    boundary = set()
    transitions = []
    for col, state in enumerate(p_states):
        for face in FACE_DATA:
            for target, sign, square in state_face_transitions(state,face):
                if target not in p_set:
                    boundary.add(target)
                transitions.append((col,target,sign,square))
    states = p_states+sorted(boundary)
    lookup = {state:i for i,state in enumerate(states)}
    retained = len(p_states)
    electric = electric_diagonal(states,kappa)
    diagonal = electric.copy()
    diagonal[:retained] += 6*nu
    rows, cols, values = [], [], []
    for col,target,sign,square in transitions:
        row = lookup[target]
        value = -nu*sign*math.sqrt(float(square))
        rows.append(row); cols.append(col); values.append(value)
        if row >= retained:
            rows.append(col); cols.append(row); values.append(value)
    rows=np.array(rows,dtype=np.int64)
    cols=np.array(cols,dtype=np.int64)
    values=np.array(values)
    lower = SymmetricOperator(diagonal,rows,cols,values)
    keep=(rows<retained)&(cols<retained)
    upper = SymmetricOperator(diagonal[:retained],rows[keep],cols[keep],values[keep])

    def solve(operator):
        if operator.size <= 1100:
            matrix=operator.dense()
            vals,vecs=np.linalg.eigh(matrix)
            residuals=[float(np.linalg.norm(matrix@vecs[:,j]-vals[j]*vecs[:,j]))
                       for j in range(2)]
            return vals[:2], {
                "method":"NumPy eigh", "residual_norms":residuals,
                "converged":max(residuals)<1e-9,
            }
        vals,_,audit=smallest_ritz(operator,levels=2,max_iterations=300,
                                 tolerance=1e-10)
        second,_,audit2=smallest_ritz(operator,levels=2,max_iterations=300,
                                    tolerance=1e-10,seed=7313)
        audit["second_seed_max_energy_difference"]=float(np.max(abs(vals-second)))
        audit["second_seed_converged"]=audit2["converged"]
        return vals,audit

    lower_values,lower_audit=solve(lower)
    upper_values,upper_audit=solve(upper)
    threshold=float(outside_electric_threshold(max_twice,Fraction(str(kappa))))
    low=np.minimum(lower_values,threshold)
    return {
        "jmax":max_twice/2, "kappa":kappa, "nu":nu,
        "retained_dimension":retained,
        "reachable_omitted_dimension":len(boundary),
        "auxiliary_dimension":len(states),
        "safe_remaining_electric_threshold":threshold,
        "auxiliary_low_eigenvalues":lower_values.tolist(),
        "projection_low_eigenvalues":upper_values.tolist(),
        "energy_intervals_evaluated_numerically":[
            [float(low[j]),float(upper_values[j])] for j in range(2)
        ],
        "gap_interval_evaluated_numerically":[
            float(low[1]-upper_values[0]),float(upper_values[1]-low[0])
        ],
        "lower_solver":lower_audit,"upper_solver":upper_audit,
        "machine_certified":False,
        "limitation":(
            "The operator comparison is exact, but these eigenvalue endpoints "
            "are numerical estimates. Residuals and multiple seeds do not "
            "certify eigenvalue ordering or omitted floating-point errors."
        ),
        "elapsed_seconds":time.perf_counter()-started,
    }


def main():
    thresholds=[]
    for cutoff in [0,1,2,3,4,5]:
        start=time.perf_counter()
        row=dict(outside_threshold_data(cutoff))
        row["elapsed_seconds"]=time.perf_counter()-start
        thresholds.append(row)
        print(json.dumps(row),flush=True)
    numerical=[]
    for cutoff in [1,2,3]:
        row=numerical_boundary_comparison(cutoff)
        numerical.append(row)
        print(json.dumps(row),flush=True)
    output={
        "title":"Omitted-state comparison for the fixed open SU(2) cube",
        "edges":EDGES,
        "hamiltonian":"H=kappa sum_links j(j+1)+nu(6-sum_faces Tr U_f/2).",
        "scope":"Untruncated gauge-invariant Hilbert space of the single open cube.",
        "exact_electric_thresholds":thresholds,
        "numerical_auxiliary_bounds":numerical,
    }
    path=ROOT/"outputs"/"SU2_Cube_Bounds_Results.json"
    path.write_text(json.dumps(output,indent=2,allow_nan=False)+"\n",encoding="utf-8")
    print(path,flush=True)


if __name__=="__main__":
    main()
