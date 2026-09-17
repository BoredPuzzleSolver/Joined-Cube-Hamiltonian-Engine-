"""Exact bases and Hamiltonian assembly for the original SU(2) cube models.

State order, doubled-spin labels, recoupling phases, and rational squared
amplitudes are inherited from the original research scripts.
"""
from __future__ import annotations

from array import array
from dataclasses import dataclass
from fractions import Fraction
from functools import cached_property, lru_cache
import itertools
import math

import numpy as np

from .geometry import EDGES, VERTICES, FACES, FACE_DATA, graph_face_data, joined_graph
from .recoupling import admissible, local_vertex_factor, singlet_exists, coupling_channels
from .recoupling import oriented_local_vertex_factor as joined_local_vertex_factor

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


def cube_energy_basis(energy_cutoff: Fraction):
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


def cube_electric_diagonal(states,kappa=1.0):
    return np.array([kappa*sum(j*(j+2) for j in s)/4 for s in states],dtype=float)


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


def cube_state_face_transitions(initial,face_data,max_twice=None):
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
            for target,sign,square in cube_state_face_transitions(state,data):
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


@lru_cache(maxsize=None)
def minimum_vertex_completion_cost(values):
    """Exact minimum unassigned physical cost sum q(q+2), without a spin cap.

    With fixed labels of sum S and maximum M, missing labels need total
    Q>=max(0,2M-S), with Q+S even. The least such Q is optimal. Convexity
    distributes it as evenly as possible among missing labels. This satisfies
    the remaining polygon inequalities; any finite spin cap can only raise
    the minimum, so it is always a safe pruning lower bound.
    """
    known=[q for q in values if q>=0]
    missing=len(values)-len(known)
    if not missing:return 0 if singlet_exists(known) else math.inf
    total=sum(known)
    required=max(0,2*max(known,default=0)-total)
    if (required+total)%2:required+=1
    base,remainder=divmod(required,missing)
    return (missing-remainder)*base*(base+2)+remainder*(base+1)*(base+3)


def _basis(model,max_physical_twice,max_four_energy=None):
    """Enumerate physical labels, then expand every independent intertwiner."""
    p=model.physical_edge_count
    assigned=[-1]*p
    answers=[]
    costs=[q*(q+2) for q in range(max_physical_twice+1)]
    def candidates(e):
        low,high,parity=0,max_physical_twice,None
        for v in model.physical_edges[e]:
            other=[assigned[k] for k in model.physical_incident[v] if k!=e]
            if min(other)>=0:
                total=sum(other)
                low=max(low,2*max(other)-total)
                high=min(high,total)
                wanted=total%2
                if parity is not None and parity!=wanted:return ()
                parity=wanted
        if parity is None:return range(low,high+1)
        if low%2!=parity:low+=1
        return range(low,high+1,2)
    def append_states():
        channels=[]
        for split in model.splits:
            a,b=(assigned[e] for e in split.first_pair)
            c,d=(assigned[e] for e in split.second_pair)
            allowed=coupling_channels(a,b,c,d)
            assert allowed
            channels.append(allowed)
        physical=tuple(assigned)
        for aux in itertools.product(*channels):answers.append(physical+aux)
    def recurse(left,cost):
        if max_four_energy is not None:
            if cost>max_four_energy:return
            bound=sum(minimum_vertex_completion_cost(tuple(sorted(assigned[e] for e in es)))
                      for es in model.physical_incident)
            # Unassigned physical links appear at exactly two real endpoints.
            if 2*cost+bound>2*max_four_energy:return
        if not left:
            append_states()
            return
        chosen=None
        for e in left:
            values=candidates(e)
            if not values:return
            known=sum(assigned[k]>=0 for v in model.physical_edges[e]
                      for k in model.physical_incident[v] if k!=e)
            score=(len(values),-known,e)
            if chosen is None or score<chosen[0]:chosen=(score,e,values)
        _,e,values=chosen
        remaining=[i for i in left if i!=e]
        for q in values:
            assigned[e]=q
            recurse(remaining,cost+costs[q])
        assigned[e]=-1
    recurse(list(range(p)),0)
    return sorted(answers)


def spin_basis(model,max_physical_twice: int):
    if int(max_physical_twice)!=max_physical_twice or max_physical_twice<0:
        raise ValueError('physical cutoff must be a nonnegative doubled-spin integer')
    return _basis(model,int(max_physical_twice))


def energy_basis(model,energy_cutoff: Fraction):
    cap=4*Fraction(energy_cutoff)
    if cap<0:return []
    max_four_energy=cap.numerator//cap.denominator
    max_physical_twice=math.isqrt(max_four_energy+1)-1
    return _basis(model,max_physical_twice,max_four_energy)


def electric_exact(model,state):
    return Fraction(sum(q*(q+2) for q in state[:model.physical_edge_count]),4)


def electric_diagonal(model,states,kappa=1.0):
    return np.array([float(kappa)*sum(q*(q+2) for q in s[:model.physical_edge_count])/4 for s in states])


@lru_cache(maxsize=16384)
def _cycle_transitions(old,spectators,order_parities):
    """Local transfer enumeration avoids blindly trying all 2^length flips."""
    size=len(old)
    choices=[tuple(q for q in (a-1,a+1) if q>=0) for a in old]
    answers=[]
    new=[0]*size
    def recurse(k,sign,square):
        if k==size:
            phase,factor=joined_local_vertex_factor(spectators[0],old[0],old[-1],new[0],new[-1])
            if order_parities[0]:phase*=(-1)**((new[0]+new[-1]-old[0]-old[-1])//2)
            if factor:answers.append((tuple(new),sign*phase,square*factor))
            return
        for q in choices[k]:
            phase,factor=joined_local_vertex_factor(spectators[k],old[k],old[k-1],q,new[k-1])
            if order_parities[k]:phase*=(-1)**((q+new[k-1]-old[k]-old[k-1])//2)
            if factor:
                new[k]=q
                recurse(k+1,sign*phase,square*factor)
    for q in choices[0]:
        new[0]=q
        recurse(1,1,Fraction(1,4))
    return tuple(answers)


def state_face_transitions(model,initial,face_index,max_physical_twice=None):
    cycle,vertices=model.face_data[face_index]
    old=tuple(initial[e] for e in cycle)
    spectators=tuple(initial[x] if x is not None else 0 for x,_,_ in vertices)
    for new,sign,square in _cycle_transitions(old,spectators,model.face_order_parities[face_index]):
        if max_physical_twice is not None and any(q>max_physical_twice
                for e,q in zip(cycle,new) if e<model.physical_edge_count):continue
        target=list(initial)
        for e,q in zip(cycle,new):target[e]=q
        yield tuple(target),sign*model.face_orientation_signs[face_index],square


def wilson_transitions(model,states):
    """Per-face sparse entries (row,col,sign,Fraction squared magnitude)."""
    index={s:i for i,s in enumerate(states)}
    output=[]
    for face in range(model.face_count):
        entries=[]
        for col,state in enumerate(states):
            for target,sign,square in state_face_transitions(model,state,face):
                row=index.get(target)
                if row is not None:entries.append((row,col,sign,square))
        output.append(entries)
    return output


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



def _integer(value, name, minimum=1):
    if isinstance(value, (bool, np.bool_)):
        raise ValueError(f"{name} must be an integer >= {minimum}")
    try:
        converted = int(value)
        valid = converted == value and converted >= minimum
    except (ValueError, TypeError, OverflowError):
        valid = False
    if not valid:
        raise ValueError(f"{name} must be an integer >= {minimum}")
    return converted


def _finite_float(value, name, strictly_positive=False):
    try:
        number = float(value)
    except (ValueError, TypeError, OverflowError) as exc:
        raise ValueError(f"{name} must be a finite real number") from exc
    if not math.isfinite(number) or (number <= 0 if strictly_positive else number < 0):
        bound = "positive" if strictly_positive else "nonnegative"
        raise ValueError(f"{name} must be finite and {bound}")
    return number


def validate_electric_cutoff(value):
    """Return the exact, dimensionless cutoff on physical sum j(j+1)."""
    try:
        # Decimal strings avoid an unintended binary-rational cutoff.
        cap = value if isinstance(value, Fraction) else Fraction(str(value))
    except (ValueError, TypeError, ZeroDivisionError, OverflowError) as exc:
        raise ValueError("electric_cutoff must be finite and nonnegative") from exc
    if cap < 0:
        raise ValueError("electric_cutoff must be finite and nonnegative")
    return cap


@dataclass(frozen=True)
class JoinedCubesModel:
    """Open N-by-1-by-1 SU(2) cube chain, with no external charges.

    H = kappa * sum_physical j(j+1) + nu * sum_faces (1 - Tr(U_face)/2).
    The auxiliary channels resolve Gauss-law invariant tensors and carry
    zero electric energy. Geometry and basis indexing match the source model.
    """

    num_cubes: int = 1
    kappa: float = 1.0
    nu: float = 1.0
    pairing: str = "longitudinal"

    def __post_init__(self):
        object.__setattr__(self, "num_cubes", _integer(self.num_cubes, "num_cubes"))
        object.__setattr__(self, "kappa", _finite_float(self.kappa, "kappa", True))
        object.__setattr__(self, "nu", _finite_float(self.nu, "nu"))
        if self.pairing not in ("longitudinal", "crossed", "mixed"):
            raise ValueError("pairing must be longitudinal, crossed, or mixed")

    @cached_property
    def graph(self):
        return joined_graph(self.num_cubes, pairing=self.pairing)


@dataclass(frozen=True)
class HamiltonianData:
    """Sparse entries for a projected Hamiltonian, in the original basis order."""

    states: tuple
    diagonal: np.ndarray
    rows: np.ndarray
    columns: np.ndarray
    entries: np.ndarray
    hermiticity_error: float

    @property
    def dimension(self):
        return len(self.states)

    def to_scipy(self):
        """Return CSR, summing duplicate entries contributed by different faces."""
        from scipy.sparse import coo_matrix, diags

        magnetic = coo_matrix(
            (self.entries, (self.rows, self.columns)),
            shape=(self.dimension, self.dimension),
        ).tocsr()
        return magnetic + diags(self.diagonal, format="csr")


def build_hamiltonian(model, electric_cutoff=16):
    """Build the complete physical-electric-cutoff compression P H P.

    The cutoff is independent of kappa and applies only to physical links.
    No spin cutoff is imposed on auxiliary intertwiner channels.
    """
    if not isinstance(model, JoinedCubesModel):
        raise TypeError("model must be a JoinedCubesModel")
    cap = validate_electric_cutoff(electric_cutoff)
    graph = model.graph
    states = tuple(energy_basis(graph, cap))
    electric = electric_diagonal(graph, states, model.kappa)
    if model.nu == 0:
        rows = columns = np.empty(0, dtype=np.int64)
        entries = np.empty(0, dtype=float)
        error = 0.0
    else:
        rows, columns, wilson = build_sparse(graph, states)
        error = numerical_hermiticity(rows, columns, wilson, len(states))
        entries = -model.nu * wilson
    diagonal = electric + graph.face_count * model.nu
    if not np.all(np.isfinite(diagonal)) or not np.all(np.isfinite(entries)):
        raise ValueError("Hamiltonian entries exceed floating-point range; rescale kappa and nu")
    return HamiltonianData(
        states=states,
        diagonal=diagonal,
        rows=rows,
        columns=columns,
        entries=entries,
        hermiticity_error=error,
    )
