"""Open N-by-1-by-1 SU(2) cube chains with complete gauge intertwiners.

All state labels are doubled spins. Physical links come first, followed by
one virtual coupling label for each original degree-four vertex. Virtual
links have ZERO electric energy and are not subject to a physical spin cut.

Splitting a four-valent vertex into two trivalent vertices and an auxiliary
tree edge resolves its invariant tensor space. For paired spins (a,b),(c,d),
k ranges over both pair-coupling intervals. Each common k labels one
orthonormal singlet intertwiner. Four spin-half links have k=0 and k=1.
The auxiliary group variable can be gauge-fixed to identity on that tree;
lifting the physical loops through the tree then gives the original Wilson
operators. The splitting changes E and V equally, creates no new loop, and
adds no propagating gauge-field degree of freedom. Its representation label
retains the original intertwiner multiplicity rather than an extra field.

H=kappa*sum_PHYSICAL_links j(j+1)+nu*sum_PHYSICAL_faces(1-Tr_fund(U_face)/2).
An N-by-1-by-1 chain increases extent only in one spatial direction; no bulk
three-dimensional or continuum limit is claimed by this finite model.

PHASE CONVENTION: orient resolved edges from their lower endpoint to their
higher endpoint and use K_j(U)=D_j(U)epsilon_j. Every trivalent vertex uses
ascending incident edge indices as its fixed, all-outgoing 3j order. Local
invariant/CG contractions are summed with exact rational radical arithmetic
in su2_oriented_vertex.py. The local order (spectator,forward,backward) is
converted to the fixed order with the old/new 3j permutation phases. The
whole face also has one minus sign per edge traversed against its orientation.
Both factors matter: a local ordering correction alone is insufficient on
some resolution trees with odd-length lifted faces. The square of each final
entry remains an exact Fraction and no floating sign decision is made.
"""
from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction
from functools import lru_cache
import itertools
import json
import math
from pathlib import Path
import time

import numpy as np

from su2_recoupling import admissible
from su2_cube_model import graph_face_data
from su2_oriented_vertex import oriented_local_vertex_factor as local_vertex_factor


@dataclass(frozen=True)
class VertexSplit:
    original_vertex: int
    first_vertex: int
    second_vertex: int
    first_pair: tuple
    second_pair: tuple
    auxiliary_edge: int


@dataclass
class JoinedCubeGraph:
    num_cubes: int
    pairing: str
    real_vertices: tuple
    real_coordinates: tuple
    physical_edges: tuple
    physical_faces: tuple
    physical_incident: tuple
    vertices: tuple
    edges: tuple
    faces: tuple
    splits: tuple
    face_data: tuple
    face_order_parities: tuple
    face_orientation_signs: tuple

    @property
    def physical_edge_count(self):return len(self.physical_edges)
    @property
    def virtual_edges(self):return tuple(s.auxiliary_edge for s in self.splits)
    @property
    def real_edges(self):return self.physical_edges
    @property
    def real_faces(self):return self.physical_faces
    @property
    def face_count(self):return len(self.physical_faces)


def joined_graph(num_cubes: int,pairing='longitudinal'):
    """Canonical chain with coordinates x+(N+1)y+2(N+1)z.

    Physical edge order is lexicographic in endpoint IDs. Faces are XY
    (z then x), XZ (y then x), then YZ (x). The default splitting pairs
    the two X links against Y/Z. 'crossed' pairs X-minus/Y against X-plus/Z,
    providing a second equivalent intertwiner resolution for validation.
    'mixed' pairs X-minus with the first transverse edge by edge index;
    its odd-length lifted faces test the orientation convention separately.
    """
    if int(num_cubes)!=num_cubes or num_cubes<1:
        raise ValueError('num_cubes must be a positive integer')
    if pairing not in ('longitudinal','crossed','mixed'):
        raise ValueError('pairing must be longitudinal, crossed, or mixed')
    n=int(num_cubes)
    vid=lambda x,y,z:x+(n+1)*y+2*(n+1)*z
    coords=tuple((v%(n+1),(v//(n+1))%2,v//(2*(n+1))) for v in range(4*(n+1)))
    vertices=tuple(range(len(coords)))
    edges=tuple((u,v) for u in vertices for v in vertices if u<v and
                sum(abs(a-b) for a,b in zip(coords[u],coords[v]))==1)
    faces=[]
    for z in (0,1):
        for x in range(n):faces.append((vid(x,0,z),vid(x+1,0,z),vid(x+1,1,z),vid(x,1,z)))
    for y in (0,1):
        for x in range(n):faces.append((vid(x,y,0),vid(x+1,y,0),vid(x+1,y,1),vid(x,y,1)))
    for x in range(n+1):faces.append((vid(x,0,0),vid(x,1,0),vid(x,1,1),vid(x,0,1)))
    incident=tuple(tuple(i for i,e in enumerate(edges) if v in e) for v in vertices)
    splits=[]
    branch={}
    for v,es in enumerate(incident):
        if len(es)==3:
            for e in es:branch[v,e]=v
            continue
        assert len(es)==4
        directional={}
        for e in es:
            u=edges[e][1] if edges[e][0]==v else edges[e][0]
            delta=tuple(a-b for a,b in zip(coords[u],coords[v]))
            directional[delta]=e
        minus=directional[(-1,0,0)];plus=directional[(1,0,0)]
        yedge=next(e for delta,e in directional.items() if delta[1])
        zedge=next(e for delta,e in directional.items() if delta[2])
        if pairing=='longitudinal':first,second=(minus,plus),(yedge,zedge)
        else:
            selected=yedge if pairing=='crossed' else min(yedge,zedge)
            other_transverse=zedge if selected==yedge else yedge
            first,second=(minus,selected),(plus,other_transverse)
        other=len(vertices)+len(splits)
        split=VertexSplit(v,v,other,tuple(first),tuple(second),len(edges)+len(splits))
        splits.append(split)
        for e in first:branch[v,e]=v
        for e in second:branch[v,e]=other
    resolved_edges=[tuple(sorted((branch[u,e],branch[v,e]))) for e,(u,v) in enumerate(edges)]
    resolved_edges += [(s.first_vertex,s.second_vertex) for s in splits]
    resolved_vertices=tuple(range(len(vertices)+len(splits)))
    original_index={e:i for i,e in enumerate(edges)}
    resolved_faces=[]
    for face in faces:
        lifted=[]
        for i,v in enumerate(face):
            incoming=original_index[tuple(sorted((face[i-1],v)))]
            outgoing=original_index[tuple(sorted((v,face[(i+1)%len(face)])))]
            before,after=branch[v,incoming],branch[v,outgoing]
            lifted.append(before)
            if after!=before:lifted.append(after)
        resolved_faces.append(tuple(lifted))
    resolved_edges=tuple(resolved_edges)
    resolved_faces=tuple(resolved_faces)
    data=tuple(graph_face_data(resolved_edges,resolved_vertices,f) for f in resolved_faces)
    parities=[]
    for _,local_vertices in data:
        local=[]
        for spectator,forward,backward in local_vertices:
            order=(spectator,forward,backward)
            assert spectator is not None
            local.append(sum(order[i]>order[j] for i in range(3) for j in range(i+1,3))%2)
        parities.append(tuple(local))
    assert len(edges)==8*n+4 and len(vertices)==4*n+4 and len(faces)==5*n+1
    assert len(splits)==4*(n-1)
    assert all(sum(v in e for e in resolved_edges)==3 for v in resolved_vertices)
    assert len(set(resolved_edges))==len(resolved_edges)
    assert len(resolved_edges)-len(resolved_vertices)==len(edges)-len(vertices)
    orientation_signs=tuple((-1)**sum(v>f[(i+1)%len(f)] for i,v in enumerate(f))
                            for f in resolved_faces)
    return JoinedCubeGraph(n,pairing,vertices,coords,edges,tuple(faces),incident,
                           resolved_vertices,resolved_edges,resolved_faces,tuple(splits),data,
                           tuple(parities),orientation_signs)


def singlet_exists(spins):
    """SU(2) tensor product contains spin zero iff parity/polygon tests pass."""
    total=sum(spins)
    return total%2==0 and 2*max(spins,default=0)<=total


def coupling_channels(a,b,c,d):
    """Doubled k values common to a⊗b and c⊗d; no extra truncation."""
    if (a+b-c-d)%2:return ()
    low=max(abs(a-b),abs(c-d));high=min(a+b,c+d)
    if (low-a-b)%2:low+=1
    return tuple(range(low,high+1,2))


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
            phase,factor=local_vertex_factor(spectators[0],old[0],old[-1],new[0],new[-1])
            if order_parities[0]:phase*=(-1)**((new[0]+new[-1]-old[0]-old[-1])//2)
            if factor:answers.append((tuple(new),sign*phase,square*factor))
            return
        for q in choices[k]:
            phase,factor=local_vertex_factor(spectators[k],old[k],old[k-1],q,new[k-1])
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


def graph_metadata(model):
    return {'num_cubes':model.num_cubes,'shape':[model.num_cubes,1,1],
            'pairing':model.pairing,'physical_vertices':len(model.real_vertices),
            'physical_edges':model.physical_edge_count,'physical_faces':model.face_count,
            'resolved_vertices':len(model.vertices),'resolved_edges':len(model.edges),
            'auxiliary_coupling_edges':len(model.splits),
            'independent_graph_cycles':model.physical_edge_count-len(model.real_vertices)+1,
            'lifted_face_lengths':[len(f) for f in model.faces],
            'electric_energy':'Sum of j(j+1) on physical links only',
            'intertwiners':'All common pair-coupling channels retained; no physical-spin cutoff is applied to auxiliary labels',
            'fixed_3j_order':'Ascending resolved incident edge indices; exact old/new permutation phases',
            'oriented_edge_convention':'K_j(U)=D_j(U)epsilon_j, lower endpoint to higher; minus sign per reversed fundamental edge',
            'magnetic_sign_arithmetic':'Exact 3j finite sums and spin-half CG contractions; exact rational squared magnitudes',
            'sources':['https://dlmf.nist.gov/34.2#E4','https://dlmf.nist.gov/34.3#ii',
                       'https://www.nature.com/articles/s42005-024-01697-4']}


if __name__=='__main__':
    for n in (1,2,3):
        model=joined_graph(n)
        for cap in (8,12):
            start=time.perf_counter();states=energy_basis(model,Fraction(cap))
            print(json.dumps({'num_cubes':n,'energy_cutoff':cap,'dimension':len(states),
                              'seconds':time.perf_counter()-start}),flush=True)
