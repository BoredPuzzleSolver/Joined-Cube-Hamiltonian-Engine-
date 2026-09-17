"""Independent SU(2) Haar multiplication via local invariant/CG tensors.

Uses angular-momentum Casimir eigenspaces and lowering operators, no 6j.
All edge matrices are D(U) epsilon; local invariants have all legs outgoing.
"""
from functools import lru_cache
import itertools
import math
import numpy as np


@lru_cache(None)
def angular(q):
    j=q/2;m=np.arange(j,-j-1,-1)
    z=np.diag(m)
    lower=np.zeros((q+1,q+1))
    for col in range(q):lower[col+1,col]=math.sqrt((j+m[col])*(j-m[col]+1))
    return z,lower


@lru_cache(None)
def clebsch(q,r):
    """C(j,m;1/2,a | J,M), J=r/2."""
    z,l=angular(q);zh,lh=angular(1)
    totalz=np.kron(z,np.eye(2))+np.kron(np.eye(q+1),zh)
    lower=np.kron(l,np.eye(2))+np.kron(np.eye(q+1),lh)
    casimir=totalz@totalz+(lower.T@lower+lower@lower.T)/2
    highest=np.where(abs(np.diag(totalz)-r/2)<1e-12)[0]
    eigenvalues,eigenvectors=np.linalg.eigh(casimir[np.ix_(highest,highest)])
    selected=int(np.argmin(abs(eigenvalues-r*(r+2)/4)))
    assert abs(eigenvalues[selected]-r*(r+2)/4)<1e-10
    vector=np.zeros(2*(q+1));vector[highest]=eigenvectors[:,selected]
    first=np.where(abs(vector)>1e-10)[0][0]
    if vector[first]<0:vector=-vector
    output=np.zeros((q+1,2,r+1))
    for col in range(r+1):
        output[:,:,col]=vector.reshape(q+1,2)
        if col<r:vector=lower@vector/math.sqrt((r-col)*(col+1))
    assert np.max(abs(output.reshape(-1,r+1).T@output.reshape(-1,r+1)-np.eye(r+1)))<1e-10
    return output


@lru_cache(None)
def invariant(labels):
    dims=[q+1 for q in labels];dim=math.prod(dims)
    totalz=np.zeros((dim,dim));lower=np.zeros((dim,dim))
    for site,q in enumerate(labels):
        z,l=angular(q)
        for value,out in ((z,totalz),(l,lower)):
            factors=[np.eye(d) for d in dims];factors[site]=value
            product=factors[0]
            for factor in factors[1:]:product=np.kron(product,factor)
            out+=product
    casimir=totalz@totalz+(lower.T@lower+lower@lower.T)/2
    zero=np.where(abs(np.diag(totalz))<1e-12)[0]
    eigenvalues,eigenvectors=np.linalg.eigh(casimir[np.ix_(zero,zero)])
    singlets=np.where(abs(eigenvalues)<1e-9)[0]
    assert len(singlets)==1,(labels,eigenvalues)
    vector=np.zeros(dim);vector[zero]=eigenvectors[:,singlets[0]]
    first=np.where(abs(vector)>1e-10)[0][0]
    if vector[first]<0:vector=-vector
    assert np.max(abs(casimir@vector))<1e-9
    return vector.reshape(dims)


@lru_cache(None)
def local_contraction(old,new,forward_position,backward_position):
    spectator=next(i for i in range(3) if i not in (forward_position,backward_position))
    order=(forward_position,backward_position,spectator)
    initial=invariant(old).transpose(order)
    final=invariant(new).transpose(order)
    cforward=clebsch(old[forward_position],new[forward_position])
    cbackward=clebsch(old[backward_position],new[backward_position])
    epsilon=np.array([[0.,1.],[-1.,0.]])
    return float(np.einsum('ijk,IJk,iaI,jbJ,ab->',initial,final,cforward,cbackward,epsilon))


def tensor_face_matrix(states,edges,vertices,face,patterns=None):
    lookup={tuple(sorted(e)):i for i,e in enumerate(edges)}
    incident={v:tuple(sorted(i for i,e in enumerate(edges) if v in e)) for v in vertices}
    ordered_edges=[lookup[tuple(sorted((v,face[(i+1)%len(face)])))] for i,v in enumerate(face)]
    reverse_count=sum(v>face[(i+1)%len(face)] for i,v in enumerate(face))
    # Independent complete retained support: fundamental multiplication shifts
    # every visited doubled spin by +/-1. No recoupling-symbol zero test is used.
    index={s:i for i,s in enumerate(states)}
    candidates=[]
    for col,state in enumerate(states):
        choices=[tuple(q for q in (state[e]-1,state[e]+1) if q>=0) for e in ordered_edges]
        for labels in itertools.product(*choices):
            target=list(state)
            for e,q in zip(ordered_edges,labels):target[e]=q
            row=index.get(tuple(target))
            if row is not None:candidates.append((row,col))
    matrix=np.zeros((len(states),len(states)))
    for row,col in candidates:
        initial,final=states[col],states[row]
        coefficient=.5*(-1)**reverse_count
        coefficient*=math.prod(math.sqrt((initial[e]+1)/(final[e]+1)) for e in ordered_edges)
        for i,v in enumerate(face):
            es=incident[v]
            old=tuple(initial[e] for e in es);new=tuple(final[e] for e in es)
            coefficient*=local_contraction(old,new,es.index(ordered_edges[i]),es.index(ordered_edges[i-1]))
        matrix[row,col]=coefficient
    if patterns is not None:
        expected=np.zeros_like(matrix)
        for row,col,_,square in patterns:expected[row,col]=math.sqrt(float(square))
        assert np.max(abs(abs(matrix)-expected))<1e-10
    assert np.max(abs(matrix-matrix.T))<1e-10
    return matrix


def probe():
    from su2_joined_cubes_independent_check import independently_resolved_graph
    from su2_cube_model import graph_wilson_transitions
    output=[]
    for pairing in ('longitudinal','crossed','mixed'):
        states,edges,vertices,faces,p=independently_resolved_graph(2,pairing)
        sparse=graph_wilson_transitions(states,edges,vertices,faces)
        matrix=np.diag([sum(q*(q+2) for q in s[:p])/4+len(faces) for s in states])
        for face,entries in zip(faces,sparse):matrix-=tensor_face_matrix(states,edges,vertices,face,entries)
        spectrum=np.linalg.eigvalsh(matrix)
        output.append(spectrum)
        print(pairing,spectrum[:6],spectrum[1]-spectrum[0],flush=True)
    print('max spectrum difference',max(np.max(abs(output[0]-s)) for s in output),flush=True)


if __name__=='__main__':probe()
