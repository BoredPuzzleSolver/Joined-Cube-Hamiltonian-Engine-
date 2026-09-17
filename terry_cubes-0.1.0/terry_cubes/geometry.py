"""Original single-cube and joined-cube graph topologies.

Edge, face, vertex, orientation, and auxiliary-channel indexing are retained
from the research scripts. Spin-network enumeration lives in hamiltonian.py.
"""


from __future__ import annotations


from dataclasses import dataclass


VERTICES = tuple(range(8))


EDGES = tuple((u,v) for u in VERTICES for v in VERTICES if u<v and (u^v) in (1,2,4))


FACES = ((0,1,3,2),(4,5,7,6),(0,1,5,4),(2,3,7,6),(0,2,6,4),(1,3,7,5))


EDGE_INDEX = {e:i for i,e in enumerate(EDGES)}


INCIDENT = tuple(tuple(i for i,e in enumerate(EDGES) if v in e) for v in VERTICES)


FACE_EDGES = tuple(tuple(EDGE_INDEX[tuple(sorted((cyc[k],cyc[(k+1)%len(cyc)])))]
                              for k in range(len(cyc))) for cyc in FACES)


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
