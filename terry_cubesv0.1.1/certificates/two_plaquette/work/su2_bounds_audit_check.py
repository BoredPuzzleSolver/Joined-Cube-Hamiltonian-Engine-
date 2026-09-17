"""Independent exact Fraction Schur-complement audit of saved certificates.

This audit uses ordinary rational LDL elimination, not integer Bareiss updates.
"""
from fractions import Fraction as F
from pathlib import Path
import json,time
import numpy as np
import su2_recoupling as rec
import su2_cutoff_bounds as bounds

ROOT=Path(__file__).resolve().parent.parent


def independent_inertia(matrix):
    a=[{j:F(x) for j,x in enumerate(row) if x} for row in matrix]
    active=list(range(len(a)))
    positive=negative=zero=0
    while active:
        pivot_index=next((i for i in active if a[i].get(i,0)),None)
        if pivot_index is None:
            assert not any(a[i].get(j,0) for i in active for j in active), 'Nontrivial zero-diagonal block requires a two-by-two pivot.'
            zero+=len(active)
            break
        k=pivot_index
        pivot=a[k][k]
        positive+=int(pivot>0)
        negative+=int(pivot<0)
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
        active.remove(k)
        a[k]={}
    return {'negative':negative,'zero':zero,'positive':positive}


saved=json.loads((ROOT/'outputs/Shared_SU2_Cutoff_Certificate_Checks.json').read_text())
records=saved['benchmarks'] if isinstance(saved,dict) else saved
assert records, 'No certificates supplied for audit.'
checks=[]
for record in records:
    start=time.perf_counter()
    retained_spin=record.get('retained_spin_cutoff',record.get('J'))
    cutoff=int(2*retained_spin)
    nu=F(record.get('nu_exact',record.get('nu')))
    kappa=F(record.get('kappa_exact',record.get('kappa',1)))
    states,congruence=rec.exact_congruence(cutoff+1,kappa,nu)
    weights=[(a+1)*(b+1)*(c+1) for a,b,c in states]
    electric=[kappa*F(3*a*(a+2)+3*b*(b+2)+c*(c+2),4) for a,b,c in states]
    p=[i for i,s in enumerate(states) if max(s)<=cutoff]
    q=[i for i,s in enumerate(states) if max(s)>cutoff]
    # Build both auxiliary rational matrices ourselves, without calling
    # rational_auxiliary_matrices or shifted_congruence from the audited file.
    lower=[row.copy() for row in congruence]
    for i in q:
        for j in q:
            lower[i][j]=electric[i]/weights[i] if i==j else F(0)
    retained=[[congruence[i][j] for j in p] for i in p]
    retained_weights=[weights[i] for i in p]
    probe_checks=[]
    for row in record['certificate']['energy_intervals']:
        for kind,matrix,metric,probe,saved in [
            ('lower',lower,weights,F(row['lower_exact']),row['lower_auxiliary_inertia']),
            ('upper',retained,retained_weights,F(row['upper_exact']),row['upper_projection_inertia'])]:
            shifted=[line.copy() for line in matrix]
            for i,weight in enumerate(metric):
                shifted[i][i]-=probe/weight
            independent=independent_inertia(shifted)
            bareiss=bounds.rational_inertia(shifted)
            assert independent==saved==bareiss
            probe_checks.append({'level':row['index'],'kind':kind,'probe':str(probe),'inertia':independent})
    # Check exact congruence reproduces the finite Hamiltonian supplied by the
    # validated recoupling matrix construction at floating-point precision.
    _,h=rec.hamiltonian(cutoff+1,float(kappa),float(nu))
    scale=np.sqrt(weights)
    recovered=np.asarray([[float(v) for v in row] for row in congruence])*scale[:,None]*scale[None,:]
    discrepancy=float(np.max(abs(h-recovered)))
    assert discrepancy<1e-12
    energies=record['certificate']['energy_intervals']
    gap=[F(energies[1]['lower_exact'])-F(energies[0]['upper_exact']),
         F(energies[1]['upper_exact'])-F(energies[0]['lower_exact'])]
    assert [str(value) for value in gap]==record['certificate']['gap_interval_exact']
    checks.append({'J':retained_spin,'nu_exact':str(nu),'kappa_exact':str(kappa),'congruence_reconstruction_error':discrepancy,
                   'probe_checks':probe_checks,'gap_interval_exact':[str(value) for value in gap],
                   'elapsed_seconds':time.perf_counter()-start})
    print(json.dumps(checks[-1]),flush=True)

# Exhaustively inspect all possible one-action targets from retained bases.
for cutoff in range(1,13):
    states=rec.basis_twice(cutoff)
    layer=set(rec.basis_twice(cutoff+1))
    for a,b,c in states:
        for sign in (-1,1):
            for middle_sign in (-1,1):
                for target in [(a+sign,b,c+middle_sign),(a,b+sign,c+middle_sign)]:
                    if rec.admissible(*target):
                        assert target in layer
probe_count=sum(len(row['probe_checks']) for row in checks)
assert probe_count==4*len(records)
print(f'All {probe_count} exact probes across {len(records)} benchmarks, matrix identifications, gap subtractions and boundary-layer completeness checks passed.',flush=True)
