"""Independent technical audit of fixed-point interval inertia.

Checks arithmetic enclosure, difficult pivots, and all four actual cap-8
cube probes against exact Bareiss inertia. Does not repeat larger probes.
"""
from fractions import Fraction as F
from pathlib import Path
import itertools
import json
import random
import time
from su2_interval_inertia import (_product_quotient,interval_inertia,
                                  adaptive_interval_inertia,UncertifiedInertia)
from su2_cube_model import energy_basis
from su2_cube_rational_certificate import comparison_matrices
from su2_cutoff_bounds import rational_inertia

ROOT=Path(__file__).resolve().parent.parent


def arithmetic_checks():
    rng=random.Random(539172)
    samples=0
    for _ in range(120):
        a=tuple(sorted([rng.randint(-1000,1000),rng.randint(-1000,1000)]))
        b=tuple(sorted([rng.randint(-1000,1000),rng.randint(-1000,1000)]))
        sign=rng.choice([-1,1])
        p=tuple(sorted([sign*rng.randint(1,1000),sign*rng.randint(1,1000)]))
        lo,hi=_product_quotient(a,b,p)
        def sample(interval):
            x,y=map(F,interval)
            return [x,y,(2*x+y)/3,(x+2*y)/3]
        for x,y,z in itertools.product(sample(a),sample(b),sample(p)):
            assert lo<=x*y/z<=hi
            samples+=1
    return samples


def pivot_checks():
    cases=[([[0,1],[1,0]],(1,1)),
           ([[0,1,2],[1,0,3],[2,3,0]],(1,2)),
           ([[0,F(-2,7),1],[F(-2,7),0,3],[1,3,0]],(2,1)),
           ([[3,F(-1,2)],[F(-1,2),-2]],(1,1))]
    results=[]
    for matrix,expected in cases:
        interval=interval_inertia(matrix,bits=64)
        exact=rational_inertia(matrix)
        assert (interval['positive'],interval['negative'])==expected
        assert all(interval[key]==exact[key] for key in ['positive','negative','zero'])
        results.append({'dimension':len(matrix),'positive':interval['positive'],
                        'negative':interval['negative'],'two_by_two_pivots':interval['two_by_two_pivots']})
    tiny=F(1,1<<200)
    for value in [tiny,-tiny]:
        answer=adaptive_interval_inertia([[value,0],[0,1]],bits=64,max_bits=256)
        assert answer['bits_attempted']==[64,128,256]
        assert answer['negative']==int(value<0)
    for matrix in [[[0]],[[1,0],[0,0]],[[1,1],[1,1]]]:
        try:
            interval_inertia(matrix,bits=128)
        except UncertifiedInertia:
            pass
        else:
            raise AssertionError('Exact singular matrix must not be certified nonsingular.')
    return results


def main():
    samples=arithmetic_checks()
    pivots=pivot_checks()
    saved=json.loads((ROOT/'outputs/SU2_Cube_Energy_Certificates.json').read_text())
    row=next(item for item in saved['certificates'] if item['electric_energy_cap_over_kappa']==8)
    retained=energy_basis(F(8))
    states,diagonal,matrix,eta_a,eta_l=comparison_matrices(retained)
    p=len(retained)
    assert p==86 and len(states)==508
    assert str(eta_a)==row['A_error_bound_exact']
    assert str(eta_l)==row['L_error_bound_exact']
    probes=[]
    for energy in row['energy_intervals']:
        for kind,probe,previous in [
            ('lower',F(energy['lower_shift_exact']),energy['lower_inertia']),
            ('upper',F(energy['upper_shift_exact']),energy['upper_inertia'])]:
            assert probe==(F(energy['lower_exact'])+eta_l if kind=='lower' else F(energy['upper_exact'])-eta_a)
            schur=[[matrix[i].get(j,F(0)) if i!=j else diagonal[i]-probe for j in range(p)] for i in range(p)]
            if kind=='lower':
                # Independently assemble the exact Schur complement, updating
                # all ordered neighbor pairs instead of a mirrored triangle.
                for q in range(p,len(states)):
                    denominator=diagonal[q]-probe
                    assert denominator>0
                    for i,qi in matrix[q].items():
                        assert i<p
                        for j,qj in matrix[q].items():
                            assert j<p
                            schur[i][j]-=qi*qj/denominator
            start=time.perf_counter()
            exact=rational_inertia(schur)
            exact_seconds=time.perf_counter()-start
            start=time.perf_counter()
            interval=interval_inertia(schur,bits=128)
            interval_seconds=time.perf_counter()-start
            assert all(interval[key]==exact[key] for key in ['positive','negative','zero'])
            full=exact.copy()
            if kind=='lower':
                full['positive']+=len(states)-p
            assert all(full[key]==previous[key] for key in ['positive','negative','zero'])
            result={'level':energy['index'],'kind':kind,'probe_exact':str(probe),
                    'schur_dimension':p,'schur_exact_inertia':exact,'full_auxiliary_inertia':full,
                    'interval_counts_match_Bareiss':True,'counts_match_saved_certificate':True,
                    'Bareiss_seconds':exact_seconds,'interval_seconds':interval_seconds,
                    'bits':interval['bits']}
            probes.append(result)
            print(json.dumps(result),flush=True)
    output={'title':'Independent technical audit of interval inertia and actual cap8 cube certificates',
            'arithmetic_exact_rational_samples':samples,
            'signed_1x1_and_coupled_2x2_pivot_cases':pivots,
            'adaptive_precision_64_128_256_checks_passed':True,
            'three_exact_singular_inputs_rejected':True,
            'actual_cap8_four_probe_comparisons':probes,
            'all_checks_passed':True,
            'scope':'Technical certification of the supplied rational matrices. The physical operator comparison, radical error bounds, and complete-basis arguments are separate audited steps.'}
    path=ROOT/'outputs/SU2_Interval_Inertia_Independent_Audit.json'
    path.write_text(json.dumps(output,indent=2)+'\n',encoding='utf-8')
    print(path,flush=True)


if __name__=='__main__':
    main()
