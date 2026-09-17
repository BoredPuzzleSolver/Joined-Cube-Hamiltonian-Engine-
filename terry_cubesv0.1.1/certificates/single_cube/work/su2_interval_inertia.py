"""Rigorous symmetric-matrix inertia with fixed-point interval elimination.

Each interval [lo,hi] represents [lo/2**bits,hi/2**bits]. All endpoints and
arithmetic operations use Python integers, with mathematical floor/ceiling
rounding. Fraction input is enclosed exactly. No floating-point estimate is
used to decide signs. At each elimination step, intervals enclose every
entry of the exact current Schur complement. A pivot sign is counted only
when its interval excludes zero; congruence then preserves the remaining
inertia. A certified indefinite 2x2 pivot handles zero diagonal pivots.

An interval containing zero is uncertainty, not a certified zero eigenvalue.
If no pivot can be certified, UncertifiedInertia is raised. Retry with more
bits; singular matrices can remain undecidable and must not be certified.
"""
from __future__ import annotations

from fractions import Fraction
import json
import math
import random
import time
from pathlib import Path


class UncertifiedInertia(ArithmeticError):
    def __init__(self, message, metadata):
        super().__init__(message)
        self.metadata = metadata


def _ceildiv(a,b):
    return -((-a)//b)


def _mul(a,b,scale):
    products = (a[0]*b[0],a[0]*b[1],a[1]*b[0],a[1]*b[1])
    return min(products)//scale,_ceildiv(max(products),scale)


def _add(a,b):
    return a[0]+b[0],a[1]+b[1]


def _sub(a,b):
    return a[0]-b[1],a[1]-b[0]


def _div(a,b,scale):
    if b[0]<=0<=b[1]:
        raise ZeroDivisionError('interval denominator contains zero')
    numerators = (a[0]*scale,a[1]*scale)
    return (min(n//d for n in numerators for d in b),
            max(_ceildiv(n,d) for n in numerators for d in b))


def _product_quotient(a,b,pivot):
    """Outward rounding of a*b/pivot in fixed-point units, in one step."""
    if pivot[0]<=0<=pivot[1]:
        raise ZeroDivisionError('interval pivot contains zero')
    if a==(0,0) or b==(0,0):
        return 0,0
    products = (a[0]*b[0],a[0]*b[1],a[1]*b[0],a[1]*b[1])
    low,high = min(products),max(products)
    # Three fixed-point factors cancel to leave integers products/pivot.
    return (min(n//d for n in (low,high) for d in pivot),
            max(_ceildiv(n,d) for n in (low,high) for d in pivot))


def _enclose(value,scale):
    value = Fraction(value)
    numerator = value.numerator*scale
    return numerator//value.denominator,_ceildiv(numerator,value.denominator)


def interval_inertia(matrix,bits=128,*,pivoting=True):
    """Return certified positive/negative/zero counts and method metadata.

    Input must be exactly symmetric, preferably with Fraction/int entries.
    If float entries are supplied their exact binary rational values are used;
    this does not certify some different intended real matrix. For a rigorous
    enclosure supplied directly, use interval_inertia_from_bounds instead.
    """
    if not isinstance(bits,int) or bits<8:
        raise ValueError('bits must be an integer at least 8')
    n = len(matrix)
    if any(len(row)!=n for row in matrix):
        raise ValueError('matrix must be square')
    scale = 1<<bits
    lower = [[0]*n for _ in range(n)]
    upper = [[0]*n for _ in range(n)]
    for i in range(n):
        for j in range(i+1):
            if matrix[i][j]!=matrix[j][i]:
                raise ValueError('matrix must be exactly symmetric')
            lo,hi = _enclose(matrix[i][j],scale)
            lower[i][j]=lower[j][i]=lo
            upper[i][j]=upper[j][i]=hi
    return interval_inertia_from_bounds(lower,upper,bits,pivoting=pivoting)


def interval_inertia_from_bounds(lower,upper,bits=128,*,pivoting=True):
    """Certify inertia of every symmetric matrix in supplied integer bounds.

    lower/upper entries are integers scaled by 2**bits. Inputs are copied.
    Symmetric copies describe the same matrix entry, not independent entries;
    interval propagation can overestimate but never omit its possible values.
    """
    start = time.perf_counter()
    n = len(lower)
    if len(upper)!=n or any(len(r)!=n for r in lower+upper):
        raise ValueError('bounds must be square matrices of the same shape')
    if not isinstance(bits,int) or bits<8:
        raise ValueError('bits must be an integer at least 8')
    lo,hi = [list(row) for row in lower],[list(row) for row in upper]
    for i in range(n):
        for j in range(n):
            if not isinstance(lo[i][j],int) or not isinstance(hi[i][j],int):
                raise TypeError('scaled interval bounds must be Python integers')
            if lo[i][j]>hi[i][j] or lo[i][j]!=lo[j][i] or hi[i][j]!=hi[j][i]:
                raise ValueError('invalid or asymmetric interval bounds')
    scale = 1<<bits
    permutation = list(range(n))
    positive=negative=two_by_two=0
    separations=[]
    def swap(a,b):
        if a==b:return
        lo[a],lo[b]=lo[b],lo[a]
        hi[a],hi[b]=hi[b],hi[a]
        for row in lo:row[a],row[b]=row[b],row[a]
        for row in hi:row[a],row[b]=row[b],row[a]
        permutation[a],permutation[b]=permutation[b],permutation[a]
    k=0
    while k<n:
        candidates=[]
        for i in range(k,n) if pivoting else [k]:
            separation = lo[i][i] if lo[i][i]>0 else -hi[i][i] if hi[i][i]<0 else 0
            if separation:candidates.append((separation,i))
        if candidates:
            _,best=max(candidates)
            swap(k,best)
            pivot=(lo[k][k],hi[k][k])
            if pivot[0]>0:positive+=1;separations.append(pivot[0])
            else:negative+=1;separations.append(-pivot[1])
            for i in range(k+1,n):
                a=(lo[i][k],hi[i][k])
                if a==(0,0):continue
                for j in range(k+1,i+1):
                    b=(lo[j][k],hi[j][k])
                    if b==(0,0):continue
                    correction=_product_quotient(a,b,pivot)
                    newlo=lo[i][j]-correction[1]
                    newhi=hi[i][j]-correction[0]
                    lo[i][j]=lo[j][i]=newlo
                    hi[i][j]=hi[j][i]=newhi
            k+=1
            continue
        # A real symmetric 2x2 block with negative determinant has inertia(1,1).
        selected=None
        for i in range(k,n):
            for j in range(i+1,n):
                aa=(lo[i][i],hi[i][i]);bb=(lo[j][j],hi[j][j]);ab=(lo[i][j],hi[i][j])
                det=_sub(_mul(aa,bb,scale),_mul(ab,ab,scale))
                if det[1]<0:
                    candidate=(-det[1],i,j)
                    if selected is None or candidate>selected:selected=candidate
        if selected is None:
            metadata={'certified':False,'bits':bits,'dimension':n,'eliminated':k,
                      'positive_confirmed':positive,'negative_confirmed':negative,
                      'unresolved_dimension':n-k}
            raise UncertifiedInertia('No pivot sign can be certified; increase bits or handle an exact nullspace separately.',metadata)
        _,i,j=selected
        swap(k,i)
        if j==k:j=i
        swap(k+1,j)
        aa=(lo[k][k],hi[k][k]);bb=(lo[k+1][k+1],hi[k+1][k+1]);ab=(lo[k][k+1],hi[k][k+1])
        det=_sub(_mul(aa,bb,scale),_mul(ab,ab,scale))
        assert det[1]<0
        for i in range(k+2,n):
            ai=(lo[i][k],hi[i][k]);bi=(lo[i][k+1],hi[i][k+1])
            for j in range(k+2,i+1):
                aj=(lo[j][k],hi[j][k]);bj=(lo[j][k+1],hi[j][k+1])
                numerator=_add(_mul(bb,_mul(ai,aj,scale),scale),_mul(aa,_mul(bi,bj,scale),scale))
                numerator=_sub(numerator,_mul(ab,_add(_mul(ai,bj,scale),_mul(bi,aj,scale)),scale))
                correction=_div(numerator,det,scale)
                lo[i][j]=lo[j][i]=lo[i][j]-correction[1]
                hi[i][j]=hi[j][i]=hi[i][j]-correction[0]
        positive+=1;negative+=1;two_by_two+=1
        k+=2
    return {'positive':positive,'negative':negative,'zero':0,'certified':True,
            'method':'Outward-rounded integer fixed-point interval Schur elimination with certified 1x1/2x2 pivots',
            'bits':bits,'dimension':n,'two_by_two_pivots':two_by_two,
            'permutation':permutation,
            'minimum_1x1_pivot_separation_exact':str(Fraction(min(separations),scale)) if separations else None,
            'elapsed_seconds':time.perf_counter()-start,
            'proof':'Exact input lies inside each entry interval. Interval arithmetic encloses every subsequent exact Schur complement. Certified pivot inertia adds under congruence; no unresolved pivot is accepted.'}


def adaptive_interval_inertia(matrix,bits=64,max_bits=1024,*,pivoting=True):
    attempts=[]
    while bits<=max_bits:
        try:
            result=interval_inertia(matrix,bits,pivoting=pivoting)
            result['bits_attempted']=attempts+[bits]
            return result
        except UncertifiedInertia:
            attempts.append(bits)
            bits*=2
    raise UncertifiedInertia('All requested interval precisions were inconclusive.',{'bits_attempted':attempts,'certified':False})


def self_checks():
    """Independent exact-congruence, Bareiss and actual cube-Schur checks."""
    from su2_cutoff_bounds import rational_inertia
    rng=random.Random(81723)
    arithmetic_cases=0
    for _ in range(500):
        a=tuple(sorted([rng.randint(-10000,10000),rng.randint(-10000,10000)]))
        b=tuple(sorted([rng.randint(-10000,10000),rng.randint(-10000,10000)]))
        sign=rng.choice([-1,1])
        p=tuple(sorted([sign*rng.randint(1,10000),sign*rng.randint(1,10000)]))
        ml,mh=_mul(a,b,256)
        assert all(ml<=Fraction(x*y,256)<=mh for x in a for y in b)
        ql,qh=_div(a,p,256)
        assert all(ql<=Fraction(x*256,y)<=qh for x in a for y in p)
        pl,ph=_product_quotient(a,b,p)
        assert all(pl<=Fraction(x*y,z)<=ph for x in a for y in b for z in p)
        arithmetic_cases+=1
    random_cases=0
    for n in range(1,15):
        for _ in range(5):
            r=[[Fraction(int(i==j) if i>=j else rng.randint(-3,3)) for j in range(n)] for i in range(n)]
            signs=[rng.choice([-1,1]) for _ in range(n)]
            d=[Fraction(signs[i]*(i+1),i+2) for i in range(n)]
            a=[[sum(r[k][i]*d[k]*r[k][j] for k in range(n)) for j in range(n)] for i in range(n)]
            got=interval_inertia(a)
            exact=rational_inertia(a)
            assert got['positive']==signs.count(1)==exact['positive']
            assert got['negative']==signs.count(-1)==exact['negative']
            random_cases+=1
    block_tests=[]
    for n in [2,4,8]:
        a=[[Fraction(0) for _ in range(n)] for _ in range(n)]
        for i in range(0,n,2):a[i][i+1]=a[i+1][i]=Fraction(i+1,3)
        got=interval_inertia(a)
        assert got['positive']==got['negative']==n//2
        block_tests.append({'dimension':n,'positive':got['positive'],'negative':got['negative']})
    coupled=[[0,1,2],[1,0,3],[2,3,0]]
    got=interval_inertia(coupled)
    exact=rational_inertia(coupled)
    assert all(got[key]==exact[key] for key in ['positive','negative','zero'])
    try:interval_inertia([[Fraction(0)]])
    except UncertifiedInertia:pass
    else:raise AssertionError('singular input was accepted')
    from su2_cube_model import cube_basis
    from su2_cube_rational_certificate import comparison_matrices
    states,diagonal,matrix,eta_a,eta_l=comparison_matrices(cube_basis(1))
    p=32
    probes=[]
    specifications=[('L0',Fraction('5.30400030')+eta_l,True),
                    ('L1',Fraction('6.20980338')+eta_l,True),
                    ('A0',Fraction('5.54178379')-eta_a,False),
                    ('A1',Fraction('8.77237429')-eta_a,False)]
    for name,probe,is_lower in specifications:
        a=[[matrix[i].get(j,Fraction(0)) if i!=j else diagonal[i]-probe for j in range(p)] for i in range(p)]
        if is_lower:
            for q in range(p,len(diagonal)):
                denominator=diagonal[q]-probe
                assert denominator>0
                neighbors=list(matrix[q])
                for offset,i in enumerate(neighbors):
                    for j in neighbors[offset:]:
                        correction=matrix[q][i]*matrix[q][j]/denominator
                        a[i][j]-=correction
                        if i!=j:a[j][i]-=correction
        start=time.perf_counter();got=interval_inertia(a);interval_seconds=time.perf_counter()-start
        start=time.perf_counter();exact=rational_inertia(a);bareiss_seconds=time.perf_counter()-start
        assert all(got[key]==exact[key] for key in ['positive','negative','zero'])
        probes.append({'name':name,'probe_exact':str(probe),'positive':got['positive'],'negative':got['negative'],
                       'interval_seconds':interval_seconds,'Bareiss_seconds':bareiss_seconds,'counts_agree':True})
    return {'method':'Rigorous fixed-point interval inertia with Python integers',
            'interval_endpoint_property_cases':arithmetic_cases,
            'exact_congruence_and_Bareiss_comparisons':random_cases,
            'random_matrix_dimensions':list(range(1,15)),
            'zero_diagonal_block_checks':block_tests,'coupled_2x2_pivot_check':True,
            'singular_matrix_rejection':True,'actual_cube_P32_checks':probes,
            'all_passed':True,
            'scope':'This verifies inertia of the supplied exact rational or interval-enclosed symmetric matrix. A field-theory certificate additionally needs its separate operator comparison, basis completeness and matrix-entry error bounds.'}


if __name__=='__main__':
    result=self_checks()
    output=Path(__file__).resolve().parents[1]/'outputs'/'SU2_Interval_Inertia_Validation.json'
    output.write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8')
    print(json.dumps(result,indent=2))
