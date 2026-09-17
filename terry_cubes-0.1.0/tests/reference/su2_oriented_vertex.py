"""Exact all-outgoing SU(2) invariant contractions for a fundamental loop.

Every amplitude is represented by (sign, Fraction squared magnitude). The
3j finite sum is DLMF 34.2.4, https://dlmf.nist.gov/34.2#E4. Clebsch--Gordan
coefficients couple j first and spin 1/2 second in Condon--Shortley phase.
No floating sign decisions are used. All nonzero summands in the local
contraction have rational ratios; this is checked by exact integer square
roots before addition, and an unsupported radical sum fails explicitly.

The normalized real invariant is the standard 3j tensor in a fixed order.
For loop legs f,b, contract old and new invariant tensors against
C(j_f,1/2;J_f), C(j_b,1/2;J_b), and epsilon=[[0,1],[-1,0]]. A Wilson loop
coefficient is one half, times (-1) per edge traversed opposite its chosen
orientation, times the product of these vertex contractions and
sqrt((2j+1)/(2J+1)) on each active edge. This follows by multiplying
K_j(U)=D_j(U) epsilon_j and recoupling at its two endpoints.
"""
from fractions import Fraction
from functools import lru_cache
import math

from su2_recoupling import admissible, triangle_delta_squared, wigner_6j_squared


def parity(exponent):
    return -1 if exponent % 2 else 1


@lru_cache(maxsize=None)
def three_j_signed_square(a,b,c,ma,mb,mc):
    """All six arguments are doubled angular momenta / magnetic labels."""
    if not admissible(a,b,c) or ma+mb+mc:
        return 0,Fraction(0)
    if any(abs(m)>q or (q-m)%2 for q,m in ((a,ma),(b,mb),(c,mc))):
        return 0,Fraction(0)
    fact=math.factorial
    prefactor=triangle_delta_squared(a,b,c)
    for q,m in ((a,ma),(b,mb),(c,mc)):
        prefactor*=fact((q+m)//2)*fact((q-m)//2)
    upper=((a+b-c)//2,(a-ma)//2,(b+mb)//2)
    offset=((c-b+ma)//2,(c-a-mb)//2)
    total=Fraction(0)
    for z in range(max(0,-min(offset)),min(upper)+1):
        denominator=fact(z)
        denominator*=math.prod(fact(u-z) for u in upper)
        denominator*=math.prod(fact(v+z) for v in offset)
        total+=Fraction(parity(z),denominator)
    if not total:return 0,Fraction(0)
    sign=parity((a-b-mc)//2)*(1 if total>0 else -1)
    return sign,prefactor*total*total


@lru_cache(maxsize=None)
def half_cg_signed_square(q,m,h,r):
    """C(q/2,m/2; 1/2,h/2 | r/2,(m+h)/2), h=+/-1."""
    if h not in (-1,1) or abs(m)>q or (q-m)%2:
        return 0,Fraction(0)
    if abs(m+h)>r or r<0:return 0,Fraction(0)
    if r==q+1:
        numerator=q+h*m+2
        return 1,Fraction(numerator,2*(q+1))
    if r==q-1:
        numerator=q-h*m
        if not numerator:return 0,Fraction(0)
        return (-1 if h==1 else 1),Fraction(numerator,2*(q+1))
    return 0,Fraction(0)


def rational_sqrt(value):
    """Return the exact nonnegative rational root, or reject other radicals."""
    n=math.isqrt(value.numerator);d=math.isqrt(value.denominator)
    if n*n!=value.numerator or d*d!=value.denominator:
        raise ArithmeticError('Local invariant sum contains non-proportional radicals')
    return Fraction(n,d)


@lru_cache(maxsize=None)
def invariant_loop_contraction(spectator,forward,backward,new_forward,new_backward):
    """Exact tensor contraction in invariant order (spectator,forward,backward).

    Its squared magnitude is (2J_f+1)(2J_b+1) times the relevant 6j square.
    The explicit finite tensor sum fixes the sign for all edge orientations;
    it does not assume that the resolved graph is bipartite.
    """
    e,f,b,F,B=spectator,forward,backward,new_forward,new_backward
    six_square=wigner_6j_squared(e,f,b,1,B,F)
    if not six_square:return 0,Fraction(0)
    common=None;total=Fraction(0)
    for me in range(e,-e-1,-2):
        for mf in range(f,-f-1,-2):
            mb=-me-mf
            so,qo=three_j_signed_square(e,f,b,me,mf,mb)
            if not qo:continue
            for h in (1,-1):
                sn,qn=three_j_signed_square(e,F,B,me,mf+h,mb-h)
                sf,qf=half_cg_signed_square(f,mf,h,F)
                sb,qb=half_cg_signed_square(b,mb,-h,B)
                square=qo*qn*qf*qb
                if not square:continue
                # epsilon_(+,-)=+1; epsilon_(-,+)=-1.
                sign=so*sn*sf*sb*h
                if common is None:common=square
                total+=sign*rational_sqrt(square/common)
    result=(1 if total>0 else -1 if total<0 else 0),total*total*(common or 0)
    assert result[1]==(F+1)*(B+1)*six_square
    return result


@lru_cache(maxsize=None)
def oriented_local_vertex_factor(spectator,forward,backward,new_forward,new_backward):
    """Assign each forward edge's dimension ratio once to its initial vertex."""
    sign,square=invariant_loop_contraction(spectator,forward,backward,new_forward,new_backward)
    return sign,square*Fraction(forward+1,new_forward+1)


def self_checks(max_twice=6):
    from su2_cube_model import local_vertex_factor as reference_vertex
    cases=0
    for e in range(max_twice+1):
        for f in range(max_twice+1):
            for b in range(max_twice+1):
                if not admissible(e,f,b):continue
                # The normalized 3j invariant has exactly unit norm.
                norm=sum((three_j_signed_square(e,f,b,me,mf,-me-mf)[1]
                          for me in range(e,-e-1,-2)
                          for mf in range(f,-f-1,-2)),Fraction(0))
                assert norm==1
                for F in (f-1,f+1):
                    for B in (b-1,b+1):
                        if min(F,B)<0 or not admissible(e,F,B):continue
                        sign,square=invariant_loop_contraction(e,f,b,F,B)
                        old_sign,old_square=reference_vertex(e,f,b,F,B)
                        assert sign==old_sign*(F-f)*(B-b)
                        assert square*Fraction(f+1,F+1)==old_square*Fraction((f+1)*(B+1),(F+1)*(b+1))
                        cases+=1
    return {'maximum_doubled_old_spin':max_twice,'exact_vertex_contractions':cases,
            'unit_3j_norms':True,'squared_6j_identity':True,
            'all_radical_ratios_exact_rational_squares':True,
            'exact_local_sign_relation_to_reference_6j_formula':True,
            'exact_local_dimension_ratio_to_reference_6j_formula':True}


if __name__=='__main__':
    import json
    print(json.dumps(self_checks(),indent=2))
