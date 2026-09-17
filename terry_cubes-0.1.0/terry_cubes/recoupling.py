"""Exact recoupling coefficients for the open two-square SU(2) ladder.

Public state labels are DOUBLED spins (integers), ordered (2*jL,2*jR,2*jM).
The physical basis consists of triangle-admissible triples, one state each.

Basis derivation:
  f_(L,R,M)(A,B) = sqrt(dL*dR/dM) Tr[P_M (D^L(A) tensor D^R(B))].
Schur orthogonality makes these functions orthonormal. P_M is the orthogonal
projector onto total spin M in V_L tensor V_R; d_j=2j+1. For L=R=1/2,M=0,
the basis function is Tr(AB^{-1}). The orientation change B -> B^{-1}
turns this into Tr(AB), preserving Tr B and product Haar measure.

Multiply f by Tr(D^{1/2}(A))/2. In V_{1/2} tensor V_L tensor V_R,
recouple ((1/2,L)->L', R)->M' versus (1/2,(L,R)->M)->M'. The squared
recoupling overlap is dM*dL' * {L R M; M' 1/2 L'}^2. Dividing by the
normalizations of f and f' gives
  <L',R,M'|W_A|L,R,M> =
  (1/2) sqrt(dL*dM*dL'*dM') * {L R M; M' 1/2 L'}^2.
All coefficients are nonnegative in this coupled-projector phase convention.
W_B follows by swapping L and R. Selection rules require L',M'=L,M +/-1/2.

Squared 6j symbols are computed with exact Fraction arithmetic using the
finite Racah sum. Only the final dimension square root uses floating point.
This module generates the orthogonal compression P_J H P_J of the original
gauge-invariant Hamiltonian, not a modified cyclic cutoff or quantum group.

Primary context: Burgio et al., hep-lat/9906036 and hep-lat/9911019 describe
orthonormal physical spin-network bases and plaquette recoupling. The compact
theta-graph specialization and projector derivation above are our derivation.
"""


from __future__ import annotations


from functools import lru_cache


from fractions import Fraction


import itertools


import json


import math


import numpy as np


SOURCES = [
    'https://arxiv.org/abs/hep-lat/9906036',
    'https://arxiv.org/abs/hep-lat/9911019',
    'https://arxiv.org/abs/2307.11829',
]


def admissible(a: int, b: int, c: int) -> bool:
    """SU(2) triangle rule for doubled nonnegative integer spins."""
    return min(a, b, c) >= 0 and abs(a-b) <= c <= a+b and (a+b+c) % 2 == 0


def basis_twice(max_twice: int) -> list[tuple[int, int, int]]:
    if int(max_twice) != max_twice or max_twice < 0:
        raise ValueError('max_twice must be a nonnegative integer')
    return [s for s in itertools.product(range(int(max_twice)+1), repeat=3)
            if admissible(*s)]


def electric_diagonal(states, kappa: float = 1.0) -> np.ndarray:
    """Three links in each outer path and one link in the shared middle."""
    return kappa*np.array([(3*a*(a+2)+3*b*(b+2)+c*(c+2))/4
                           for a, b, c in states], dtype=float)


@lru_cache(maxsize=None)
def triangle_delta_squared(a: int, b: int, c: int) -> Fraction:
    if not admissible(a, b, c):
        return Fraction(0)
    f = math.factorial
    return Fraction(f((a+b-c)//2)*f((a-b+c)//2)*f((-a+b+c)//2),
                    f((a+b+c)//2+1))


@lru_cache(maxsize=None)
def wigner_6j_squared(a: int, b: int, c: int,
                      d: int, e: int, f: int) -> Fraction:
    """Return {a/2 b/2 c/2; d/2 e/2 f/2} squared, exactly."""
    deltas = [triangle_delta_squared(*t)
              for t in [(a,b,c), (a,e,f), (d,b,f), (d,e,c)]]
    if not all(deltas):
        return Fraction(0)
    lower = [(a+b+c)//2, (a+e+f)//2, (d+b+f)//2, (d+e+c)//2]
    upper = [(a+b+d+e)//2, (a+c+d+f)//2, (b+c+e+f)//2]
    racah_sum = Fraction(0)
    for z in range(max(lower), min(upper)+1):
        denominator = math.prod(math.factorial(z-x) for x in lower)
        denominator *= math.prod(math.factorial(x-z) for x in upper)
        racah_sum += Fraction((-1)**z*math.factorial(z+1), denominator)
    return math.prod(deltas)*racah_sum*racah_sum


def wilson_element(initial, final, which: str = 'A') -> float:
    """Normalized fundamental trace W=Tr(U)/2 on the chosen plaquette."""
    a, b, c = initial
    ap, bp, cp = final
    if which.upper() == 'B':
        a, b, ap, bp = b, a, bp, ap
    elif which.upper() != 'A':
        raise ValueError('which must be A or B')
    if b != bp or abs(ap-a) != 1 or abs(cp-c) != 1:
        return 0.0
    if not admissible(a,b,c) or not admissible(ap,bp,cp):
        return 0.0
    square = wigner_6j_squared(a,b,c,cp,1,ap)
    return 0.5*math.sqrt((a+1)*(c+1)*(ap+1)*(cp+1))*float(square)


def wilson_matrices(max_twice: int):
    """Return (lexicographic doubled-spin states, W_A, W_B)."""
    states = basis_twice(max_twice)
    lookup = {s: i for i, s in enumerate(states)}
    matrices = []
    for which in ('A', 'B'):
        w = np.zeros((len(states), len(states)))
        for i, s in enumerate(states):
            active = 0 if which == 'A' else 1
            for delta, delta_middle in itertools.product((-1,1), repeat=2):
                target = list(s)
                target[active] += delta
                target[2] += delta_middle
                target = tuple(target)
                j = lookup.get(target)
                if j is not None:
                    w[j, i] = wilson_element(s, target, which)
        matrices.append(w)
    return states, *matrices


def hamiltonian(max_twice: int, kappa: float = 1.0, nu: float = 1.0):
    """Return states,H for H=kappa sum E^2+nu(2-W_A-W_B)."""
    states, wa, wb = wilson_matrices(max_twice)
    h = np.diag(electric_diagonal(states, kappa)+2*nu)-nu*(wa+wb)
    return states, h


def exact_congruence(max_twice: int, kappa=Fraction(1),
                     nu=Fraction(1), shift=Fraction(0)):
    """Return states and exact rational D^-1 (H-shift*I) D^-1.

    D_ii=sqrt(K_i), K_i=dL*dR*dM. Congruence preserves inertia, not
    eigenvalues. For unchanged R, WA_ij/(sqrt(K_i*K_j)) equals
    (6j)^2/(2*dR), and similarly with L for WB. Integer, string or Fraction
    parameters avoid accidentally specifying a binary floating-point rational.
    """
    kappa, nu, shift = map(Fraction, (kappa,nu,shift))
    states = basis_twice(max_twice)
    lookup = {s:i for i,s in enumerate(states)}
    n = len(states)
    out = [[Fraction(0) for _ in range(n)] for _ in range(n)]
    for i,(a,b,c) in enumerate(states):
        kinetic = kappa*Fraction(3*a*(a+2)+3*b*(b+2)+c*(c+2),4)
        out[i][i] = (kinetic+2*nu-shift)/((a+1)*(b+1)*(c+1))
        for da,dc in itertools.product((-1,1),repeat=2):
            target = (a+da,b,c+dc)
            j = lookup.get(target)
            if j is not None:
                out[j][i] -= nu*wigner_6j_squared(a,b,c,c+dc,1,a+da)/(2*(b+1))
        for db,dc in itertools.product((-1,1),repeat=2):
            target = (a,b+db,c+dc)
            j = lookup.get(target)
            if j is not None:
                out[j][i] -= nu*wigner_6j_squared(b,a,c,c+dc,1,b+db)/(2*(a+1))
    return states,out


@lru_cache(maxsize=None)
def spin_operators(twice_j: int):
    j = twice_j/2
    m = np.arange(-j, j+1)
    raising = np.zeros((twice_j+1, twice_j+1), dtype=complex)
    for k, mk in enumerate(m[:-1]):
        raising[k+1,k] = math.sqrt((j-mk)*(j+mk+1))
    return ((raising+raising.T)/2, (raising-raising.T)/(2j), np.diag(m))


def projector(casimir: np.ndarray, twice_j: int):
    values, vectors = np.linalg.eigh(casimir)
    j = twice_j/2
    chosen = vectors[:, abs(values-j*(j+1)) < 1e-8]
    return chosen @ chosen.conj().T


@lru_cache(maxsize=None)
def triple_projectors(a: int, b: int):
    """Independent angular momentum matrices on V_half tensor V_a tensor V_b."""
    dims = [2, a+1, b+1]
    site_j = []
    for site, twice in enumerate([1,a,b]):
        ops = []
        for op in spin_operators(twice):
            mats = [np.eye(d) for d in dims]
            mats[site] = op
            ops.append(np.kron(np.kron(mats[0],mats[1]),mats[2]))
        site_j.append(ops)
    def casimir(sites):
        out = np.zeros((math.prod(dims),)*2, dtype=complex)
        for k in range(3):
            total = sum(site_j[site][k] for site in sites)
            out += total @ total
        return out
    hl = casimir([0,1])
    lr = casimir([1,2])
    all_spins = casimir([0,1,2])
    p_hl = {ap: projector(hl,ap) for ap in range(abs(a-1),a+2,2)}
    p_lr = {c: projector(lr,c) for c in range(abs(a-b),a+b+1,2)}
    p_total = {cp: projector(all_spins,cp) for cp in range((1+a+b)%2,1+a+b+1,2)}
    return p_hl,p_lr,p_total


def independent_projector_element(initial, final):
    """No 6j symbols: trace of three angular-momentum spectral projectors."""
    a,b,c = initial
    ap,bp,cp = final
    if b != bp or abs(ap-a) != 1 or abs(cp-c) != 1:
        return 0.0
    p_hl,p_lr,p_total = triple_projectors(a,b)
    overlap = float(np.trace(p_hl[ap] @ p_lr[c] @ p_total[cp]).real)/(cp+1)
    norm_ratio = math.sqrt((a+1)*(cp+1)/((c+1)*(ap+1)))
    return 0.5*norm_ratio*overlap


def self_checks() -> dict:
    assert wigner_6j_squared(0,0,0,1,1,1) == Fraction(1,2)
    assert wigner_6j_squared(0,1,1,0,1,1) == Fraction(1,4)
    states,wa,wb = wilson_matrices(1)
    four_order = [(0,0,0),(1,0,1),(0,1,1),(1,1,0)]
    ids = [states.index(s) for s in four_order]
    expected_a = np.array([[0,.5,0,0],[.5,0,0,0],[0,0,0,.25],[0,0,.25,0]])
    expected_b = np.array([[0,0,.5,0],[0,0,0,.25],[.5,0,0,0],[0,.25,0,0]])
    four_error = max(float(np.max(abs(wa[np.ix_(ids,ids)]-expected_a))),
                     float(np.max(abs(wb[np.ix_(ids,ids)]-expected_b))))
    assert four_error < 1e-15
    states,wa,wb = wilson_matrices(3)
    projector_errors = []
    for i,s in enumerate(states):
        for j,t in enumerate(states):
            if s[1] == t[1] and abs(s[0]-t[0]) == 1 and abs(s[2]-t[2]) == 1:
                projector_errors.append(abs(wa[j,i]-independent_projector_element(s,t)))
    assert max(projector_errors) < 1e-12
    rational_states,rational_matrix = exact_congruence(3,nu=Fraction(3,2),shift=Fraction(1,7))
    hs,h = hamiltonian(3,nu=1.5)
    assert hs == rational_states
    d = np.sqrt([math.prod(x+1 for x in s) for s in hs])
    expected_congruence = (h-np.eye(len(h))/7)/np.outer(d,d)
    actual_congruence = np.array(rational_matrix,dtype=float)
    congruence_error = float(np.max(abs(expected_congruence-actual_congruence)))
    assert congruence_error < 1e-14
    summaries = []
    for max_twice in [1,2,3,4,6,8,10,12]:
        states,wa,wb = wilson_matrices(max_twice)
        symmetry = max(float(np.max(abs(wa-wa.T))),float(np.max(abs(wb-wb.T))))
        assert symmetry < 1e-14
        swap = [states.index((b,a,c)) for a,b,c in states]
        swap_error = float(np.max(abs(wa[np.ix_(swap,swap)]-wb)))
        assert swap_error < 1e-14
        r_zero = [i for i,(a,b,c) in enumerate(states) if b == 0]
        a_zero = [i for i,(a,b,c) in enumerate(states) if a == 0]
        character_expected = (np.diag(np.full(max_twice,.5),1)+
                              np.diag(np.full(max_twice,.5),-1))
        assert np.max(abs(wa[np.ix_(r_zero,r_zero)]-character_expected)) < 1e-14
        assert np.max(abs(wb[np.ix_(a_zero,a_zero)]-character_expected)) < 1e-14
        row = {'jmax':max_twice/2,'dimension':len(states),'WA_nonzero':int(np.count_nonzero(wa)),
               'WB_nonzero':int(np.count_nonzero(wb)),'Hermiticity_residual':symmetry,
               'A_B_exchange_residual':swap_error}
        if max_twice <= 6:
            norms = [float(np.max(abs(np.linalg.eigvalsh(w)))) for w in (wa,wb)]
            assert max(norms) <= 1+1e-12
            row['Wilson_operator_norms'] = norms
        summaries.append(row)
    return {'sources':SOURCES,'four_state_matrix_residual':four_error,
            'independent_projector_entries_checked':len(projector_errors),
            'independent_projector_max_residual':max(projector_errors),
            'rational_congruence_vs_float_residual':congruence_error,
            'cutoff_matrix_checks':summaries}


# Original signed trivalent recoupling convention from su2_cube_model.py.


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


# Original orientation-aware invariant contractions from su2_oriented_vertex.py.


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


def oriented_self_checks(max_twice=6):
    from .recoupling import local_vertex_factor as reference_vertex
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
