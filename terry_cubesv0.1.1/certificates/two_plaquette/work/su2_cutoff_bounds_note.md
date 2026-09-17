# Certified spin-cutoff bounds for the shared-edge SU(2) patch

## Result and scope

The two-square graph supports a rigorous comparison between its **untruncated gauge-invariant Hamiltonian** and two finite matrices. This bounds both lowest energies individually and therefore bounds their difference. Agreement between successive numerical gaps alone would not establish those bounds.

For \(\kappa=1\), the exact rational checks in verify_su2_cutoff_certificates.py establish:

| Magnetic coupling \(\nu\) | Retained cutoff \(J\) | Retained / auxiliary dimensions | Full-patch gap interval |
|---:|---:|---:|---:|
| 1 | 2 | 42 / 69 | \([3.09195898,\;3.09195907]\) |
| 4 | 2 | 42 / 69 | \([4.20913941,\;4.20995910]\) |
| 4 | 3 | 106 / 154 | \([4.20942742,\;4.20942757]\) |

These concern the fixed open graph with two plaquettes, seven links, no external charges, and **all** allowed link spins in the target Hilbert space. They do not concern an infinite spatial lattice or a continuum limit. Energies use the stated Hamiltonian normalization.

## 1. Hamiltonian and omitted electric threshold

The physical spin-network labels are \((j_L,j_R,j_M)\), obeying the triangle inequalities and integer total spin. The three outer links on each side share their spin. Write

\[
H=T+V,\qquad
T=\kappa\{3j_L(j_L+1)+3j_R(j_R+1)+j_M(j_M+1)\},
\qquad V=\nu(2-W_A-W_B).
\]

Here \(W_A=\operatorname{Tr}U_A/2\) and similarly for \(B\). Each is multiplication by a real function in \([-1,1]\), so \(0\leq V\leq4\nu I\) when \(\nu\geq0\). The electric operator has finite-dimensional sublevel spaces and tends to infinity with the largest spin. Consequently \(H\) is self-adjoint with compact resolvent for \(\kappa>0\), and its eigenvalues can be ordered with multiplicity.

Let \(P\) retain all three spins at most \(J\). Set \(K=2J\). The exact electric minimum outside \(P\) is

\[
t_Q=\frac{\kappa}{4}
\{3a(a+2)+3b(b+2)+m(m+2)\},
\quad m=K+1,\quad
a=\lfloor m/2\rfloor,\quad b=\lceil m/2\rceil.
\]

To prove the formula, use doubled integer labels. For any admissible triple, putting its largest label on the shared edge cannot increase the energy: that edge has coefficient 1 and the others have coefficient 3. At fixed largest label \(m\), the remaining sum is at least \(m\) and has the same parity. Convexity minimizes the two remaining Casimirs at the smallest admissible sum, \(a+b=m\), split as evenly as possible. This minimum increases strictly with \(m\). Therefore the first omitted largest label is sufficient; no farther state can have smaller electric energy.

For \(J=1/2,1,3/2,2\), this gives \(t_Q/\kappa=6.5,12,18,26\).

## 2. A finite lower comparison retaining the actual boundary couplings

Decompose \(Q=I-P\) into its complete first spin layer \(Q_1\), with largest spin \(J+1/2\), and the remaining states \(R\). A fundamental Wilson loop changes each affected spin by \(1/2\). Thus \(PHR=0\): every coupling out of \(P\) is contained in \(Q_1\).

With \(A=PHP\) and \(B=PHQ_1\), form the finite auxiliary matrix

\[
L=
\begin{pmatrix}
A&B\\
B^\dagger&T_{Q_1}
\end{pmatrix}.
\]

This is obtained from the exact Hamiltonian at cutoff \(J+1/2\) by replacing the entire \(Q_1\)-to-\(Q_1\) block by its electric diagonal. Retain all \(P\)-to-\(Q_1\) entries.

In the full Hilbert space,

\[
H-(L\oplus T_R)=0_P\oplus QVQ\ \geq\ 0.
\]

In particular, this comparison drops a **positive compressed operator**, not selected magnetic entries in isolation. Dropping arbitrary off-diagonal entries would not justify an operator inequality.

Let \(a_j\) be the ordered eigenvalues of \(A\), \(\ell_j\) those of \(L\), and \(t_R\) the exact electric threshold outside cutoff \(J+1/2\). Min–max gives, for the retained levels under discussion,

\[
\min(\ell_j,t_R)\ \leq\ E_j(H)\ \leq\ a_j.
\]

The minimum accounts for every omitted state beyond the auxiliary matrix. The benchmarks are far below \(t_R\), which is \(34.5\) at \(J=2\) and \(56\) at \(J=3\).

If \(L_j\leq E_j\leq U_j\) are individual level bounds, then

\[
L_1-U_0\ \leq\ E_1-E_0\ \leq\ U_1-L_0.
\]

The difference of two variational upper levels is not used as a gap bound.

### Relation to the Schur-complement approach

For a trial energy \(z<t_Q\), the exact effective operator on \(P\), now using \(B=PHQ\), is

\[
A-z-B(QHQ-z)^{-1}B^\dagger
\ \geq\
A-z-\frac{BB^\dagger}{t_Q-z}.
\]

Replacing the omitted block by \(t_QI\) therefore gives a valid, weaker lower comparison. Keeping its individual electric energies as above improves this without repeated nonlinear eigenvalue solves. The coarse norm-only alternative is

\[
E_j\geq
\frac{a_j+t_Q-\sqrt{(t_Q-a_j)^2+4\beta^2}}{2},
\qquad \beta\geq\|PHQ\|,\quad \beta=2\nu\ \text{is sufficient}.
\]

The auxiliary matrix is usually much tighter because it retains which boundary states actually couple to the low-energy states.

## 3. Exact rational certificates, without a floating-point allowance

The recoupling implementation provides an exact rational congruence. For a state with doubled labels \((a,b,c)\), define \(w=(a+1)(b+1)(c+1)\) and \(D=\operatorname{diag}(\sqrt w)\). Although the Hamiltonian contains square roots, \(D^{-1}(H-zI)D^{-1}\) has rational entries at rational \(\kappa,\nu,z\). Its diagonal is \((T+2\nu-z)/w\). For a left-loop transition, its off-diagonal entry is

\[
-\frac{\nu}{2(b+1)}
\left\{\begin{matrix}
j_L&j_R&j_M\\
j'_M&1/2&j'_L
\end{matrix}\right\}^{\!2},
\]

with the analogous exchange of left and right labels for the right loop. The squared \(6j\) symbol is computed as an exact fraction. For \(L\), the omitted diagonal becomes \(T/w\), and all off-diagonal entries inside that omitted block become zero.

Congruence preserves inertia: the counts of negative, zero and positive eigenvalues. After a positive common denominator is cleared, symmetric Bareiss elimination uses only integers and exact divisions. Its successive LDL pivot signs determine those counts. Symmetric permutations and integer congruences handle vanishing diagonal pivots; a zero residual block gives the exact nullity.

For each level \(j\), a rational lower probe is certified when \(L-zI\) has at most \(j\) negative eigenvalues and \(z\leq t_R\). A rational upper probe is certified when \(A-zI\) has at least \(j+1\) negative eigenvalues. Numerical eigensolves are useful for selecting probes, but do not enter these logical tests.

The fixed probe intervals are:

| \(\nu\) | \(J\) | \(E_0\) interval | \(E_1\) interval |
|---:|---:|---:|---:|
| 1 | 2 | \([1.83601099,\;1.83601102]\) | \([4.92797000,\;4.92797006]\) |
| 4 | 2 | \([5.79687356,\;5.79692505]\) | \([10.00606446,\;10.00683266]\) |
| 4 | 3 | \([5.79688776,\;5.79688779]\) | \([10.00631521,\;10.00631533]\) |

Every displayed endpoint is interpreted as an exact decimal rational. The JSON records each fraction and the exact inertia counts. These are certificates for the stated mathematical model, conditional on the physical basis and recoupling identities defining that model; integer arithmetic does not independently prove those identities.

## 4. Reproduction

Run:

    python work/verify_su2_cutoff_certificates.py

It imports su2_recoupling.py and su2_cutoff_bounds.py, checks the exact inertia routine on singular and indefinite examples, reconstructs both matrices, checks the four rational probes for each benchmark, and writes outputs/Shared_SU2_Cutoff_Certificate_Checks.json.

The numerical_cutoff_bounds function remains available for wider sweeps. Its output explicitly says that floating-point evaluations alone are not machine-certified. It preserves raw endpoint ordering and flags widths near numerical resolution. That diagnostic is a heuristic, not an error bound. Use certify_energy_probes when an exact certificate is required.
