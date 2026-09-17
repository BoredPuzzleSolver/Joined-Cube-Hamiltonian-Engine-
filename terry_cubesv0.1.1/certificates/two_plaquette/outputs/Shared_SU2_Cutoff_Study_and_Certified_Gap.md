# Two interacting SU(2) plaquettes: larger bases and a certified gap

**Research development note — 15 September 2026**

## Main result

For the open graph consisting of two squares sharing one link, with Gauss law at every vertex and the Hamiltonian defined below, the complete, untruncated physical Hilbert space has a gap satisfying

\[
\boxed{3.09195898\leq \Delta/\kappa\leq3.09195907}
\qquad (\nu/\kappa=1).
\]

These endpoints are exact terminating decimals established by rational arithmetic and an operator inequality controlling all omitted spins. The larger numerical calculation gives approximately **3.091959020089**. The additional numerical digits are estimates; the displayed interval is the certificate.

The earlier four-state result, 3.159602831021, changes when higher-spin states are included. At the same normalization, the certified answer is below π. This calculation supplies no universal value of the Yang–Mills mass gap or identification with a cosmological energy scale.

At a second coupling, the corresponding certificate is

\[
\boxed{4.20942742\leq \Delta/\kappa\leq4.20942757}
\qquad (\nu/\kappa=4).
\]

**Scope:** these are computer-assisted spectral bounds for a specified finite lattice graph, including its infinite tower of SU(2) representations. We have controlled the spin cutoff for these cases. We have not taken an infinite spatial volume or a continuum limit, and do not claim novelty or a solution of the Millennium problem.

## 1. The actual gauge model

The graph has seven links, six vertices and two plaquettes. Boundary vertices also obey Gauss law; there are no external charges. Each link variable belongs to SU(2). The physical Hilbert space is the gauge-invariant subspace of the product Haar-measure space on the seven links.

We use

\[
H=T+V,
\quad T=\kappa\sum_{\ell=1}^{7}E_\ell^2,
\quad V=\nu(2-W_A-W_B),
\quad W_p=\tfrac12\operatorname{Tr}_{1/2}U_p,
\]

with \(\kappa>0\), \(\nu\geq0\), and link Casimir eigenvalue \(j(j+1)\). Since \(-1\leq W_p\leq1\), \(0\leq V\leq4\nu\). The electric operator has discrete eigenvalues tending to infinity with finite multiplicities. A bounded magnetic perturbation therefore leaves a well-defined self-adjoint operator with compact resolvent on this fixed graph, permitting the ordered eigenvalues and min–max arguments used here.

All numerical tables set \(\kappa=1\). Thus energies and gaps are in units of \(\kappa\); time is in units of \(1/\kappa\), with \(\hbar=1\). Changing both coefficients by a common factor changes every energy difference by that factor.

Hamiltonian lattice gauge theory and the construction of orthonormal gauge-invariant spin-network bases are established methods. The compact two-plaquette specialization used here is derived explicitly in the accompanying code. See [Burgio et al., *The basis of the physical Hilbert space of lattice gauge theories*](https://arxiv.org/abs/hep-lat/9906036) and [D’Andrea et al., *A new basis for Hamiltonian SU(2) simulations*](https://arxiv.org/abs/2307.11829).

### Physical basis and electric energy

Gauss law identifies the spins along each of the two three-link outer paths. A physical basis state is labeled by \((j_L,j_R,j_M)\): left outer path, right outer path, and shared middle link. The labels obey

\[
|j_L-j_R|\leq j_M\leq j_L+j_R,
\qquad j_L+j_R+j_M\in\mathbb Z.
\]

For each admissible triple, the trivalent intertwiners have multiplicity one. The electric energy is

\[
T_{L,R,M}=\kappa\big[3j_L(j_L+1)+3j_R(j_R+1)+j_M(j_M+1)\big].
\]

The cutoff \(J\) retains every admissible triple with all three spins at most \(J\). Increasing \(J\) adds physical gauge states on the same seven links. It does not enlarge the graph or refine the lattice spacing.

### Magnetic interaction

Write \(d_j=2j+1\). The normalized fundamental trace on the left plaquette has matrix elements

\[
\langle j_L',j_R,j_M'|W_A|j_L,j_R,j_M\rangle
=\frac12\sqrt{d_Ld_Md_{L'}d_{M'}}
\begin{Bmatrix}j_L&j_R&j_M\\j_M'&\tfrac12&j_L'\end{Bmatrix}^{\!2}.
\]

The allowed transitions shift \(j_L,j_M\) by \(\pm\tfrac12\), leaving \(j_R\) fixed. Exchange left and right for \(W_B\). The braces are a Wigner 6j symbol, the coefficient for changing how angular momenta are coupled. The code computes its square exactly as a rational number before any final floating-point square root.

One derivation uses normalized invariant functions

\[
f_{L,R,M}(A,B)=\sqrt{d_Ld_R/d_M}\,
\operatorname{Tr}\big[P_M(D^L(A)\otimes D^R(B))\big],
\]

where \(P_M\) projects onto total spin \(M\). Recoupling an extra spin \(1/2\) produces the displayed coefficient. With this loop orientation, the \((1/2,1/2,0)\) function is \(\operatorname{Tr}(AB^{-1})\); reversing the orientation assigned to \(B\) gives the earlier \(\operatorname{Tr}(AB)\) convention and preserves Haar measure and \(\operatorname{Tr}B\).

## 2. What changes as the basis grows

At \(\nu/\kappa=1\), diagonalizing each orthogonal compression gives:

| Maximum link spin J | Physical states | Vacuum energy E₀/κ | First excited energy E₁/κ | Projected gap/κ |
|---:|---:|---:|---:|---:|
| 0.5 | 4 | 1.840397169 | 5.000000000 | 3.159602831 |
| 1 | 11 | 1.836046029 | 4.930560435 | 3.094514407 |
| 1.5 | 23 | 1.836011113 | 4.927983118 | 3.091972005 |
| 2 | 42 | 1.836011003 | 4.927970050 | 3.091959047 |
| 3 | 106 | 1.836011002 | 4.927970022 | 3.091959020 |
| 4 | 215 | 1.836011002 | 4.927970022 | 3.091959020 |
| 5 | 381 | 1.836011002 | 4.927970022 | 3.091959020 |
| 6 | 616 | 1.836011002 | 4.927970022 | 3.091959020 |

Each projected energy is an upper bound on its corresponding exact ordered energy. Their difference alone is not a guaranteed upper or lower bound on the exact gap. The next section supplies the missing opposite bounds.

The sweep covered nine coefficient ratios and eight cutoffs, for 72 diagonalizations. Representative results:

| ν/κ | Four-state projected gap/κ | J=6 projected gap/κ |
|---:|---:|---:|
| 0 | 3.000000000 | 3.000000000 |
| 0.1 | 3.001665896 | 3.000928491 |
| 1 | 3.159602831 | 3.091959020 |
| 4 | 4.789221847 | 4.209427462 |
| 40 | 32.740443235 | 14.819841043 |

Larger magnetic-to-electric ratios require larger spin bases. At ratio 40, the gap estimate changes from 14.848634897 at J=3 to 14.819841043 at J=6. The analytic enclosure evaluated in floating point at J=6 is approximately [14.819840657, 14.819841111]. That case has not received the separate exact-arithmetic certificate. The expected difficulty of electric-representation truncations toward weak gauge coupling is also discussed by [D’Andrea et al.](https://arxiv.org/abs/2307.11829).

## 3. How all omitted states are controlled

Let \(P=P_J\) project onto the retained states and \(Q=1-P\). Define \(A=PHP\). The fundamental magnetic operator only changes affected spins by one half, so every coupling from \(P\) to an omitted state lies within the complete cutoff \(J+1/2\).

Construct a finite auxiliary matrix \(L\) on that extended basis:

1. Keep the retained block \(A\).
2. Keep all magnetic couplings from retained states to the first omitted layer.
3. Replace the entire omitted-to-omitted block by its electric diagonal.

On the remaining infinite space keep the electric operator alone, denoted \(T_{\rm rest}\). The exact difference is

\[
H-(L\oplus T_{\rm rest})=QVQ\geq0.
\]

Positivity holds because the entire compressed positive operator \(QVQ\) is removed. Removing individual positive-looking matrix entries would not justify the inequality.

The minimum electric energy outside a cutoff has an explicit formula. Put \(K=2J\), \(m=K+1\), \(a=\lfloor m/2\rfloor\), \(b=\lceil m/2\rceil\). Then

\[
t_Q(J)=\frac\kappa4\big[3a(a+2)+3b(b+2)+m(m+2)\big].
\]

To minimize, place the largest spin on the shared link, which has electric weight one instead of three. Triangle admissibility requires the two outer doubled spins to sum to at least the largest. The minimizing choice has sum exactly \(m\) and is as balanced as possible. The resulting minimum increases with \(m\), so the first omitted value suffices. This gives, for example, \(t_Q(2)=26\kappa\) and \(t_Q(2.5)=34.5\kappa\).

Min–max comparison now gives, for levels \(j=0,1\),

\[
\ell_j:=\min\{\lambda_j(L),t_Q(J+1/2)\}
\leq E_j(H)\leq u_j:=\lambda_j(A).
\]

Consequently

\[
\boxed{\ell_1-u_0\leq E_1-E_0\leq u_1-\ell_0.}
\]

This controls the infinite omitted spin tower while leaving the spatial graph fixed.

### Exact arithmetic certificate

Floating-point diagonalization chooses convenient candidate endpoints. It does not certify them. For state \(i\), define \(K_i=d_{L_i}d_{R_i}d_{M_i}\) and \(D_{ii}=\sqrt{K_i}\). For rational parameters and a rational trial energy \(q\),

\[
D^{-1}(H-qI)D^{-1}
\]

has entirely rational entries. This is a congruence: it preserves the numbers of positive, negative and zero eigenvalues, rather than their values. After clearing denominators, integer fraction-free elimination counts these signs exactly. Four such checks, two for \(L\) and two for \(A\), certify the two energy intervals.

For \(\kappa=\nu=1\), retained J=2 and extended J=2.5 give:

| Level | Certified lower endpoint | Certified upper endpoint | Negative eigenvalues at lower / upper probe |
|---|---:|---:|---:|
| E₀ | 1.83601099 | 1.83601102 | 0 in L / 1 in A |
| E₁ | 4.92797000 | 4.92797006 | 1 in L / 2 in A |

The omitted electric remainder begins at 34.5, above these probes. Subtracting opposite endpoints gives the main gap certificate. Every rational probe and inertia result is recorded in the certificate JSON.

For \(\nu/\kappa=4\), retained J=3 and extended J=3.5 give E₀/κ in [5.79688776, 5.79688779] and E₁/κ in [10.00631521, 10.00631533]. The remaining electric spectrum begins at 56κ. These establish the second gap certificate.

At very large cutoffs for modest coupling, estimated enclosure widths reach floating-point precision and can even become slightly negative through rounding. Those raw numerical widths are not certificates. The exact rational probes above avoid any assumed eigensolver error allowance.

## 4. Independent checks and what a correlator can miss

The magnetic matrices were checked in three ways:

- Exact agreement with the earlier four-state matrix.
- Forty-eight independent angular-momentum projector contractions, with maximum difference 4.44×10⁻¹⁶.
- Direct SU(2) Haar integration of invariant basis functions, constructing projectors from Casimir eigenspaces without 6j coefficients. For J=1 and J=1.5, the largest matrix difference from the recoupling formula was 5.25×10⁻¹⁶. Increasing product quadrature from 8³ to 12³ points changed entries by at most 8.33×10⁻¹⁶.

Hermiticity, left–right exchange symmetry, the fundamental vacuum transition, character-sector reductions and nested variational energy ordering also passed. These floating-point cross-checks support the formula and implementation; the certificate separately uses the exact formula and rational arithmetic.

An independent audit also rebuilt the certificate blocks and spectral shifts and recomputed all twelve rational probes using ordinary Fraction Schur-complement elimination. Its sign counts agreed with the integer Bareiss calculation for every probe. The operator inequality, omitted-layer completeness and gap subtraction were separately reviewed.

The J=6 calculation at coefficient ratio one also evaluates connected Euclidean-time correlations for

\[
O_+=(W_A+W_B)/\sqrt2,
\qquad O_-=(W_A-W_B)/\sqrt2.
\]

The vacuum and lowest excitation are even under left–right exchange. The next excitation is odd. Thus the symmetric operator couples to the actual first excitation, whereas the antisymmetric operator misses it by symmetry.

| Operator | Lowest gap with nonzero overlap, in κ units |
|---|---:|
| Symmetric O₊ | 3.091959020089 |
| Antisymmetric O₋ | 3.113957723441 |

The connected spectral sum is

\[
C_O(t)=\sum_{n>0}|\langle n|O|0\rangle|^2e^{-(E_n-E_0)t}.
\]

Using \(\delta t=0.25\), the effective gap \(\log[C(t)/C(t+\delta t)]/\delta t\) of the symmetric operator falls from 3.095553626 at t=1 to 3.091959080 at t=8. The antisymmetric operator approaches the larger odd-sector gap. Tiny computed weights below 10⁻¹⁸ are excluded explicitly; total discarded weight is recorded in the results.

This is a practical lesson for larger calculations: a clean exponential decay in one observable need not reveal the lowest excitation. These correlators are constructed from the same Hamiltonian spectrum, so they check operator overlap rather than constitute independent physical measurements.

## 5. What remains for Yang–Mills

The official problem requires a nontrivial quantum Yang–Mills theory on four-dimensional spacetime, for any compact simple gauge group, with the stated axiomatic properties and a positive mass gap. The finite two-square system does not meet that target. See [Jaffe and Witten, *Quantum Yang–Mills Theory*, especially sections 4–6](https://www.claymath.org/wp-content/uploads/2022/06/yangmills.pdf).

The next mathematical tasks are distinct:

1. **Three spatial directions.** Build an actual three-dimensional lattice, starting with a cube. Higher-degree vertices require additional intertwiner labels; the unique-triple simplification used here no longer covers the basis.
2. **Growing volume.** Control the excitation energy above the interacting vacuum as more cubes are added. A positive gap on every finite graph can still approach zero with graph size. One needs an estimate that remains useful as volume grows.
3. **Continuum refinement.** Send the lattice spacing to zero while adjusting the bare coupling and fixing a physical reference scale. A dimensionless lattice gap can tend to zero while the physical mass remains positive; demanding a fixed gap in arbitrary lattice units would test the wrong condition.
4. **Existence and axioms.** Establish the limiting theory and its required properties, then carry the gap estimate into that limit. Extending SU(2) evidence alone also does not fulfill the statement for every compact simple group.

The completed step is a controlled interacting SU(2) benchmark and a reproducible method for bounding its omitted states. The difficult open step is obtaining estimates that survive spatial growth and continuum refinement.

## 6. Files and reproduction

- `Shared_SU2_Cutoff_Study_Results.json`: 72 cutoff/coupling results, numerical enclosures, matrix checks and correlator data.
- `Shared_SU2_Cutoff_Certificate_Checks.json`: exact rational endpoints and inertia certificates.
- `SU2_Direct_Haar_Validation.json`: independent integration matrices and residuals.
- `shared_su2_cutoff_study.py`, `su2_recoupling.py`, `su2_cutoff_bounds.py`: model, spectra and bound implementation.
- `verify_su2_cutoff_certificates.py`: repeatable exact certificate calculations.
- `su2_direct_haar_validation.py`: independent magnetic-matrix validation.
- `su2_bounds_audit_check.py`: second exact-arithmetic verification of all twelve probes.

The reproduction archive includes the scripts, results, proof notes, a README and a SHA-256 manifest. It requires Python with NumPy; it needs no symbolic algebra package or network connection to run. The user's original research document is unchanged.
