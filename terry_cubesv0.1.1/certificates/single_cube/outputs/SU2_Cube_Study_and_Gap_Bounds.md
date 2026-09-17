# SU(2) on a cube: spectra, omitted states, and a positive gap

**Research calculation — 15 September 2026**

## Outcome

We extended the interacting gauge calculation from two adjacent squares to all six faces of an open cube. The model now has twelve links, eight vertices and square loops in the xy, xz and yz planes. Gauss law is imposed at every vertex.

At the same coefficient ratio used in the earlier benchmark, ν/κ=1, the largest numerical calculation gives

\[
\Delta/\kappa\approx 2.95104.
\]

This estimate comes from a complete total-electric-energy cutoff containing **34,680 physical states**. Its underlying projected gap is 2.951042834527. The additional digits describe that finite projection; they are not claimed as certified digits of the complete cube gap.

A separate calculation with rigorously bounded integer arithmetic proves a positive gap on the **untruncated physical Hilbert space of this fixed cube**. The strongest completed certificate is

\[
\boxed{2.941733\leq\Delta/\kappa\leq2.958693}
\qquad(\nu/\kappa=1).
\]

These endpoints are exact terminating decimals. The certificate uses 398 retained states and 1,644 directly coupled omitted states; an analytic comparison covers every remaining state. The numerical estimate and the certified interval have different precision and are reported separately.

This is a finite-lattice benchmark. It includes arbitrarily large link spins in the target operator, but it has no infinite-volume or continuum limit. It is not a solution of the Millennium problem or a claim of a new result in continuum quantum field theory.

## 1. What was actually built

The Hamiltonian is

\[
H=T+V,
\qquad T=\kappa\sum_{e=1}^{12}j_e(j_e+1),
\qquad V=\nu\left(6-\sum_{f=1}^{6}W_f\right),
\qquad W_f=\tfrac12\operatorname{Tr}_{1/2}U_f.
\]

We use κ>0 and ν≥0. Each normalized SU(2) trace lies between −1 and 1, so 0≤V≤12ν. All quoted energies are in units of κ. No value in electronvolts, cosmological density or universal geometric constant is inferred.

The graph has five independent cycles, although the magnetic Hamiltonian contains all six face operators. The sixth face was included, with its actual shared-link dependence; the faces were not modeled as six independent oscillators.

At each vertex, the three incident spins satisfy triangle inequalities and have an integer sum. These are exactly the conditions for a singlet in the product of three SU(2) representations. Its multiplicity is one, so the link-spin configuration specifies a physical basis state without additional local labels.

**A clarification to the previous roadmap:** the vertices of a single open cube still have degree three. Additional coupling labels become necessary as cubes are joined and vertices meet more than three links. For example, four incident spin-½ representations can combine into a singlet in more than one way. This basis issue is discussed in [Kavaki and Lewis, *From square plaquettes to triamond lattices for SU(2) gauge theory*](https://www.nature.com/articles/s42005-024-01697-4).

All eight vertices here are boundary vertices. The single cube supplies three face orientations, but it does not yet provide a three-dimensional bulk with interior vertices.

## 2. Magnetic interactions and independent validation

The magnetic matrix elements use the signed trivalent recoupling formula in Eq. 22 of [Kavaki and Lewis](https://www.nature.com/articles/s42005-024-01697-4). At a visited vertex, let e denote the spectator link, f the forward face link, and b the backward face link. Lowercase spins are initial labels and uppercase spins are final labels. The local trace factor is

\[
(-1)^{j_e+j_f+J_b+1/2}
\sqrt{(2J_f+1)(2j_b+1)}
\begin{Bmatrix}j_e&j_f&j_b\\1/2&J_b&J_f\end{Bmatrix}.
\]

Multiplying the four vertex factors and dividing by two gives the normalized Wilson trace matrix element. Each affected link changes by ±½. The code retains the sign and the exact rational square of the coefficient's magnitude.

The checks address both normalization and signs:

- An independent exact Haar-integration construction evaluated all **192 directed matrix entries** at maximum link spin J=½, without using 6j symbols.
- The independent loop-character basis and recoupling basis agree after one common diagonal change of state signs across all six faces. Twelve of the 32 state signs differ between the conventions. This is a change of basis, not a modification of the physics.
- In the loop-character convention, the J=½ amplitudes are exactly ½, ¼ and ⅛, depending on the surrounding flux.
- The general formula reproduces both earlier two-square matrices through J=3/2 after their corresponding common basis rephasing.
- Exact Hermiticity, the vacuum-to-fundamental-face amplitude ½, and gauge-admissibility checks pass.

The independent cube Haar calculation validates the J=½ retained block. It does not independently integrate every higher-spin coefficient. Those use the sourced general recoupling formula, its exact arithmetic implementation, and the additional consistency tests. Taking absolute values of the signed higher-spin matrices is not part of this calculation.

## 3. Two ways of increasing the physical basis

### A. Increase the maximum spin on every link

At ν/κ=1:

| Maximum link spin J | Physical states | E₀/κ | E₁/κ | Projected gap/κ |
|---:|---:|---:|---:|---:|
| ½ | 32 | 5.541783784 | 8.772374285 | 3.230590500 |
| 1 | 1,013 | 5.506884796 | 8.466665143 | 2.959780347 |
| 3/2 | 14,879 | 5.506498362 | 8.457590259 | 2.951091897 |

The J=2 basis contains 148,678 states. It was counted exactly, but its full magnetic matrices were not generated. A different complete cutoff reaches higher spins with fewer low-energy states to diagonalize.

### B. Retain every state below a total electric-energy cap

Here P contains **all** physical states with T/κ≤C, without an additional chosen spin cutoff. Since each term is nonnegative, a doubled link label q obeys q(q+2)≤4C. This gives a sufficient finite label bound, and exact energy pruning enumerates the full subspace.

At ν/κ=1:

| Electric-energy cap C | Physical states | Largest link spin present | Projected gap/κ |
|---:|---:|---:|---:|
| 12 | 398 | 3/2 | 2.958505193 |
| 18 | 2,920 | 2 | 2.951119238 |
| 24 | 11,409 | 2 | 2.951043688 |
| 30 | 34,680 | 5/2 | 2.951042835 |

The cap-30 basis includes 1,464 states beyond J=2. This is not a relabeling of a fixed maximum-spin calculation. Agreement between the two cutoff routes is a useful check on the approximate low spectrum.

Both E₀ and E₁ decrease in each nested variational sequence. Each is a variational upper bound on its corresponding ordered energy. Their difference alone is not necessarily an upper or lower bound on the true gap.

### Larger magnetic coefficient

The study also tested ν/κ=¼ and 4. At ratio four, cutoff effects remain more pronounced:

| Electric-energy cap | Projected gap/κ at ν/κ=4 |
|---:|---:|
| 12 | 4.308033188 |
| 18 | 3.817699983 |
| 24 | 3.703289394 |
| 30 | 3.681572293 |

The last value is still changing visibly with the cutoff and is not reported as a converged full-cube gap. The exact positive-gap certificates in this package concern ν/κ=1.

## 4. Controlling every omitted state

For a retained projection P, include every omitted state reached by one face operator in a finite set Q₁. Let R contain all remaining states and Q=Q₁⊕R. Construct

\[
A=PHP,
\qquad
L=\begin{pmatrix}A&PHQ_1\\Q_1HP&T_{Q_1}\end{pmatrix}.
\]

Every P-to-Q coupling is retained. The exact difference is

\[
H-(L\oplus T_R)=QVQ\geq0.
\]

All omitted magnetic terms are removed as one compressed positive operator. This includes interactions between Q₁ and R. Positivity would not follow from dropping arbitrary individual matrix entries.

If t_R is a proven electric lower bound on R, ordered eigenvalues satisfy

\[
\min\{\lambda_j(L),t_R\}\leq E_j(H)\leq\lambda_j(A).
\]

Subtracting opposite bounds for E₀ and E₁ then encloses the gap.

### An exact omitted electric threshold

Outside J=½, the least possible electric energy is exactly **6.5κ**. The proof uses the cube's cycles and the even-sum Gauss constraint, with a separate bound excluding every doubled spin at least three. A state attaining 6.5κ has a spin-one edge and two adjacent three-edge return paths carrying spin ½.

The separate bound note gives the analytic proof and a terminating integer search for higher cutoffs. Their minima are:

| Retained J | Exact minimum omitted T/κ |
|---:|---:|
| ½ | 6.5 |
| 1 | 12 |
| 3/2 | 18 |
| 2 | 26 |
| 5/2 | 34 |

The last value differs from the two-square formula, which gives 34.5. Additional cube paths change the energy minimization; the earlier graph's formula cannot simply be reused.

### Exact arithmetic, including the square roots

For the first certificate, P has 32 states and Q₁ has 180, so L has 212 states. Its magnetic coefficients contain square roots of rational numbers. Integer square-root brackets enclose each coefficient by rationals. The largest sum of entry-error bounds in a row bounds the operator norm of the entire symmetric error matrix.

With rounding denominator 10¹², the retained matrix is exact and the auxiliary matrix error is at most 8×10⁻¹². Lower probes are shifted upward by this error; upper probes are shifted downward by their respective error.

Below the omitted electric diagonal, an exact Schur complement reduces the sign count of L−qI to a 32-dimensional matrix. Rational elimination then certifies the number of eigenvalues below each probe. Floating-point eigensolves may suggest candidate probes, but the exact sign checks decide whether they are valid.

This gives

| Level | Certified lower endpoint | Certified upper endpoint |
|---|---:|---:|
| E₀/κ | 5.30400030 | 5.54178379 |
| E₁/κ | 6.20980338 | 8.77237429 |

The resulting gap interval is [0.66801959, 3.46837399]. An independent audit rebuilt the auxiliary space by explicit face-label shifts, checked radical brackets by integer binary search, and reproduced all four sign counts using a different exact elimination method.

Complete electric-energy subspaces strengthen this certificate. Their remaining states have T greater than the chosen cap. A quarter-unit lower bound follows immediately from the Casimir formula; on this bipartite graph, parity in fact gives half-integer total electric energies. The certificates use the conservative quarter-unit bound.

| Retained electric-energy cap | Retained states | Auxiliary states | Certified full-cube gap interval, in κ units |
|---:|---:|---:|---:|
| 8 | 86 | 508 | [2.236257, 3.068193] |
| 10 | 212 | 1,117 | [2.889649, 2.979921] |
| 12 | 398 | 2,042 | **[2.941733, 2.958693]** |

For the last row, the level bounds are E₀/κ in [5.506419, 5.506607] and E₁/κ in [8.448340, 8.465112]. Subtracting opposite endpoints gives the displayed gap interval. The omitted electric remainder is at least 12.25κ, safely above the probes.

These larger sign counts use an additional rigorous arithmetic method to keep the computation practical. Exact rational matrix entries are enclosed by integer endpoints scaled by 2¹²⁸. Every Schur update rounds its endpoints outward, and a pivot sign is accepted only if its entire interval excludes zero. An unresolved sign causes failure rather than an assumed answer. This is integer interval arithmetic, not a presumed error tolerance for floating-point eigenvalues.

The interval sign counter was tested against known exact inertias and integer Bareiss elimination, including indefinite blocks and the four original cube probes. It also underwent a separate algorithm review and comparison on the cap-eight probes. The independent audit files identify precisely which larger calculations were rerun by a second method.

### Much tighter numerical evaluations

The same operator inequality was also evaluated on larger auxiliary matrices:

| Link cutoff J | Retained states | Auxiliary states | Numerically evaluated gap endpoints |
|---:|---:|---:|---:|
| ½ | 32 | 212 | 0.6680196053, 3.4683739797 |
| 1 | 1,013 | 5,927 | 2.9386864844, 2.9604884366 |
| 3/2 | 14,879 | 69,377 | 2.9510181075, 2.9510937766 |

These narrow final endpoints are numerical evaluations, not exact certificates. Direct residual checks and two starting seeds support the eigensolver results but do not themselves certify eigenvalue ordering or all arithmetic errors. They are kept separate from the exact endpoints above.

## 5. What follows for the research argument

Three conclusions are supported:

1. The interacting SU(2) construction extends consistently to all faces of a cube.
2. Small spin bases can substantially shift the predicted gap. Increasing them and controlling their complement are both necessary.
3. The full fixed-cube Hamiltonian has a positive gap at the certified coefficient ratio, with explicit checkable endpoints.

The next spatial extension requires gluing cubes and retaining the additional coupling labels at vertices of degree four or higher. Simply keeping the trivalent link-label basis would omit physical states there. The resulting bounds would then need to remain useful as the graph grows.

Separately, a continuum construction requires controlling lattice refinement with appropriately adjusted coupling and a fixed physical scale. The official Yang–Mills problem requires a nontrivial theory on four-dimensional spacetime, its stated axiomatic properties, and a positive gap for any compact simple gauge group. The present finite SU(2) benchmark does not establish those requirements. [Jaffe and Witten, *Quantum Yang–Mills Theory*](https://www.claymath.org/wp-content/uploads/2022/06/yangmills.pdf)

Our next unknown is therefore precise: **can the excitation bound be controlled when cubes are joined, without losing it as the volume grows?** The current calculation does not extrapolate a positive infinite-volume limit from a few finite graphs.

## Reproduction

The accompanying archive contains the model, both cutoff studies, exact certificates, independent checks, mathematical notes, JSON results, and a README. It uses Python and NumPy, with no SciPy or symbolic algebra dependency. The package's reproduction command reruns the calculations and checks the reported outputs. Exact endpoint fractions and sign counts are separate from numerical residuals and tolerances.

The earlier two-square results and the original research document remain unchanged.
