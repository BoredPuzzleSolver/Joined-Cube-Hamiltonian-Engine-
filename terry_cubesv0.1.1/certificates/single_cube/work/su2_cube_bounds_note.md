# Omitted-state bounds for the untruncated open SU(2) cube

## Scope and result

The target is the gauge-invariant Hilbert space of one open cube with twelve links, eight trivalent vertices, six square faces and no external charges. Each link spin may be arbitrarily large. A finite retained basis is a comparison space, not a restriction on the target Hilbert space.

The operator is

\[
H=T+V,\qquad
T=\kappa\sum_{e=1}^{12}j_e(j_e+1),\qquad
V=\nu\left(6-\sum_{f=1}^{6}\frac{\operatorname{Tr}U_f}{2}\right).
\]

For \(\kappa>0,\nu\geq0\), \(0\leq V\leq12\nu I\). Each electric sublevel space is finite dimensional; \(H\) consequently has compact resolvent on this fixed graph. Ordered low energies and their difference are well-defined.

At cutoff \(J=1/2\), the exact minimum electric energy among all omitted states is \(13\kappa/2\). A positive omitted-block comparison then permits a certificate for the full graph using a retained basis of only 32 states. Larger comparisons sharpen the numerical bounds.

These statements remove the spin cutoff for this fixed graph. They do not take infinite spatial volume or a continuum limit.

## 1. Short analytic proof of the \(J=1/2\) electric threshold

Use doubled integer labels \(a_e=2j_e\). At each vertex the three labels obey triangle inequalities and have even sum. The energy is \((\kappa/4)\sum_e a_e(a_e+2)\).

First suppose an edge carries a label \(q\geq3\). At either endpoint, its other two labels sum to at least \(q\). Convexity minimizes their Casimir sum at \(\lfloor q/2\rfloor,\lceil q/2\rceil\). The four neighboring edges are distinct on the cube, so the total energy obeys

\[
T\geq\frac{\kappa}{4}
\left[q(q+2)+2\left\{
\lfloor q/2\rfloor(\lfloor q/2\rfloor+2)
+\lceil q/2\rceil(\lceil q/2\rceil+2)\right\}\right]
\geq\frac{37\kappa}{4}.
\]

This increasing bound excludes every larger label from a putative state below \(6.5\kappa\).

It remains to consider labels \(0,1,2\), with at least one 2. Parity at every vertex makes the label-1 edges a disjoint union of cycles. The cube has eight vertices and girth four, so the number of such edges is 0, 4, 6 or 8.

- With at least six label-1 edges, one label-2 edge already gives \(T\geq(6\cdot3+8)\kappa/4=6.5\kappa\).
- With no label-1 edges, the support of label-2 edges has minimum degree at least two, by the triangle inequality. It contains a cycle of at least four edges, costing at least \(8\kappa\).
- With four label-1 edges, they form a chordless square. A label-2 edge cannot join two vertices of that square without using one of its already label-1 edges. It therefore has an endpoint outside the square. That endpoint has no label-1 edge, and the triangle inequality requires a second label-2 edge. The cost is at least \((4\cdot3+2\cdot8)\kappa/4=7\kappa\).

Finally, one label-2 edge together with two adjacent three-edge return paths, each labelled 1, is admissible and costs \(6.5\kappa\). Thus

\[
\min_{Q(J=1/2)}T=\frac{13\kappa}{2}.
\]

This proof is independent of numerical diagonalization and of the threshold search below.

## 2. Exact larger-cutoff thresholds

The module su2_cube_bounds.py also implements a terminating exact integer minimization. Cube edge transitivity allows a largest label \(q\) to be placed on the edge \((0,1)\). All remaining labels then lie in \(0,\ldots,q\). Triangle/parity constraints and positive energy bounds prune an otherwise exhaustive search.

An explicit two-face configuration supplies an initial upper candidate. The increasing marked-edge bound displayed above excludes every sufficiently large \(q\), so the search does not silently assume that the first omitted label is the largest relevant one.

At an intermediate assignment, the algorithm minimizes the incident energy separately at each vertex. Summing these minima and dividing by two is a valid lower bound, because each edge occurs at two vertices. All comparisons use integer units of \(\kappa/4\); no eigenvalue or floating-point tolerance enters this threshold computation.

| Retained cutoff \(J\) | Exact omitted electric minimum divided by \(\kappa\) |
|---:|---:|
| 0 | 3 |
| \(1/2\) | \(13/2\) |
| 1 | 12 |
| \(3/2\) | 18 |
| 2 | 26 |
| \(5/2\) | 34 |

The final row is less than the value \(34.5\) supplied by the two-face configuration. Additional cube paths can lower the electric cost. The two-plaquette formula should therefore not be extrapolated to all cube cutoffs.

The output records a minimizing admissible state, search-node counts, the last tested largest label, and the bound excluding all larger labels.

## 3. Finite auxiliary lower comparison

Let \(P\) be any finite retained set of physical basis states. Let \(Q_1\) contain every omitted state reached from \(P\) by one fundamental plaquette multiplication. The remaining states form \(R\). A square Wilson loop shifts its four affected spins by \(1/2\), so this neighbor set is finite and \(PHR=0\).

Set

\[
A=PHP,\qquad B=PHQ_1,\qquad
L=\begin{pmatrix}A&B\\B^\dagger&T_{Q_1}\end{pmatrix}.
\]

Keeping the complete \(P\)-to-\(Q_1\) coupling and replacing the whole omitted block by its electric part gives

\[
H-(L\oplus T_R)=0_P\oplus QVQ\geq0.
\]

If \(a_j,\ell_j\) are the ordered finite eigenvalues and \(t_R\) is any proven electric lower bound on the remainder, then

\[
\min(\ell_j,t_R)\leq E_j(H)\leq a_j.
\]

For the reachable-neighbor construction, the minimum outside \(P\) is always a safe \(t_R\). A threshold outside the next complete spin cutoff would require actually including that complete layer or separately excluding every omitted lower-energy state.

An electric-energy sublevel basis is also useful. On the bipartite cube, \(T/\kappa\) is a half-integer. Indeed, sum the even incident-label sums over either vertex color class: each edge appears once, so \(\sum_e a_e\) is even. Since \(a(a+2)\equiv a\pmod2\), the total doubled-label Casimir is even. Thus retaining every state with \(T\leq8\kappa\) or \(T\leq10\kappa\) gives omitted thresholds of at least \(8.5\kappa\) or \(10.5\kappa\), respectively. This conclusion requires a complete energy-sublevel basis, without a hidden per-link spin cutoff.

Individual bounds \(L_j\leq E_j\leq U_j\) imply

\[
L_1-U_0\leq\Delta\leq U_1-L_0.
\]

As a coarse alternative, min–max on \(T\leq H\) and the constant electric-vacuum trial state give \(\Delta\geq3\kappa-6\nu\). A nonpositive value of that coarse estimate does not imply a zero gap.

## 4. Numerical refinement and certification

At \(\kappa=\nu=1\), the following are numerical evaluations of the exact auxiliary comparison:

| \(J\) | Retained states | Reachable omitted states | Evaluated gap endpoints |
|---:|---:|---:|---:|
| \(1/2\) | 32 | 180 | \(0.6680196053,\;3.4683739797\) |
| 1 | 1,013 | 4,914 | \(2.9386864844,\;2.9604884366\) |
| \(3/2\) | 14,879 | 54,498 | \(2.9510181075,\;2.9510937766\) |

The largest auxiliary matrix has 69,377 states. Sparse Lanczos calculations used full reorthogonalization, direct residual checks and two independent starting seeds. These checks support numerical reliability; they do not certify that every endpoint encloses the corresponding ordered eigenvalue. The JSON explicitly labels these results as not machine-certified.

For exact endpoints at small \(P\), each signed square-root matrix element can instead be enclosed by rational numbers, with an exact symmetric row-sum bound on the total operator error. At probes below every \(Q_1\) electric energy, the diagonal omitted block can be eliminated:

\[
L-zI\ \sim\
\left[A-zI-B(T_{Q_1}-zI)^{-1}B^\dagger\right]
\oplus(T_{Q_1}-zI).
\]

The second block is positive. Rational inertia therefore needs only the retained-size Schur matrix. The separate su2_cube_rational_certificate.py implements this certification path; its output, rather than raw Ritz residuals, supplies certified endpoints.

## Reproduction

Run:

    python work/su2_cube_bounds.py

It writes outputs/SU2_Cube_Bounds_Results.json. Dependencies are the existing cube recoupling and sparse-spectrum modules, Python and NumPy. It changes none of the earlier two-plaquette certificate files.
