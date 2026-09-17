# Joined-cube SU(2) Hamiltonian: basis and magnetic conventions

## Physical model

The graph is an open N×1×1 chain of unit cubes, with 4N+4 vertices,
8N+4 physical links, and 5N+1 distinct physical square faces. Every face,
including an interior shared face, appears once. There are no external
charges. The dimensionless normalization is

H = κ Σ_physical links j(j+1) + ν Σ_physical faces [1 − Tr_fund(U_face)/2].

Both κ and ν are specified independently in the finite study. Restoring
lattice units requires the usual coupling and spacing factors; these finite
calculations establish neither an infinite-volume nor a continuum gap.
The chain grows in one direction, while its two transverse widths remain
one cube. It is not a sequence approaching a three-dimensional bulk limit.

## Complete four-valent intertwiners

Each interior physical vertex has four incident links. Split it into two
trivalent vertices joined by an auxiliary tree edge. Pair physical spins
(a,b) on one branch and (c,d) on the other. The auxiliary spin k runs over
the entire intersection of the SU(2) coupling ranges a⊗b and c⊗d, with the
usual integer/half-integer parity condition. Each common k labels one
orthonormal invariant tensor. In particular, four spin-half links have two
singlets, k=0 and k=1. Keeping only k=0 loses a physical state.

The auxiliary edge has zero electric cost and no independent spin cutoff.
Its allowed labels are bounded by the incident physical spins. Gauge-fixing
its group element to the identity identifies the resolved graph with the
original graph. The split increases both the number of links and the number
of vertices by one, creates no new graph cycle, and resolves the original
intertwiner multiplicity. It does not add a propagating gauge field.

The implementation labels all spins by integers q=2j. State tuples list
physical link labels first and auxiliary labels last. A physical cutoff
q≤q_max applies only to the first group. At a complete electric cutoff T≤C,
the sufficient physical bound is q≤floor(sqrt(4C+1)−1). The enumeration
then retains every admissible physical assignment and every auxiliary
coupling channel. Exact integer energy and vertex-completion lower bounds
prune the search without adding another cutoff.

## Wilson operators with a fixed orientation convention

Orient each resolved edge from the lower vertex index to the higher index.
Write K_j(U)=D_j(U) ε_j, where
ε_j(m,n)=(-1)^(j−m) δ_(m,−n). Contract K matrices with normalized real 3j
invariant tensors at their endpoints. Every invariant uses ascending
incident edge indices as its fixed order.

Multiplication by a fundamental Wilson loop shifts each visited spin by
±1/2. It leaves every unvisited spin unchanged. At a visited vertex let f
and b denote the forward and backward loop legs, and e the spectator. Form
the contraction of the initial and final 3j tensors, the two Clebsch–Gordan
coefficients C(j_f,1/2;J_f), C(j_b,1/2;J_b), and ε_(1/2). The complete
matrix element is

  (1/2) (-1)^r ∏_visited edges sqrt((2j+1)/(2J+1)) ∏_visited vertices C_v,

where r counts loop edges traversed opposite their chosen orientation.
This follows from the product decomposition of D_j(U)ε_j and its
fundamental counterpart into irreducible matrices at both endpoints.
The factor 1/2 is the requested normalized trace, not a vertex factor.

The local contraction is first evaluated in order (e,f,b). An odd
permutation to the fixed incident-edge order multiplies it by the product
of the old and new 3j permutation phases. In doubled labels this is
(-1)^[(q'_f+q'_b−q_f−q_b)/2]. The global (-1)^r is also necessary.
Omitting it can produce different spectra for valid resolution trees
whose lifted loops have odd length.

`su2_oriented_vertex.py` evaluates the 3j finite sum and elementary
spin-half Clebsch–Gordan coefficients with exact Fraction arithmetic.
Each term is a signed square root of a rational. Before adding terms it
checks by integer square roots that their ratios are rational. The sign of
the exact sum determines the sign of the matrix element. There is no
floating-point sign test. Its squared magnitude also agrees exactly with
the corresponding 6j identity. Sparse entries remain
(row, column, sign, exact rational squared magnitude).

For comparison with the earlier fixed-order trivalent formula, the local
sign ratio is (q'_f−q_f)(q'_b−q_b). Those factors multiply to +1 around a
closed loop. The corresponding dimension ratios telescope to one. Thus the
whole-face difference is precisely (-1)^r. For the default longitudinal
resolution this is one common basis transformation

  |q> → (-1)^[Σ q_aux at original vertices with y≠z] |q>.

It leaves the electric diagonal and every spin or energy cutoff invariant.
The old and final default Hamiltonians therefore have the same spectrum
at every κ and ν. The face-incidence identity was checked exactly for
N=1,2,3,4, and the local exact identities are separately regression-tested.
The root study is rerun with the final code as a reproducibility check.

## Validation and feasible sizes

Complete physical-spin cutoff j≤1/2, with all auxiliary channels:

| Cubes | Physical states |
|---|---:|
| 1 | 32 |
| 2 | 868 |
| 3 | 25,676 |

Independent real-link enumeration reproduces all three counts. The
two-cube count would be only 512 if one incorrectly assigned a single
intertwiner to each admissible real-link pattern.

Complete electric-energy cutoffs:

| Cubes | T≤8 | T≤12 | T≤18 |
|---|---:|---:|---:|
| 1 | 86 | 398 | 2,920 |
| 2 | 622 | 8,289 | 201,248 |
| 3 | 1,555 | 40,038 | 2,848,545 |

An independent implementation constructs invariant tensors from total
angular-momentum eigenspaces and Clebsch–Gordan tensors from lowering
operators; it uses no 6j symbols. On all three two-cube resolution trees
(longitudinal, crossed, and mixed), all 868-state face matrices agree with
the final model after one common diagonal basis phase. Their complete
spectra agree within 1.0×10^-13. At κ=ν=1 the j≤1/2 two-cube gap is
3.22259354077756. This is a cutoff value for this finite lattice.

Higher-spin checks also pass: N2 T≤8 (622 states, 5,584 directed entries),
N1 T≤12 (398 states, 4,260 entries), and N3 T≤6 (307 states, 1,936 entries).
The largest entry-magnitude discrepancy is below 9×10^-16. The earlier
single-cube results were checked separately by exact matrix rephasing and
full spectra at physical j≤1 and at T≤12; their spectra are unchanged.

See `outputs/SU2_Joined_Cubes_Independent_Validation.json` for the independent
checks and the root study outputs for the larger cutoff and volume runs.

## Reproducibility and sources

Main API: `work/su2_joined_cubes_model.py`. Required local modules:
`work/su2_oriented_vertex.py`, `work/su2_cube_model.py`, and
`work/su2_recoupling.py`. These use the Python standard library and NumPy.

Independent validation: `work/su2_joined_cubes_independent_check.py` and
`work/su2_joined_cubes_tensor_probe.py`. No earlier result file was edited.

- [NIST DLMF, exact 3j finite sum, Eq.34.2.4](https://dlmf.nist.gov/34.2#E4).
- [NIST DLMF, 3j permutation symmetry](https://dlmf.nist.gov/34.3#ii).
- [From square plaquettes to triamond lattices for SU(2) gauge theory, Eq.22](https://www.nature.com/articles/s42005-024-01697-4),
  for the trivalent plaquette recoupling context. The all-outgoing
  orientation bookkeeping and the chain implementation above are the
  present derivation and are checked independently by tensor contraction.
