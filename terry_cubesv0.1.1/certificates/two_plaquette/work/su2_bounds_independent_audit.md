# Independent audit of the fixed-graph SU(2) cutoff certificates

**Outcome: the analytic enclosure and all three supplied certificates pass.** No blocking mathematical error was found. This conclusion concerns the stated open two-plaquette SU(2) Hamiltonian with kappa = 1 and nu = 1 or 4, including its unbounded spin labels. It does not concern an infinite lattice or the four-dimensional continuum limit.

## 1. Operator lower bound

Write `H = T + V`, where `V = nu(2 - W_A - W_B)`. In the exact holonomy representation, `W_A = Tr(A)/2` and `W_B = Tr(B)/2` are multiplication by real numbers in [-1,1]. Consequently `V >= 0` for nu >= 0.

For the retained projection P and Q = I - P, define

`S = PHP + PHQ + QHP + QTQ`.

Then `H - S = QVQ >= 0` as a quadratic-form inequality. This keeps the entire P-to-Q coupling; it does not discard individual off-diagonal entries under an unjustified positivity assumption.

The fundamental Wilson operator changes the active doubled outer spin and doubled middle spin by plus or minus one, preserving the other outer spin. Therefore every direct coupling from the cutoff `max(2jL,2jR,2jM) <= K` lies in the complete basis with cutoff K+1. Calling that first omitted layer Q1 and the remaining space R, the exact auxiliary operator splits as `S = L direct-sum T_R`, with L precisely the finite matrix used by the code. I also exhaustively checked these target-selection rules for retained cutoffs K = 1 through 12.

## 2. Omitted electric threshold

The doubled-spin electric energy is

`T(a,b,c) = kappa [3a(a+2) + 3b(b+2) + c(c+2)]/4`.

For an omitted admissible triple, permuting its largest label into c cannot increase this expression, because c has the smallest coefficient and the triangle/parity conditions are symmetric. At fixed c, the triangle rule requires a+b >= c. Monotonicity and convexity give the minimum at a+b = c with a and b as balanced as possible. This minimum increases with c. Thus the exact omitted threshold occurs at `c = K+1`, `a = floor(c/2)`, `b = ceil(c/2)`, as implemented. The remainder beyond the extended layer uses the same formula with K replaced by K+1.

## 3. Min-max and probe logic

For kappa > 0 the electric eigenvalues tend to infinity with the spin labels; V is bounded. The fixed-graph Hamiltonian consequently has discrete eigenvalues to which min-max applies.

With eigenvalues ordered from index zero, the valid inequalities are

`min(lambda_j(L), t_R) <= E_j(H) <= lambda_j(PHP)`.

The lower rational probe is valid when L has at most j eigenvalues strictly below it and it does not exceed t_R. The upper rational probe is valid when PHP has at least j+1 eigenvalues strictly below it. These are exactly the implemented inertia conditions. The gap interval correctly subtracts the upper E0 bound from the lower E1 bound, and the lower E0 bound from the upper E1 bound.

## 4. Independent exact-arithmetic verification

I rebuilt the auxiliary blocks and spectral shifts from the rational recoupling congruence, without calling the audited block/shift helpers. I recomputed all twelve actual probe inertias using ordinary exact Fraction Schur-complement elimination, independently of the integer Bareiss routine. Every inertia equals both the saved certificate and the Bareiss result. No zero pivot requiring a two-by-two step occurred in these actual probes.

The rational congruence's diagonal and off-diagonal factors agree algebraically with `D^-1 H D^-1`, for `D_ii = sqrt(dL dR dM)`: a W_A coefficient becomes `(6j)^2/(2dR)`, and W_B similarly uses dL. Reconstructing the floating-point Hamiltonian from this congruence agrees to at most 2.85e-14. All saved exact gap subtractions were independently reproduced.

| Retained J | nu | Certified full-spin gap interval |
|---:|---:|---:|
| 2 | 1 | [3.09195898, 3.09195907] |
| 2 | 4 | [4.20913941, 4.20995910] |
| 3 | 4 | [4.20942742, 4.20942757] |

Audit reproduction: `work/su2_bounds_audit_check.py`, using `work/su2_recoupling.py`, `work/su2_cutoff_bounds.py`, and `outputs/Shared_SU2_Cutoff_Certificate_Checks.json`.

The final `work/verify_su2_cutoff_certificates.py` generator was rerun, followed by the independent checker. Its object schema with a `benchmarks` array retains all three cases above. All twelve exact probes passed again with identical certified intervals. The independent checker accepts both this final schema and the earlier list schema.

## Applicability requirements

The general helper functions do not enforce every theorem hypothesis. Their use requires kappa > 0, nu >= 0, the stated exact physical Hamiltonian, and the complete first omitted layer. The provided three certificates meet those requirements. Arbitrary matrix input, negative nu, or a partial boundary layer must not inherit this certification claim. This audit did not modify either mathematical implementation.
