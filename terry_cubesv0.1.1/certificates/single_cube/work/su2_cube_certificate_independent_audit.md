# Independent audit of the complete-cube gap certificate

**The 32-state retained-basis certificate passes.** I found no error in the operator comparison, signed radical enclosure, norm allowance, shifted probes, or exact inertia calculations. Its full-spin fixed-cube gap interval is

`66801959/100000000 <= Delta <= 346837399/100000000`,

or **[0.66801959, 3.46837399]**, for kappa = nu = 1.

## Operator comparison and omitted couplings

For `H = T + V`, the potential is `V = nu(6 - sum_f W_f) >= 0`, since each normalized SU(2) trace lies in [-1,1]. Let P be the retained physical basis, Q1 every omitted state directly coupled to it by a face, and R the remaining states. The auxiliary matrix L keeps PHP and every P-to-Q1 entry and uses only the diagonal electric energy on Q1. Then

`H - (L direct-sum T_R) = 0_P direct-sum QVQ >= 0`, where `Q = Q1 + R`.

Missing Q1-to-Q1, Q1-to-R, and R-to-R magnetic couplings do **not** invalidate the inequality: they are all part of the single positive compressed operator QVQ being removed. This would fail if some nonzero P-to-R coupling were missed.

I independently enumerated all sixteen sign choices on each face's four doubled-spin labels, applied the triangle/parity constraints at all eight vertices, and recovered precisely the same 180 reachable omitted states. Thus the auxiliary matrix has 32 + 180 = 212 states. Distinct faces cannot contribute duplicate amplitudes to the same ordered pair: every Wilson transition changes all four edges of its particular face. The builder's deduplication identifies the two Hermitian directions of a retained-to-retained pair; it does not accidentally drop a sum over different faces.

I also reviewed the analytic `13/2` omitted electric threshold in `su2_cube_bounds_note.md`. Its marked-edge argument excludes labels at least three, and its cycle/triangle analysis for labels 0,1,2 is valid. The threshold applies to every state outside P and therefore safely bounds the subset R.

## Exact rounding and probe directions

For every signed radical amplitude, integer binary search independently reproduced the proposed rational grid point and its exact square inequalities. This check does not call the certificate's integer-square-root routine. The resulting symmetric error bounds are

- `eta_A = 0`;
- `eta_L = 1/125000000000 = 8e-12`.

The maximum absolute row-error sum bounds the spectral norm because the error matrix is symmetric. Hence `L >= Lhat - eta_L I` and `A <= Ahat + eta_A I` in operator order. The correct lower probe is consequently `lo + eta_L`; the correct upper probe is `up - eta_A`. The implementation uses these signs.

Every lower shifted probe is below every omitted diagonal electric energy. Eliminating that positive block preserves inertia and leaves the retained-size Schur matrix. I independently eliminated the complete rational 212-by-212 matrix using ordinary Fraction Schur steps, ordering omitted indices first, rather than calling either certificate inertia helper or Bareiss. All four saved probe inertias matched exactly. The upper E1 probe has three negative eigenvalues, which is sufficient: the required condition is at least two, not exactly two.

The energy-to-gap subtraction is correct. Floating-point eigensolver results are not needed to justify these saved rational endpoints.

## Generalization to complete electric-energy sublevel bases

The same construction is valid for any finite physical retained set, provided all one-face omitted targets are included. The new complete sublevel bases also pass independent verification:

| Retain every state with T at most | Physical states |
|---:|---:|
| 8 | 86 |
| 10 | 212 |
| 12 | 398 |

A separate fixed-edge-order enumeration, using only the remaining exact integer energy budget and completed-vertex triangle/parity tests, reproduces all three complete sorted basis lists. It does not use the candidate's vertex-completion pruning.

The per-link bound follows from `q(q+2) <= floor(4 cap)`, equivalent to `(q+1)^2 <= floor(4 cap)+1`. The vertex pruning is also sound: each local completion minimum underestimates the corresponding true remaining local cost; summing and dividing by two compensates for each unknown edge's two endpoints.

Because every T is a quarter-integer, the proposed remainder threshold `cap + 1/4` for integer caps is valid. The cube's parity constraints allow an even stronger half-integer statement, but that improvement is not required. The shifted-probe and Schur logic extends unchanged; its explicit check that the lower shifted probe lies below each omitted diagonal must remain.

The complete energy-basis certificates now use the separately reviewed outward-rounded integer interval-inertia engine. Its independent audit, `work/su2_interval_inertia_independent_audit.md`, recomputes all four actual cap-8 probe counts against Bareiss, with exact agreement. The completed saved cap-10 and cap-12 intervals are [2.889649, 2.979921] and [2.941733, 2.958693]. Those larger Bareiss runs were not independently repeated. The independently verified full basis lists include the strongest certificate's cap-12 space of 398 states.

## Reproduction and limits

- `work/su2_cube_certificate_audit_check.py` reproduces the independent checks.
- `outputs/SU2_Cube_Certificate_Independent_Audit.json` records them.
- Audited implementations: `work/su2_cube_rational_certificate.py`, `work/su2_cube_energy_certificate.py`, and the energy-basis construction in `work/su2_cube_model.py`.

The certification remains conditional on the stated physical spin-network basis and recoupling identities. The independent loop-character check separately verifies all retained J = 1/2 Wilson blocks; it does not independently rederive the higher-spin transition formula. The generic builder assumes kappa = 1 and nu >= 0; a negative nu cannot be inserted into its squared-amplitude construction without changing its sign logic. The actual certificates use kappa = nu = 1. No audited implementation was modified.
