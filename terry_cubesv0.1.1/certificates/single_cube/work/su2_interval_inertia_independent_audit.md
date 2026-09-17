# Independent audit of the interval-inertia certificate engine

**The algorithm and all four actual cap-8 cube probe comparisons pass.** Its interval arithmetic supplies rigorous pivot-sign decisions; it does not use floating-point eigenvalue estimates to certify inertia. No mathematical blocker was found.

## Arithmetic enclosure

Let S = 2^bits. An integer interval [a0,a1] represents [a0/S,a1/S]. For a Schur correction a*b/p, the corresponding fixed-point numerator is `A*B/P`: the scale factors cancel exactly. `_product_quotient` first encloses the product by the minimum and maximum of its four endpoint products, then divides those endpoints by both denominator endpoints. When the denominator interval excludes zero, these corner ratios enclose the full quotient range. Integer floor and ceiling round outward, including when the denominator is negative. The early return for an exactly zero factor is valid.

The separate multiplication and division helpers have the required scale factors. Their floor/ceiling rules remain valid for signed numerators and denominators. Treating repeated appearances of the same uncertain entry as independent can widen an interval, but cannot exclude the true value.

The independent checker tested 7,680 exact rational endpoint/interior combinations, including positive and negative pivot intervals, against `_product_quotient`.

## Elimination invariant

The active symmetric interval block encloses the exact current Schur complement. A diagonal interval strictly above or below zero certifies the corresponding 1-by-1 pivot's sign. The update `A_ij - A_ik*A_jk/pivot`, with outward-rounded correction, maintains the enclosure. Symmetric permutations preserve inertia and the invariant. Exactly zero coupling intervals can safely skip an update.

If no 1-by-1 pivot is certifiable, the algorithm can use a symmetric block `B = [[a,c],[c,b]]` whose determinant interval is strictly negative. Every enclosed B then has one positive and one negative eigenvalue. For outside coupling vectors `(x_i,y_i)`, the Schur correction is

`[b*x_i*x_j + a*y_i*y_j - c*(x_i*y_j + y_i*x_j)] / (a*b - c*c)`.

The implemented numerator, determinant, and division enclose this expression with the correct fixed-point scaling. Both pivot permutations retain the selected block. Counting one sign of each kind and proceeding with the enclosed complement is therefore valid by congruence.

If no pivot can be certified, the routine raises `UncertifiedInertia`. It never turns an interval containing zero into a claimed zero eigenvalue. A successful result consequently certifies a nonsingular matrix and reports zero nullity. This is appropriate for spectral probes separated from eigenvalues.

Additional independent checks exercised signed pivots, coupled indefinite 2-by-2 pivots, rejection of three exact singular inputs, and adaptive precision at 64, 128 and 256 bits for both signs of a diagonal entry of magnitude 2^-200.

## Actual cap-8 probe comparison

I reconstructed all four rational probe matrices from the saved cap-8 certificate. For lower probes, I assembled the exact Schur complement independently by updating every ordered neighbor pair. I then evaluated each retained 86-by-86 matrix with both the interval engine and exact integer Bareiss inertia.

| Probe | Negative Schur eigenvalues | Positive Schur eigenvalues | Comparison |
|---|---:|---:|---|
| E0 lower | 0 | 86 | Exact agreement |
| E0 upper | 1 | 85 | Exact agreement |
| E1 lower | 1 | 85 | Exact agreement |
| E1 upper | 2 | 84 | Exact agreement |

The 422 eliminated electric directions are positive at each lower probe. Adding them reproduces the saved 508-dimensional auxiliary inertia exactly. All counts also match the final certificate JSON. These checks reproduce the cap-8 gap interval **[2.236257, 3.068193]**. The larger cap-10 and cap-12 Bareiss calculations were not repeated.

## Reproduction and scope

- Script: `work/su2_interval_inertia_independent_audit.py`.
- Results: `outputs/SU2_Interval_Inertia_Independent_Audit.json`.
- Input certificates: `outputs/SU2_Cube_Energy_Certificates.json`.

The separate cube-certificate audit verifies radical error allowances, probe shifts, complete omitted couplings and physical operator ordering. Its independent energy-basis enumeration now covers caps 8, 10 and 12, giving 86, 212 and 398 states respectively. Those steps are necessary in addition to this matrix-inertia audit.

The saved cap-12 interval **[2.941733, 2.958693]** uses this same reviewed interval algorithm with the completed certified probes. It is a fixed-graph statement, not a continuum or infinite-volume theorem.
