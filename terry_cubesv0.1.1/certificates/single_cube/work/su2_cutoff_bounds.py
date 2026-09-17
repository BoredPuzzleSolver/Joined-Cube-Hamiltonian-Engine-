"""Truncation bounds for the actual two-plaquette SU(2) Hamiltonian.

NumPy values are estimates. Optional rational congruence/inertia probes
certify bounds without assuming an eigensolver rounding allowance.
"""
from fractions import Fraction
import math
import numpy as np


def as_fraction(value):
    return value if isinstance(value, Fraction) else Fraction(str(value))


def outside_electric_threshold(max_twice, kappa=Fraction(1)):
    """Exact minimum T outside max(2jL,2jR,2jM) <= max_twice."""
    if int(max_twice) != max_twice or max_twice < 0:
        raise ValueError("max_twice must be a nonnegative integer")
    m = int(max_twice)+1
    a, b = m//2, (m+1)//2
    return as_fraction(kappa)*Fraction(
        3*a*(a+2)+3*b*(b+2)+m*(m+2), 4
    )


def auxiliary_lower_matrix(extended_h, electric, p_mask):
    """P blocks remain exact; replace the entire omitted block by T."""
    h = np.asarray(extended_h, dtype=float)
    electric = np.asarray(electric, dtype=float)
    p = np.asarray(p_mask, dtype=bool)
    if h.shape != (len(p), len(p)) or electric.shape != (len(p),):
        raise ValueError("Matrix, electric diagonal and mask sizes disagree")
    if not np.isfinite(h).all() or not np.isfinite(electric).all():
        raise ValueError("Inputs must be finite")
    q = np.flatnonzero(~p)
    lower = h.copy()
    lower[np.ix_(q, q)] = np.diag(electric[q])
    return lower


def numerical_cutoff_bounds(extended_h, electric, p_mask, max_twice,
                           kappa=1, levels=2):
    """Require the complete J+1/2 basis, including all first-boundary states.

    The result evaluates analytic bounds using floating point. Exact rational
    probes are separate; no fixed epsilon is represented as a certificate.
    """
    h = np.asarray(extended_h, dtype=float)
    p = np.asarray(p_mask, dtype=bool)
    a = h[np.ix_(p, p)]
    if len(a) < levels:
        raise ValueError("Insufficient retained states")
    lower = auxiliary_lower_matrix(h, electric, p)
    hermitian_error = float(np.max(abs(h-h.T)))
    if hermitian_error > 1e-9:
        raise ValueError("Input matrix is not numerically symmetric")
    # This averaging only defines the floating-point matrix being evaluated.
    upper_values = np.linalg.eigvalsh((a+a.T)/2)[:levels]
    lower_values = np.linalg.eigvalsh((lower+lower.T)/2)[:levels]
    t_q = outside_electric_threshold(max_twice, kappa)
    t_rest = outside_electric_threshold(max_twice+1, kappa)
    lower_values = np.minimum(lower_values, float(t_rest))
    rows = [
        {"index": j, "lower_estimate": float(lo), "upper_estimate": float(up),
         "width_estimate": float(up-lo),
         "numerical_ordering_consistent": bool(lo <= up)}
        for j, (lo, up) in enumerate(zip(lower_values, upper_values))
    ]
    # Diagnostic only: this heuristic identifies resolution limits; it is
    # expressly not a rigorous eigenvalue error or truncation allowance.
    resolution_indicator = float(
        64*np.finfo(float).eps*len(h)*max(1.0, float(np.linalg.norm(h, np.inf)))
    )
    for row in rows:
        row["roundoff_limited"] = bool(
            abs(row["width_estimate"]) <= resolution_indicator
        )
    result = {
        "retained_dimension": len(a),
        "extended_dimension": len(h),
        "first_omitted_electric_threshold_exact": str(t_q),
        "remaining_electric_threshold_exact": str(t_rest),
        "hermiticity_error": hermitian_error,
        "energy_intervals_evaluated_numerically": rows,
        "numerical_ordering_consistent": all(
            row["numerical_ordering_consistent"] for row in rows
        ),
        "roundoff_limited": any(row["roundoff_limited"] for row in rows),
        "roundoff_diagnostic_scale": resolution_indicator,
        "roundoff_diagnostic_is_error_bound": False,
        "machine_certified": False,
        "reason": "Eigenvalue evaluations and matrix inputs use floating point.",
    }
    if levels >= 2:
        result["gap_interval_evaluated_numerically"] = [
            float(lower_values[1]-upper_values[0]),
            float(upper_values[1]-lower_values[0]),
        ]
    return result


def rational_auxiliary_matrices(extended_congruence, electric, weights, p_mask):
    """Return rational D^-1 A D^-1 and D^-1 L D^-1, and A weights.

    weights[i] = (2jL+1)(2jR+1)(2jM+1), so D_ii=sqrt(weights[i]).
    extended_congruence must have zero spectral shift.
    """
    c = [[as_fraction(x) for x in row] for row in extended_congruence]
    n = len(c)
    if any(len(row) != n for row in c):
        raise ValueError("Congruence must be square")
    p = [i for i in range(n) if p_mask[i]]
    q = [i for i in range(n) if not p_mask[i]]
    lower = [row.copy() for row in c]
    for i in q:
        for j in q:
            lower[i][j] = (as_fraction(electric[i])/weights[i]
                           if i == j else Fraction(0))
    a = [[c[i][j] for j in p] for i in p]
    return a, lower, [weights[i] for i in p]


def rational_inertia(matrix):
    """Exact inertia using integer symmetric fraction-free elimination.

    Clear one positive common denominator. Bareiss pivots encode LDL pivot
    signs via current_pivot / previous_pivot. Symmetric permutations and,
    for a zero-diagonal block, one integer congruence avoid zero pivots.
    """
    n = len(matrix)
    if any(len(row) != n for row in matrix):
        raise ValueError("Matrix must be square")
    fractions = [[as_fraction(x) for x in row] for row in matrix]
    if any(fractions[i][j] != fractions[j][i]
           for i in range(n) for j in range(i)):
        raise ValueError("Exact matrix must be symmetric")
    denominator = 1
    for row in fractions:
        for x in row:
            denominator = math.lcm(denominator, x.denominator)
    a = [[x.numerator*(denominator//x.denominator) for x in row]
         for row in fractions]
    positive = negative = zero = 0
    previous = 1

    def symmetric_swap(i, j):
        if i == j:
            return
        a[i], a[j] = a[j], a[i]
        for row in a:
            row[i], row[j] = row[j], row[i]

    for k in range(n):
        if a[k][k] == 0:
            candidate = next((i for i in range(k+1, n) if a[i][i]), None)
            if candidate is not None:
                symmetric_swap(k, candidate)
            else:
                pair = next(((i, j) for i in range(k, n)
                             for j in range(i+1, n) if a[i][j]), None)
                if pair is None:
                    zero += n-k
                    break
                i, j = pair
                symmetric_swap(k, i)
                # Add row j to row k and column j to column k. This is
                # an invertible integer congruence; its new diagonal is
                # twice a previously nonzero off-diagonal element.
                for col in range(k, n):
                    a[k][col] += a[j][col]
                for row in range(k, n):
                    a[row][k] += a[row][j]
        pivot = a[k][k]
        if pivot == 0:
            raise ArithmeticError("Internal zero-pivot failure")
        if pivot*previous > 0:
            positive += 1
        else:
            negative += 1
        for i in range(k+1, n):
            for j in range(i, n):
                numerator = pivot*a[i][j]-a[i][k]*a[k][j]
                value, remainder = divmod(numerator, previous)
                if remainder:
                    raise ArithmeticError("Bareiss division was not exact")
                a[i][j] = a[j][i] = value
        for i in range(k+1, n):
            a[i][k] = a[k][i] = 0
        previous = pivot
    return {"negative": negative, "zero": zero, "positive": positive}


def shifted_congruence(base_congruence, weights, probe):
    probe = as_fraction(probe)
    result = [[as_fraction(x) for x in row] for row in base_congruence]
    if len(weights) != len(result):
        raise ValueError("Metric weights and matrix sizes disagree")
    for i, weight in enumerate(weights):
        if weight <= 0:
            raise ValueError("Congruence weights must be positive")
        result[i][i] -= probe/weight
    return result


def certify_energy_probes(a_congruence, lower_congruence, a_weights,
                         lower_weights, lower_probes, upper_probes,
                         remainder_threshold):
    """Certify full-space energy intervals via exact inertia at rational probes.

    Returns proof checks, not rounded eigenvalue error guesses. A lower
    probe ell_j is valid if L has at most j eigenvalues strictly below it
    and the omitted electric remainder has none below it. An upper probe
    u_j is valid if A has at least j+1 eigenvalues below it.
    """
    rows = []
    threshold = as_fraction(remainder_threshold)
    for j, (lo, up) in enumerate(zip(lower_probes, upper_probes)):
        lo, up = as_fraction(lo), as_fraction(up)
        low_inertia = rational_inertia(shifted_congruence(
            lower_congruence, lower_weights, lo
        ))
        up_inertia = rational_inertia(shifted_congruence(
            a_congruence, a_weights, up
        ))
        low_pass = low_inertia["negative"] <= j and lo <= threshold
        upper_pass = up_inertia["negative"] >= j+1
        rows.append({
            "index": j, "lower_exact": str(lo), "upper_exact": str(up),
            "lower_decimal": float(lo), "upper_decimal": float(up),
            "lower_auxiliary_inertia": low_inertia,
            "upper_projection_inertia": up_inertia,
            "lower_pass": low_pass, "upper_pass": upper_pass,
            "certified": low_pass and upper_pass and lo <= up,
        })
    result = {
        "method": "Exact rational congruence plus integer Bareiss inertia.",
        "energy_intervals": rows,
        "all_certified": all(row["certified"] for row in rows),
        "scope": (
            "Certificate assumes the supplied exact congruence is the stated "
            "Hamiltonian and that the complete first omitted layer contains "
            "all P-to-Q magnetic couplings."
        ),
    }
    if len(rows) >= 2:
        gap_lo = as_fraction(rows[1]["lower_exact"])-as_fraction(rows[0]["upper_exact"])
        gap_up = as_fraction(rows[1]["upper_exact"])-as_fraction(rows[0]["lower_exact"])
        result["gap_interval_exact"] = [str(gap_lo), str(gap_up)]
        result["gap_interval_decimal"] = [float(gap_lo), float(gap_up)]
    return result


def outward_probes(numerical_bounds, decimal_digits=8):
    """Choose convenient rational probes; exact inertia decides validity."""
    denominator = 10**decimal_digits
    rows = numerical_bounds["energy_intervals_evaluated_numerically"]
    lower = [Fraction(math.floor(r["lower_estimate"]*denominator)-1, denominator)
             for r in rows]
    upper = [Fraction(math.ceil(r["upper_estimate"]*denominator)+1, denominator)
             for r in rows]
    return lower, upper


def self_test():
    cases = [
        ([[0, 1], [1, 0]], (1, 0, 1)),
        ([[0, 0, 0], [0, 2, 0], [0, 0, -3]], (1, 1, 1)),
        ([[1, 2, 3], [2, 4, 6], [3, 6, 9]], (0, 2, 1)),
        ([[Fraction(1, 3), Fraction(2, 7)],
          [Fraction(2, 7), Fraction(-4, 5)]], (1, 0, 1)),
    ]
    for matrix, expected in cases:
        got = rational_inertia(matrix)
        assert (got["negative"], got["zero"], got["positive"]) == expected
    rng = np.random.default_rng(73191)
    for n in range(2, 11):
        for _ in range(4):
            matrix = rng.integers(-6, 7, size=(n, n))
            matrix = matrix+matrix.T
            eigenvalues = np.linalg.eigvalsh(matrix.astype(float))
            got = rational_inertia(matrix.tolist())
            assert got["negative"] == int(np.sum(eigenvalues < -1e-8))
            assert got["positive"] == int(np.sum(eigenvalues > 1e-8))
    for cutoff in range(7):
        candidates = [
            Fraction(3*a*(a+2)+3*b*(b+2)+c*(c+2), 4)
            for a in range(cutoff+5) for b in range(cutoff+5)
            for c in range(cutoff+5)
            if max(a,b,c)>cutoff and abs(a-b)<=c<=a+b
            and (a+b+c)%2 == 0
        ]
        assert min(candidates) == outside_electric_threshold(cutoff)
    print("Exact inertia and omitted-electric-threshold checks passed.")


if __name__ == "__main__":
    self_test()
