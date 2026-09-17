"""Reproduce exact certificates for the untruncated two-square SU(2) patch.

Python + NumPy; run this script from any directory. Sibling modules
su2_recoupling.py and su2_cutoff_bounds.py are required. Floating-point
eigenvalues are reported for comparison only. The fixed rational probes
below, checked by integer inertia, establish the certificates.
"""
from fractions import Fraction
from pathlib import Path
import json
import time

from su2_recoupling import hamiltonian, exact_congruence
from su2_cutoff_bounds import (
    outside_electric_threshold, numerical_cutoff_bounds,
    rational_auxiliary_matrices, certify_energy_probes, self_test,
)

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "outputs" / "Shared_SU2_Cutoff_Certificate_Checks.json"

BENCHMARKS = [
    {
        "max_twice": 4,
        "kappa": Fraction(1),
        "nu": Fraction(1),
        "lower_probes": ["1.83601099", "4.92797000"],
        "upper_probes": ["1.83601102", "4.92797006"],
    },
    {
        "max_twice": 4,
        "kappa": Fraction(1),
        "nu": Fraction(4),
        "lower_probes": ["5.79687356", "10.00606446"],
        "upper_probes": ["5.79692505", "10.00683266"],
    },
    {
        "max_twice": 6,
        "kappa": Fraction(1),
        "nu": Fraction(4),
        "lower_probes": ["5.79688776", "10.00631521"],
        "upper_probes": ["5.79688779", "10.00631533"],
    },
]


def main():
    self_test()
    rows = []
    for benchmark in BENCHMARKS:
        start = time.perf_counter()
        cutoff = benchmark["max_twice"]
        kappa, nu = benchmark["kappa"], benchmark["nu"]
        states, h = hamiltonian(cutoff+1, float(kappa), float(nu))
        exact_states, congruence = exact_congruence(cutoff+1, kappa, nu)
        assert states == exact_states
        p_mask = [max(state) <= cutoff for state in states]
        weights = [(a+1)*(b+1)*(c+1) for a, b, c in states]
        electric = [
            kappa*Fraction(3*a*(a+2)+3*b*(b+2)+c*(c+2), 4)
            for a, b, c in states
        ]
        a, lower, a_weights = rational_auxiliary_matrices(
            congruence, electric, weights, p_mask
        )
        probes_low = [Fraction(x) for x in benchmark["lower_probes"]]
        probes_high = [Fraction(x) for x in benchmark["upper_probes"]]
        certificate = certify_energy_probes(
            a, lower, a_weights, weights, probes_low, probes_high,
            outside_electric_threshold(cutoff+1, kappa),
        )
        if not certificate["all_certified"]:
            raise AssertionError("An exact inertia probe failed")
        rows.append({
            "kappa_exact": str(kappa),
            "nu_exact": str(nu),
            "retained_spin_cutoff": cutoff/2,
            "boundary_spin_cutoff": (cutoff+1)/2,
            "retained_dimension": sum(p_mask),
            "auxiliary_dimension": len(states),
            "first_omitted_electric_threshold_exact": str(
                outside_electric_threshold(cutoff, kappa)
            ),
            "remaining_electric_threshold_exact": str(
                outside_electric_threshold(cutoff+1, kappa)
            ),
            "certificate": certificate,
            "floating_point_comparison": numerical_cutoff_bounds(
                h, [float(t) for t in electric], p_mask, cutoff, float(kappa)
            ),
            "elapsed_seconds": time.perf_counter()-start,
        })
    output = {
        "title": "Exact spin-cutoff certificates for a fixed shared-edge SU(2) patch",
        "scope": (
            "The untruncated gauge-invariant Hilbert space of two square "
            "plaquettes sharing one link, with open boundary and no external "
            "charges. All link spins are allowed in the bounded operator "
            "comparison. The two cutoffs define finite comparison matrices."
        ),
        "hamiltonian": (
            "H = kappa[3jL(jL+1)+3jR(jR+1)+jM(jM+1)] "
            "+ nu(2-WA-WB), WA=Tr U_A/2, WB=Tr U_B/2."
        ),
        "assumptions": [
            "kappa>0 and nu>=0, so the compressed omitted magnetic potential is positive.",
            "The admissible triangle/parity spin-network basis is complete for this graph.",
            "The exact recoupling matrix elements are those in su2_recoupling.py.",
            "A fundamental Wilson loop changes affected spins only by plus or minus one half, so the complete next cutoff contains every P-to-Q coupling.",
            "No infinite spatial volume, lattice-spacing limit or four-dimensional continuum Yang-Mills claim is made.",
        ],
        "method": (
            "Variational upper levels and finite auxiliary lower levels; "
            "exact rational congruence and integer Bareiss inertia establish "
            "all displayed rational endpoints. Floating point selects or "
            "compares values and is not used as an error certificate."
        ),
        "benchmarks": rows,
    }
    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(json.dumps(output, indent=2, allow_nan=False)+"\n", encoding="utf-8")
    print(json.dumps({
        "output": str(OUT),
        "benchmarks": [
            {
                "nu": row["nu_exact"],
                "J": row["retained_spin_cutoff"],
                "gap_interval": row["certificate"]["gap_interval_decimal"],
                "all_certified": row["certificate"]["all_certified"],
                "elapsed_seconds": row["elapsed_seconds"],
            } for row in rows
        ],
    }, indent=2))


if __name__ == "__main__":
    main()
