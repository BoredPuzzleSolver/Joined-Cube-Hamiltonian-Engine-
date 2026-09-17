"""Isolated verification worker. Only fixed, packaged reference inputs are used."""
from fractions import Fraction as F
from hashlib import sha256
from pathlib import Path
import contextlib
import io
import json
import sys
import time


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def check_manifest(root):
    count = 0
    for line in (root / "SHA256SUMS.txt").read_text().splitlines():
        if line.strip():
            expected, name = line.split(None, 1)
            path = root / name.strip()
            require(path.is_file() and sha256(path.read_bytes()).hexdigest() == expected,
                    "Certificate source manifest mismatch: " + name)
            count += 1
    return count


def single_cube(root, cutoff):
    import su2_cube_model as reference
    from su2_cube_rational_certificate import comparison_matrices, electric, inertia_a, inertia_lower
    from su2_interval_inertia import interval_inertia
    from terry_cubes.hamiltonian import JoinedCubesModel, electric_exact, energy_basis, state_face_transitions

    saved = json.loads((root / "reference_results/SU2_Cube_Energy_Certificates.json").read_text())
    record = next(r for r in saved["certificates"] if r["electric_energy_cap_over_kappa"] == cutoff)
    retained = reference.energy_basis(F(cutoff))
    graph = JoinedCubesModel(1).graph
    require(retained == energy_basis(graph, F(cutoff)), "Package and certificate basis order differ")
    require(graph.physical_edge_count == 12 and graph.face_count == 6 and not graph.splits,
            "Certificate bridge requires the original single cube")
    checked = 0
    phase_constraints = {}
    for state in retained:
        for face_index, face in enumerate(reference.FACE_DATA):
            old = {s: (phase, square) for s, phase, square in reference.state_face_transitions(state, face)}
            new = {s: (phase, square) for s, phase, square in state_face_transitions(graph, state, face_index)}
            require(old.keys() == new.keys(), "Package and certificate transition support differ")
            for target, (old_sign, old_square) in old.items():
                new_sign, new_square = new[target]
                require(old_square == new_square, "Package and certificate exact amplitudes differ")
                ratio = old_sign * new_sign
                require(ratio in (-1, 1), "Unexpected nonzero transition phase")
                phase_constraints.setdefault(state, []).append((target, ratio))
                phase_constraints.setdefault(target, []).append((state, ratio))
            checked += len(old)
    # The original cube and joined-graph basis conventions may use different
    # state signs. One common diagonal unitary must identify ALL face blocks,
    # including every P-to-Q target; separate per-face phases are insufficient.
    phases = {}
    for start in phase_constraints:
        if start in phases:
            continue
        phases[start] = 1
        queue = [start]
        for state in queue:
            for target, ratio in phase_constraints[state]:
                wanted = phases[state] * ratio
                if target in phases:
                    require(phases[target] == wanted, "No common package-to-certificate basis phase exists")
                else:
                    phases[target] = wanted
                    queue.append(target)
    states, diagonal, matrix, eta_a, eta_l = comparison_matrices(retained)
    p = len(retained)
    require(all(electric_exact(graph, state) == electric(state) for state in states),
            "Package and certificate exact electric energies differ")
    require(all(value == electric_exact(graph, state) + (6 if i < p else 0)
                for i, (state, value) in enumerate(zip(states, diagonal))),
            "Certificate comparison diagonal differs from the preserved construction")
    require(p == record["retained_dimension"] and len(states) == record["auxiliary_dimension"],
            "Certificate dimensions differ from the published construction")
    require(eta_a == F(record["A_error_bound_exact"]) and eta_l == F(record["L_error_bound_exact"]),
            "Radical-rounding error bounds differ")
    threshold = F(4*cutoff+1, 4)
    require(threshold == F(record["remaining_electric_threshold_exact"]), "Remainder threshold differs")
    results = []
    for j, row in enumerate(record["energy_intervals"]):
        lower, upper = F(row["lower_exact"]), F(row["upper_exact"])
        low = inertia_lower(diagonal, matrix, p, lower+eta_l, counter=interval_inertia)
        high = inertia_a(diagonal, matrix, p, upper-eta_a, counter=interval_inertia)
        require(low["negative"] <= j and lower <= threshold and high["negative"] >= j+1,
                "A fixed energy endpoint failed its exact inertia test")
        for key in ("negative", "positive", "zero"):
            require(low[key] == row["lower_inertia"][key] and high[key] == row["upper_inertia"][key],
                    "An inertia count differs from the reference")
        results.append({"index": j, "lower_exact": str(lower), "upper_exact": str(upper),
                        "lower_decimal": float(lower), "upper_decimal": float(upper),
                        "lower_inertia": low, "upper_inertia": high})
    gap = [F(results[1]["lower_exact"])-F(results[0]["upper_exact"]),
           F(results[1]["upper_exact"])-F(results[0]["lower_exact"])]
    require(gap == [F(x) for x in record["gap_interval_exact"]], "Gap subtraction differs")
    return {"model_id": "single_cube", "kappa_exact": "1", "nu_exact": "1",
            "electric_cutoff": cutoff, "retained_dimension": p, "comparison_dimension": len(states),
            "remaining_electric_threshold_exact": str(threshold),
            "A_error_bound_exact": str(eta_a), "L_error_bound_exact": str(eta_l),
            "energy_intervals": results, "gap_interval_exact": [str(x) for x in gap],
            "gap_interval_decimal": [float(x) for x in gap],
            "package_transitions_checked_including_omitted_targets": checked,
            "package_electric_energies_checked": len(states),
            "package_bridge": "Exact support and squared amplitudes, with one common diagonal state-sign transformation across every face and retained-to-omitted coupling; exact electric energies and comparison diagonal checked",
            "package_bridge_verified": True, "certified": True,
            "scope": "Full spin space of one fixed open SU(2) cube at kappa=nu=1, under the preserved basis, recoupling, positivity and omitted-state assumptions. No multi-cube, volume or continuum certificate."}


def two_plaquette(root):
    import verify_su2_cutoff_certificates as reference
    from terry_cubes.recoupling import exact_congruence

    for spec in reference.BENCHMARKS:
        args = spec["max_twice"]+1, spec["kappa"], spec["nu"]
        states, matrix = reference.exact_congruence(*args)
        current_states, current_matrix = exact_congruence(*args)
        require(states == current_states and matrix == current_matrix,
                "Package and two-plaquette exact congruence differ")
    with contextlib.redirect_stdout(io.StringIO()):
        reference.main()
    result = json.loads(reference.OUT.read_text())
    saved = json.loads((root / "reference_results/Shared_SU2_Cutoff_Certificate_Checks.json").read_text())
    require(len(result["benchmarks"]) == len(saved["benchmarks"]), "Missing two-plaquette benchmark")
    for row, old in zip(result["benchmarks"], saved["benchmarks"]):
        require(row["certificate"]["all_certified"], "An exact two-plaquette probe failed")
        require(row["certificate"]["gap_interval_exact"] == old["certificate"]["gap_interval_exact"],
                "Two-plaquette gap endpoints differ")
    result.update(model_id="two_plaquette", package_bridge_verified=True, certified=True)
    return result


def main():
    require(sys.flags.optimize == 0, "Reference assertions must be enabled")
    kind, directory, cap = sys.argv[1:]
    root = Path(directory)
    verified_files = check_manifest(root)
    # The preserved scripts resolve only inside their temporary reference bundle.
    sys.path.insert(0, str(root / "work"))
    sys.path.insert(1, str(Path(__file__).resolve().parent.parent))
    from terry_cubes._version import __version__
    start = time.perf_counter()
    if kind == "single_cube":
        result = single_cube(root, int(cap))
    elif kind == "two_plaquette":
        result = two_plaquette(root)
    else:
        raise ValueError("Unsupported certificate model")
    result.update(package_version=__version__, source_manifest_files_verified=verified_files,
                  elapsed_seconds=time.perf_counter()-start,
                  verification="Fresh matrix reconstruction, package equivalence checks and exact endpoint probes")
    print(json.dumps(result, allow_nan=False))


if __name__ == "__main__":
    main()
