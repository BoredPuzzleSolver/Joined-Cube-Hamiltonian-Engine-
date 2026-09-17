"""Cutoff comparisons on the unchanged joined-cube Hamiltonian."""
import numpy as np

from ._version import __version__
from .hamiltonian import JoinedCubesModel, validate_electric_cutoff, _integer
from .solver import SpectrumSolver


def cutoff_study(model, cutoffs, *, backend="arpack", tolerance=1e-10,
                 max_iterations=1000, seed=1097, second_seed=8917,
                 embedding_cutoff=None):
    """Compare E0/E1 across complete electric cutoffs at one coupling.

    Results include direct residuals and residuals in the largest test
    basis. Neither diagnostic bounds the full omitted space. Cutoffs are
    sorted and deduplicated exactly; their meaning is independent of kappa.
    Only two low vectors are retained for each solve.
    """
    if not isinstance(model, JoinedCubesModel):
        raise TypeError("model must be a JoinedCubesModel")
    caps = sorted({validate_electric_cutoff(c) for c in cutoffs})
    if not caps or caps[0] < 3:
        raise ValueError("cutoffs must contain values >= 3")
    largest = caps[-1] if embedding_cutoff is None else validate_electric_cutoff(embedding_cutoff)
    if largest < caps[-1]:
        raise ValueError("embedding_cutoff must be at least the largest requested cutoff")
    if second_seed is not None:
        second_seed = _integer(second_seed, "second_seed", 0)
        if second_seed == seed:
            raise ValueError("second_seed must differ from seed")
    settings = dict(backend=backend, tolerance=tolerance, max_iterations=max_iterations,
                    seed=seed, levels=2)
    largest_solver = SpectrumSolver(model, largest, **settings)
    outer = largest_solver.build_hamiltonian()
    outer_matrix = outer.to_scipy()
    lookup = {state: i for i, state in enumerate(outer.states)}
    previous = None
    rows = []
    for cap in caps:
        solver = largest_solver if cap == largest else SpectrumSolver(model, cap, **settings)
        eigen = solver.compute_eigensystem()
        result = eigen.spectrum
        row = result.to_dict()
        energies = np.array(result.eigenvalues[:2])
        changes = None if previous is None else energies - previous
        row["change_in_e0"] = None if changes is None else float(changes[0])
        row["change_in_e1"] = None if changes is None else float(changes[1])
        allowance = max(10 * tolerance, 1e-9 * max(model.kappa, model.nu))
        row["nested_energy_ordering_consistent"] = bool(changes is None or np.all(changes <= allowance))
        if cap < largest:
            embedded = np.zeros((outer.dimension, 2))
            indices = [lookup[s] for s in eigen.states]
            embedded[indices] = eigen.eigenvectors[:, :2]
            residuals = np.linalg.norm(outer_matrix @ embedded - embedded * energies, axis=0)
            row["residual_in_largest_test_basis"] = residuals.tolist()
        if cap == caps[-1] and second_seed is not None:
            # Reuse the largest requested matrix rather than enumerate it again.
            solver.seed = second_seed
            repeated = solver.compute_mass_gap()
            row["second_seed"] = second_seed
            row["second_seed_max_energy_difference"] = float(
                np.max(np.abs(energies - np.array(repeated.eigenvalues[:2]))))
            row["second_seed_residual_norms"] = list(repeated.residual_norms)
        rows.append(row)
        previous = energies
    return {"schema_version": 1, "package_version": __version__, "model_id": "joined_cubes",
            "num_cubes": model.num_cubes, "embedding_electric_cutoff": str(largest),
            "embedding_dimension": outer.dimension, "results": rows, "certified": False,
            "scope": "Numerical cutoff comparison at fixed finite geometry. Larger-basis residuals and changes in gaps are not untruncated error bounds."}
