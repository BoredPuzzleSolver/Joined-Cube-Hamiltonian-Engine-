"""Complete spectrum of a small retained basis, not of the full spin space."""
import json
from terry_cubes import JoinedCubesModel, SpectrumSolver, plaquette_operator, connected_correlator

if __name__ == "__main__":
    model = JoinedCubesModel(1)
    solver = SpectrumSolver(model, 6, backend="arpack")
    solver.levels = solver.build_hamiltonian().dimension
    eigen = solver.compute_eigensystem()
    observable = plaquette_operator(model, eigen.states, face_index=0)
    result = connected_correlator(eigen.spectrum.eigenvalues, eigen.eigenvectors,
                                 observable, [0, .5, 1, 2, 4],
                                 complete_spectrum=eigen.complete_retained_spectrum)
    print(json.dumps(result, indent=2, allow_nan=False))
