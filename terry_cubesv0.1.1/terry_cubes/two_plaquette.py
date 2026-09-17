"""Public interface to the preserved seven-link two-square benchmark."""
from dataclasses import dataclass
import numpy as np

from ._version import __version__
from .hamiltonian import _finite_float, _integer
from .recoupling import wilson_matrices, electric_diagonal
from .observables import connected_correlator
from .solver import SpectrumConvergenceError


@dataclass(frozen=True)
class TwoPlaquetteSpectrum:
    model: "TwoPlaquetteModel"
    max_twice: int
    states: tuple
    eigenvalues: np.ndarray
    eigenvectors: np.ndarray
    residual_norms: np.ndarray
    wilson_a: np.ndarray
    wilson_b: np.ndarray

    @property
    def gap(self):
        return float(self.eigenvalues[1] - self.eigenvalues[0])

    def correlators(self, times=(0.25, 0.5, 1, 2, 4, 8), *, dt=0.25, weight_threshold=1e-18):
        times = tuple(times)
        return {name: connected_correlator(
            self.eigenvalues, self.eigenvectors, operator, times, dt=dt,
            weight_threshold=weight_threshold, complete_spectrum=True)
            for name, operator in (("symmetric", (self.wilson_a+self.wilson_b)/np.sqrt(2)),
                                   ("antisymmetric", (self.wilson_a-self.wilson_b)/np.sqrt(2)))}

    def to_dict(self):
        return {"schema_version": 1, "package_version": __version__, "model_id": "two_plaquette",
                "geometry": {"physical_links": 7, "physical_vertices": 6, "plaquettes": 2},
                "kappa": self.model.kappa, "nu": self.model.nu,
                "max_twice": self.max_twice, "max_spin": self.max_twice/2,
                "cutoff_definition": "Each physical spin <= max_twice/2; this is not an electric-energy cutoff",
                "state_labels": "(2*jL, 2*jR, 2*jM), lexicographic order",
                "dimension": len(self.states), "e0": float(self.eigenvalues[0]),
                "e1": float(self.eigenvalues[1]), "gap": self.gap,
                "eigenvalues": self.eigenvalues.tolist(), "residual_norms": self.residual_norms.tolist(),
                "method": "dense NumPy eigh", "certified": False,
                "scope": "Numerical finite-spin projection of two squares sharing one link. This is not a joined-cube graph or an untruncated certificate."}


@dataclass(frozen=True)
class TwoPlaquetteModel:
    """Two squares sharing one link; the original (L,R,M) labeling is retained."""
    kappa: float = 1.0
    nu: float = 1.0

    def __post_init__(self):
        object.__setattr__(self, "kappa", _finite_float(self.kappa, "kappa", True))
        object.__setattr__(self, "nu", _finite_float(self.nu, "nu"))

    def compute_spectrum(self, max_twice=12, *, tolerance=1e-10):
        max_twice = _integer(max_twice, "max_twice", 1)
        tolerance = _finite_float(tolerance, "tolerance", True)
        states, a, b = wilson_matrices(max_twice)
        matrix = np.diag(electric_diagonal(states, self.kappa)+2*self.nu)-self.nu*(a+b)
        if not np.all(np.isfinite(matrix)):
            raise ValueError("Hamiltonian entries exceed floating-point range; rescale kappa and nu")
        values, vectors = np.linalg.eigh(matrix)
        residuals = np.linalg.norm(matrix @ vectors-vectors*values, axis=0)
        if not np.all(np.isfinite(values)) or np.max(residuals) > tolerance:
            raise SpectrumConvergenceError("The retained two-plaquette residual did not meet tolerance")
        return TwoPlaquetteSpectrum(self, max_twice, tuple(states), values, vectors, residuals, a, b)
