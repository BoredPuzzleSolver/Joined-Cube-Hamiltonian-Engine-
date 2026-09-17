"""NumPy-only sparse symmetric spectra with explicit Ritz residual checks.

These are numerical approximations, not certified eigenvalue enclosures.
The COO arrays contain both directions of every off-diagonal entry.
"""
import numpy as np


class SymmetricOperator:
    def __init__(self, diagonal, rows, columns, entries):
        self.diagonal = np.asarray(diagonal, dtype=float)
        self.rows = np.asarray(rows, dtype=np.int64)
        self.columns = np.asarray(columns, dtype=np.int64)
        self.entries = np.asarray(entries, dtype=float)
        self.size = len(self.diagonal)
        assert self.rows.shape == self.columns.shape == self.entries.shape
        assert not len(self.rows) or (min(self.rows.min(), self.columns.min()) >= 0
                                    and max(self.rows.max(), self.columns.max()) < self.size)

    def matvec(self, vector):
        return self.diagonal*vector + np.bincount(
            self.rows, weights=self.entries*vector[self.columns], minlength=self.size)

    def dense(self):
        result = np.diag(self.diagonal)
        np.add.at(result, (self.rows, self.columns), self.entries)
        return result


def smallest_ritz(operator, levels=3, max_iterations=240, tolerance=1e-10, seed=1097):
    """Fully reorthogonalized symmetric Lanczos; report direct residual norms.

    A small residual does not prove that no unobserved eigenvalue lies below
    a computed Ritz value. Exact enclosure verification is a separate step.
    """
    n = operator.size
    if n <= 150:
        values, vectors = np.linalg.eigh(operator.dense())
        residuals = [float(np.linalg.norm(operator.matvec(vectors[:,i])-values[i]*vectors[:,i]))
                     for i in range(levels)]
        return values[:levels], vectors[:,:levels], {
            'method':'dense NumPy eigh', 'iterations':n,
            'residual_norms':residuals, 'converged':max(residuals) <= tolerance,
        }
    count = min(n, max_iterations)
    basis = np.empty((n, count+1), dtype=float)
    rng = np.random.default_rng(seed)
    basis[:,0] = rng.standard_normal(n)
    basis[:,0] /= np.linalg.norm(basis[:,0])
    alpha, beta = [], []
    converged = False
    for k in range(count):
        w = operator.matvec(basis[:,k])
        alpha.append(float(basis[:,k] @ w))
        w -= alpha[-1]*basis[:,k]
        if k:
            w -= beta[-1]*basis[:,k-1]
        for _ in range(2):
            q = basis[:,:k+1]
            w -= q @ (q.T @ w)
        norm = float(np.linalg.norm(w))
        should_check = k+1 >= levels and ((k+1)%12 == 0 or k+1 == count or norm < 1e-13)
        if should_check:
            tridiagonal = np.diag(alpha)
            if k:
                tridiagonal += np.diag(beta,1)+np.diag(beta,-1)
            values, small_vectors = np.linalg.eigh(tridiagonal)
            vectors = basis[:,:k+1] @ small_vectors[:,:levels]
            # Rayleigh-Ritz on the actual Krylov basis removes accumulated
            # recurrence roundoff from the returned low subspace.
            images = np.column_stack([operator.matvec(vectors[:,i]) for i in range(levels)])
            reduced = vectors.T @ images
            values, rotation = np.linalg.eigh((reduced+reduced.T)/2)
            vectors = vectors @ rotation
            images = images @ rotation
            residuals = np.linalg.norm(images-vectors*values, axis=0)
            converged = bool(np.max(residuals) <= tolerance)
            if converged or norm < 1e-13:
                break
        if k+1 < count:
            beta.append(norm)
            basis[:,k+1] = w/norm
    return values[:levels], vectors[:,:levels], {
        'method':'fully reorthogonalized Lanczos with direct Ritz residuals',
        'iterations':k+1, 'residual_norms':residuals.tolist(),
        'converged':converged, 'seed':seed,
        'scope':'Numerical Ritz values; a residual check does not certify spectral ordering',
    }


def self_test():
    # A connected tridiagonal with known analytic spectrum checks the sparse
    # path against eigenvalues independent of our Lanczos implementation.
    n = 240
    rows = np.concatenate([np.arange(n-1),np.arange(1,n)])
    cols = np.concatenate([np.arange(1,n),np.arange(n-1)])
    op = SymmetricOperator(np.full(n,2.),rows,cols,np.full(len(rows),-1.))
    values, _, check = smallest_ritz(op,max_iterations=n,tolerance=1e-10)
    expected = 2-2*np.cos(np.arange(1,4)*np.pi/(n+1))
    assert np.max(abs(values-expected)) < 1e-11
    assert check['converged']
    print('Sparse Lanczos analytic-spectrum check passed.',flush=True)



# The routines above are preserved from su2_sparse_spectrum.py.
# The high-level interface below validates inputs and surfaces convergence.
from dataclasses import dataclass
from ._version import __version__

from .hamiltonian import (
    JoinedCubesModel, _integer, _finite_float,
    validate_electric_cutoff, build_hamiltonian,
)


@dataclass(frozen=True)
class SpectrumResult:
    """Numerical Ritz spectrum in a specified finite basis, not a certificate."""

    e0: float
    e1: float
    gap: float
    eigenvalues: tuple
    residual_norms: tuple
    converged: bool
    dimension: int
    num_cubes: int
    kappa: float
    nu: float
    electric_cutoff: str
    pairing: str
    backend: str
    method: str
    iterations: int
    seed: int
    tolerance: float
    max_iterations: int
    hermiticity_error: float

    def to_dict(self):
        """Return JSON-safe values and explicit numerical scope."""
        from dataclasses import asdict

        result = asdict(self)
        result["eigenvalues"] = list(self.eigenvalues)
        result["residual_norms"] = list(self.residual_norms)
        result["schema_version"] = 1
        result["package_version"] = __version__
        result["model_id"] = "joined_cubes"
        ratio = self.nu / self.kappa
        normalized_gap = self.gap / self.kappa
        result["nu_over_kappa"] = ratio if np.isfinite(ratio) else None
        result["gap_over_kappa"] = normalized_gap if np.isfinite(normalized_gap) else None
        if not np.isfinite(ratio) or not np.isfinite(normalized_gap):
            result["normalization_note"] = (
                "A normalized ratio exceeds floating-point range and is null; "
                "unnormalized energy values are retained."
            )
        result["certified"] = False
        result["scope"] = (
            "Numerical finite-basis spectrum. Residuals do not certify spectral "
            "ordering, omitted states, infinite volume, or the continuum limit. "
            "Single-vector Lanczos may omit copies of degenerate eigenvalues."
        )
        result["cutoff_definition"] = "sum_physical_links j(j+1) <= electric_cutoff"
        return result


@dataclass(frozen=True)
class Eigensystem:
    """Numerical vectors in the exact row order of ``states``.

    Partial spectra give partial observable spectral sums. Single-vector
    Lanczos does not guarantee multiplicities or spectral completeness.
    """

    spectrum: SpectrumResult
    states: tuple
    eigenvectors: np.ndarray
    complete_retained_spectrum: bool


class SpectrumConvergenceError(RuntimeError):
    """Raised when the direct retained-matrix residual check fails."""

    def __init__(self, message, result=None):
        super().__init__(message)
        self.result = result


class SpectrumSolver:
    """Compute low eigenvalues using the original Lanczos or SciPy ARPACK.

    Parameters
    ----------
    model : JoinedCubesModel
        Geometry and physical coefficients.
    electric_cutoff : int, float, str, or Fraction
        Complete cutoff on sum_physical j(j+1), independent of kappa.
    backend : {"lanczos", "arpack"}
        The default preserves the original fully reorthogonalized routine.
    tolerance : float
        Absolute direct-residual tolerance in the energy units of H.
    max_iterations : int
        Original Lanczos basis limit or ARPACK iteration limit.
    levels : int
        Requested low Ritz values, at least two. Single-vector Lanczos
        can omit repeated copies of degenerate eigenvalues.
    seed : int
        Nonnegative seed for the starting vector.
    """

    def __init__(self, model, electric_cutoff=16, *, backend="lanczos",
                 tolerance=1e-10, max_iterations=240, levels=3, seed=1097):
        if not isinstance(model, JoinedCubesModel):
            raise TypeError("model must be a JoinedCubesModel")
        if backend not in ("lanczos", "arpack"):
            raise ValueError("backend must be lanczos or arpack")
        self.model = model
        self.electric_cutoff = validate_electric_cutoff(electric_cutoff)
        self.backend = backend
        self.tolerance = _finite_float(tolerance, "tolerance", True)
        self.levels = _integer(levels, "levels", 2)
        self.max_iterations = _integer(max_iterations, "max_iterations")
        if self.max_iterations < self.levels:
            raise ValueError("max_iterations must be at least levels")
        self.seed = _integer(seed, "seed", 0)
        self._hamiltonian = None
        self._hamiltonian_key = None

    def build_hamiltonian(self):
        """Build once and retain the exact ordering used for eigenvectors."""
        if not isinstance(self.model, JoinedCubesModel):
            raise TypeError("model must be a JoinedCubesModel")
        cap = validate_electric_cutoff(self.electric_cutoff)
        key = (self.model, cap)
        if self._hamiltonian is None or key != self._hamiltonian_key:
            self._hamiltonian = build_hamiltonian(self.model, cap)
            self._hamiltonian_key = key
        return self._hamiltonian

    def compute_mass_gap(self):
        """Return E0, E1 and E1-E0; raise on insufficient states or nonconvergence."""
        return self._compute(retain_vectors=False)

    def compute_eigensystem(self):
        """Return the same spectrum and its vectors for observable calculations."""
        return self._compute(retain_vectors=True)

    def _compute(self, retain_vectors):
        # Revalidate public settings so notebook parameter changes are safe.
        if self.backend not in ("lanczos", "arpack"):
            raise ValueError("backend must be lanczos or arpack")
        self.tolerance = _finite_float(self.tolerance, "tolerance", True)
        self.levels = _integer(self.levels, "levels", 2)
        self.max_iterations = _integer(self.max_iterations, "max_iterations")
        self.seed = _integer(self.seed, "seed", 0)
        if self.max_iterations < self.levels:
            raise ValueError("max_iterations must be at least levels")
        self.electric_cutoff = validate_electric_cutoff(self.electric_cutoff)
        data = self.build_hamiltonian()
        if data.dimension < 2:
            raise ValueError(
                "electric_cutoff retains only the vacuum; increase it to at "
                "least 3 to include an excited gauge-invariant state"
            )
        levels = min(self.levels, data.dimension)
        operator = SymmetricOperator(data.diagonal, data.rows, data.columns, data.entries)
        if self.model.nu == 0:
            # Single-vector Lanczos can omit exact multiplicities. A diagonal
            # electric Hamiltonian has its full spectrum available directly.
            values = np.sort(data.diagonal)[:levels]
            if retain_vectors:
                order = np.argsort(data.diagonal, kind="stable")[:levels]
                vectors = np.zeros((data.dimension, levels))
                vectors[order, np.arange(levels)] = 1.0
            info = {
                "method": "exact diagonal electric spectrum",
                "iterations": 0,
                "residual_norms": [0.0] * levels,
                "converged": True,
            }
        elif self.backend == "lanczos":
            values, vectors, info = smallest_ritz(
                operator, levels=levels, max_iterations=self.max_iterations,
                tolerance=self.tolerance, seed=self.seed,
            )
        else:
            from scipy.sparse.linalg import eigsh, ArpackNoConvergence

            matrix = data.to_scipy()
            try:
                if data.dimension <= 150 or levels >= data.dimension:
                    values, vectors = np.linalg.eigh(operator.dense())
                    values, vectors = values[:levels], vectors[:, :levels]
                    method = "dense NumPy eigh (small ARPACK input)"
                else:
                    start = np.random.default_rng(self.seed).standard_normal(data.dimension)
                    values, vectors = eigsh(
                        matrix, k=levels, which="SA", v0=start, tol=0.0,
                        maxiter=self.max_iterations,
                    )
                    order = np.argsort(values)
                    values, vectors = values[order], vectors[:, order]
                    method = "SciPy ARPACK eigsh"
            except ArpackNoConvergence as exc:
                raise SpectrumConvergenceError(
                    "ARPACK did not converge; increase max_iterations"
                ) from exc
            residuals = np.linalg.norm(matrix @ vectors - vectors * values, axis=0)
            info = {
                "method": method,
                # SciPy does not expose the number of iterations used.
                "iterations": -1,
                "residual_norms": residuals.tolist(),
                "converged": bool(np.max(residuals) <= self.tolerance),
            }
        result = SpectrumResult(
            e0=float(values[0]), e1=float(values[1]), gap=float(values[1] - values[0]),
            eigenvalues=tuple(float(v) for v in values),
            residual_norms=tuple(float(r) for r in info["residual_norms"]),
            converged=bool(info["converged"]), dimension=data.dimension,
            num_cubes=self.model.num_cubes, kappa=self.model.kappa, nu=self.model.nu,
            electric_cutoff=str(self.electric_cutoff), pairing=self.model.pairing,
            backend=self.backend, method=info["method"], iterations=info["iterations"],
            seed=self.seed, tolerance=self.tolerance, max_iterations=self.max_iterations,
            hermiticity_error=data.hermiticity_error,
        )
        if not result.converged or not all(np.isfinite(result.eigenvalues)):
            raise SpectrumConvergenceError(
                "The retained-matrix residual did not meet tolerance; "
                "increase max_iterations or inspect the requested tolerance.",
                result=result,
            )
        if retain_vectors:
            complete = levels == data.dimension and (
                self.model.nu == 0 or info["method"].startswith("dense NumPy")
            )
            return Eigensystem(result, data.states, vectors, complete)
        return result
