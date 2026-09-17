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


if __name__ == '__main__':
    self_test()
