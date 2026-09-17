# Terry Cubes (`terry_cubes`)

**A Python research library for gauge-invariant SU(2) Hamiltonian calculations on finite open chains of joined cubes.**

`terry_cubes` constructs spin-network bases, evaluates local SU(2) recoupling coefficients, assembles sparse Hamiltonians, and computes numerical ground-state energies, low-lying excitation energies, and spectral-gap estimates. It provides a Python API, a command-line interface, automated tests, and reference datasets for one, two, and three cubes.

The package is intended for finite-lattice studies, numerical benchmarking, cutoff comparisons, and educational exploration of Hamiltonian lattice gauge theory. It implements a finite-graph SU(2) model based on the Kogut–Susskind Hamiltonian framework.

**Status:** research prototype, version 0.1.0. The published spectra are finite-basis numerical results. An infinite-volume gap, a three-dimensional bulk bound, and a continuum Yang–Mills mass gap have not been established by this package.

## Implemented features

- Open `N × 1 × 1` chains of cubes, with shared physical links and faces.
- Gauge-invariant state enumeration using SU(2) singlet constraints.
- Complete sets of compatible intertwiner channels at four-valent joining vertices.
- Exact signed local recoupling data, with floating-point conversion for numerical matrix construction.
- Sparse electric and magnetic Hamiltonian terms.
- A fully reorthogonalized Lanczos solver and an optional SciPy ARPACK backend.
- Eigenpair residual checks and explicit convergence reporting.
- Python and command-line interfaces, with JSON output.
- Reference spectra and regression tests for `N = 1, 2, 3`.

The current implementation supports **SU(2)** and the specified chain geometry. General three-dimensional lattices, six-valent bulk vertices, and SU(3) are not implemented.

## Model and truncation

The Hamiltonian is

$$
H = \kappa\sum_{e\in E_{\mathrm{physical}}}j_e(j_e+1)
  + \nu\sum_{p\in F}\left(1-\frac{1}{2}\operatorname{Tr}_{\mathrm{fund}}U_p\right),
\qquad \kappa>0,\quad \nu\geq0.
$$

Here, `kappa` is the electric coefficient, `nu` is the magnetic plaquette coefficient, and `U_p` is the ordered product of link matrices around a face.

The basis includes every compatible physical state satisfying the **total electric-energy cutoff**

$$
\sum_{e\in E_{\mathrm{physical}}}j_e(j_e+1)\leq C.
$$

This is the meaning of `electric_cutoff`. It is independent of `kappa` and is different from imposing the same fixed maximum spin on every link. All compatible auxiliary coupling channels are retained; these labels resolve vertex intertwiners and carry no additional electric energy.

For the largest supplied calculation, `N = 3` and `C = 16` produce **771,425 basis states**. The largest physical spin present in that basis is `j = 3/2`; auxiliary coupling spins reach `2`. This basis size counts quantum states, not spatial voxels.

The reported gap is the finite-basis estimate

$$
\Delta = E_1-E_0.
$$

At fixed cutoff and `nu / kappa`, a common rescaling of `kappa` and `nu` rescales all energies. Comparison with the supplied `kappa = 1` datasets should therefore use `E0 / kappa`, `E1 / kappa`, and `gap / kappa`. Conversion to physical units requires a separate scale-setting prescription.

## Installation

Requires **Python 3.10 or later**, NumPy, SciPy, and SymPy.

The installable package currently lives inside the repository's `terry_cubes-0.1.0` directory:

```bash
git clone https://github.com/BoredPuzzleSolver/Joined-Cube-Hamiltonian-Engine-.git
cd Joined-Cube-Hamiltonian-Engine-/terry_cubes-0.1.0
python -m venv .venv
```

Activate the virtual environment using the command for your platform.

**Windows PowerShell:**

```powershell
.\.venv\Scripts\Activate.ps1
```

**macOS / Linux:**

```bash
source .venv/bin/activate
```

Then install:

```bash
python -m pip install -e .
```

This installs the Python library and the `terry-cubes` command. These instructions install from the repository; they do not assume a PyPI release.

## Python quick start

A small example:

```python
from terry_cubes import JoinedCubesModel, SpectrumSolver

model = JoinedCubesModel(num_cubes=1, kappa=1.0, nu=1.0)
solver = SpectrumSolver(model, electric_cutoff=12)
results = solver.compute_mass_gap()

print("Ground-state energy E0:", results.e0)
print("First-excitation energy E1:", results.e1)
print("Finite-basis gap:", results.gap)
print("Basis dimension:", results.dimension)
print("Residual norms:", results.residual_norms)
print("Residual tolerance satisfied:", results.converged)
```

Change `num_cubes`, `electric_cutoff`, and the coupling ratio to study other finite systems. Increasing the cutoff can change the result even when the eigensolver residual is very small.

To use SciPy's sparse eigensolver:

```python
solver = SpectrumSolver(
    model,
    electric_cutoff=12,
    backend="arpack",
    max_iterations=1000,
)
results = solver.compute_mass_gap()
```

## Command line

```bash
terry-cubes run --cubes 2 --kappa 1.0 --nu 0.25 --cutoff 16 --output results.json
```

Additional options include `--backend`, `--pairing`, `--tolerance`, `--max-iterations`, `--levels`, and `--seed`.

```bash
terry-cubes run --help
```

## Reference calculations

The bundled data contain 40 finite-basis reference cases across `N = 1, 2, 3`, several electric cutoffs, and the ratios `nu / kappa = 0.25` and `1.0`.

Selected equal-coupling results (`kappa = nu = 1`):

| Cubes | Electric cutoff | Basis dimension | Numerical gap |
|---:|---:|---:|---:|
| 1 | 18 | 2,920 | 2.95111924 |
| 2 | 18 | 201,248 | 2.91087079 |
| 3 | 16 | 771,425 | 2.90154382 |

These are eigenvalue differences of the specified finite matrices. Their displayed digits do not provide an error bound for the untruncated theory. The different cutoffs and small number of volumes do not establish exponential saturation or a positive infinite-volume limit.

The data are in [`terry_cubes-0.1.0/data/`](terry_cubes-0.1.0/data/). See the [validation record](terry_cubes-0.1.0/VALIDATION.md) for the documented numerical reproduction checks and their scope.

### Numerical convergence and certification

A small eigenpair residual measures how accurately the computed vector satisfies the retained matrix equation. It does not bound the effect of omitted states, establish spectral multiplicities, or prove a bulk or continuum gap. Solver output is explicitly marked `certified: false`.

The default single-vector Lanczos backend can omit repeated copies of degenerate eigenvalues. Researchers studying degeneracies or spectral counting should use additional solver and symmetry checks. Failed residual convergence raises `SpectrumConvergenceError`.

A separate research calculation gives the fixed single-cube interval `[2.941733, 2.958693] * kappa` at `nu / kappa = 1`, with control of omitted spin states. The package tests check numerical consistency with that reported interval. **The separate certificate verifier and derivation are not included in this release**, and interval membership is not a certification procedure. No corresponding untruncated multi-cube certificate is supplied.

## Intended uses

### Current applications

- Calculate low-energy spectra for the implemented finite SU(2) cube chains.
- Investigate sensitivity to the electric cutoff and magnetic/electric ratio.
- Compare equivalent intertwiner coupling trees and numerical eigensolvers.
- Supply finite-system reference matrices and spectra for testing other numerical methods.
- Demonstrate spin networks, SU(2) recoupling, and local Gauss-law constraints.

### Exploratory extensions

The generated Hamiltonians and spectra may serve as inputs to future work on observables, finite-system time evolution, quantum-algorithm benchmarks, or visual demonstrations. Those applications require additional implementation and validation.

For a game or visualization, a custom adapter could map a computed value to an interaction parameter. That mapping would be part of the demonstration's chosen model. The current package supplies no game-engine integration or validated model of fluids, elasticity, friction, fracture, or material strength. A computed spectral gap alone does not determine those properties.

## Limits and computational cost

- Increasing `num_cubes` extends only one spatial direction. The transverse dimensions remain one cube wide.
- The package does not perform physical scale setting, continuum extrapolation, or a two-loop scaling test.
- Cutoff convergence, volume dependence, and continuum scaling are separate questions.
- The implementation constructs a global finite Hamiltonian. It does not implement independent per-cube solves, a linear-time parallel physics algorithm, GPU acceleration, or exact voxel level-of-detail decimation.
- Large bases can require several gigabytes of memory. The default Lanczos implementation stores its Krylov vectors explicitly, in addition to state labels and sparse matrix arrays.

## Tests

From the `terry_cubes-0.1.0` directory:

```bash
python -m pip install -e ".[dev]"
python -m pytest
```

To recompute all 40 recorded baseline cases:

```bash
python -m pytest --run-slow -m slow tests/test_spectrum.py
```

The full reproduction takes more time and memory than the routine suite. Passing a reference comparison establishes numerical agreement for the specified finite model; it does not supply a new rigorous enclosure of an untruncated spectrum.

## Repository layout

```text
Joined-Cube-Hamiltonian-Engine-/
├── README.md
├── LICENSE
├── terry_cubes-0.1.0/       # Installable package and its documentation
│   ├── terry_cubes/        # Current implementation and public API
│   ├── data/              # Reference spectra
│   ├── tests/             # Tests and preserved reference implementations
│   ├── pyproject.toml
│   ├── README.md
│   └── VALIDATION.md
└── research_scripts/      # Archived standalone research scripts
```

Use the installable `terry_cubes` package for new calculations. The archived scripts document the earlier implementation.

## Citation

When reporting results, include the package version or commit, geometry, cutoff, coupling ratio, solver settings, and any additional error analysis.

```bibtex
@software{mullee2026terrycubes,
  author = {Mullee, Terrance},
  title = {Terry Cubes: Finite SU(2) Hamiltonian Spectra on Joined Cubes},
  version = {0.1.0},
  year = {2026},
  url = {https://github.com/BoredPuzzleSolver/Joined-Cube-Hamiltonian-Engine-}
}
```

## License

MIT. See [LICENSE](LICENSE).

Original research implementation: Terrance Mullee / BoredPuzzleSolver.





