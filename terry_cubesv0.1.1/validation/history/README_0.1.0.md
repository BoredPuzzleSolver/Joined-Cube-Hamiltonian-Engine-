# terry_cubes

A Python library for finite SU(2) Hamiltonian lattice gauge calculations on
open chains of joined cubes. It provides gauge-invariant state enumeration,
exact local recoupling coefficients, sparse Hamiltonian construction, and
numerical low-energy spectra through a Python API and the `terry-cubes` command.

**Scope:** the included JSON files contain finite-basis numerical spectra for
one, two, and three cubes. They are reference results for reproducibility.
The reported single-cube interval
`[2.941733, 2.958693] * kappa` is checked for consistency with a numerical
single-cube result at `nu / kappa = 1`; the solver does not produce that
untruncated certificate. This package does not establish an infinite-volume
or continuum Yang–Mills mass gap.

## Installation

Requires Python 3.10 or later, NumPy, SciPy, and SymPy.

Clone this repository and install the package from its root:

```bash
git clone https://github.com/BoredPuzzleSolver/Joined-Cube-Hamiltonian-Engine-.git
cd Joined-Cube-Hamiltonian-Engine-
python -m venv .venv
```

Activate the environment on macOS or Linux:

```bash
source .venv/bin/activate
```

Or activate it in Windows PowerShell:

```powershell
.\.venv\Scripts\Activate.ps1
```

Then install:

```bash
python -m pip install -e .
```

This installs the library and the `terry-cubes` terminal command. The dependency
list is also available in `requirements.txt`; installation with pip resolves it
automatically. These are source-install instructions; no PyPI release is assumed.

## Python quick start

The requested two-cube calculation is:

```python
from terry_cubes import JoinedCubesModel, SpectrumSolver

model = JoinedCubesModel(num_cubes=2, kappa=1.0, nu=1.0)
solver = SpectrumSolver(model, electric_cutoff=16)
results = solver.compute_mass_gap()

print("Ground state E0:", results.e0)
print("First excited state E1:", results.e1)
print("Mass gap Delta:", results.gap)
print("Basis dimension:", results.dimension)
print("Residual norms:", results.residual_norms)
print("Residual tolerance satisfied:", results.converged)
```

For a smaller first run, use `num_cubes=1` and `electric_cutoff=6`.
This changes the finite basis and is not a converged-cutoff calculation.

The result also provides `eigenvalues` and `to_dict()` for serialization:

```python
import json
from pathlib import Path

Path("results.json").write_text(
    json.dumps(results.to_dict(), indent=2) + "\n",
    encoding="utf-8",
)
```

### Model and solver options

| Argument | Default | Meaning |
| --- | --- | --- |
| `num_cubes` | `1` | Number of cubes in an open N-by-1-by-1 chain. |
| `kappa` | `1.0` | Electric energy coefficient. |
| `nu` | `1.0` | Magnetic plaquette coefficient. |
| `pairing` | `"longitudinal"` | Four-valent intertwiner resolution; also `"crossed"` and `"mixed"`. |
| `electric_cutoff` | `16` | Maximum sum of physical-link Casimirs, before multiplying by kappa. |
| `backend` | `"lanczos"` | Original NumPy Lanczos solver, or optional SciPy `"arpack"`. |
| `tolerance` | `1e-10` | Direct eigenvector residual tolerance. |
| `max_iterations` | `240` | Iteration budget; its interpretation follows the selected solver. |
| `levels` | `3` | Number of low-energy levels requested, including the ground state. |
| `seed` | `1097` | Seed for the initial iterative-solver vector. |

`JoinedCubesModel` accepts the first four arguments; `SpectrumSolver` accepts
the model and the remaining arguments. A mass-gap calculation requires at
least two states and two requested levels.

The default backend preserves the original fully reorthogonalized NumPy
Lanczos routine, including its dense diagonalization path for small matrices.
The ARPACK backend uses the same Hamiltonian with SciPy's sparse eigensolver:

```python
solver = SpectrumSolver(
    model,
    electric_cutoff=16,
    backend="arpack",
    tolerance=1e-10,
    max_iterations=1000,
)
results = solver.compute_mass_gap()
```

The API raises `SpectrumConvergenceError` if the requested residual tolerance
is not met. Where available, the exception exposes a partial result as `.result`.
Check `results.converged` and the reported residuals before using a result.
A small residual measures how well an eigenpair satisfies the finite matrix
equation; it does not by itself certify that every lower eigenvalue has been
found or control the omitted states above the cutoff. The original
single-vector Lanczos routine can omit repeated copies of degenerate
eigenvalues, so its returned `eigenvalues` are observed Ritz values and do
not establish multiplicities. The purely electric case (`nu=0`) is handled
directly from the diagonal, preserving its exact degeneracies. The public
API raises `SpectrumConvergenceError` if residual convergence fails.

## Command line

```bash
terry-cubes run --cubes 2 --kappa 1.0 --nu 0.25 --cutoff 16 --output results.json
```

Use `--backend arpack`, `--tolerance`, `--max-iterations`, `--levels`,
`--seed`, and `--pairing` to configure a run. View the complete options with:

```bash
terry-cubes run --help
```

## Hamiltonian, units, and basis conventions

The original Hamiltonian and its additive plaquette constant are preserved:

```text
H = kappa * sum_physical_links j(j + 1)
    + nu * sum_physical_faces (1 - Tr_fund(U_face) / 2)
```

The cutoff is imposed on

```text
sum_physical_links j(j + 1) <= electric_cutoff
```

It is independent of `kappa`. At `kappa = 1`, returned energies agree with
the units used in the bundled records. To compare a run with another positive
kappa to those records, match `nu / kappa` and the cutoff, then compare
`e0 / kappa`, `e1 / kappa`, and `gap / kappa`.

The refactor retains the original conventions:

- All representation labels are doubled spins, so a label of 1 denotes
  physical spin 1/2.
- Physical links precede virtual coupling labels in each state. Virtual
  labels resolve four-valent invariant tensor spaces, carry zero electric
  energy, and retain every compatible pair-coupling channel.
- Physical edge and face order, state enumeration order, and matrix indexing
  remain as defined by the original scripts.
- Wigner 6j and oriented 3j/Clebsch–Gordan contractions use the original exact
  rational arithmetic and phase rules. Conversion to floating-point matrix
  entries occurs at the same stage as before.
- The oriented-edge convention is `K_j(U) = D_j(U) epsilon_j`, with a fixed
  ascending incident-edge order at each resolved vertex and the original
  signs for reversed face traversals.

Increasing N grows only one spatial direction. Cutoff convergence and chain
length dependence are separate questions; neither is a three-dimensional
bulk or continuum limit.

## Reference data

The three original JSON files live in the repository's `data/` directory.
They are installed as resources in `terry_cubes.data`, so loading works from
both editable installations and wheels without relying on the working directory.

```python
from terry_cubes.data import load_spectra

reference = load_spectra(2)
print(reference["status"])
for run in reference["spectra"]:
    print(run["energy_cap"], run["nu_over_kappa"], run["E0"], run["E1"])
```

Each record preserves the original model description, cutoff, eigenvalues,
and solver diagnostics. A comparison must match the number of cubes,
electric cutoff, intertwiner pairing, and magnetic/electric ratio.
Values labeled `residual_in_largest_test_basis` in the original records
are distinct from residuals in the diagonalized finite matrix.

## Tests and reproducibility

Install the development dependencies and run the routine suite:

```bash
python -m pip install -e ".[dev]"
python -m pytest
```

The spectrum tests load the original JSON files. Routine baseline comparisons
cover N=1 at cutoff 12, N=2 at cutoff 8, and N=3 at cutoff 6, each at
`nu / kappa = 0.25` and `1.0`. The single-cube consistency test checks the
reported interval without treating interval membership as a certificate.
Additional tests check recoupling, geometry, preserved matrix construction,
solver behavior, and the API/CLI.

To recompute all 40 recorded baseline cases:

```bash
python -m pytest --run-slow -m slow tests/test_spectrum.py
```

These full reproductions can require substantial time and several gigabytes
of memory. The largest supplied basis has 771,425 states. The original
Lanczos implementation stores its Krylov vectors explicitly; its basis
allocation alone is approximately
`8 * dimension * (min(dimension, max_iterations) + 1)` bytes. Sparse
Hamiltonian arrays, state labels, and work arrays require additional memory.

Reference eigenvalues are compared with numerical tolerances, not bitwise
identity: floating-point results may vary across numerical-library versions
and platforms. Reproducing a reference spectrum is a regression check for
the specified finite model, not a new rigorous spectral enclosure.

The untouched original five scripts are retained under `tests/reference/`
for comparison. [SOURCE_PROVENANCE.json](SOURCE_PROVENANCE.json) records the
source commit, file hashes, and the origin of the sparse-builder routines.
The original scripts are included in the source distribution and excluded
from the installed library. Use `terry_cubes` for new calculations.

## Package layout

```text
terry_cubes/
    __init__.py       # Public model and solver API
    __main__.py       # python -m terry_cubes entry point
    geometry.py       # Single-cube and joined-cube graph topology
    recoupling.py     # Exact SU(2) recoupling and oriented intertwiners
    hamiltonian.py    # State basis, electric terms, plaquette matrix terms
    solver.py         # Sparse operators, eigensolvers, result objects
    cli.py            # argparse command-line entry point
data/
    __init__.py       # Installed as terry_cubes.data
    SU2_Joined_Cubes_N1_Spectra.json
    SU2_Joined_Cubes_N2_Spectra.json
    SU2_Joined_Cubes_N3_Spectra.json
tests/
    test_spectrum.py
    reference/        # Original mathematical implementations
pyproject.toml        # Authoritative package metadata and CLI entry point
setup.py              # Minimal setuptools compatibility entry point
requirements.txt
README.md
LICENSE
```

Build a source archive and an installable wheel with:

```bash
python -m build
```

The files are written to `dist/`. CI runs the routine suite on Linux and
Windows, builds both distributions, and checks the installed wheel's data
resources and terminal entry point. The full large-basis runs are opt-in.

## Attribution and license

Original research implementation: Terrance Mullee / BoredPuzzleSolver,
September 2026. The existing MIT license and copyright notice are preserved
in [LICENSE](LICENSE).
