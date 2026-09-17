# terry_cubes 0.1.0 — validation record

## Source and scope

Refactored from commit 71c3190363791c0f4464ee5e24606876e9131636 of
https://github.com/BoredPuzzleSolver/Joined-Cube-Hamiltonian-Engine-.

The original five Python scripts are preserved in tests/reference/.
The three spectra JSON files are preserved in data/. Git recognizes all
eight as 100% unchanged renames. SOURCE_PROVENANCE.json records original
checkout checksums and LF-normalized checksums for cross-platform comparison.
The original MIT license is unchanged.

This validates a software refactor and numerical reproduction of the
specified finite-basis Hamiltonians. It does not rerun the separate
single-cube exact-arithmetic certificate or establish an infinite-volume
or continuum Yang–Mills mass gap.

## Tests executed

| Check | Result |
| --- | --- |
| Routine source-checkout suite, python -m pytest -q | 68 passed; 40 large cases skipped by default |
| All published spectra, python -m pytest --run-slow -m slow tests/test_spectrum.py -q | 40 passed; 30 routine spectrum tests deselected; 614.22 seconds |
| Routine suite against an installed wheel outside the checkout | 68 passed; 40 large cases skipped |
| Editable installation and actual console executable | Passed |
| Source distribution and wheel built from that source distribution | Passed |
| All three datasets loaded from an isolated wheel installation | Passed; original JSON checksums matched |
| Original five scripts in source archive and absence of test code from wheel | Passed |
| License retained in wheel | Passed |
| Git whitespace check | Passed |

The routine and full-baseline test sets contain **108 distinct passing test
cases**. The wheel run additionally checks installation independently of
imports from the working source directory.

### Mathematical preservation checks

The 26 dedicated preservation tests compare the package with the original
implementations. They cover:

- Graph topology, physical edge and face ordering, orientation phases, and
  metadata for one, two, and three cubes with all three vertex pairings.
- Exact ordered states and signed rational-squared Wilson amplitudes.
- All compatible auxiliary intertwiner channels, including auxiliary spin
  one with physical links restricted to spin one-half.
- Zero electric energy on auxiliary channels.
- Original cube bases, sparse entry order, and contributions from multiple
  faces to the same matrix entry.
- Independent SymPy Wigner 6j checks and 62 exact oriented-vertex contractions.

### Numerical baseline coverage

Every row in the three original spectra files was recomputed, using absolute
eigenvalue/gap tolerance 1e-9 and zero relative tolerance:

- N=1: cutoffs 6, 8, 10, 12, 14, 16, 18, at both coupling ratios.
- N=2: cutoffs 6, 8, 10, 12, 14, 16, 18, at both coupling ratios.
- N=3: cutoffs 6, 8, 10, 12, 14, 16, at both coupling ratios.
- Coupling ratios: nu/kappa = 0.25 and 1.0.
- Largest retained basis: 771,425 states.

The N=1 cutoff-12 estimate at kappa=nu=1 lies within the previously reported
single-cube interval [2.941733, 2.958693]. This is an interval-consistency
regression, not a replacement for that certificate's proof machinery.

## Exact requested CLI example

Command:

    terry-cubes run --cubes 2 --kappa 1.0 --nu 0.25 --cutoff 16 --output results.json

Observed output:

| Quantity | Value |
| --- | --- |
| Basis dimension | 76,326 |
| E0 | 2.6927571061973903 |
| E1 | 5.69136685244009 |
| Gap | 2.9986097462427 |
| Maximum retained-matrix residual | 2.13e-12 |
| Original Lanczos iterations | 132 |

Absolute differences from the matching original JSON row were approximately
2.66e-15 for E0, 8.88e-16 for E1, and 1.78e-15 for the gap.
The saved run is validation/two_cubes_cutoff16.json.

## API and operational checks

Tests also cover parameter validation, insufficient cutoff, direct
pure-electric degeneracies, common rescaling of the couplings, ARPACK
agreement, an independently known analytic sparse spectrum, strict JSON
serialization, and failure without overwriting a previous output file.

Changing a solver's model or cutoff rebuilds its cached Hamiltonian.
Nonconverged solves raise SpectrumConvergenceError. The original
single-vector Lanczos algorithm remains unchanged and can omit repeated
copies of degenerate eigenvalues; its converged Ritz values alone do not
certify spectral ordering or multiplicities.

## Validation environment

- Windows 11, build 26200
- Python 3.12.14
- NumPy 2.3.5
- SciPy 1.18.1
- SymPy 1.14.0
- pytest 9.1.1
- setuptools 84.0.0

Runs used OPENBLAS_NUM_THREADS=1 and OMP_NUM_THREADS=1.
Timings and residual roundoff will differ across machines.

CI is configured for Linux with Python 3.10 and 3.12, and Windows with
Python 3.12. Those remote CI jobs have not been run during this local refactor.
