[README.md](https://github.com/user-attachments/files/32347692/README.md)
# terry_cubes 0.1.1

A Python research library for finite SU(2) Hamiltonian lattice gauge
calculations on open chains of joined cubes. It provides gauge-invariant
state enumeration, exact local recoupling, sparse Hamiltonians, numerical
spectra, cutoff comparisons, and finite-basis plaquette correlations.

Version 0.1.1 incorporates supporting research omitted from 0.1.0. The
original geometry, Hamiltonian, recoupling formulas, state indexing, and
NumPy Lanczos routine are preserved. The original 0.1.0 package remains a
separate release. See [CHANGELOG.md](CHANGELOG.md) and
[RELEASE_PROVENANCE.json](RELEASE_PROVENANCE.json).

## What the results establish

| Calculation | Supported interpretation |
| --- | --- |
| Single open cube, `nu/kappa = 1` | Numerical gap near `2.951043 * kappa`; separately verified interval `[2.941733, 2.958693] * kappa` for the full spin space of this fixed graph. |
| Two and three joined cubes | Finite-basis numerical spectra with recorded cutoffs and residuals. |
| Two squares sharing one link | A separate seven-link benchmark with its own spin cutoffs, correlators, and exact certificates. It is not two joined cubes. |
| Each fixed finite cube chain | The archived operator argument supports a qualitative positive gap at fixed positive kappa and finite nonnegative nu, under its stated assumptions; it supplies no bound uniform in chain length. |

Infinite-volume and continuum claims, physical glueball masses, and gameplay
applications remain unestablished. Earlier speculative ideas in the research
chronology are hypotheses. Increasing an N-by-1-by-1 chain extends only one
spatial direction. Neither that sequence nor spin-cutoff convergence defines
a three-dimensional continuum limit.

Numerical outputs use `certified: false`. Reading a saved certificate does
not reverify it. Only a successful call to the separate fixed-case certificate
verifier returns `certified: true`, with its model, parameters, and scope.

## Install this release

Requires Python 3.10 or later, NumPy, SciPy, and SymPy. Open a terminal in
the extracted **terry_cubes-0.1.1** folder, which contains `pyproject.toml`:

```bash
python -m venv .venv
```

Activate on Windows PowerShell:

```powershell
.\.venv\Scripts\Activate.ps1
```

Or on macOS/Linux:

```bash
source .venv/bin/activate
```

Then install:

```bash
python -m pip install .
terry-cubes --version
```

Alternatively, install `dist/terry_cubes-0.1.1-py3-none-any.whl` with pip.
Use `python -m pip install -e ".[dev]"` for development. These instructions
work from this versioned package folder; no PyPI publication or repository
root layout is assumed.

## Joined-cube spectra

```python
from terry_cubes import JoinedCubesModel, SpectrumSolver

model = JoinedCubesModel(num_cubes=1, kappa=1.0, nu=1.0)
result = SpectrumSolver(model, electric_cutoff=12, backend="arpack",
                        max_iterations=1000).compute_mass_gap()
print(result.e0, result.e1, result.gap)
print(result.dimension, result.residual_norms, result.converged)
```

```bash
terry-cubes run --cubes 1 --cutoff 12 --backend arpack --max-iterations 1000
terry-cubes run --cubes 2 --kappa 1 --nu 0.25 --cutoff 16 --output results.json
```

The original `JoinedCubesModel`, `SpectrumSolver`, result properties, and
`run` command remain available. `to_dict()` serializes a numerical result.
`compute_eigensystem()` additionally returns `.spectrum`, `.states`, and
`.eigenvectors` (one column per returned energy in the original state order).

| Argument | Default | Meaning |
| --- | --- | --- |
| `num_cubes` | `1` | Open N-by-1-by-1 chain. |
| `kappa`, `nu` | `1.0`, `1.0` | Positive electric and nonnegative magnetic coefficients. |
| `pairing` | `"longitudinal"` | Also `"crossed"` and `"mixed"`. All compatible intertwiner channels remain present. |
| `electric_cutoff` | `16` | Bound on physical-link Casimir sum, independent of kappa. |
| `backend` | `"lanczos"` | Preserved NumPy implementation; `"arpack"` uses the same matrix with SciPy. |
| `tolerance` | `1e-10` | Absolute residual tolerance in the retained matrix. |
| `max_iterations` | `240` | Solver iteration budget; interpretation depends on the backend. |
| `levels` | `3` | Requested low energies, at least two. |
| `seed` | `1097` | Initial iterative vector seed. |

The first three rows configure the model; the remaining rows configure the
solver. Failed residual convergence raises `SpectrumConvergenceError`, with
a partial `.result` where available. A small residual does not certify that
every lower state was found or control omitted states. Single-vector Lanczos
can omit repeated copies of degenerate eigenvalues. Use a dense calculation
where feasible when multiplicities matter. The purely electric case retains
its exact diagonal degeneracies.

## Cutoff comparisons

```python
from terry_cubes import JoinedCubesModel, cutoff_study

study = cutoff_study(JoinedCubesModel(1), [6, 8, 10, 12],
                     embedding_cutoff=18)
for row in study["results"]:
    print(row["electric_cutoff"], row["gap"],
          row.get("residual_in_largest_test_basis"))
```

```bash
terry-cubes study --cubes 1 --cutoffs 6 8 10 12 --embedding-cutoff 18 --output study.json
```

This compares E0 and E1 in nested complete electric cutoffs, checks numerical
consistency with their nonincreasing variational ordering, and repeats the
largest requested solve with a second seed. It reports both direct residuals
and, for smaller bases, residuals after embedding in a larger test basis.
The gap difference itself has no monotonicity guarantee. Larger-basis
residuals and cutoff changes are diagnostics, not full-space error bounds.
The study defaults to ARPACK, 1,000 iterations, and seeds 1097 and 8917.

## Fresh certificate verification

```bash
terry-cubes certify single-cube --output cube_certificate.json
terry-cubes certify two-plaquette --output square_certificates.json
```

```python
from terry_cubes import verify_certificate

certificate = verify_certificate("single_cube", cutoff=12)
print(certificate["gap_interval_exact"])
# ['2941733/1000000', '2958693/1000000']
```

The single-cube verifier supports the original complete electric cutoffs
8, 10, and 12, only at `kappa = nu = 1`. Cutoff 12 reconstructs the 398-state
retained block and the 2,042-state comparison construction, including 1,644
directly coupled omitted states. It checks fixed rational energy endpoints
with exact interval inertia, radical-rounding bounds, and the remaining-state
electric threshold. Scaling both coefficients by the same positive kappa
scales the resulting interval by kappa; a different ratio needs its own work.

Before reporting success, the verifier checks the archived source hashes,
matches the retained basis to the actual package, checks all retained-to-retained
and retained-to-omitted transition amplitudes, and checks exact electric energies.
The original standalone cube and joined-graph conventions are related by a
single diagonal state-sign transformation, verified consistently across all
faces. This does not change state indexing or any eigenvalue.

The two-plaquette verifier checks all three original benchmark constructions
at ratios 1 and 4, including their exact agreement with the package's preserved
two-square recoupling implementation. Its strongest saved intervals are
`[3.09195898, 3.09195907] * kappa` at ratio 1 and
`[4.20942742, 4.20942757] * kappa` at ratio 4. These belong to the two-square graph.

The preserved programs execute in a temporary directory, with assertions
enabled. The installed sources and saved records remain unchanged. Full source
bundles and their proof notes are included as wheel resources and as readable
source folders under `certificates/`. Verification may take considerably longer
than a small numerical solve. `--timeout SECONDS` fails without a certificate
if the time limit is exceeded. There is no general multi-cube certificate API.

## Two-plaquette benchmark and observables

```python
from terry_cubes import TwoPlaquetteModel

result = TwoPlaquetteModel(kappa=1, nu=1).compute_spectrum(max_twice=12)
channels = result.correlators()
print(result.gap)
print(channels["symmetric"]["lowest_detected_gap"])
print(channels["antisymmetric"]["lowest_detected_gap"])
```

```bash
terry-cubes two-plaquette --max-twice 12 --correlators --output two_squares.json
```

Here `max_twice=12` means physical spins at most 6. This is a per-link spin
cutoff, not the joined-cube electric cutoff. States retain the original
`(2*jL, 2*jR, 2*jM)` labels and lexicographic order. Full dense diagonalization
is practical for the reference 616-state case; larger cutoffs grow its cost.

The symmetric channel detects the approximately 3.091959 gap at ratio 1.
The antisymmetric channel first detects approximately 3.113958 because of
symmetry selection. A correlator's first detected excitation need not be the
lowest Hamiltonian excitation.

For joined cubes, `plaquette_operator(model, states, face_index)` builds
`P (Tr_fund U_face / 2) P` using the preserved face and state order.
`connected_correlator(energies, eigenvectors, operator, times)` computes

```text
C_connected(t) = sum_{n>0} |<n|O|0>|^2 exp[-(En-E0)t]
m_effective(t) = log[C_connected(t)/C_connected(t+dt)] / dt
```

Time is in inverse Hamiltonian energy units with hbar=1. Partial eigenvectors
give a partial spectral sum. Output records unrepresented and threshold-discarded
weight. Even a complete retained spectrum leaves physical states above the
cutoff unaccounted for. Examples are in [examples/](examples/).

## Hamiltonian, units, and indexing

The original formula and additive constant are preserved:

```text
H = kappa * sum_physical_links j(j + 1)
    + nu * sum_physical_faces (1 - Tr_fund(U_face) / 2)

sum_physical_links j(j + 1) <= electric_cutoff
```

All representation labels are doubled spins. Physical links precede virtual
coupling labels. Virtual labels resolve four-valent invariant tensor spaces,
carry zero electric energy, and retain every compatible pair-coupling channel.
Physical edge and face order, state enumeration, and matrix indexing are
unchanged. Exact rational-squared Wigner 6j and oriented 3j/Clebsch-Gordan
contractions keep the original phase rules and conversion to floating point.
The oriented-edge convention remains `K_j(U) = D_j(U) epsilon_j`, with ascending
incident-edge order and the original reversed-face traversal signs.

At fixed positive kappa, compare `E0/kappa`, `E1/kappa`, and `gap/kappa` with
reference data at the same `nu/kappa`, graph, pairing, and cutoff. There is no
implemented calibration from these units to physical glueball masses.

## Research records and source separation

```python
from terry_cubes.data import load_spectra, list_research_datasets, load_research_dataset

joined = load_spectra(2)  # Original N=1,2,3 records remain available unchanged.
catalog = list_research_datasets()
extended = load_research_dataset("single_cube_energy")
assert extended["verified_in_this_call"] is False
```

```bash
terry-cubes data
terry-cubes data --name single_cube_energy_certificates
```

Seven additional JSON datasets include the single-cube spin and electric
studies through electric cutoff 30, the strong and earlier weaker certificates,
cutoff bounds, and two-square spectra, correlators, and certificates. Every
catalog entry identifies its graph, cutoff type, and evidence status.

| Location | Role |
| --- | --- |
| `terry_cubes/` | Actual runtime library and public interfaces. |
| `data/` | Unchanged records and preserved certificate ZIPs, installed as `terry_cubes.data`. |
| `tests/reference/` | Five original mathematical scripts used by preservation tests. |
| `certificates/` | Separate, preserved fixed-graph verification programs and proof notes. |
| `validation/reference_joined/` | Independent tensor, cycle-space, orientation, and volume-study research archive. |
| `docs/theory/` | Preserved derivations and limits of the finite-graph conclusions. |
| `docs/research/` | Chronology and an index mapping research evidence into this release. |
| `validation/` | This release's fresh results; `history/` explicitly retains older reports. |

The chronology is a guide to the research sequence; its narrative claims are
not automatically calculator features. The research scripts and certificate
programs retain separate roles from the public numerical solver. Read
[Research evidence index](docs/research/RESEARCH_INDEX.md) for the mapping.

## Tests and building

```bash
python -m pip install -e ".[dev]"
python -m pytest
python -m pytest --run-research -m research
python -m pytest --run-slow -m slow tests/test_spectrum.py
python -m build
```

The routine suite covers the original API and preservation checks plus the
new interfaces and a small fresh certificate. The opt-in research suite checks
3,368 exact oriented contractions, independent tensor matrices and cycle-space
counts, all three pairing choices, extended single-cube results at cutoffs 24
and 30, and the stronger fresh certificates. The separate 40-case slow suite
recomputes every original joined-cube numerical record.

The slow suite can require several gigabytes of memory. Its largest supplied
basis has 771,425 states. The original Lanczos Krylov allocation alone is
approximately `8 * dimension * (min(dimension, max_iterations) + 1)` bytes,
in addition to sparse matrices, labels, and work arrays.

The included CI workflow is configured for Linux and Windows, source and wheel
tests, and an optional manual research run. It becomes active when this folder
is the repository root; a nested package needs the workflow at the repository's
top-level `.github/workflows/` and appropriate working-directory paths. No remote
CI run or publication is implied by this local release. See
[VALIDATION.md](VALIDATION.md) for checks actually executed.

## Attribution and license

Original research implementation: Terrance Mullee / BoredPuzzleSolver,
September 2026. The existing MIT license and copyright notice are preserved
in [LICENSE](LICENSE). [SOURCE_PROVENANCE.json](SOURCE_PROVENANCE.json) records
the original refactor's source commit and hashes; the release provenance
records this update's preserved files and added source bundles.
