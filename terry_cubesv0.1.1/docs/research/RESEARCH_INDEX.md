# Research evidence incorporated in 0.1.1

The supplied chronology was read as historical context, then compared with the
actual `terry_cubes-0.1.0` package and the supporting files in `outputs/`.
The mathematical core already contained the corrected recoupling and joined-cube
construction. The main omissions were accessible evidence, independent validation,
certificate execution, and reusable research interfaces.

| Evidence in the supplied outputs | Place and use in 0.1.1 |
| --- | --- |
| `SU2_Cube_Reproduction` | Preserved as `certificates/single_cube/` and an installed ZIP resource; fresh cube verifier uses its original comparison matrices and inertia routines. |
| `Shared_SU2_Cutoff_Reproduction` | Preserved as `certificates/two_plaquette/` and an installed ZIP resource; separate two-square certificate interface. |
| `SU2_Joined_Cubes_Reproduction` | Preserved as `validation/reference_joined/`; independent tensor and cycle-space checks, orientation evidence, numerical study, and volume-bound audit. |
| `SU2_Cube_Energy_Study_Results.json` | Dataset `single_cube_energy`; regression at cutoffs 24 and 30 for ratios 1 and 4. |
| `SU2_Cube_Study_Results.json` | Dataset `single_cube_spin`; original physical-spin cutoff sweep, distinct from electric cutoffs. |
| `SU2_Cube_Energy_Certificates.json` | Dataset `single_cube_energy_certificates`; complete electric-cutoff constructions at 8, 10, 12, including the stronger 398-state certificate. |
| `SU2_Cube_Gap_Certificate.json` | Dataset `single_cube_early_certificate`; earlier weaker 32-state construction, retained as history. |
| `SU2_Cube_Bounds_Results.json` | Dataset `single_cube_bounds`; exact omitted-state thresholds with numerical endpoint estimates, explicitly identified as mixed evidence. |
| `Shared_SU2_Cutoff_Study_Results.json` | Dataset `two_plaquette_study`; public numerical model and symmetric/antisymmetric correlator regressions. |
| `Shared_SU2_Cutoff_Certificate_Checks.json` | Dataset `two_plaquette_certificates`; reference for fresh fixed-case arithmetic checks. |
| `su2_joined_cubes_derivation.md`, `joined_cubes_volume_bound_audit.md` | Unchanged copies in `docs/theory/`, with their finite-graph assumptions and missing uniform-volume bounds. |
| `Yang_Mills_Research_Chronology.docx` | Unchanged historical context in this directory. Its speculative entries are hypotheses, not supported calculator capabilities. |

## Boundaries between models and evidence

The public joined-cube solver builds the same Hamiltonian as 0.1.0. It was not
replaced with standalone research code. The separate two-square benchmark
exposes low-level functions already present in the original recoupling module;
its seven-link graph does not represent two cubes.

A finite-matrix residual tests the computed eigenpair. A residual in a larger
test basis also samples omitted couplings within that larger basis. Neither
alone encloses the untruncated gap. The certificate programs add the comparison
construction, omitted-state control, radical rounding, and exact inertia probes.

The single-cube source and joined-graph source can choose different signs for
basis vectors. The verifier checks one consistent diagonal sign transformation
across all face transitions and all directly coupled omitted targets, preserving
indexing and exact spectra. It checks exact electric energies as well.

The independent higher-spin tensor tests use the external archived constructor,
not the package's recoupling helper. They compare every face on N=1 cutoff 12,
N=2 cutoff 8 with longitudinal/crossed/mixed pairings, and N=3 cutoff 6. The
independent cycle-space counts are 32, 868, and 25,676 at physical spin 1/2.
The full oriented-vertex sweep contains 3,368 exact contractions through doubled
spin 14. These are checks of the finite-model implementation; they do not
establish a continuum limit.

## Reading archived notes

Archived files are deliberately preserved, including their original relative
links, dates, exploratory language, and output-writing commands. Their `work/`
and `reference_results/` paths refer to the surrounding original reproduction
bundle. Execute the supported `terry-cubes certify` interface for fresh
verification without overwriting preserved records.

The current README and result metadata define the supported release scope.
The fixed-graph positive-gap argument does not give a chain-length-independent
constant. There is no package implementation of a continuum limit, physical
glueball calibration, Standard Model derivation, or gameplay dynamics.

## Provenance

The original three joined-cube JSON datasets and five reference scripts remain
unchanged. Added datasets and certificate ZIP hashes are recorded in
`RELEASE_PROVENANCE.json`; the ZIPs contain the original internal manifests.
`validation/original_package_sha256.json` snapshots the pre-update package.
Historical 0.1.0 validation is in `validation/history/`; fresh 0.1.1 results
are identified in `VALIDATION.md` and the accompanying JSON/XML files.
