# Changelog

## 0.1.1 — 2026-09-17

Research-completeness update to the preserved 0.1.0 finite SU(2) library.

- Added explicit access to eigenvectors in the original state order, nested
  electric-cutoff studies, larger-test-basis residuals, and a second-seed check.
- Added normalized plaquette operators and connected imaginary-time spectral
  sums with partial-spectrum and symmetry-selection metadata.
- Exposed the existing seven-link two-square benchmark as `TwoPlaquetteModel`,
  keeping its graph and spin-cutoff conventions distinct from joined cubes.
- Added fresh verification interfaces for the preserved single-cube and
  two-plaquette certificates. The verifier checks source integrity and agreement
  with the actual package before checking exact endpoints. The single-cube
  cutoff-12 result is `[2.941733, 2.958693] * kappa` at `nu/kappa = 1`.
- Bundled seven additional unchanged research datasets and complete original
  certificate source bundles in installable wheels. Added the independent joined
  archive, theory notes, chronology, evidence index, and provenance to source releases.
- Added routine interface regressions and opt-in independent tensor, cycle-space,
  full orientation, higher-cutoff, and exact-certificate checks.
- Updated installation, CI, evidence-status, units, and scope documentation.
  Infinite-volume, continuum, physical-mass, and gameplay claims remain unestablished.

`geometry.py`, `hamiltonian.py`, and `recoupling.py` are byte-for-byte unchanged
from 0.1.0. The mathematical bodies of `SymmetricOperator` and `smallest_ritz`
are unchanged. Existing model/solver entry points, state indexing, baseline
data, reference scripts, license, and default numerical backend are preserved.
Output metadata gains explicit model identifiers; new features are additive.

## 0.1.0

Original packaged refactor. Its validation report and README are retained under
`validation/history/`; its original package folder was not modified.
