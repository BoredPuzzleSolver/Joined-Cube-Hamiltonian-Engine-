# terry_cubes 0.1.1 — validation record

Validation date: 17 September 2026. This release was built from a separate
copy of 0.1.0. All 37 files in the original package, including its existing
distribution artifacts, were rehashed and found unchanged; its file set also
matched the pre-update snapshot.

## Fresh checks executed for this release

| Check | Result |
| --- | --- |
| Routine source suite | 94 passed; 17 opt-in research and 40 opt-in large cases skipped; 4.77 seconds. |
| Independent research suite | 17 passed; 134 other cases deselected; 58.41 seconds. |
| Combined routine and research suite against an installed wheel outside the source checkout | 111 passed; 40 large cases skipped; 62.92 seconds. |
| Actual installed console entry point | Version 0.1.1; fresh single-cube certificate commands at cutoffs 8 and 10 passed. |
| Strongest single-cube certificate | Fresh cutoff-12 construction and exact endpoint checks passed against source and installed package. |
| Two-plaquette certificates | All three original fixed cases passed, including exact agreement with the package's existing recoupling implementation. |
| Numerical examples | Cutoff comparison, complete small retained cube spectrum/correlator, and two-square symmetry correlators ran successfully. |
| Original package preservation | All 37 original files unchanged; no files added to or removed from the original folder. |
| Mathematical core | `geometry.py`, `hamiltonian.py`, `recoupling.py` byte-identical to 0.1.0; `smallest_ritz` and `SymmetricOperator` syntax trees unchanged. |
| Reference integrity | Original three datasets and five reference scripts unchanged; all 42 cube, 19 two-square, and 22 joined-archive manifest entries verified. |
| Installed resources | All JSON resources byte-identical to the release data; both complete certificate ZIP hashes verified; license present. |

Source and wheel runs cover **111 distinct selected test cases**. The wheel
run repeats them against the installed distribution; it does not add 111 new
mathematical cases. Distribution construction and archive-completeness checks
are recorded separately in `dist/BUILD_CHECKS.json`, with artifact hashes in
`dist/SHA256SUMS.json`.

JUnit records are `validation/routine_tests.xml`, `research_tests.xml`, and
`installed_wheel_tests.xml`. Environment and preservation checks are in
`validation/preservation_fresh.json`; installed-package paths, resource checks,
and tested runtime hashes are in `validation/installed_wheel_checks.json`.

## Single-cube arithmetic result

At `kappa = nu = 1`, the complete electric-cutoff-12 certificate has 398
retained states and a 2,042-state comparison construction, including 1,644
directly coupled omitted states. The fresh result is:

| Quantity | Exact decimal endpoints |
| --- | --- |
| E0 | [5.506419, 5.506607] |
| E1 | [8.448340, 8.465112] |
| E1 - E0 | [2.941733, 2.958693] |

The four inertia probes, rational radical-rounding bounds, and conservative
remaining electric threshold `49/4` were recomputed. The wrapper separately
checked the actual package basis, transition support and exact squared
amplitudes, one common basis-sign transformation across all faces and omitted
targets, and exact electric energies on all comparison states. See
`validation/single_cube_certificate_fresh.json` for the full result.

The certificate is for the full spin space on this fixed single cube, under
the preserved representation, positivity, and omitted-state assumptions. It
is not a multi-cube, infinite-volume, or continuum certificate. The numerical
solver still returns `certified: false`.

## Independent and extended numerical checks

- All 3,368 exact oriented-vertex contractions through doubled spin 14 and
  the common orientation-phase identity for N=1,2,3,4 passed.
- An independent cycle-space enumeration reproduced the complete physical
  spin-half state lists for N=1,2,3: 32, 868, and 25,676 states.
- Independent Casimir/tensor face matrices agreed up to one common diagonal
  basis-sign transformation for N=1 cutoff 12, N=2 cutoff 8 with all three
  pairings, and N=3 cutoff 6.
- Dense diagonalization checked the low ARPACK levels for N=2 and N=3 at
  cutoff 8. The first two Lanczos energies agreed; repeated degenerate copies
  are deliberately not claimed for the single-vector algorithm.
- Four larger single-cube calculations, at cutoffs 24 and 30 and ratios 1
  and 4, reproduced the saved energies and gaps to absolute tolerance 1e-9.
  The new cutoff-30 ratio-one numerical gap is approximately
  `2.951042834527371`, in the 34,680-state basis. It is a numerical estimate,
  not an improved exact enclosure.
- The two-square symmetric and antisymmetric correlators reproduced the
  recorded symmetry selection, with first detected gaps approximately
  3.091959020 and 3.113957723 at ratio one.

Fresh numerical and certificate records live directly under `validation/`;
the three numerical example outputs are under `validation/example_runs/`.

## Environment and limits of this validation

Executed locally on Windows 11 with Python 3.12.14, NumPy 2.3.5, SciPy 1.18.1,
SymPy 1.14.0, and pytest 9.1.1, with numerical thread counts limited to one.
Floating-point comparisons use tolerances, not bitwise eigenvalue equality.

The 40 large original joined-cube baseline runs were not repeated for this
release. Their previous successful execution belongs to 0.1.0 and is retained
as historical evidence in `validation/history/VALIDATION_0.1.0.md`. The core
mathematics is unchanged and the current suites exercise its original
preservation checks, representative spectra, new interfaces, extended
single-cube calculations, and independent constructions.

The CI workflow is included and updated, but no remote Linux/Windows CI run
or publication was performed. Infinite-volume and continuum conclusions,
physical glueball masses, and gameplay applications remain unestablished.
