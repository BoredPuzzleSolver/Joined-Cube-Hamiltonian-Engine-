# Open SU(2) cube: reproducible spectra and certified gap bounds

## Result

At kappa=nu=1, the complete physical Hilbert space on a single open cube has
gap/kappa in the certified interval [2.941733, 2.958693]. A separate larger
numerical calculation estimates about 2.95104. The cube has 12 links, 8 vertices,
6 faces and Gauss law at every vertex, with no external charges.

Read outputs/SU2_Cube_Study_and_Gap_Bounds.md for the model, exact certificate,
numerical limitations, source references and next research steps. The result
concerns one finite spatial graph; no continuum or infinite-volume claim is made.

## Run

Requires Python 3.10 or later and NumPy. Tested with Python 3.12 and NumPy 2.3.5.
No network connection, SciPy or symbolic algebra package is needed.

From this directory:

```text
python reproduce.py
```

The full reproduction may take a few minutes. It checks the signed recoupling
matrices, independent exact Haar contractions, sparse eigensolver, both cutoff
studies, omitted-state bounds, rational certificates, integer interval arithmetic
and independent audits. It writes only outputs/ and prints completion per script.

reference_results/ holds the supplied baseline JSON. Generated floating-point
values are compared with 1e-9 relative/absolute tolerance for numerical-library
variation. Exact endpoint fractions, integer inertia counts and other discrete
results must match exactly. Runtime measurements are excluded from comparison.
The exact proofs do not rely on the floating-point comparison tolerance.

## Organization

- work/: complete source code and mathematical/audit notes.
- outputs/: report, calculated JSON, and reproduction logs.
- reference_results/: the supplied baseline results.
- SHA256SUMS.txt: hashes of the delivered files before any rerun.

The two spectral cutoffs are different: maximum spin per link, and total
electric energy. The latter is complete without an additional spin cap.
The most precise numerical endpoints are explicitly separate from certified ones.

The physical proof assumes the stated complete spin-network basis and signed
SU(2) recoupling identities. Exact arithmetic certifies matrix comparisons and
endpoint sign counts; it does not independently prove those physical identities.
Independent Haar checks cover the J=1/2 cube block, with further model checks
and source-derived formulas for higher spins.

The original research document and earlier two-square study were not modified.
