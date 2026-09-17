# Joined SU(2) cubes — reproduction

Start with **SU2_Joined_Cubes_Volume_Study.md** for the results and
**Space_Mathematical_Reality_and_Proof.md** for the conceptual discussion.

## Dependencies

Python 3.12 and NumPy were used (NumPy 2.3.5). Only NumPy and Python's standard
library are required. Run commands from this directory. The three-cube cap 16
run retains 771,425 states, uses several gigabytes of memory, and takes several
minutes on the original machine. Exact timing depends on hardware and BLAS.

## Independent model checks

```text
python su2_joined_cubes_independent_check.py
python su2_joined_orientation_validation.py
```

This checks complete spin-half basis counts, local invariant tensors,
three resolutions of the two-cube graph, and finite higher-spin matrix fixtures
against a separately implemented angular-momentum tensor calculation.
The independent checker writes SU2_Joined_Cubes_Independent_Validation.json.
The second command checks exact rational local tensor contractions through
doubled spin 14 and the old/final basis-phase identities for chains 1–4. Its
range covers the link and coupling labels allowed by the delivered energy caps.

## Reproduce the numerical study

```text
python su2_joined_cubes_study.py --cubes 1 --max-cap 18
python su2_joined_cubes_study.py --cubes 2 --max-cap 18
python su2_joined_cubes_study.py --cubes 3 --max-cap 16
```

Each command writes the corresponding N1/N2/N3 Spectra JSON in this directory.
The output includes all nested caps, κ=1 with ν=0.25 and 1, direct retained-matrix
residuals and a second starting-vector check on each largest retained matrix.
Compare reported energies at approximately 1e-9 tolerance; timing, residual
roundoff and iteration counts may vary between numerical libraries.

For a quick import and small-case check:

```text
python su2_joined_cubes_study.py --cubes 1 --max-cap 8 --caps 6 8
```

This quick command overwrites the local N1 JSON with the smaller run; retain a
copy of the reference outputs or extract a fresh copy of the ZIP first.

## Interpretation

The supplied numerical gaps are finite-basis Ritz estimates. They are not
certified enclosures of the untruncated joined-cube gap. The accompanying
analytic note proves qualitative gap positivity for each finite graph but
gives no positive bound uniform in volume. Increasing N in N×1×1 is a chain
limit with fixed transverse width, not a full 3D thermodynamic limit or a
continuum construction. Earlier single-cube certificates referenced in the
report are in the separate SU2_Cube_Reproduction package from the prior study.

The original research DOCX was not modified. SHA256SUMS.json records packaged
file hashes at creation; rerunning scripts changes the result files.
