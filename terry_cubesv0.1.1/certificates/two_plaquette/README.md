# Shared SU(2) plaquettes: cutoff study and exact gap bounds

## Result and scope

This package studies two square plaquettes sharing one link, with SU(2) link
variables, open boundary, and Gauss law at every vertex. It includes the
full-spin fixed-graph gap certificates [3.09195898, 3.09195907] at nu/kappa=1
and [4.20942742, 4.20942757] at nu/kappa=4. Energies are in units of kappa.

Read `outputs/Shared_SU2_Cutoff_Study_and_Certified_Gap.md` for the derivation,
results, sources, limitations, and next research steps. These certificates
do not establish the infinite-volume or continuum Yang-Mills mass gap.

## Run

Requires Python 3.10 or later and NumPy. Tested with NumPy 2.3.5.
From this directory, run:

```text
python reproduce.py
```

The studies require no network access, SciPy, or symbolic algebra package.
The command writes only the package's outputs directory. It reruns:

1. The 72-point coupling/cutoff study and connected correlators.
2. The exact rational certificates, with fixed published probe values.
3. Independent SU(2) Haar integration and comparison to magnetic recoupling.
4. A second exact rational elimination check of all twelve certificate probes.

It then compares generated JSON with the supplied reference_results directory.
Exact endpoints, counts and other discrete values must match exactly.
Floating-point values allow 1e-10 relative/absolute variation across numerical
libraries; elapsed timing fields and heuristic roundoff flags are ignored.
Those diagnostic flags can vary across platforms near machine precision; the exact certificate checks
are the authoritative endpoint verification.

Individual scripts can also be run from work/. Files use paths relative to
their own location, so no original Windows path is needed. Scripts intentionally
fail with an assertion or exception if a required check fails.

## Reproducibility and provenance

The source is included in work/. The numerical and exact results are separated
in outputs/. Mathematical notes and an independent audit are also in work/.
SHA256SUMS.txt records hashes of the delivered files. Rerunning replaces result
and log files, including timing data, so their hashes can change afterward.

The main model and normalization are explicit in su2_recoupling.py. The
analytic bound requires kappa>0, nu>=0, the physical admissible basis, and the
complete first omitted layer. Generic helper functions must not be used with
arbitrary inputs while inheriting the physical certificate claim.

This is a reproducible research calculation, not a claim of a new theorem in
continuum Yang-Mills theory. The original user document was left unchanged.
