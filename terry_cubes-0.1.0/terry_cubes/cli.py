"""Command-line interface to the same model and solver used by Python users."""
from __future__ import annotations

import argparse
import json
import os
from pathlib import Path
import sys
import tempfile

from . import __version__, JoinedCubesModel, SpectrumSolver, SpectrumConvergenceError


def build_parser():
    parser = argparse.ArgumentParser(
        prog="terry-cubes",
        description="Numerical finite-basis SU(2) cube-chain spectra.",
    )
    parser.add_argument("--version", action="version", version=f"%(prog)s {__version__}")
    commands = parser.add_subparsers(dest="command", required=True)
    run = commands.add_parser("run", help="Compute the lowest levels and mass-gap estimate")
    run.add_argument("--cubes", type=int, required=True, help="Number of cubes in an open N-by-1-by-1 chain")
    run.add_argument("--kappa", type=float, default=1.0, help="Positive electric coefficient (default: 1)")
    run.add_argument("--nu", type=float, default=1.0, help="Nonnegative magnetic coefficient (default: 1)")
    run.add_argument("--cutoff", default="16", help="Cutoff on physical sum j(j+1), e.g. 16 or 33/2")
    run.add_argument("--pairing", choices=("longitudinal", "crossed", "mixed"), default="longitudinal")
    run.add_argument("--backend", choices=("lanczos", "arpack"), default="lanczos")
    run.add_argument("--tolerance", type=float, default=1e-10, help="Absolute retained-matrix residual tolerance")
    run.add_argument("--max-iterations", type=int, default=240)
    run.add_argument("--levels", type=int, default=3, help="Number of low eigenvalues, at least 2")
    run.add_argument("--seed", type=int, default=1097)
    run.add_argument("--output", type=Path, help="Save JSON to this file; otherwise print JSON")
    return parser


def _write_json(path, text):
    """Atomically replace the requested output after a successful calculation."""
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = None
    try:
        with tempfile.NamedTemporaryFile(
            mode="w", encoding="utf-8", dir=path.parent,
            prefix=f".{path.name}.", suffix=".tmp", delete=False,
        ) as stream:
            temporary = Path(stream.name)
            stream.write(text)
        os.replace(temporary, path)
    finally:
        if temporary is not None and temporary.exists():
            temporary.unlink()


def main(argv=None):
    parser = build_parser()
    args = parser.parse_args(argv)
    try:
        model = JoinedCubesModel(
            num_cubes=args.cubes, kappa=args.kappa, nu=args.nu, pairing=args.pairing,
        )
        result = SpectrumSolver(
            model, electric_cutoff=args.cutoff, backend=args.backend,
            tolerance=args.tolerance, max_iterations=args.max_iterations,
            levels=args.levels, seed=args.seed,
        ).compute_mass_gap()
        text = json.dumps(result.to_dict(), indent=2, allow_nan=False) + "\n"
        if args.output is None:
            sys.stdout.write(text)
        else:
            _write_json(args.output, text)
            print(f"Saved {args.output}")
            print(f"E0={result.e0:.12g}  E1={result.e1:.12g}  gap={result.gap:.12g}")
        return 0
    except (ValueError, TypeError, OSError, SpectrumConvergenceError) as exc:
        print(f"terry-cubes: error: {exc}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
