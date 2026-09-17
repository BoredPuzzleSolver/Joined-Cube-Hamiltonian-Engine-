"""Command-line interface to the same model and solver used by Python users."""
from __future__ import annotations

import argparse
import json
import os
from pathlib import Path
import sys
import tempfile

from . import __version__, JoinedCubesModel, SpectrumSolver, SpectrumConvergenceError
from .studies import cutoff_study
from .two_plaquette import TwoPlaquetteModel
from .certification import verify_certificate, CertificateVerificationError
from .data import list_research_datasets, load_research_dataset


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
    run.add_argument("--levels", type=int, default=3,
                     help="Requested low Ritz values, at least 2; single-vector Lanczos can omit degenerate copies")
    run.add_argument("--seed", type=int, default=1097)
    run.add_argument("--output", type=Path, help="Save JSON to this file; otherwise print JSON")
    study = commands.add_parser("study", help="Compare complete electric cutoffs on one cube chain")
    study.add_argument("--cubes", type=int, required=True)
    study.add_argument("--kappa", type=float, default=1.0)
    study.add_argument("--nu", type=float, default=1.0)
    study.add_argument("--cutoffs", nargs="+", required=True, help="Dimensionless electric cutoffs, e.g. 6 8 10")
    study.add_argument("--embedding-cutoff", help="Larger test basis; defaults to the largest requested cutoff")
    study.add_argument("--pairing", choices=("longitudinal", "crossed", "mixed"), default="longitudinal")
    study.add_argument("--backend", choices=("lanczos", "arpack"), default="arpack")
    study.add_argument("--tolerance", type=float, default=1e-10)
    study.add_argument("--max-iterations", type=int, default=1000)
    study.add_argument("--seed", type=int, default=1097)
    study.add_argument("--second-seed", type=int, default=8917)
    study.add_argument("--no-second-seed", action="store_true")
    study.add_argument("--output", type=Path)
    pair = commands.add_parser("two-plaquette", help="Separate seven-link two-square numerical benchmark")
    pair.add_argument("--kappa", type=float, default=1.0)
    pair.add_argument("--nu", type=float, default=1.0)
    pair.add_argument("--max-twice", type=int, default=12, help="Maximum doubled physical spin; 12 means j<=6")
    pair.add_argument("--tolerance", type=float, default=1e-10)
    pair.add_argument("--correlators", action="store_true", help="Include symmetric/antisymmetric imaginary-time channels")
    pair.add_argument("--times", type=float, nargs="+", default=[.25, .5, 1, 2, 4, 8])
    pair.add_argument("--dt", type=float, default=.25)
    pair.add_argument("--output", type=Path)
    certify = commands.add_parser("certify", help="Fresh verification of supported fixed-graph certificates")
    certify.add_argument("model", choices=("single-cube", "two-plaquette"))
    certify.add_argument("--cutoff", type=int, help="Single cube only: 8, 10, or 12 (default 12), kappa=nu=1")
    certify.add_argument("--timeout", type=float, help="Optional verification time limit in seconds")
    certify.add_argument("--output", type=Path)
    data = commands.add_parser("data", help="List or read additional historical research records")
    data.add_argument("--name", help="Dataset name from the list; reading a record does not verify it")
    data.add_argument("--output", type=Path)
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
        if args.command in ("run", "study"):
            model = JoinedCubesModel(
                num_cubes=args.cubes, kappa=args.kappa, nu=args.nu, pairing=args.pairing,
            )
        if args.command == "run":
            result = SpectrumSolver(
                model, electric_cutoff=args.cutoff, backend=args.backend,
                tolerance=args.tolerance, max_iterations=args.max_iterations,
                levels=args.levels, seed=args.seed,
            ).compute_mass_gap()
            payload = result.to_dict()
        elif args.command == "study":
            payload = cutoff_study(model, args.cutoffs, backend=args.backend,
                tolerance=args.tolerance, max_iterations=args.max_iterations,
                seed=args.seed, second_seed=None if args.no_second_seed else args.second_seed,
                embedding_cutoff=args.embedding_cutoff)
        elif args.command == "two-plaquette":
            result = TwoPlaquetteModel(args.kappa, args.nu).compute_spectrum(args.max_twice, tolerance=args.tolerance)
            payload = result.to_dict()
            if args.correlators:
                payload["correlation_channels"] = result.correlators(args.times, dt=args.dt)
        elif args.command == "certify":
            cutoff = 12 if args.model == "single-cube" and args.cutoff is None else args.cutoff
            payload = verify_certificate(args.model.replace("-", "_"), cutoff=cutoff, timeout=args.timeout)
        else:
            payload = load_research_dataset(args.name) if args.name else {"datasets": list_research_datasets()}
        text = json.dumps(payload, indent=2, allow_nan=False) + "\n"
        if args.output is None:
            sys.stdout.write(text)
        else:
            _write_json(args.output, text)
            print(f"Saved {args.output}")
            if "e0" in payload:
                print(f"E0={payload['e0']:.12g}  E1={payload['e1']:.12g}  gap={payload['gap']:.12g}")
        return 0
    except (ValueError, TypeError, OSError, SpectrumConvergenceError, CertificateVerificationError) as exc:
        print(f"terry-cubes: error: {exc}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
