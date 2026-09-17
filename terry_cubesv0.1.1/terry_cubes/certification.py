"""Run the preserved finite-graph verifiers with fixed published probes.

The reference programs and their proof notes are preserved separately from
the numerical solver. Verification takes place in a temporary directory;
installed sources, reference results, and historical logs are never changed.
"""
from hashlib import sha256
from importlib.resources import files
from io import BytesIO
from pathlib import Path, PurePosixPath
import json
import math
import subprocess
import sys
import tempfile
from zipfile import ZipFile


class CertificateVerificationError(RuntimeError):
    """A source-integrity, model-bridge, or exact arithmetic check failed."""


def verify_certificate(model="single_cube", *, cutoff=None, timeout=None):
    """Rebuild and check a supported reference certificate.

    ``single_cube`` accepts complete electric cutoffs 8, 10, or 12 at
    kappa=nu=1. ``two_plaquette`` checks the three original fixed probesets
    at ratios 1 and 4; pass cutoff=None for that different graph. These
    settings are deliberately not a general-volume certification API.
    """
    if model not in ("single_cube", "two_plaquette"):
        raise ValueError("model must be single_cube or two_plaquette")
    if model == "single_cube":
        if cutoff is None:
            cutoff = 12
        if isinstance(cutoff, bool) or cutoff not in (8, 10, 12):
            raise ValueError("single_cube certificate cutoff must be 8, 10, or 12")
        cutoff = int(cutoff)
    elif cutoff is not None:
        raise ValueError("two_plaquette uses its own fixed spin cutoffs; pass cutoff=None")
    if timeout is not None and (not math.isfinite(timeout) or timeout <= 0):
        raise ValueError("timeout must be finite and positive")
    data = files("terry_cubes.data")
    archive_bytes = data.joinpath(model + "_certificate_sources.zip").read_bytes()
    expected = json.loads(data.joinpath("certificate_archive_hashes.json").read_text())[model]
    if sha256(archive_bytes).hexdigest() != expected:
        raise CertificateVerificationError("The preserved certificate archive failed its SHA-256 check")
    with tempfile.TemporaryDirectory(prefix="terry-cubes-certificate-") as temporary:
        root = Path(temporary)
        with ZipFile(BytesIO(archive_bytes)) as archive:
            for member in archive.infolist():
                relative = PurePosixPath(member.filename)
                if relative.is_absolute() or ".." in relative.parts or "\\" in member.filename or ":" in member.filename:
                    raise CertificateVerificationError("Unsafe path in the certificate archive")
                target = root.joinpath(*relative.parts)
                if member.is_dir():
                    target.mkdir(parents=True, exist_ok=True)
                else:
                    target.parent.mkdir(parents=True, exist_ok=True)
                    target.write_bytes(archive.read(member))
        worker = Path(__file__).with_name("_certificate_worker.py")
        command = [sys.executable, "-I", str(worker), model, str(root), str(cutoff)]
        try:
            process = subprocess.run(command, text=True, encoding="utf-8", capture_output=True,
                                     cwd=root, timeout=timeout)
        except subprocess.TimeoutExpired as exc:
            raise CertificateVerificationError("Certificate verification exceeded the requested timeout") from exc
        if process.returncode:
            raise CertificateVerificationError("Exact verification failed: " + process.stderr.strip())
        try:
            result = json.loads(process.stdout)
        except ValueError as exc:
            raise CertificateVerificationError("The verifier did not return a complete JSON result") from exc
        if not result.get("certified") or not result.get("package_bridge_verified"):
            raise CertificateVerificationError("The verifier did not establish all required checks")
        result["source_archive_sha256"] = expected
        return result
