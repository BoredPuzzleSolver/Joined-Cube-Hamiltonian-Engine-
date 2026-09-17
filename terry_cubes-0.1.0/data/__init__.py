"""Bundled, unchanged finite-basis spectra from the original research runs."""

from importlib.resources import files
import json
from numbers import Integral

__all__ = ["load_spectra"]


def load_spectra(num_cubes: int) -> dict:
    """Load the original JSON record for one, two, or three joined cubes.

    The record includes its model, cutoff, normalization, solver diagnostics,
    and finite-basis status. These data are reference calculations, not a
    certificate for the untruncated Hamiltonian.
    """
    if isinstance(num_cubes, bool) or not isinstance(num_cubes, Integral):
        raise ValueError("num_cubes must be one of the integers 1, 2, or 3")
    if num_cubes not in (1, 2, 3):
        raise ValueError("reference spectra are available only for 1, 2, or 3 cubes")
    name = f"SU2_Joined_Cubes_N{int(num_cubes)}_Spectra.json"
    return json.loads(files(__package__).joinpath(name).read_text(encoding="utf-8"))
