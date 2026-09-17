"""Bundled, unchanged finite-basis spectra from the original research runs."""

from importlib.resources import files
import json
from numbers import Integral

__all__ = ["load_spectra", "list_research_datasets", "load_research_dataset"]

_RESEARCH = {
    "single_cube_energy": ("SU2_Cube_Energy_Study_Results.json", "single_cube", "complete_electric_energy", "numerical"),
    "single_cube_spin": ("SU2_Cube_Study_Results.json", "single_cube", "physical_spin", "numerical"),
    "single_cube_energy_certificates": ("SU2_Cube_Energy_Certificates.json", "single_cube", "complete_electric_energy", "saved_certificate"),
    "single_cube_early_certificate": ("SU2_Cube_Gap_Certificate.json", "single_cube", "physical_spin", "saved_certificate"),
    "single_cube_bounds": ("SU2_Cube_Bounds_Results.json", "single_cube", "physical_spin", "mixed_exact_thresholds_and_numerical_endpoints"),
    "two_plaquette_study": ("Shared_SU2_Cutoff_Study_Results.json", "two_plaquette", "physical_spin", "numerical"),
    "two_plaquette_certificates": ("Shared_SU2_Cutoff_Certificate_Checks.json", "two_plaquette", "physical_spin", "saved_certificate"),
}


def list_research_datasets():
    """List separately identified graph/cutoff records; reading is not verification."""
    return [{"name": key, "filename": value[0], "model_id": value[1],
             "cutoff_type": value[2], "status": value[3]}
            for key, value in _RESEARCH.items()]


def load_research_dataset(name):
    """Load an unchanged historical record with explicit model and status metadata.

    A saved certificate record does not mean it was verified in this call.
    Use verify_certificate() for fresh arithmetic verification.
    """
    if not isinstance(name, str) or name not in _RESEARCH:
        raise ValueError("Unknown research dataset; use list_research_datasets()")
    filename, model, cutoff, status = _RESEARCH[name]
    return {"name": name, "model_id": model, "cutoff_type": cutoff, "status": status,
            "verified_in_this_call": False,
            "record": json.loads(files(__package__).joinpath(filename).read_text(encoding="utf-8"))}


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
