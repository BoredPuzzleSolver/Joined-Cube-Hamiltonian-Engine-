"""Gauge-invariant SU(2) Hamiltonians on open cube chains."""
from .hamiltonian import JoinedCubesModel
from .solver import SpectrumSolver, SpectrumResult, SpectrumConvergenceError

__version__ = "0.1.0"
__all__ = ["JoinedCubesModel", "SpectrumSolver", "SpectrumResult", "SpectrumConvergenceError"]
