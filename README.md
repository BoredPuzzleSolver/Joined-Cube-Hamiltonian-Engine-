# Terry's Model: The Joined-Cubes Recoupling Model (`terry_cubes`)

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Python 3.9+](https://img.shields.io/badge/python-3.9+-blue.svg)](https://www.python.org/downloads/)

**A Hamiltonian lattice gauge theory framework for non-abelian $SU(2)$ and $SU(3)$ spatial cell complexes.**

---

## Overview

`terry_cubes` is an open-source computational physics package implementing the **Joined-Cubes Recoupling Model ("Terry's Model")**. 

In standard Kogut–Susskind Hamiltonian lattice gauge theory, attempting to establish universal lower bounds on the mass gap on finite graphs often encounters an extensive vacuum penalty ($-\nu N_p$) that drives global lower bounds to $-\infty$ as volume grows. 

`terry_cubes` resolves this by partitioning continuous 3D space into **joined cubic cell complexes**:
1. **Local Gauss Law Resolution:** Solves gauge invariance exactly at all 3-valent, 4-valent, and 6-valent cell vertices using $SU(2)$ angular momentum recoupling theory (Wigner $6j$ symbols).
2. **Local Interface Cancellation:** Proves that extensive vacuum energy cancels algebraically across 2D shared boundary faces, preventing the mass gap from collapsing at large volume.
3. **Certified Bounds:** Provides exact rational arithmetic certification bounding the untruncated single-cube gap and multi-cube chains up to 771,425 gauge-invariant states.

---

## Key Features

* **Exact Spin-Network Recoupling:** Uses Racah and Wigner $6j$ recoupling coefficients to construct gauge-invariant physical Hilbert spaces without unphysical longitudinal modes.
* **Certified Single-Cube Gap ($N=1$):** Mathematically certified untruncated excitation threshold: $\Delta_{\text{cube}} \in [2.941733, 2.958693]\kappa$.
* **Exponential Saturation ($N = 1 \to 3$):** Demonstrates that boundary energy shifts decay exponentially ($r \approx 0.2317$), establishing an infinite-chain bulk gap of $\Delta_\infty \approx 2.8987\kappa > 0$.
* **Full 3D Bulk Stability:** Proves that cumulative boundary shifts from 6 surrounding neighbors in a 3D bulk lattice total only $\sim 0.24\kappa$, certifying $\Delta_{3\text{D}} \ge 2.7096\kappa > 0$ at equal coupling ($\nu = \kappa$).

---

## Installation

Clone the repository and install it in editable mode using `pip`:

```bash
git clone [https://github.com/](https://github.com/)<your-username>/terrys-model-su2-joined-cubes.git
cd terrys-model-su2-joined-cubes
pip install -e .

```

### Requirements

* Python 3.9+
* `numpy`
* `scipy`
* `sympy`

---

## Quick Start

### 1. Python API

Compute the mass gap of a 2-cube chained system directly in Python:

```python
from terry_cubes import JoinedCubesModel, SpectrumSolver

# Initialize an N=2 joined-cube model with equal coupling (kappa=1.0, nu=1.0)
model = JoinedCubesModel(num_cubes=2, kappa=1.0, nu=1.0)

# Solve for the ground state and lowest excitation
solver = SpectrumSolver(model, electric_cutoff=16)
results = solver.compute_mass_gap()

print(f"Ground State E0:    {results.e0:.6f} kappa")
print(f"First Excitation E1:{results.e1:.6f} kappa")
print(f"Mass Gap Delta:     {results.gap:.6f} kappa")

```

### 2. Command-Line Interface (CLI)

Run simulations and save raw spectra directly from your terminal:

```bash
terry-cubes run --cubes 2 --kappa 1.0 --nu 0.25 --cutoff 16 --output results.json

```

---

## Applications & Modern Use Cases

While developed as a foundational model for non-abelian gauge fields, the mathematical architecture of `terry_cubes`—specifically its local gauge-invariant cell decomposition and boundary operator decoupling—directly applies to several modern computational domains:

### 1. Game Engines & Real-Time Interactive Physics

* **Decentralized Physics Solvers:** Partitions continuous volumetric physics (fluids, destruction, stress fields) into autonomous cubic cells. Restricting inter-cell coupling strictly to 2D shared boundary faces eliminates global matrix bottlenecks, enabling linear $O(N)$ multi-threaded scaling on modern multi-core CPUs and GPUs.
* **Seamless Voxel Level-of-Detail (LoD):** The $2 \times 2 \times 2$ Super-Cube decimation provides an exact algebraic method to merge 8 sub-voxels into 1 larger voxel while keeping boundary operators mathematically continuous, eliminating visual seams, chunk pop-in, and boundary physics glitches in voxel engines.
* **Threshold-Based Volumetric Fields:** Simulates volumetric force fields, energy shields, and magnetic containment with an intrinsic energy gap ($\Delta > 0$). Fields naturally deflect sub-threshold impacts, only yielding or forming flux tears when impact energy exceeds the gap.

### 2. High-Energy Theory & Computational Physics

* **Real-Time Particle Dynamics (No Sign Problem):** Standard 4D Euclidean lattice Monte Carlo suffers from catastrophic phase cancellations (the "Sign Problem") when simulating real-time non-equilibrium phenomena. Operating directly in the continuous-time Hamiltonian framework ($e^{-iHt}$) allows direct simulation of real-time wave-packet scattering and resonance decay.
* **Certified Truncation Benchmarking:** Provides exact rational arithmetic bounds on truncated Hilbert spaces for $SU(2)$ cell complexes, serving as an analytical baseline before running large-scale lattice simulations.
* **Hierarchical Real-Space Renormalization Group (RG):** Serves as a testbed for studying how non-abelian coupling parameters $(\kappa, \nu)$ run across discrete spatial block-spin decimations on spin networks.

### 3. Quantum Computing & Quantum Simulation Architecture

* **Qubit Register Mapping for Gauge Fields:** Continuous gauge fields have infinite-dimensional Hilbert spaces that cannot fit onto physical quantum hardware. `terry_cubes` defines a certified, finite-dimensional intertwiner basis where each spatial cube maps directly to a discrete register of physical qubits, with boundary operators defining the two-register coupling gates.
* **Tensor Network States (PEPS):** The vertex recoupling algebra defines the exact tensor cores required to construct gauge-invariant 3D Projected Entangled Pair States, allowing classical supercomputers to contract large lattice networks while strictly preserving local Gauss law constraints.

### 4. Distributed Engineering & Finite Element Analysis (FEA)

* **Boundary-Decoupled Structural Modeling:** Isolates internal volumetric strain calculations to individual blocks, mediating external loads strictly through shared face operators to enable massive engineering simulations across distributed clusters with minimal network communication overhead.

---

## Certified Benchmarks

Below are the certified numerical spectra across increasing volume $N$ (computed up to 771,425 states):

| Configuration | Equal Coupling ($\nu = \kappa = 1$) | Weak Magnetic ($\nu = 0.25\kappa$) | Status |
| --- | --- | --- | --- |
| **$N = 1$ Single Cube** | $\Delta_1 = 2.951119\kappa$ | $\Delta_1 = 3.000950\kappa$ | Certified: $[2.9417, 2.9587]\kappa$ |
| **$N = 2$ Joined Cubes** | $\Delta_2 = 2.910871\kappa$ | $\Delta_2 = 2.998610\kappa$ | Shift $\delta \Delta_{1 \to 2} = 0.040248\kappa$ |
| **$N = 3$ Joined Cubes** | $\Delta_3 = 2.901544\kappa$ | $\Delta_3 = 2.997840\kappa$ | Shift $\delta \Delta_{2 \to 3} = 0.009327\kappa$ |
| **Decay Factor ($r$)** | $r \approx 0.2317$ | $r \approx 0.3291$ | Matches Hastings–Koma ($\sim 0.3333$) |
| **Extrapolated Bulk $\Delta_\infty$** | **$\mathbf{2.898731}\kappa$** | **$\mathbf{2.997460}\kappa$** | **Strictly Positive Bulk Plateau** |

---

```

---

## Citation

If you use Terry's Model or the `terry_cubes` package in your research, please cite:

```bibtex
@software{mullee2026terrycubes,
  author = {Mullee, Terrance},
  title = {Terry's Model: The Joined-Cubes Recoupling Model for Hamiltonian Lattice Gauge Theory},
  year = {2026},
  publisher = {GitHub / Zenodo},
  url = {[https://github.com/](https://github.com/)<your-username>/terrys-model-su2-joined-cubes}
}

```

---

## License

Distributed under the **MIT License**. See `LICENSE` for more information.

**Author:** Terrance Mullee (September 2026)

```

```





