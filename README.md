# Terry Cubes — Joined-Cube Hamiltonian Engine

**A Python research toolkit for exploring finite quantum gauge systems, reproducing spectral calculations, and developing research benchmarks.**

Terry Cubes (`terry_cubes`) calculates the allowed quantum states and low-energy behavior of an SU(2) lattice gauge model on open chains of joined cubes. It constructs states that satisfy the model's local symmetry constraints, builds the Hamiltonian—the operator describing the system's energy—and computes energy levels, excitation gaps, and selected observables.

In practical terms, it lets you investigate how a small quantum system changes when you adjust its couplings, connect more cubes, or increase the number of states retained in a calculation. Reference datasets, numerical diagnostics, and separate verification tools make those calculations reproducible and useful for comparison.

**Current release: 0.1.1.** The project supports finite-system research, numerical-method development, and education. Game applications and broader physical simulations are potential extensions described below.

## What the tool does today

- **Builds gauge-invariant quantum bases** for open `N × 1 × 1` chains of joined cubes, including the compatible coupling channels at joining vertices.
- **Constructs sparse Hamiltonians** using exact local SU(2) recoupling data, converted to floating point for numerical calculations.
- **Computes low-energy spectra and eigenvectors**, including ground-state energies, excitation energies, and spectral-gap estimates.
- **Compares electric cutoffs and solver results**, helping users examine numerical convergence and sensitivity to omitted states.
- **Calculates selected plaquette observables and imaginary-time correlations** within the retained quantum basis.
- **Provides reproducible reference calculations** for one, two, and three cubes, plus a separate two-square benchmark.
- **Includes fixed-case certificate verifiers** that establish bounds for specified finite models under their documented assumptions.
- **Offers a Python API and command-line interface**, with JSON output, examples, research records, and automated tests.

The implementation uses established Hamiltonian lattice gauge methods. Its contribution is the specific joined-cube implementation, documented conventions, reference calculations, and connection to reproducible finite-graph verification.

## Practical uses

| Application | How Terry Cubes can help |
| --- | --- |
| **Finite quantum-system research** | Explore low-energy spectra, coupling dependence, operator overlaps, and cutoff behavior on the supported graphs. |
| **Numerical-method benchmarking** | Compare eigensolvers, truncation strategies, and independent implementations against documented reference cases. |
| **Quantum-algorithm development** | Supply classical reference energies and states for a matching quantum simulation or variational calculation. Qubit encodings and hardware integration would be additional work. |
| **Computational mathematics** | Study concrete examples of constrained state spaces, angular-momentum recoupling, and finite-graph spectral bounds. |
| **Education and visualization** | Build demonstrations of spin networks, local symmetry constraints, energy levels, and how numerical approximations change a result. |
| **Testing approximate models** | Generate small-system reference data against which proposed reduced models or other approximations can be evaluated. |

Benchmark comparisons must match the geometry, boundary conditions, Hamiltonian normalization, and retained state space. The project does not currently claim measured speed or accuracy advantages over other research packages.

## Potential game and interactive applications

The existing calculations could provide a mathematical basis for educational games and deliberately fictional behavior. Possible applications include:

- **Quantum puzzles:** players choose allowed link states, satisfy local constraints, or change couplings to reach a target energy pattern.
- **Fictional materials and devices:** a game adapter maps calculated energy gaps or observables to a shield's activation threshold, a reactor's operating modes, or a fictional material's response.
- **Procedural visuals and sound:** selected observables or energy differences drive colors, patterns, animation parameters, or musical relationships.
- **Interactive simulation exhibits:** users explore how a small quantum model responds to changes in its parameters and computational cutoff.

A practical first prototype could precompute a small library of results and expose them through a game interface, while the game's existing physics system handles movement and contact. These mappings would be authored gameplay rules; their usefulness and performance would need to be tested.

Further development could add real-time quantum evolution, state preparation, and measurement rules for small interactive models. The current imaginary-time correlation calculations do not provide that complete interaction loop. New fictional interaction laws would also require new operators or an explicit adapter.

**These are proposed applications. Version 0.1.1 does not include a game-engine adapter or a general gameplay physics system.**

## Possible future versions: mechanics and materials

A broader version could introduce a separate mechanics module, or integrate with an existing simulator, to explore friction, stress, deformation, and different materials. That would require new physical models, state variables, and validation alongside the current quantum calculations.

| Possible addition | What would need to be implemented | Potential use after validation |
| --- | --- | --- |
| **Motion, collisions, and friction** | Positions, velocities, masses, contact detection, friction laws, and stable time integration. | Interactive mechanical systems and object-contact simulations. |
| **Elasticity and stress** | Displacement and strain fields, material response laws, loads, boundary conditions, and a suitable mechanical solver. | Deformation studies and comparisons of modeled structural responses. |
| **Different material behaviors** | Calibrated density, stiffness, damping, plasticity, or other relevant properties. | Comparing how specified material models respond to the same conditions. |
| **Damage and fracture** | Failure criteria, damage evolution, and methods for updating connectivity as objects break. | Destruction effects and controlled fracture simulations. |

This direction could support fictional material systems, teaching tools, and eventually selected real-world modeling tasks if the added models were calibrated and independently validated. Friction and other dissipative effects need explicit treatment of energy loss.

The present SU(2) spectral gap does not directly determine friction, stiffness, stress, or material strength. Any physical connection would need to be derived or validated; assigning a connection for a game is a design choice. Existing quantum certificates would continue to apply only to the models and parameters they verify.

These are possible development directions, with no committed release schedule.

## Reference results and scope

For one open cube at `nu / kappa = 1`, the reference calculations give a numerical gap near **`2.951043 × kappa`**. A separate verified certificate bounds the gap in the full spin space of that fixed graph within **`[2.941733, 2.958693] × kappa`**.

Version 0.1.1 includes tools to rerun the specified certificate checks. Ordinary numerical results remain marked `certified: false`; loading a saved certificate does not perform fresh verification. Multi-cube spectra are finite-basis numerical results, and the two-square benchmark concerns a different graph.

The current model is a finite SU(2) cube chain. Increasing its length extends one spatial direction, and the number of quantum states can grow rapidly. The supplied three-cube calculation with 771,425 basis states describes three spatial cubes with many quantum configurations. It does not establish real-time performance for large worlds.

These calculations do not establish an infinite-volume or continuum Yang–Mills mass gap, or calibrated predictions for real materials.

## Getting started

Use the **`terry_cubesv0.1.1`** directory for the current release. It requires Python 3.10 or later, NumPy, SciPy, and Sym
