# Joined SU(2) cubes: volume, cutoff and gap study

**Research results — 15 September 2026**

## Outcome

The positive gap survives joining one, two and three cubes in the numerical tests. Its estimated value changes with the volume and with the basis cutoff. A standard analytic argument also proves a positive gap for every fixed finite graph in this family, without truncating the spins. Neither result establishes a positive lower bound uniform in volume or a continuum Yang–Mills mass gap.

This study uses actual SU(2) gauge states on shared links. The extra invariant coupling channels at the shared vertices are included. It does not replace those channels with an independent-cube approximation.

## 1. Precisely specified model

The open lattice has shape N×1×1, with N cubes in a row. It has 4(N+1) physical vertices, 8N+4 physical links and 5N+1 square faces. Each pair of neighboring cubes shares a face, including its four physical links and vertices. A shared face occurs once in the Hamiltonian.

\[
H_N=\kappa\sum_{e\text{ physical}}j_e(j_e+1)
+\nu\sum_{p=1}^{5N+1}\left(1-\frac12\operatorname{Tr}U_p\right).
\]

Gauss law is imposed at every vertex, including the boundary. We set κ=1 and separately test ν=1 and ν=0.25. The reported energies are in these model units; no physical lattice spacing or energy in electronvolts has been inferred.

Four links meet at the interior shared vertices. Their representations can combine to a gauge singlet in more than one way. For example, four spin-½ links have two independent singlet channels. We represent these channels by resolving a four-leg vertex into two three-leg vertices with a virtual coupling label. This virtual label has zero electric energy and is not restricted by a physical-link spin cutoff. All compatible values are included.

The finite basis contains **every physical state with total electric energy T≤C**, including every compatible invariant coupling channel. This is a complete energy cutoff, not a restriction to a few selected loops. Virtual labels are automatically bounded by the physical spins they combine.

## 2. Numerical spectra at equal electric and magnetic coefficients

Here κ=ν=1. Each row uses the largest complete energy cutoff run for that volume.

| Joined cubes | Electric cap C | States retained | E₀ | E₁ | Estimated gap E₁−E₀ |
|---|---:|---:|---:|---:|---:|
| 1 | 18 | 2,920 | 5.506498182 | 8.457617419 | 2.951119238 |
| 2 | 18 | 201,248 | 10.094496597 | 13.005367386 | 2.910870790 |
| 3 | 16 | 771,425 | 14.682735497 | 17.584279317 | 2.901543820 |

The last decimal places report computed matrix eigenvalues; they are not error bars for the untruncated theory. In particular, the largest three-cube basis is still changing noticeably under refinement.

### Cutoff sensitivity

| Joined cubes | Last cutoff increase | Gap before → after | Absolute change |
|---|---|---|---:|
| 1 | 16 → 18 | 2.951370356 → 2.951119238 | 0.000251119 |
| 2 | 16 → 18 | 2.912483772 → 2.910870790 | 0.001612983 |
| 3 | 14 → 16 | 2.917872080 → 2.901543820 | 0.016328260 |

Each retained E₀ and E₁ decreases under the nested basis enlargement, as required by the variational principle. The difference of two variational upper bounds is **neither a rigorous upper nor a rigorous lower bound on the true gap**. The changes above diagnose remaining truncation effects; they are not certified error intervals.

For comparison, the earlier one-cube calculation at energy cap 30 gave a numerical gap 2.9510428345, and a separate exact-arithmetic comparison certified the untruncated one-cube interval [2.941733,2.958693]. Those single-cube certificates do not automatically apply to joined cubes.

## 3. Smaller magnetic coefficient

Here κ=1 and ν=0.25, with the same complete bases.

| Joined cubes | Electric cap C | States retained | E₀ | E₁ | Estimated gap E₁−E₀ |
|---|---:|---:|---:|---:|---:|
| 1 | 18 | 2,920 | 1.468778128 | 4.469725755 | 3.000947628 |
| 2 | 18 | 201,248 | 2.692757106 | 5.691366705 | 2.998609599 |
| 3 | 16 | 771,425 | 3.916736084 | 6.914576886 | 2.997840802 |

At ν=0 the electric gap is exactly 3κ for every chain length. A nonzero gauge-invariant spin configuration contains a physical loop; the shortest loop has four links, each costing at least 3κ/4. A fundamental square loop attains 3κ. The number of independent first electric excitations is 5N+1. Thus the interacting estimates near 3 have a concrete electric-spectrum reference. Proximity to 3 or π alone does not identify an exact interacting gap.

## 4. Validation of the shared-vertex calculation

The independent validation includes:

- Counting all spin-½ physical states through binary cycle configurations and invariant multiplicities: 32, 868 and 25,676 states for one, two and three cubes. For two cubes there are 512 binary cycle patterns, but 868 states after all coupling channels are included.
- Constructing the four-spin-½ invariant tensors directly and checking their two-dimensional singlet space.
- Computing Wilson matrix elements independently from angular-momentum eigenvectors and Clebsch–Gordan multiplication, without the Wigner 6j formulas used in the main sparse construction.
- Comparing the entire 868-state two-cube spectrum under three different resolutions of the same shared vertices, including one with odd auxiliary cycles.
- Checking the complete two-cube energy-cap 8 basis of 622 states, which includes higher physical spins, against the independent tensor calculation.
- Checking the three-cube energy-cap 6 basis of 307 states against the same independent route, covering both interior shared planes.
- Recovering the earlier 398-state one-cube energy-cap 12 matrices up to a single common diagonal change of basis signs.
- Checking 3,368 local tensor contractions in exact rational arithmetic through doubled spin 14, covering all link and auxiliary labels allowed by the reported energy caps.

A missing convention-dependent sign was caught while changing the vertex resolution. A stronger test exposed an additional orientation issue for odd auxiliary cycles. The final implementation includes exact local tensor signs and oriented-edge signs; all three resolutions pass the independent tensor comparison. Its default matrices differ from the earlier validated default only by a common change of basis signs, which preserves the spectrum. The full volume study was rerun with this final implementation. The validation artifact records the convention checks and their scope; these are numerical matrix checks, not a formal proof of every floating-point output.

Every study matrix passed a Hermiticity check, every reported lowest pair passed direct eigenvector residual checks, and the largest retained matrix at each volume and coupling was rerun with a second random starting vector.

The largest direct retained-matrix residual was 6.56e-11; the largest repeated-run energy difference was 5.68e-14. A small residual measures how well the vectors solve the retained matrix problem. It does not by itself exclude an unobserved lower state or control the entire omitted Hilbert space.


## 5. What can be proved beyond the finite matrices

### Positive gap for each finite graph

For finite N, the full link configuration space SU(2)^(8N+4) is compact and connected. The electric Laplacian has discrete spectrum with finite multiplicities; the bounded real plaquette potential preserves that property. Imaginary-time evolution turns every nonzero, nonnegative wavefunction into a strictly positive one. This positivity-improving property makes the lowest eigenstate unique and strictly positive. Gauge transformations must fix that state, so it is a physical gauge-invariant vacuum. Its next physical eigenvalue is separated by a positive gap.

This applies standard spectral facts to the present finite graph; it is not a new solution of the continuum problem. The positivity theorem used in the argument is stated in [Güneysu, Theorem 2.9 and Corollary 2.10](https://arxiv.org/pdf/1109.0151). The full graph-specific proof and assumptions are in the accompanying volume audit.

### The quantifier that remains open

\[
\underbrace{\Delta_N>0\quad\text{for each finite }N}_{\text{established for this model}}
\qquad\not\Longrightarrow\qquad
\underbrace{\exists\delta>0:\ \Delta_N\geq\delta\quad\text{for all }N}_{\text{not established here}}.
\]

For example, 1/N is positive for every finite N and tends to zero. Three positive numerical gaps cannot settle the infinite sequence.

A simple bound is Δ_N≥3κ−ν(5N+1). It becomes ineffective as N grows. This failure is a limitation of that estimate, not evidence that the physical gap closes.

The full vacuum energy itself satisfies

\[
E_0(H_N)\geq\frac{N+1}{2}\min(\nu,3\kappa/4).
\]

Consequently, a fixed bare-electric cutoff eventually cannot place the omitted-state energy threshold above the vacuum. Subtracting a scalar vacuum energy shifts both quantities equally and does not repair that comparison. The volume audit gives a possible stronger local-energy comparison, with matching changes to the whole omitted block; its numerical performance is not claimed here.

## 6. What this means for the research argument

The calculation supports continuing the SU(2) Hamiltonian route. Joining cubes does not immediately eliminate the gap. It also shows why keeping all shared-vertex channels, checking representation invariance and increasing the cutoff are necessary.

The next mathematical target is a quantitative lower bound that stays positive as all three spatial extents grow. A chain of fixed 1×1 cross-section tests growth along one direction only. Widening the lattice introduces further coupling channels and must be treated explicitly. A subsequent continuum construction must shrink the lattice spacing and adjust the bare coupling while controlling a physical energy scale.

The official problem asks for the required nontrivial quantum Yang–Mills theory on four-dimensional spacetime and a positive mass gap. Qualitative finite-volume positivity supplies only part of the much earlier lattice setup. [Jaffe and Witten, official problem statement](https://www.claymath.org/wp-content/uploads/2022/06/yangmills.pdf)

The separate note, **Space, mathematical reality, and what a Yang–Mills proof would establish**, answers the proposals about mathematical ontology, Level IV and Gödel. Those interpretations can motivate questions, but they supply no additional Hamiltonian term or gap estimate until expressed as a precise mathematical construction.
