# Volume audit for an open chain of joined SU(2) cubes

## 1. Physical graph and gauge space

An open \(N\times1\times1\) chain has

\[
|V_N|=4(N+1),\qquad |E_N|=8N+4,\qquad
n_p=5N+1.
\]

There are four degree-4 vertices on each interior transverse plane and eight degree-3 end vertices. The Hamiltonian is

\[
H_N=T_N+V_N,\qquad
T_N=\kappa\sum_{\text{physical }e}j_e(j_e+1),\qquad
V_N=\nu\sum_{p=1}^{n_p}(1-W_p),\quad
W_p=\operatorname{Tr}U_p/2.
\]

For \(\kappa>0,\nu\geq0\), \(0\leq V_N\leq2\nu n_p I\).

Splitting a degree-4 vertex into two trivalent vertices is a representation of its intertwiner space. The connecting auxiliary label is a coupling channel, with **zero electric energy**. It ranges over the intersection of the two pair-coupling intervals. At fixed physical spins, that interval is finite. With all physical spins zero, the only allowed auxiliary label is zero. Thus the representation introduces neither extra electric energy nor additional electric-vacuum states.

Retained physical-spin or energy cutoffs must include every allowed coupling channel. A cap on the physical link spins does not impose that same cap on the auxiliary spins. Statements about shortest cycles below refer to the original physical graph.

## 2. The electric gap remains exactly \(3\kappa\)

In any nonzero gauge-admissible spin configuration, the support of nonzero physical links has no degree-1 vertex: a single nontrivial SU(2) representation cannot couple to a singlet against trivial representations. A finite graph with minimum nonzero degree at least two contains a cycle.

The physical joined-cube graph has girth four. Every nonzero physical edge costs at least \((3/4)\kappa\). Therefore every non-vacuum electric state has energy at least \(3\kappa\).

Putting \(j=1/2\) around any elementary square and zero on all other links attains this value, with admissible intertwiners. The electric ground state is unique, and

\[
\Delta_N(\nu=0)=3\kappa
\]

for every chain length. Its first excited eigenspace has one basis state for each elementary square, hence dimension \(5N+1\).

For the interacting operator, min–max gives \(E_1(H_N)\geq3\kappa\). The constant electric-vacuum trial state has \(\langle W_p\rangle=0\), so \(E_0(H_N)\leq\nu(5N+1)\). Consequently

\[
\Delta_N\geq3\kappa-\nu(5N+1).
\]

This is a valid finite-volume lower estimate. Losing its positivity as \(N\) increases says that the estimate becomes ineffective; it does not establish gap closing.

### Every fixed finite chain has a positive interacting gap

There is also a structural statement valid for every finite \(N\), every \(\kappa>0\) and every \(\nu\geq0\), without a spin cutoff.

Before imposing Gauss law, the Hilbert space is \(L^2(M_N)\), where \(M_N=SU(2)^{|E_N|}\) is compact and connected. The electric Laplacian has compact resolvent; adding the bounded real plaquette potential preserves this. The ground energy is therefore an eigenvalue.

Feynman–Kac positivity improvement implies that this eigenvalue is simple and its normalized eigenfunction can be chosen strictly positive. [Güneysu, On generalized Schrödinger semigroups, Theorem 2.9 and Corollary 2.10](https://arxiv.org/pdf/1109.0151).

Gauge transformations commute with \(H_N\) and act as positivity-preserving, measure-preserving pullbacks. They must fix that unique normalized positive ground state. It belongs to the closed gauge-invariant subspace, whose restricted Hamiltonian retains compact resolvent. The physical space is infinite dimensional for these graphs with cycles, so its next distinct eigenvalue exists and is separated from the vacuum:

\[
\Delta_N(\kappa,\nu)>0
\quad\text{for every fixed finite }N,\ \kappa>0,\ \nu\geq0.
\]

This does not bound \(\inf_N\Delta_N\) away from zero. It explains why qualitative finite-volume positivity is already expected structurally; explicit certificates quantify the gap and test the implementation. The argument uses positivity in link-coordinate space, not a sign assumption on matrix entries in a chosen spin-network basis.

## 3. Absolute vacuum energy grows with volume

There is a direct lower bound, so the growth of the energy origin need not be inferred from numerical fits.

The \(N+1\) transverse \(yz\) plaquettes use disjoint physical links. Drop all other positive electric and magnetic terms:

\[
H_N\geq \sum_{r=0}^{N} h_{\square,r},\qquad
h_\square=\kappa\sum_{e\in\square}j_e(j_e+1)+\nu(1-W_\square).
\]

This operator inequality holds on the full link Hilbert space and therefore on its gauge-invariant subspace.

A useful bound for one square can be derived before imposing local Gauss law. Let \(g=3\kappa/4\), the first nonzero electric energy on the full four-link tensor-product space, and let \(b=\min(\nu,g)\). Reducing the magnetic coefficient from \(\nu\) to \(b\) lowers the operator. Split off its constant product vector \(\Omega\). The corresponding diagonal entry is \(b\), the norm of the off-diagonal block is \(b/2\), and the complementary diagonal block is at least \(gI\). These facts use

\[
\langle\Omega,W_\square\Omega\rangle=0,\qquad
\|W_\square\Omega\|=1/2,
\]

which follow from Haar orthogonality of the fundamental SU(2) character.

The smaller eigenvalue of the resulting two-by-two lower comparison is

\[
c(\kappa,\nu)=
\frac{b+g-\sqrt{(b-g)^2+b^2}}{2}
\geq \frac b2>0\qquad(\nu>0).
\]

Thus

\[
E_0(H_N)\geq(N+1)c(\kappa,\nu)
\geq \frac{N+1}{2}\min(\nu,3\kappa/4).
\]

At \(\kappa=\nu=1\), this gives \(E_0(H_N)\geq3(N+1)/8\). It is deliberately coarse, but proves linear growth. Together with the vacuum trial bound:

\[
\frac38(N+1)\leq E_0(H_N)\leq5N+1
\qquad(\kappa=\nu=1).
\]

This is a bound on the extensive vacuum energy, not on the excitation gap.

## 4. Complete energy caps and all omitted channels

Let \(P_C\) contain every physical state with \(T_N\leq C\), including all of its intertwiner channels. This is finite dimensional. A physical link has \(j(j+1)\leq C/\kappa\), and every auxiliary label is bounded by the incident physical labels.

There is stronger energy spacing than the generic quarter-integer Casimir spacing. At every physical vertex, the sum of doubled labels \(a_e=2j_e\) is even. Sum this parity condition over one color class of the bipartite physical graph: every physical edge occurs once, hence \(\sum_e a_e\) is even. Since \(a(a+2)\equiv a\pmod2\),

\[
T_N/\kappa\in\tfrac12\mathbb Z_{\geq0}.
\]

Therefore a valid threshold on the full omitted space is

\[
t_C=\frac{\kappa}{2}
\left(\left\lfloor\frac{2C}{\kappa}\right\rfloor+1\right).
\]

For a cap already on this grid, the omitted threshold is at least \(C+\kappa/2\). The weaker quarter-step bound is also safe.

Let \(Q_1\) be the complete set of omitted states reached from \(P_C\) by one plaquette operator, including changed auxiliary channels. Define \(R=I-P_C-Q_1\). The neighbor set is finite because multiplication by the fundamental loop has finite fusion alternatives. With

\[
A=P_CH_NP_C,\quad B=P_CH_NQ_1,\quad
L=\begin{pmatrix}A&B\\B^\dagger&T_{Q_1}\end{pmatrix},
\]

the same positive comparison used on the cube holds:

\[
H_N-(L\oplus T_R)=0_{P_C}\oplus QV_NQ\geq0.
\]

Thus, with ordered finite eigenvalues \(a_j,\ell_j\),

\[
\min(\ell_j,t_C)\leq E_j(H_N)\leq a_j.
\]

The rational matrix-enclosure and Schur-inertia procedure remains valid if its exact matrix-element and rounding hypotheses are checked for the degree-4 channels. In particular, symmetry, all reachable channels and the full remainder threshold belong to the certificate.

### Why a fixed cap eventually stops being useful

The comparison above cannot yield a positive gap bound once its fixed tail threshold lies below the growing true ground energy. Its lower level endpoint is at most \(t_C\), while its variational ground upper endpoint is at least \(E_0(H_N)\). The bound in section 3 proves that this happens eventually at any fixed \(C\) and fixed \(\nu>0\).

Merely subtracting a scalar vacuum-energy estimate does not repair the denominator:

\[
(t_C-c_N)-(z-c_N)=t_C-z.
\]

A cap growing with volume, or a different comparison organized around the interacting vacuum and local excitation content, is needed. This limitation concerns the method, not the actual gap.

### A stronger local-energy comparison to test

Let \(d\) be the largest number of plaquettes meeting a physical link: \(d=2\) for one cube and \(d=3\) for a chain of at least two cubes. Let \(e_\square(\lambda,\nu)\) be any certified lower bound on the full four-link square Hamiltonian with electric coefficient \(\lambda\). Allocating \((1-s)\kappa/d\) of electric energy to each incident face gives, for \(0\leq s\leq1\),

\[
H_N\geq sT_N+
n_p\,e_\square((1-s)\kappa/d,\nu)\,I.
\]

Each local face Hamiltonian is bounded below by its scalar ground-energy bound. Summing these operator inequalities is valid even though the faces overlap. No independence or commuting assumption is required.

Consequently a scalar Schur bound on the entire omitted block can use

\[
\sup_{0\leq s\leq1}
\left\{s\,t_C+
n_p\,e_\square((1-s)\kappa/d,\nu)\right\}.
\]

To retain the stronger finite-boundary construction, choose a fixed \(s\) and replace its whole omitted diagonal by \(sT_{Q_1}+c_N(s)I\), with the matching remainder threshold \(s t_C+c_N(s)\). Here \(c_N(s)=n_p e_\square((1-s)\kappa/d,\nu)\). Keeping the old \(T_{Q_1}\) auxiliary unchanged while raising only the remainder threshold would not follow from this inequality.

The elementary bound from section 3 supplies an immediate valid local estimate; sharper certified one-face ground energies may improve it. This proposal is an exact comparison framework. Its numerical usefulness for the current joined graphs has not been evaluated here, and the local lower bound must hold on the full four-link space, not merely on a selected local charge sector.

## 5. Why small certified patches cannot simply be joined

A uniform thermodynamic claim would require a number \(\delta>0\) such that \(\Delta_N\geq\delta\) for every sufficiently large \(N\). Certifying several positive finite-volume gaps does not establish that quantifier.

For independent systems, a tensor-sum Hamiltonian has gap equal to the smaller factor gap. Joined cubes share physical links, local electric terms and Gauss-law constraints; their interacting vacuum is not supplied by a tensor product of isolated cube vacua. No independent-system gap formula applies.

Published finite-size criteria such as Lemm–Mozgunov's require frustration-free hypotheses and quantitative finite-volume thresholds, with appropriate boundary information. A positive local gap alone is not the criterion. [Lemm and Mozgunov, Spectral gaps of frustration-free spin systems with boundary](https://arxiv.org/abs/1801.08915).

Our natural decomposition into positive electric and plaquette terms is not frustration-free at \(\nu>0\). A state annihilated by every electric term would be the constant product vector, which has positive magnetic expectation. Subtracting the global ground energy does not construct local positive terms with a common kernel. Applying a frustration-free criterion would require a separately justified transformation or extension.

## 6. A constructive strong-coupling route

Yarotsky proves uniform spectral-gap stability near a class of gapped product-state models. His setup permits infinite-dimensional local Hilbert spaces and unbounded classical terms, with sufficiently small local perturbations. The inspected theorem is stated for translation-invariant interactions and periodic finite volumes. [Yarotsky, Ground states in relatively bounded quantum perturbations of classical lattice systems, Theorem 1](https://arxiv.org/pdf/math-ph/0412040).

This is a plausible route to audit: use the unprojected link tensor product with electric product vacuum, group links into finite cells, and regard bounded plaquette multiplication as a local perturbation. Each perturbation has size controlled by \(\nu/\kappa\), independently of the total volume. Gauge invariance and a unique invariant ground state can then relate an established full-space gap to the physical sector.

This mapping is a proposal for checking the theorem's assumptions, not an application completed here. It still needs the precise cell normalization, boundary treatment and allowed smallness constant. In particular, the paper does not certify the current ratio \(\nu/\kappa=1\). Such a small-perturbation result would concern an existing strong-coupling regime at fixed lattice spacing, not establish a continuum Yang–Mills mass gap.

## 7. What the joined-cube study can establish next

- Exact Gauss-law and recoupling checks with all degree-4 channels.
- Individual energy certificates at increasing complete caps for each fixed chain.
- The behavior of the gap and excitation observables as the length increases, with error intervals kept separate from fitted trends.
- A clearly stated theorem linking local estimates to a uniform bound, if its assumptions can be proved.

Finally, \(N\times1\times1\) growth keeps the transverse dimensions fixed. Even a uniform gap along this sequence would first be a result for a strip with fixed cross-section. A three-dimensional thermodynamic limit requires increasing the transverse sizes as well, followed by the separate continuum analysis.
