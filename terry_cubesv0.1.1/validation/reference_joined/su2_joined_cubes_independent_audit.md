# Independent audit: joined SU(2) cubes

## Result and scope

The complete physical spin-half basis on two joined cubes has **868 states**, including both invariant channels at every vertex carrying four half-spin links. Independent angular-momentum tensor contractions give the same full 868-level spectrum for three different auxiliary pairing trees, to less than `1e-12`. At `kappa = nu = 1`, the spin-half truncated gap is approximately **3.22259354077756**.

The final longitudinal, crossed and mixed model matrices agree with those independent tensor matrices, allowing a single common diagonal change of basis signs across every face. The higher-spin two-cube electric-energy cutoff `T <= 8` also passes this check on 622 states and 5,584 directed face entries. A three-cube `T <= 6` fixture passes on 307 states and 1,936 directed face entries, checking both interior planes. These are finite retained-matrix validations. They do not prove convergence under removal of a cutoff, an infinite-volume mass gap, or the continuum Yang–Mills problem.

Reproducer: `su2_joined_cubes_independent_check.py`, using the independent tensor module `su2_joined_cubes_tensor_probe.py`. Machine-readable results: `SU2_Joined_Cubes_Independent_Validation.json`.

## Independent physical basis count

The real chain graph has `4N+4` vertices, `8N+4` edges and `5N+1` physical faces. There are `4(N-1)` vertices of degree four. The independent enumeration first constructs a spanning tree and its fundamental cycles, then enumerates the entire cycle space over the two-element field. This does not use the model's spin-assignment recursion.

With physical spins restricted to zero and one half, Gauss law requires an even number of occupied edges at each vertex. At a degree-four vertex with four occupied edges, the singlet multiplicity is two. All other allowed local patterns have multiplicity one. Thus each physical cycle mask with `r` four-half vertices contributes `2^r` states.

| Cubes | Physical cycle masks | Full spin-half dimension |
|---|---:|---:|
| 1 | 32 | 32 |
| 2 | 512 | 868 |
| 3 | 8,192 | 25,676 |

For two cubes the mask counts by `r=0,1,2,4` are `324,128,56,4`, respectively. The weighted count is `324 + 2*128 + 4*56 + 16*4 = 868`. The complete independently generated two-cube state tuples agree with the model for all three supported pairing choices.

## Why the auxiliary channel must remain complete

Direct Pauli-matrix Casimir calculations in the 16-dimensional tensor product of four fundamental representations find a two-dimensional total-spin-zero space. The Casimir of either chosen pair has eigenvalues zero and two in that space, corresponding to intermediate spins `k=0,1`.

The absolute overlap matrix between two pairing bases is

```
1/2       sqrt(3)/2
sqrt(3)/2 1/2
```

Both channels contribute to contractions across an alternative pairing: the squared normalized contractions are `1/4` and `3/4`. A generic invariant operator gives the same two eigenvalues `6 ± sqrt(3)` in all three pairing bases. Removing the `k=1` channel is therefore a physical truncation, not a harmless convention choice.

Splitting one degree-four vertex into two trivalent vertices joined by a tree edge increases edges and vertices equally. The auxiliary group variable can be gauge-fixed on that tree, and it creates no new cycle. Its representation label resolves the original invariant-tensor multiplicity. It has zero physical electric energy, and the physical link spin cutoff must not be applied to it. The allowed auxiliary values are all common pair-coupling channels.

## Independent Wilson construction

The independent route uses no Wigner 6j symbols. It builds angular-momentum generators directly. Each normalized three-leg invariant is the one-dimensional zero eigenspace of the total Casimir. Clebsch–Gordan tensors for multiplying a representation by spin one half are obtained from the appropriate highest-weight Casimir eigenvector and repeated total-spin lowering. This preserves the relative phases between magnetic components.

Write each oriented edge matrix as `K_j(U)=D_j(U) epsilon_j`, so every vertex index transforms as an outgoing index. Tensor-product Clebsch–Gordan decomposition gives

```
K_j(m,n) K_h(a,b)
  = sum_(J,M,N) C(jm,ha|JM) C(jn,hb|JN) K_J(M,N).
```

Haar orthogonality and normalized spin-network edge factors then supply `sqrt(d_old/d_new)` on each visited edge. At a face vertex, contract the initial and final normalized invariant tensors, the two endpoint Clebsch–Gordan tensors and the fundamental antisymmetric tensor. With fundamental indices ordered as forward and backward, call this real contraction `c_v`. For a face traversing `r` edges opposite their fixed orientation, the normalized trace coefficient is

```
(1/2) (-1)^r product_face_edges sqrt(d_old/d_new) product_face_vertices c_v.
```

The code enumerates every retained target obtained by shifting each visited doubled spin by plus or minus one. It does not use a 6j zero test to determine support. The complete independent matrices, including expected zeros, are compared with the model. Their signs are allowed to differ only by one state-diagonal basis phase common to all face operators. Independently changing signs face by face would not pass this test.

## Pairing-tree and higher-spin checks

Three independently routed trees are checked on the same physical two-cube graph:

1. Pair both longitudinal links together.
2. Pair the negative longitudinal link with the transverse Y link.
3. Pair that longitudinal link with the first transverse edge in the independent edge ordering. This mixes Y and Z choices between vertices.

All three independent tensor Hamiltonians use zero auxiliary electric weight and have the same complete 868-level spectrum. The independent one-cube routing also reduces exactly to the original cube's face set and matrices before basis phases.

Additional retained-matrix checks:

| Fixture | Dimension | Directed face entries | Largest physical spin | Maximum magnitude error |
|---|---:|---:|---:|---:|
| Prior single cube, `T <= 12` | 398 | 4,260 | 3/2 | `6.2e-16` |
| Two cubes, `T <= 8` | 622 | 5,584 | 1 | `8.4e-16` |
| Three cubes, `T <= 6` | 307 | 1,936 | 1/2 | `8.9e-16` |

Every face in each fixture agrees under one common basis-sign change. The first check confirms that the previously used single-cube retained matrix at the strongest audited cutoff has the correct spectrum under the independent tensor convention. It checks retained entries; it is not a rerun of all omitted-sector certificate blocks.

## Phase-convention issue caught during validation

The first direct extension of a face-local trivalent formula gave different spectra for equivalent pairing trees. One missing factor came from expressing each local `(spectator, forward, backward)` 3j order in a fixed incident-edge order. An odd permutation multiplies a 3j symbol by `(-1)^(j1+j2+j3)`; the initial and final intertwiner phases must both be included. This property is documented in [NIST DLMF, section 34.3(ii)](https://dlmf.nist.gov/34.3.ii).

That correction alone made two uniform pairing choices agree, but a third mixed choice exposed an additional orientation-convention failure. The independent all-outgoing tensor construction resolves the orientation directly and supplies the regression reference. Consequently agreement of only two specially chosen resolutions is insufficient as a general validation argument. The local trivalent magnitude formula being audited is from [From square plaquettes to triamond lattices for SU(2) gauge theory](https://www.nature.com/articles/s42005-024-01697-4), equation 22; adapting local formulas to a new graph requires consistent global intertwiner and edge conventions.

The final model now uses orientation-safe exact signed contractions in `su2_oriented_vertex.py`; it does not rely on the incomplete permutation correction. Its 3j finite sums and half-spin Clebsch–Gordan formulas are distinct from the independent Casimir-eigenvector implementation. All three final model resolutions agree with the independent tensor matrices under a single common basis phase per resolution. The JSON and all three additional coefficient fixtures were regenerated successfully against this final implementation.
