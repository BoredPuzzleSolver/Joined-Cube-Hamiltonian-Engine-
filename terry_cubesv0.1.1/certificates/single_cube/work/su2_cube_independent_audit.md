# Independent SU(2) cube Wilson-matrix audit

**Result: all six Wilson matrices at link-spin cutoff J = 1/2 pass an exact, independent Haar-contraction check.** The signed recoupling matrices agree after one common change of state-basis signs. No absolute-value replacement of recoupling coefficients is needed or justified.

## Independent construction

The graph has eight vertices, twelve edges, and six square faces, with Gauss law at every vertex and no external charges. At this cutoff each vertex has zero or two occupied fundamental edges. The 32 allowed edge masks therefore describe disjoint closed loops. For a configuration C, take its wavefunction to be the product of fundamental holonomy traces around those loops. Each such wavefunction has Haar norm one; different masks are orthogonal because some edge then carries a single fundamental factor.

For `W_f = Tr_fund(U_f)/2`, a nonzero matrix element toggles precisely the four face edges. Every used edge in the resulting integral carries two fundamental factors. The independent calculation uses only

`integral U_ab conjugate(U_cd) dU = delta_ac delta_bd / 2`,

`integral U_ab U_cd dU = epsilon_ac epsilon_bd / 2`,

with `epsilon_01 = 1`, `epsilon_10 = -1`. It retains epsilon signs and counts all consistent binary index assignments using parity constraints. It calls no 6j or Clebsch-Gordan routine to construct these integrals.

The resulting 192 directed nonzero entries are exactly 1/2, 1/4, or 1/8 in this loop-character basis. Their magnitudes are `2^(-1-s/2)`, where s is the number of occupied spectator edges at the four face vertices. Here s is 0, 2, or 4. All six matrices are exactly symmetric.

## Recoupling comparison

Only after constructing the independent matrices, the checker loads `su2_cube_model.py` and compares its signed rational squared coefficients. One diagonal matrix D, with entries plus or minus one, satisfies

`W_recoupling,f = D W_loop,f D`

for every face simultaneously. Twelve of the 32 states need a minus sign. The complete phase assignment is saved with the results. Every entry agrees exactly as a Fraction; this comparison uses no floating-point tolerance.

A separate exhaustive enumeration of all `3^12` assignments with doubled spins 0, 1, or 2 gives **1,013 physical states at J = 1**, agreeing with the recoupling module's backtracking count.

## Reproduction and scope

- Script: `work/su2_cube_independent_check.py`.
- Exact matrices, transitions, and phase assignment: `outputs/SU2_Cube_Independent_Validation.json`.
- The optional comparison needs sibling `su2_cube_model.py` and its `su2_recoupling.py` dependency.

This audit independently checks the J = 1/2 retained Wilson blocks and counts the J = 1 basis. It does not independently construct Wilson transitions involving spin-one links, prove an untruncated cube gap, or address increasing lattice volume and the continuum limit.
