# Unit screen: E8 and the paper's triples are one series, and a bound on every compatible unit (no κ claim)

Prepared by G. V. Chudnovsky with Anthropic Claude assistance. Apache-2.0. Base: main `32890bc5d`.

**No κ claim, no word, no supplier.** This note scores line configurations ("units") on the repository's own unit
yardstick. That yardstick is the certified B_K-format ceiling `skel.X` of `research/format-limits-corrected`, with the
window, label-rank and copy-cost conventions of `research/label-field-census` and
`research/e8-decoder-floors-and-cold14-notes`. The ceiling belongs to the older B_K format. It ranks units against
each other; it does not predict κ under the current productive-helper interface.

## Findings

1. **E8 and the upstream triples belong to one series.** Under A8, the E8 roots are the A8 roots plus the 84 triples
   of a 9-element set. In general, let T_n be the A_{n−1} roots together with the paper's triples, viewed as norm-2
   vectors of Z^{n,1}. These are the roots of levels 0 and 1 of the E_n lattice. T_9 is E8, T_10 lies in E10
   (II_{9,1}), and larger n is Lorentzian. The paper uses only the level-1 part.
2. **Every unit compatible with the current label rule has at most C(h+2,3) lines at rank h,** in any signature.
   E8 attains this exactly (it is the tight spherical 7-design), and Leech attains its own quintic analogue exactly.
   T_n fills (n−1)/(n+2) of the bound. For units whose label rank is at least h, as for E8, the triples and the series,
   no compatible unit can exceed about 1.07·10⁻² on the yardstick. That is about twice E8's ceiling.
3. **The series beats E8 only at the shortest schedule.** At s = 3, T_18–T_22 reach 8.1–8.2·10⁻³ with window h, and
   7.7–7.8·10⁻³ with the complex window c. E8's best is 6.0·10⁻³ (bit, P111) and 5.19·10⁻³ (complex, P112). At
   s = 9, where today's words run, E8 is better: 4.61·10⁻³ against 3.08·10⁻³ for T_20.
4. **Helper cover numbers do not separate them.** On the cover measure of format-limits §D, every unit here needs
   2–3 hyperplane sections per side. The earlier objection to the triples was the raw count of orthogonal neighbours
   (630 against 63); it does not show up in this measure.
5. **The special lattices in dimensions 16 and 24 do not work as units.** Leech fails the label rule in every integral
   scaling and has at least 124 pairwise non-orthogonal lines. BW16 and the Golay octads are far below their bounds.
   Incidence units from designs or curves embed in the series or have too few inner-product values (§4).

## The yardstick and its calibration (`calib.py`)

A unit is v lines with an integer Gram matrix G. The quantities are:

- **window h:** rank_Q G;
- **label rank c:** rank_F2 of B = (J + G) mod 2, which needs a unit diagonal. This is the complex window (9 for E8);
- **copy cost x:** the exact cost per port of the coordinate-total decoder, Σ_f rank_F2 B[S_f, S_f] / v, where S_f
  is the set of lines not orthogonal to frame vector f. The repository's proxy c·h/v is printed alongside;
- **ceiling:** X(s, h, 1, x), taken at its best odd s from 3 to 13.

The rule B = (J + G) mod 2 applies when the norm is 2 and every off-diagonal entry is 0 or ±1. B must vanish where
the decoder acts (non-orthogonal pairs) and be 1 on the diagonal.

| check | screener | published |
|---|---|---|
| bit P111 (x = 63/111, h = 8) | 6.0005·10⁻³ at s = 5 | 6.0005·10⁻³ (format-limits A.3) |
| complex P112 (x = 64/112, h = 9) | 5.1863·10⁻³ at s = 5 | 5.1863·10⁻³ |
| E8, Hamming frame clique bound 8/15 | 8.0093·10⁻³ at s = 3 | 8.009·10⁻³ (label-field census) |
| E7 (x = 7/9) | 1.8468·10⁻³ | 1.847·10⁻³ |
| E8 decoder cost | 72 on 120 lines, c = 9, 63 orthogonal per line | 72 (#446 Theorem C value) |
| paper triples, h = 24 | 7.388·10⁻³, 630 orthogonal per line | ≈ 7.4·10⁻³, 630 (cold14 notes, float screen) |

## 1. One series: E8, E10 and the triples (`en_series.py`)

Let T_n be the following norm-2 vectors of Z^{n,1}, with form x·y − m m′:

- the A_{n−1} roots (e_i − e_j; 0);
- the paper's triples (t_S; 1).

All of them are orthogonal to ρ = (1^n; 3), and the pairwise inner products are 0 or ±1.

- **n = 9.** ρ is all-odd (so ρ^⊥ is even), primitive and isotropic, and ρ^⊥/ρ is even unimodular and positive
  definite of rank 8, that is, E8. The Gram matrix of T_9 is positive semidefinite of rank 8, and its 120 lines are
  distinct. They lie among the roots of a rank-8 root lattice with at least 240 roots, which forces E8: the 36 + 84
  lines are all of E8. This is the decomposition 240 = 72 + 84 + 84 of E8 under A8, and it gives the "Lorentzian E8
  chart" already recorded in `research/c5-dyadic-decoder-gate`.
- **n = 10.** ρ has norm 1, so ρ^⊥ is even unimodular of signature (9,1): the E10 lattice II_{9,1}. T_10 is exactly its
  level-0 and level-1 roots.
- **n ≤ 12.** No root of height ≤ 2 in Z^{n,1} is compatible with all of T_n, so T_n is maximal there. An unseeded
  clique search (`lorentz.py`, about five minutes at n = 10 over 17,370 roots of heights ≤ 2) reached 165 lines, the
  size of T_10, within a second and found nothing larger.

In the table, "orthogonal" is the largest number of lines orthogonal to one line, and the covers are bit/complex. The
E8 row takes x and the ceilings from the Hamming frame. In the Z^{9,1} chart the coordinate-total decoder costs
x = 0.750 and gives 1.985·10⁻³ (window h) and 1.743·10⁻³ (window c); c, the orthogonal count and the covers are the
same in both charts.

| n | lines v | rank h | c | fill of C(h+2,3) | x | ceiling, window h | ceiling, window c | orthogonal | covers |
|---|---|---|---|---|---|---|---|---|---|
| 9 (E8) | 120 | 8 | 9 | 1.000 | 0.600 | 5.160·10⁻³ | 4.536·10⁻³ | 63 | 3/2 |
| 10 (E10) | 165 | 10 | 11 | 0.750 | 0.667 | 2.695·10⁻³ | 2.432·10⁻³ | 92 | 3/2 |
| 12 | 286 | 12 | 13 | 0.786 | 0.545 | 4.760·10⁻³ | 4.382·10⁻³ | 175 | 3/2 |
| 16 | 680 | 16 | 17 | 0.833 | 0.400 | 7.800·10⁻³ | 7.328·10⁻³ | 469 | 3/2 |
| 18 | 969 | 18 | 19 | 0.850 | 0.353 | 8.136·10⁻³ | 7.697·10⁻³ | 696 | 3/2 |
| 20 | 1330 | 20 | 21 | 0.864 | 0.316 | **8.171·10⁻³** | **7.773·10⁻³** | 987 | 3/3 |
| 24 | 2300 | 24 | 25 | 0.885 | 0.261 | 7.849·10⁻³ | 7.529·10⁻³ | 1793 | 3/3 |

For comparison, the paper's triples alone at h = 22 (1540 lines, x = 0.314) reach 7.444·10⁻³. Their covers are 2/3
(`covers.py`, h = 16, 20, 24). E8's complex cover is 2 with the Z^{9,1} candidates, as format-limits §D finds for
P112; the Hamming-chart candidates in `covers.py` give 3.

**The gain sits at s = 3 only** (certified upper ends, window h, except the complex rows):

| unit | s = 3 | s = 5 | s = 7 | s = 9 | s = 11 | s = 13 |
|---|---|---|---|---|---|---|
| E8 bit P111 | 5.947·10⁻³ | **6.000·10⁻³** | 5.209·10⁻³ | 4.607·10⁻³ | 4.149·10⁻³ | 3.788·10⁻³ |
| E8 complex P112 (window 9) | 5.047·10⁻³ | **5.186·10⁻³** | 4.513·10⁻³ | 3.995·10⁻³ | 3.599·10⁻³ | 3.286·10⁻³ |
| triples h = 22 | **7.444·10⁻³** | 4.293·10⁻³ | 3.324·10⁻³ | 2.796·10⁻³ | 2.448·10⁻³ | 2.195·10⁻³ |
| T_20 | **8.171·10⁻³** | 4.724·10⁻³ | 3.660·10⁻³ | 3.080·10⁻³ | 2.697·10⁻³ | 2.419·10⁻³ |
| T_20 complex (window 21) | **7.773·10⁻³** | 4.491·10⁻³ | 3.479·10⁻³ | 2.927·10⁻³ | 2.563·10⁻³ | 2.298·10⁻³ |

A wide window pays idle mass that grows with s·(h − 1). The low copy cost of the series wins only when the schedule is
short. The frontier words run at s = 9–11, where format-limits §A.5 attributes 87–93% of the gap below the ceiling to
helpers. So the open question for the series is whether a short-schedule word (s = 3) can carry its helpers. Its
deficit per port at s = 3 is only 2 − 3x ≈ 1.05, against 8 − 9·(63/111) ≈ 2.89 for E8 at s = 9.

## 2. A bound for every compatible unit (`bound.py`)

**Claim.** Let G be a real symmetric matrix of rank h with diagonal 2 and off-diagonal entries in {−1, 0, 1}, of any
signature. Then v ≤ C(h+2, 3).

**Proof.** Write G = R Q R^T with R a real v×h matrix and Q = diag(±1): take the spectral decomposition and keep the h
nonzero eigenvalues. Put f(t) = t³ − t, which vanishes on {−1, 0, 1} and has f(2) = 6. Since r_iᵀ Q r_i = 2 for every
i, the matrix F_ij = G_ij³ − G_ij G_ii G_jj / 4 equals f(G_ij) entrywise. So F = 6I, which has rank v. Each entry of F
is a form of degree 3 in r_i and of degree 3 in r_j, so rank F ≤ C(h+2, 3).

The same homogenization works for any odd (or any even) polynomial vanishing on the off-diagonal values of a
constant-norm configuration. A polynomial with terms of both parities adds C(h+s−2, s−1) to C(h+s−1, s). This gives
the bounds quoted below for BW16, the Golay octads and Leech.

How full is the bound?

| unit | fill |
|---|---|
| E8 (120 = C(10,3)) | 1.000 (tight spherical 7-design) |
| Leech, 98,280 lines against C(28,5) = 98,280 | 1.000 (tight 11-design) |
| T_n | (n − 1)/(n + 2): 0.750 at n = 10, rising to 0.885 at n = 24 |
| paper triples | 0.545 (h = 10) → 0.778 (h = 24) |
| D_n, E8+E8 | ≤ 0.47 |

E8 and Leech are the tight objects of their respective classes, which is the precise sense in which they are
special.

**What the bound allows.** Take a unit exactly at the bound with label rank c = h. Its ceiling is at most
1.071·10⁻² at h = 13 (the maximum), 1.034·10⁻² at h = 16 and 8.61·10⁻³ at h = 24. This is a proxy cap, not a theorem
about copy floors. It assumes c ≥ h, which holds for E8, the triples and the series but not for every unit (D_n has
c = h − 1, and the Steiner systems less). It also uses the coordinate-total decoder cost c·h.

**Concrete targets.** To beat T_20's 8.17·10⁻³ (at s = 3, with c = h + 1), a unit needs:

| rank h | lines needed | bound C(h+2,3) | fill needed | T_h has |
|---|---|---|---|---|
| 10 | none suffices | 220 | impossible: even 220 lines give 7.92·10⁻³ | 165 |
| 12 | ≥ 340 | 364 | 0.934 | 286 |
| 14 | ≥ 497 | 560 | 0.887 | 455 |
| 16 | ≥ 703 | 816 | 0.862 | 680 |

Something "E8-like" at rank 12–16 would mean norm-2 vectors with inner products 0 or ±1 filling 0.86–0.93 of
C(h+2,3), in Lorentzian signature, since positive-definite such systems are root systems. That is the precise opening
for a new algebraic construction.

## 3. Dimensions 16 and 24 (`families.py`, `leech.py`)

In this table, "non-orthogonal set" is the largest set of pairwise non-orthogonal lines found greedily, a lower bound.

| configuration | lines | rank | inner products | fill of its bound | orthogonal per line | non-orthogonal set | verdict |
|---|---|---|---|---|---|---|---|
| D16 | 240 | 16 | 0, ±1 | 0.294 | 183 | ≥ 15 | x = 1, no saving |
| E8 + E8 | 240 | 16 | 0, ±1 | 0.294 | 183 | ≥ 8 | x = 17/15, no saving |
| Barnes–Wall BW16 (Gram/2, norm 4) | 2160 | 16 | 0, ±1, ±2 | 0.139 | 855 | ≥ 41 (72 by #895) | fails the label rule (±2); complex window ≥ 72 |
| Golay octads, (G − 4J) | 759 | 24 | 0, −2, −4 | 0.262 | 280 | ≥ 9 | dyadic decoder only; 759 lines at rank 24 against 2300 for T_24 |
| Leech (norm 32 scaling) | 98,280 | 24 | 0, ±8, ±16 | 1.000 (C(28,5)) | 46,575 | ≥ 124 | fails the label rule in every integral scaling; complex window ≥ 124; bit cover ≤ 9 sections (greedy, 3000 sampled hyperplanes) against E8's 3 |

On the label rule for Leech: at G/8 the entries ±2 are even and nonzero, so B = 1 on pairs where it must vanish. Odd
multiples of that scaling keep those entries even, and even multiples make the ±1 entries even.

Leech is E8's analogue in being tight, but every unit-level obstruction above applies to it. The 24-element object
that does work on this yardstick is T_24, the triples of a 24-set together with A23.

## 4. Set systems and curves

For 0/1 vectors (incidence vectors of k-sets) under a form I − λJ, norm 2 forces λ = (k − 2)/k². Inner products
0 or ±1 then force |S ∩ T| ∈ {k − 3, k − 2, k − 1}. The families of this kind found here embed isometrically in the
series. For example, the 4-sets {S : |S ∩ D| ≥ 2} with |D| = 3 embed in T_n by sending D ∪ {x} to the A-roots.

Designs with only two Gram values satisfy v ≤ C(h+1, 2) + h by the same argument, so v/h² stays below about 1/2. On
the proxy x = c·h/v, a saving then needs a label rank c below v/h. The Steiner triple systems checked have x = 3.0
(PG(3,2)) and 4.2 (PG(4,2)). Unital blocks with q ≥ 4 and pairs of ovals have intersections outside the allowed range,
so they do not fit the norm-2 rule at all.

Curves could still enter through signed or integral constructions, the way E8 comes from the Hamming code, but not as
incidence units.

## Scope and caveats

- The ceilings are certified upper ends of `skel.X` in the B_K format. They do not include helper, hosting, compiler
  or analytic costs. Realised B_K words sat 2.5–2.7× below them (format-limits A.5). The current record uses a
  different interface.
- Copy costs are those of the coordinate-total decoder: exact for that decoder and an upper bound on the copy floor.
  A lower floor would only raise the corresponding ceiling.
- Cover numbers are taken over the stated candidate hyperplanes, so they are upper bounds. The value reported is the
  greedy count, or 2 when two candidate sections already cover every reader (`two_max` equal to the reader count,
  as for the paper's triples on the bit side). For E8, format-limits §D proves the exact bit value 3, and two
  candidates cover at most 52 of 63 here, matching it.
- Cliques are greedy, so they are lower bounds. Ranks are exact (Bareiss) for v ≤ 1400 and floating point above.
- The label rule used is B = (J + G) mod 2. A unit whose B is chosen freely on orthogonal pairs (a minrank question)
  may have smaller c. The Leech, BW16 and octad rows use only rule-independent quantities.
- Realising a new unit would need a native call, decoders, helpers and hosting for that window. None is attempted here.

## Reproduce

From the repository root, with numpy installed (the engine scripts are those of `research/format-limits-corrected`):

```sh
python3 -B research/unit-screen-en-series/calib.py       # calibration against published numbers   (~1 min)
python3 -B research/unit-screen-en-series/bound.py       # the bound, cap and targets               (~7 min)
python3 -B research/unit-screen-en-series/families.py    # ADE, triples, BW16, octads, STS          (~1 min)
python3 -B research/unit-screen-en-series/leech.py       # Leech                                    (< 1 min)
python3 -B research/unit-screen-en-series/covers.py      # cover numbers, E8 and triples            (< 1 min)
python3 -B research/unit-screen-en-series/en_series.py   # the series T_n, n = 9..24                (~5 min)
python3 -B research/unit-screen-en-series/lorentz.py 10 012 400   # optional unseeded clique search
```

Each `.log` file here is the output of the script of the same name. `families.json` and `en_series.json` hold the
same rows in machine-readable form. Only the full default runs refresh those two files: `en_series.py` with arguments
prints its rows without overwriting the receipt, and `lorentz.py` writes its result to the current directory.

## Credits

The ceiling engine, cover model and helper floors are from `research/format-limits-corrected` (DreamingOfClouds, with
#446). The label-field census is Rohan Guliani's. The ADE/extension screens and the triple side finding are in Joel
Pulikkan's `research/e8-decoder-floors-and-cold14-notes`. The Lorentzian E8 chart is in
`research/c5-dyadic-decoder-gate`, and Douglas Colkitt's subset-core screens are in `research/kappa-nine`. The BW16
72-clique is from DaysSky's #895, and the Leech/Turyn and Niemeier screens are Dominik Scholz's #598/#606. The triple
construction itself is the upstream paper's, after Alon.
