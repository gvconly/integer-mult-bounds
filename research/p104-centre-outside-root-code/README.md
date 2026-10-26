# P104 centre contents outside the root code: a Collapse theorem and where a cheaper centre must live (no κ claim)

The P104 word pays **56** per FULL for its centre: eight rank-7 owned descents. [research/p104-multigen-centre-floor](../p104-multigen-centre-floor/README.md) (open PR, same author) bounds that centre when its contents are root-code vectors and descents may re-use landed contents. This package drops the root-code hypothesis: centre contents are **arbitrary vectors** on the 104 port coordinates.

## Results

| Statement | Status |
|---|---|
| **Lemma S.** Every delivered content pays, within its own ancestry of owned steps, at least the Model C price of that one delivery. Mixing is allowed. | Proven |
| **Collapse theorem.** A centre whose owned steps never *mix* (each step carries fresh first-generation content, or exactly one earlier landed content) costs at least the Model C cost of its collapsed decoder. Contents are arbitrary and re-descents are allowed. | Proven |
| **The unit.** The word's 104 live and 16 removed lines are isometric to #968's P104 = E8 − F0 − H4. | Checked |
| **Floors.** A mixing-free centre costs ≥ 56 if its collapsed decoder has ≤ 8 terms, or 9 terms with a rank-8 realisation, or a rank-8 realisation with every factor inside its column/row space. | Proven, from the cited Model C theorems |
| **One-sided, every m.** A Model C decoder with root-code columns and arbitrary rows (case (b)) costs ≥ 56 in the deep profiles (R_8 ≥ 7, or R_8 = 6 with R_7 = 8). By transposition the same holds for arbitrary contents with root-code coefficients. | Proven, via the predecessor's Theorems A–C |
| **Mixing gap.** A two-content mixing fragment pays 15. The best natural Model C decoders giving every port the same span pay 16, and term-by-term pricing pays 23. So the Collapse theorem genuinely needs mixing-freeness. | Checked |

So **a P104 centre of cost ≤ 55, whatever its contents, must have one of:**
- a step that **mixes** landed contents;
- at least **nine** pairwise non-proportional paid first-generation contents whose collapsed decoder is in an open Model C class: case (c) at m = 9 (rank-9 square), or a factor outside the column/row space at m ≥ 10.

By transposition, the one-sided version of "contents outside the root code" (arbitrary contents, root-code coefficients) is exactly Model C case (b). For P104 it is closed at m ≤ 9, and at every m in the deep profiles. It remains open only in the shallow profiles.

**Evidence** (exact, not proofs):
- **Truncation.** Truncating a frame-pair content to a 7-dimensional subframe saves one unit. All 1,016 such truncations leave all 52 users of that content without a clean vector, even with every combination re-solved.
- **Overcomplete removed-root directions.** With the 16 removed roots as directions and any minimal spanning set of size ≤ 4 per port, the frame-pair is a strict local optimum in that neighbourhood: all 524 single-port changes cost ≥ 57.
- **Mixing near-misses.** In the all-planes family (R_8 = 0), one mixed rotation e_j ± e_k can replace two pure rotations for ten direction pairs, except for exactly one sign-mismatched port each.
  - Four perfect matchings give centres of cost 52 that leave four ports unserved.
  - Standalone repairs cost at least one unit per port, so these designs come back to exactly 56.
  - An exhaustive C search (RECORDED) finds no mixed rotation with zero failures.
- **Frame counting.** No frame-counting proof of 56 exists. A cost-6 landing frame inside a removed hyperplane meets all 759 live frames, where counting would allow 650.

## How it is proved

1. **Model** (as in the predecessor). Owned MOVEs split into unit steps. A step's content is first-generation material at its start frame plus contents landed at subframes. World entry (Lemma W) needs no root code.
2. **Lemma S.** Adding the start frames of an ancestry from latest to earliest raises the dimension by at most one per step. So the ancestry size bounds the label span of the content, hence its Model C price.
3. **Union lemma.** A content re-descended alone through several steps pays at least dim F_root − dim(F_root ∩ all its landings). This covers every user it reaches.
4. **Collapse.** In a mixing-free centre the steps form trees, one per fresh content. Each tree becomes one Model C term whose rank is at most the tree's size.
5. **Consequences.** These come from #2338's Model C results for #968's unit: m ≤ 8 ⇒ ≥ 56 whatever the factors; Theorem P; the m = 9 trichotomy; and case (b) at m = 9 (research/nine-term-outside-colspace, claim C2).

The full argument is in [PROOF.md](PROOF.md).

## Scope

- **Covered:** words with the P104 mains, arbitrary centre contents, any number of owned lifetimes, and multi-generation re-descents, as long as no step mixes.
- **Not covered:**
  - mixing steps;
  - Model C case (c) at m = 9, and factors outside the column/row space at m ≥ 10;
  - the r^t moment (the bounds are on mass).
- y- and helper coordinates are dropped. This is a relaxation, so it does not weaken the lower bounds.
- **No supplier, κ or selected-result change** is claimed.

## Reproduce

Python 3.12, standard library only, from the repository root. About 90 s.

```sh
python3.12 -I -B research/p104-centre-outside-root-code/verify.py
python3.12 -B scripts/research_package.py --package research/p104-centre-outside-root-code
```

verify.py reads the word and #2338's `lib/units.py` with its data, and nothing else. It recomputes every number, compares the result with RESULT.json, and refuses to run under `-O`, because its checks are assertions.

## Open next questions

1. **Mixing.** Can an owned step that combines landed contents ever lower the centre? This is the only multi-generation mechanism Model C does not price. The pair-mix designs of §7 reach 52 + 4 unserved ports. First step: decide whether **chained** repairs, which re-use a repair's landed content, can serve those four ports for 3 units. That would give a 55 centre.
2. **Case (c) at m = 9** for #968's P104: a square decoder of a rank-9 realisation. First step: classify rank-9 realisations of the pattern with the frame-module theorem (#2338 §3.3) and price the dead sets of their dual bases.
3. **Case (b) in the shallow profiles** (R_8 ≤ 5, or R_8 = 6 with R_7 ≤ 7). This is the shallow covering problem of the predecessor, now with single-generation prices: every shallow delivery pays its full level. The frame-touch witnesses show that the label-span structure, not frame incidence, must carry the proof.
4. **Two-sided factors outside the column/row space at m ≥ 10** (case (d)).

## Credit

- **P104 word, controller and invoice:** research/p104-productive-controller (Ali Romman), with lineage PR2305/2152/2219 (David Arias, Cinna Davis, Gareth Davidson).
- **Model C, Theorem P, the m ≤ 8 theorem, the m = 9 trichotomy and the frame-module theorem:** decoder-floors-frontier-format (#2338) and its lineage (#1180, #1224, #1212, #2228, #2277).
- **Case (b) at m = 9:** research/nine-term-outside-colspace (claim C2, PROVEN from CHECKED there, by exact enumeration not re-run here).
- **The P104 unit and frame-pair decoder:** #968 (themilkman216).
- **New here:**
  - Lemma S, the Union lemma and the Collapse theorem;
  - the isometry check;
  - the mixing gap;
  - the truncation, overcomplete-direction and frame-touch evidence.

Prepared by G. V. Chudnovsky with Anthropic Claude assistance. Apache-2.0.
