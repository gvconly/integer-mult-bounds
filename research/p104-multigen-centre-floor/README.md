# The P104 centre with multi-generation descents: floor 56 except in multi-shallow profiles (no κ claim)

The P104 word pays **56** per FULL for its centre: eight rank-7 owned descents. This is the only place where port data can leave its starting frame. Single-generation decoder floors (Theorem 4 of [ideal-decoder-floors](../ideal-decoder-floors/README.md), Theorem P of [decoder-floors-frontier-format](../decoder-floors-frontier-format/README.md)) show that the 56 cannot be lowered when every owned descent carries first-generation material.

The word model also allows **multi-generation** descents. An owned register may copy content that earlier owned descents landed, and descend it further. Model C does not price this: two contents landed at a common line can share their next descent.

This package proves lower bounds that include multi-generation, for root-code centre contents (the column-factor half of Theorem 4's hypothesis). Row factors and the number of lifetimes are free.

## Results

| Profile of the centre | Floor | Status |
|---|---:|---|
| Contents landed at lines and 0 span Q^8 (R_8 = 8) | **56** | Proven (Theorem A) |
| Exactly one direction never reaches a line (R_8 = 7) | **56** | Proven (Theorem B) |
| Two such directions, both reaching planes (R_8 = 6, R_7 = 8) | **56** | Proven (Theorem C) |
| ≥ 3 such directions, or two with one ending above the planes | — | **Open**: reduced to an explicit shallow covering problem |

So a root-code P104 centre below 56 needs **at least two independent content directions that never reach a line**. If there are only two, at least one must stop at a frame of dimension ≥ 3.

## How it is proved

1. **Entry lemma.** Rows inside a port's world can only come from owned landings there. Every port is "good": its orthogonal live roots span r_p^⊥. Hence each port p needs r_p in the span of contents landed in its world.
2. **Sources.** Root-code first-generation content exists only at the full frame and at the 16 removed hyperplanes.
3. **Feed inequality, by levels.** Each unit of descent into a level can add at most one dimension to the content span below it. Content that does not propagate further is charged as fresh dimension. Summing gives the master inequality:
   cost ≥ 7·R_8 + Σ_{d=3..7}(R_d − R_8) + R_9 + (shallow fresh count).
   The selected word meets it with equality.
4. **Shallow coverage.** A direction that never reaches a line must be delivered to its ports, at least 49 of them for any direction, from larger frames.
   - A frame at level d holds at most m_d = 1, 2, 4, 8, 12, 20, 36, 60 ports (exhaustive over all subspaces of F2^9).
   - With one shallow direction this forces the missing cost back, at every level where it can stop.
   - With two directions on planes, a two-plane cover forces three removed roots into a common 2-plane. The two removed frames contain no A2, so that cannot happen.

All finite facts are recomputed exactly from the published word (integers and Fractions):
- **Realisation.** The integral E8 realisation: the eight centre contents form a Z-basis.
- **Lemmas.** The capacities, the |T_z| ≥ 49 census (13,320 lattice vectors and 3,360 A2 continua), the two-plane lemma (4,016 candidate sets against 43,435 planes), and the A2 check.
- **Controls.** All three designs below satisfy all 104 ports, so 56 is attained in very different profiles.
  - The selected word is tight in the master inequality.
  - A genuinely multi-generation centre of cost 56 (four coupled pairs) has R_8 = 8, as Theorem A requires.
  - A centre of cost 56 in which **no content ever reaches a line** (R_8 = 0) also exists. Each content lands on a plane inside its copy frame (rank 6), then rotates through a 3-dimensional frame onto a second plane meeting the first in its frame-pair line (rank 1).

The full argument is in [PROOF.md](PROOF.md).

## Scope

- **Covered:** words with the P104 mains, root-code owned contents, any number of owned lifetimes, any per-port combinations, and multi-generation descents.
- **Not covered:**
  - non-root-code contents (rank-9 realisations; the open decoder cases (b)–(d) of the frontier-format package);
  - the residual multi-shallow profiles;
  - r^t moment effects: the floor is on mass, and multi-generation designs at 56 use more, smaller jobs, so their moment is slightly higher.
- **No supplier, κ or selected-result change** is claimed.

## Reproduce

Python 3.12, standard library only, from the repository root. About 90 s.

```sh
python3.12 -I -B research/p104-multigen-centre-floor/verify.py
python3.12 -B scripts/research_package.py --package research/p104-multigen-centre-floor
```

The verifier:
- reads only `research/p104-productive-controller/WORD.json.gz`;
- recomputes every number;
- compares the result with RESULT.json;
- refuses to run under `-O`, because its checks are assertions.

## Open next questions

**The shallow covering problem** (PROOF.md §4). Let Z = (V'_8)^⊥, with dim Z = s. Choose pairs (β_j ∈ Z, frame F_j of level ≤ 7) such that proj_Z(r_p) ∈ span{β_j : p ∈ world(F_j)} for every port.
- Show the number of pairs is ≥ Σ_i (9 − max(δ_i, 2)) over the shallow directions i.
- At every δ_i = 7 this means at least 2s pairs. The all-planes control is a solution with s = 8 and exactly 16 pairs, so the question there is precisely whether **15 plane-level deliveries** can serve all ports.
- A proof closes the root-code centre at 56. A counterexample, if it survives the availability of its contents, is a sub-56 centre.
- The Singleton-type argument gives only s + 1. A plane-level analogue of Theorem 4 seems to be needed.

## Credit

- The P104 word, controller and invoice: research/p104-productive-controller (Ali Romman) and its lineage (PR2305/2152/2219; David Arias, Cinna Davis, Gareth Davidson).
- The frame-pair decoder and the single-generation 56 floors: #968 (themilkman216), #1180, #1212, #2228 and the frontier-format and ideal-decoder-floors packages (#2338, #2639).
- The independent-tree centre bound: ctoth (Discussion 1516) and its audit package.
- New here:
  - the multi-generation model and the entry and source lemmas;
  - the level feed inequality and the master inequality;
  - the capacity and shallow-set censuses, the two-plane lemma, and Theorems A–C;
  - the explicit multi-generation 56 control.

Prepared by G. V. Chudnovsky with Anthropic Claude assistance. Apache-2.0.
