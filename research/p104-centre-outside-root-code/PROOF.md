# Proof: P104 centre contents outside the root code

`verify.py` recomputes every finite fact below from the published word `research/p104-productive-controller/WORD.json.gz`, with integers, Fractions and F2 bit arithmetic. RESULT.json records the numbers.

## 1. Model

**Word model.** This is the P104 local word model of research/p104-productive-controller (`local.py`), as in research/p104-multigen-centre-floor:
- public registers only climb along nested frames, and owned registers only descend;
- every MOVE is paid |Δdim|, and ADD/COPY need a common frame;
- the mains are fixed: x_p goes from [q_p] to F2^9, and y_p goes from 0 to q_p^⊥, ending with e_{y_p} + e_{x_p}.

The **centre cost** is the total rank of the owned descents. The word pays 56: eight rank-7 lifetimes.

**Contents, with no root-code hypothesis.**
- The *content* of an owned register is its row restricted to the 104 port coordinates (the x-coordinates).
- y- and helper coordinates are dropped. They only add constraints, so this is a relaxation and lower bounds are unaffected.
- A content is an arbitrary vector of Q^104.
- Pub(F) is the set of vectors supported on {s : q_s ∈ F}: the combinations of x-registers that can sit at frame F.

**Unit steps** (as in research/p104-multigen-centre-floor).
- Replace each owned MOVE of rank r by r codimension-1 steps, with the content landing at every intermediate frame. Cost is unchanged and availability only grows.
- Step δ starts at F_δ, lands at G_δ (dim G_δ = dim F_δ − 1), and carries

  u_δ = v_δ + Σ_{δ' ∈ E(δ)} a_{δδ'} u_{δ'},  with v_δ ∈ Pub(F_δ) and E(δ) = {δ' earlier : G_{δ'} ⊆ F_δ}.

- v_δ is the *first-generation* part.
- The u_{δ'} are contents landed at subframes of F_δ: in public registers that climbed to F_δ, or the owned register's own content after its previous step.
- Copying a public y- or helper register is covered by this form. Its x-part is first-generation at F_δ plus contents landed at frames ⊆ F_δ.

**Lemma W (world entry, any contents).** For every port t there are coefficients β_{tδ}, on the steps δ with G_δ ⊆ q_t^⊥, such that w_t := Σ_δ β_{tδ} u_δ satisfies w_t(t) ≠ 0 and w_t(s) = 0 for every s ⊥ t.

*Proof.* These are steps 1–3 of Lemma 1 of research/p104-multigen-centre-floor, which do not use the root code.
1. Rows of registers inside the world of t come from registers starting there and from owned landings there.
2. Registers starting inside the world vanish on M_t = {t} ∪ {s : s ⊥ t}.
3. y_t ends with e_{x_t}. ∎

**Model C.** A decoder is N = Σ_j p_j c_jᵀ such that N vanishes on orthogonal pairs and N_tt ≠ 0, with entries on non-orthogonal pairs free (#2338, *decoder-floors-frontier-format*).
- The row factor p_j holds the per-port coefficients of content j. The column factor c_j is the content.
- The cost is Σ_j rank_F2 B[R_j, A_j], where R_j = supp p_j, A_j = supp c_j and B[t,s] = ⟨q_t, q_s⟩ is the dense table.
- Write U_X for the span of the labels of X. The exact rank formula (nine-term-outside-colspace §2) is

  rank B[R, A] = dim U_A − dim(U_A ∩ U_R^⊥).

  This is standard: B[R, A] is the matrix of the pairing U_R × U_A → F2, and its rank is dim U_A minus the dimension of the part of U_A that pairs to zero with U_R.
- verify.py cross-checks the formula against direct F2 elimination on every block it prices: the collapse controls, the random tests, the mixing gap and the direction prices of §7.

## 2. Lemma S: every content pays for itself

The ancestry anc(δ) of a step is the smallest set containing δ and closed under δ'' ∈ anc, a_{δ''δ'} ≠ 0 ⇒ δ' ∈ anc. Let σ(u) be the span of the labels of the support of u.

**Lemma S (mixing allowed).** For every step δ and every set R of ports seeing G_δ:

  |anc(δ)| ≥ dim(Σ_{δ'∈anc(δ)} F_{δ'}) − dim G_δ ≥ dim(σ(u_δ) + G_δ) − dim G_δ ≥ rank B[R, supp u_δ].

*Proof.*
1. List anc(δ) from latest to earliest, starting with δ. Add the start frames one at a time.
2. dim F_δ = dim G_δ + 1.
3. Every later-listed δ' feeds an earlier-listed member δ'' with G_{δ'} ⊆ F_{δ''}. So the sum so far already contains G_{δ'}, and adding F_{δ'} raises the dimension by at most dim F_{δ'} − dim G_{δ'} = 1. This gives the first inequality.
4. By the recursion, σ(u_δ) ⊆ Σ F_{δ'}, and G_δ ⊆ F_δ. This gives the second.
5. U_R^⊥ ⊇ G_δ, so dim σ − dim(σ ∩ U_R^⊥) ≤ dim(σ + G_δ) − dim G_δ. ∎

So an owned descent can never deliver a content more cheaply than Model C prices that one delivery. If multi-generation gains at all, it gains because one ancestry serves several deliveries.

## 3. Pure steps and the Union lemma

The definitions are on the content, not on its representation (which need not be unique).
- A step is **carried** if u_δ = a u_{δ'} for some δ' ∈ E(δ) and a ≠ 0. Fix one such δ' as its parent.
- Otherwise it is **fresh** if u_δ ∈ Pub(F_δ).
- Otherwise it **mixes**: no representation uses only first-generation material or only one landed content.

A centre is **mixing-free** if no step mixes.
- In a mixing-free centre, the parent links make the steps a forest rooted at fresh steps.
- A *class* C is one tree. All its contents are nonzero multiples of the root's v_C. Zero contents are dropped; they only cost.
- An owned lifetime that copies x-registers once and then only descends is a class. Re-descending one landed content alone extends its class.

**Union lemma.** For a class C with root δ_0:

  |C| ≥ dim F_{δ_0} − dim(F_{δ_0} ∩ ⋂_{δ∈C} G_δ).

*Proof.*
1. Process C in time order, with X := F_{δ_0} and X ← X ∩ G_δ at each step.
2. When δ is processed, X ⊆ F_δ. For δ_0 this holds by definition. For a carried δ with parent δ', X ⊆ G_{δ'} ⊆ F_δ.
3. G_δ is a hyperplane of F_δ, so each step lowers dim X by at most one. ∎

## 4. The Collapse theorem

**Theorem (Collapse).** Let Ω be a mixing-free centre, with arbitrary contents and multi-generation re-descents allowed. Then cost(Ω) is at least the Model C cost of the decoder

  N_Ω = Σ_C p_C v_Cᵀ,  p_C(t) = Σ_{δ∈C, G_δ ⊆ q_t^⊥} β_{tδ} a_δ (u_δ = a_δ v_C).

N_Ω is a realisation of the P104 pattern.

*Proof.*
1. By Lemma W, row t of N_Ω is w_t. It vanishes on orthogonal pairs and N_tt ≠ 0.
2. Every user t ∈ supp p_C sees some landing of C. Hence U_{supp p_C}^⊥ ⊇ ⋂_{δ∈C} G_δ, and U_{supp v_C} ⊆ F_{δ_0}.
3. By the rank formula and submodularity:

   rank B[supp p_C, supp v_C] ≤ dim F_{δ_0} − dim(F_{δ_0} ∩ ⋂_C G_δ),

   and the Union lemma bounds the right side by |C|.
4. The classes partition the steps, so summing over C gives the claim. ∎

Model C's normalisations do not raise the Model C cost, so the bound survives them (#2338 §3.3, items 1.1 and 1.2):
- delete free terms (B[R_j, A_j] = 0); they change N only on non-orthogonal pairs, so the result is again a realisation, possibly of another rank;
- merge terms with proportional column factors (#1224 Lemma D), and likewise for rows; merging does not change N.

Write N'_Ω for the normalised decoder and m for its number of terms. m is at most the number of pairwise non-proportional paid first-generation contents. **All ranks below are ranks of N'_Ω**, as in #2338 item 1.4, where ρ is the rank of the paid part.

**Controls (CHECKED).**
- The selected centre has eight classes. Each collapses to a Model C term of cost 7, which equals its Union bound and its seven steps (total 56 = 56).
  The x_j form a basis, so each port's expansion is unique and the users of x_j are exactly its 52 needers. Each content's support lies in its copy frame h_j^⊥.
- The R_8 = 0 centre of research/p104-multigen-centre-floor collapses to eight cost-7 terms (56 = 56). It re-descends each content once: 6 steps to a plane P_j ⊆ h_j^⊥, then 1 step to a plane P'_j with P_j ∩ P'_j = [g_j].
  - The second descent needs a public carrier for the landed content up to P_j + P'_j, followed by a second owned lifetime.
  - It is pure when the carrier holds only that content, under the predecessor's convention of free helper endpoints.
- Random tests, exact:
  - 2,584 random pure classes satisfy the Union lemma;
  - 2,827 random mixing ancestries satisfy all three inequalities of Lemma S;
  - the minimum slack is 0 in both;
  - users are all ports seeing any landing (the worst case), and landings at 0 may feed later steps;
  - of 3,000 attempts each, those with empty support are skipped.

## 5. Consequences for P104

**The unit (CHECKED).** A linear isometry carries the word's x-basis onto removed roots of #2338's E8, and the word's 104 live and 16 removed lines onto #968's P104 = E8 − F0 − H4 (`decoder-floors-frontier-format/lib/units.py`). So the Model C floors proved for #968's unit apply to the word as stated.

Let Ω be mixing-free, N'_Ω its normalised collapsed decoder, and m its number of terms.

**Corollary 1.** If m ≤ 8, then cost(Ω) ≥ 56. This holds whatever the contents, row factors or realisation.
- #2338 §3.3, item 1.3: rank N'_Ω ≥ 8 forces m = 8 = rank and square factors, and Theorem P gives ≥ 56 on #968's unit.

**Corollary 2.** If m = 9 and rank N'_Ω = 8, then cost(Ω) ≥ 56.
- By Sylvester (#2338 1.4), s + t ≤ 1, where s and t are the excess ranks of the two factor matrices over 8. At most one factor matrix has rank 9.
- With s = t = 0, this is Theorem P.
- With s + t = 1, it is case (b) at m = 9, which research/nine-term-outside-colspace (claim C2) states for #968's P104. Its status there is "PROVEN from CHECKED", by an exact class-level enumeration not re-run here.

**Corollary 3.** If rank N'_Ω = 8 and all factors lie in col N'_Ω / row N'_Ω, then cost(Ω) ≥ 56 for every m. This is Theorem P.

**Corollary 4 (one-sided, every m).** Let N = Σ_j p_j c_jᵀ be a Model C decoder of #968's P104 with column factors in a gauged root code, c_j(s) = d_s⟨z_j, r_s⟩, and arbitrary row factors (case (b)). Define:
- the depth of term j: λ_j := dim(U_{R_j} + U_{A_j}^⊥), the level of its cheapest landing frame U_{A_j} ∩ U_{R_j}^⊥;
- R_d := dim span{z_j : λ_j ≥ d}.

Then cost(N) ≥ 56 whenever R_8 ≥ 7, or R_8 = 6 and R_7 = 8.

*Proof.* N is a single-generation root-code centre of the abstract model of research/p104-multigen-centre-floor:
1. Term j is one owned descent. It copies c_j at σ(c_j) = U_{A_j} (F2^9, or ℓ_ρ^⊥ when z_j ∝ a removed root ρ, by that package's Lemma 2) and lands at U_{A_j} ∩ U_{R_j}^⊥. Its cost dim U_A − dim(U_A ∩ U_R^⊥) is the Model C rank by the rank formula.
2. Every port t of R_j sees the landing, and the entry condition r_t ∈ span{z_j : t ∈ R_j} holds. Row t of N is (d_s⟨Σ_j p_j(t) z_j, r_s⟩)_s. It vanishes on the live roots orthogonal to r_t, which span r_t^⊥, so Σ_j p_j(t) z_j ∝ r_t.
3. The master inequality and Theorems A–C of that package use only the entry condition, Lemma 2 and the unit-step levels, and its R_d is the R_d above. ∎

Transposition gives the same floor for decoders with **arbitrary contents and root-code coefficient vectors**, with λ_j defined by exchanging R and A. This is the one-sided form of "contents outside the root code".

**Where a cheaper centre must live.** A P104 centre of cost ≤ 55 must have at least one of:
1. a step that **mixes** landed contents (or a landed content with fresh material);
2. at least **nine** pairwise non-proportional paid first-generation contents, with N'_Ω outside the closed classes:
   - m = 9 with rank N'_Ω = 9 (case (c), a square decoder of a rank-9 realisation), or
   - m ≥ 10 outside the scope of Theorem P: rank N'_Ω ≥ 9, or rank 8 with a factor outside col N'_Ω / row N'_Ω;
   - if one side is root code, also in a shallow profile of Corollary 4.

This holds for any contents, inside or outside the root code.

**Remark (transposition).** N ↦ Nᵀ exchanges contents and coefficient vectors and preserves the pattern and the cost, since B is symmetric. Hence:
- the one-sided version of "contents outside the root code" (contents arbitrary, coefficient vectors in the root code) is exactly case (b) for Nᵀ;
- for P104 that case is closed for m ≤ 9 (Corollaries 1–2), and for every m in the deep profiles of Corollary 4;
- it stays open for m ≥ 10 in the shallow profiles.

## 6. The mixing gap (CHECKED)

The Collapse theorem needs mixing-freeness. Here is the smallest mixing pattern, priced exactly:
1. The frame-pair contents x_0 and x_1 are copied at their frames ℓ^⊥ and descend to a common line g ⊆ ℓ_{x_0}^⊥ ∩ ℓ_{x_1}^⊥: 7 + 7 steps.
2. A third owned step copies their sum at g and descends to 0: 1 step. This step mixes.

The centre pays 15. The 56 ports of world(g) get span{x_0, x_1}, and every other port gets x_0 + x_1.

Model C prices:
- the three deliveries term by term: 7 + 7 + 9 = 23. The sum's support labels span F2^9;
- the two natural Model C decoders giving every port at least the same span: 16 each. These are (x_0, world(g)) + (x_0 + x_1, all), at 7 + 9, and (x_0, all) + (x_1, all), at 8 + 8.

Over all 28 pairs and every common line, the cheapest of these alternatives is 16. Optimality over all Model C decoders is not claimed.

So in this fragment mixing pays one unit less than the best Model C alternative found. Lemma S still holds for the sum: 15 ≥ 9. Whether such savings survive in a complete P104 centre is the open question this package isolates.

## 7. Evidence about contents outside the root code (CHECKED, not proofs)

**Truncation.** One natural non-root-code content is a frame-pair content truncated to the ports whose labels lie in a 7-dimensional subframe h^⊥ ∩ k^⊥ that still contains its landing line. Its first descent costs 6 instead of 7.
- 1,016 distinct truncations were tested. Every port re-solves its combination exactly.
- Every truncation leaves all 52 users of that content without a clean vector from the remaining contents.
- So a single truncation must come with repair deliveries reaching all 52 users. The census does not price such repairs; one repair term could serve many users.

**Overcomplete directions.** Directions are drawn from the 16 removed roots, each with the −1 label bonus.
- Port t takes any minimal set D_t of at most four removed roots spanning r_t, or its frame-pair support. The smallest such sets have sizes 2, 3, 4 for 32, 48, 24 ports.
- Larger minimal sets exist (port 0 has 12, 76, 240 and 1,401 of sizes 5–8) and are not explored.
- Direction ρ is delivered to its users at the cheapest landing frame. The price uses U_{supp ρ} = ℓ_ρ^⊥, which is asserted for all 16, and is cross-checked by both rank computations.
- Within this neighbourhood the frame-pair (cost 56) is a strict local optimum: all 524 single-port changes cost ≥ 57.

**Why frame counting cannot prove 56.** Every live 8-frame (759 of them) must meet the worlds of at least eight terms, since their directions span Q^8. A counting proof would need every term of cost c to meet at most 759c/7 frames. Landing frames inside removed hyperplanes break this:

| Level | Cost | World | Frames met | Counting allowance |
|---:|---:|---:|---:|---:|
| 7 | 6 | 32 | 759 | 650.6 |
| 6 | 5 | 19 | 699 | 542.1 |
| 5 | 4 | 11 | 542 | 433.7 |

The three witnesses are fixed in `checks.FRAME_TOUCH_WITNESSES` and checked exactly. A scratch census in C over all 8.3 M subspaces of F2^9 found them; it is not shipped and not needed for the claim. Any proof for one-sided or general contents must use the span structure, not frame incidence alone.

## 8. Open

- **Mixing**, with any contents. This is the one multi-generation mechanism Model C does not price.
- **Case (c) at m = 9** for #968's P104: the rank-9 square decoder. The frame-module theorem of #2338 §3.3 is the available tool.
- **m ≥ 10 with a factor outside the column/row space.**
  - The one-sided case (b) remains open only in the shallow profiles of Corollary 4: R_8 ≤ 5, or R_8 = 6 with R_7 ≤ 7. This is the shallow covering problem of research/p104-multigen-centre-floor, now with single-generation prices.
  - The two-sided case (d) is fully open.
- For comparison, P112's one-sided case is closed only in part: in the Gram realisation, with at most 2, 3 and then 4 non-full corrections (one-sided-decoders-all-m, one-sided-six-multisets, one-sided-excess-three).

No supplier, κ or selected-result change is claimed. The bounds are on descent mass, not on the r^t moment.
