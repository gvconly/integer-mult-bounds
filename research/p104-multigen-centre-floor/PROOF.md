# Proof: the P104 centre floor with multi-generation descents

All finite facts are recomputed by `verify.py` from the published word `research/p104-productive-controller/WORD.json.gz`. RESULT.json records them.

## 1. Model

**The word model.** We use the P104 local word model of that package (`local.py`):
- Public registers only climb, along nested frames. Owned registers only descend.
- Every MOVE is paid its rank |Δdim|. ADD and COPY need a common frame.
- The 208 mains are fixed. x_p goes from [q_p] to F2^9, and y_p goes from 0 to q_p^⊥, ending with e_{y_p} + e_{x_p}.

The **centre cost** of a word is the total rank of its owned descents. The published word pays 56: eight rank-7 lifetimes.

**Root realisation** (`common.Geometry`, recomputed). Write the port content of owned lifetime j as u_j.
- Then u_j(s) = ⟨x_j, r_s⟩/4.
- The vectors x_j form a Z-basis of an E8 lattice: the Gram determinant is 1.
- The r_s are roots, unique up to scale: the symmetric solutions of r_sᵀ Y r_t = 0 on F2-non-orthogonal label pairs form a 1-dimensional space.
- The labels satisfy ⟨q_p, q_s⟩ = 1 iff r_p ⊥ r_s or p = s.
- The 16 non-live lines (removed roots) are two orthogonal 8-frames. Each has an F2 label ℓ_ρ, characterised by ⟨ℓ_ρ, q_p⟩ = [ρ ⊥ r_p].

**The selected centre** is the frame-pair decoder:
- The content of lifetime j is the removed root x_j.
- It is copied at the frame ℓ_{x_j}^⊥ and lands on the line [ℓ_{y_j}], where (y_j) is the dual basis (also removed roots).
- For every port, supp(r_p) in the x-basis = {j : ⟨q_p, ℓ_{y_j}⟩ = 0}.

**Hypothesis RC.** The port projection of every owned content is a root-code vector, c_x = (⟨x, r_s⟩)_s. Per-port combinations (row factors) and the number of lifetimes are unrestricted.

**Multi-generation** is allowed: an owned descent may carry content landed by earlier owned descents. This is the case not priced by Model C (cost Σ rank B[R_j, A_j]), which charges each term from first-generation material.

## 2. Lemmas

**Lemma 1 (world entry).** In any valid word, for every port p:

  r_p ∈ span{ x : an owned content c_x landed at a frame ⊆ q_p^⊥ }.

*Proof.*
1. Call frames ⊆ q_p^⊥ the world of p. Public frames only grow, and ADD/COPY need equal frames. So rows of registers in the world are combinations of initial rows of registers that start in the world, and of owned contents that land there.
2. Project onto the port coordinates M_p = {p} ∪ {s : r_s ⊥ r_p}. Initial rows starting in the world vanish there: x_s with s ≠ p starts in the world iff ⟨q_s, q_p⟩ = 0, i.e. s ∉ M_p.
3. y_p ends with the coordinate e_{x_p}. So some combination y = Σ λ x satisfies ⟨y, r_s⟩ = 0 for s ∈ M_p \ {p}, and ⟨y, r_p⟩ ≠ 0.
4. Every port is good: the live roots orthogonal to r_p span r_p^⊥ (checked for all 104). Hence y ∝ r_p. ∎

**Lemma 2 (first-generation sources).** A root-code vector made at frame F from first-generation material is supported on {s : q_s ∈ F}. So it is c_x with x ∈ K(F) := (span{r_s : q_s ∉ F})^⊥, and:
- K(F2^9) = Q^8;
- K(ℓ_ρ^⊥) = [ρ] for each removed root ρ;
- K(F) = 0 for every other frame.

This is checked for all 511 hyperplanes. Smaller frames lie in at least three hyperplanes, at most one of which can be removed-type with a given ρ. Under RC the first-generation part of any owned content is a difference of root-code vectors, hence root-code.

**Unit-step relaxation.**
1. Replace each owned MOVE F → G of rank r by r codimension-1 steps through any chain, with the content landing at every intermediate frame. Cost is unchanged and availability only grows, so lower bounds for the relaxed process hold for the word.
2. Each step δ has a start frame F_δ and a landing frame G_δ with dim G_δ = dim F_δ − 1.
3. Its content satisfies x_δ ∈ K(F_δ) + span{x_{δ'} : δ' earlier, G_{δ'} ⊆ F_δ}.

**Notation.**
- level(G) = 9 − dim G, and n_d = #steps landing at level d. The cost is Σ_d n_d.
- C_G = the landed contents at G, and A(G) = span{x_δ : G_δ ⊆ G}.
- S_F = Σ_{G⊊F} A(G), the fresh dimension φ(F) = dim A(F) − dim S_F, and f_d = Σ_{level F = d} φ(F).
- V'_d = Σ_{level G ≥ d} A(G), with R_d = dim V'_d.
- F_sh = f_2 + … + f_7 is the shallow fresh count.

**Lemma 3 (feed inequality).** For 2 ≤ d ≤ 8: n_d ≥ R_{d+1} + f_d. Also n_9 ≥ R_9 and n_1 ≥ f_1.

*Proof.*
1. Take a step landing at level ≥ d+1 whose start frame has level ≥ d+1. That frame has dim ≤ 7, so K = 0 there and the step adds nothing to V'_{d+1}.
2. Take a step from a frame F of level d. Again K(F) = 0, since dim F ≤ 7. So x_δ = c_δ + v_δ, with c_δ ∈ span C_F and v_δ ∈ V'_{d+1}(earlier).
3. Put K_F = span{c_δ}. Then V'_{d+1} = Σ_F K_F, and K_F ⊆ span C_F ∩ S_F, because x_δ and v_δ are contents at proper subframes of F.
4. Therefore n_d = Σ_F |C_F| ≥ Σ_F [dim K_F + dim(span C_F/(span C_F ∩ S_F))] ≥ R_{d+1} + Σ_F φ(F). ∎

**Lemma 4 (fresh basis).** Choose B_F ⊂ A(F) of size φ(F) completing S_F to A(F). Then:
- A(G) = span ⋃_{F⊆G} B_F and V'_d = span ⋃_{level F ≥ d} B_F;
- by Lemma 1, r_p ∈ span ⋃_{F⊆q_p^⊥} B_F for every port.

**Master inequality.** Sum Lemma 3 over all levels and use f_8 ≥ R_8 − R_9 (Lemma 4) and R_d ≥ R_8 for d ≤ 8:

  cost ≥ 7·R_8 + Σ_{d=3}^{7} (R_d − R_8) + R_9 + F_sh + f_1.

**Lemma 5 (capacities).**
- A frame F at level d serves the ports whose labels lie in the d-dimensional space F^⊥, and m_d is the maximum number of port labels in such a space.
- Exhaustively, m_d = 1, 2, 4, 8, 12, 20, 36, 60 for d = 1..8.
- Every subspace W ⊂ F2^9 was enumerated: 511, 43,435, 788,035 and 3,309,747 per dimension.

**Lemma 6 (shallow port sets).** Put T_z = {p : ⟨z, r_p⟩ ≠ 0} for z ≠ 0.
1. **|T_z| ≥ 49.** Let S be the span of the live roots in z^⊥.
   - The roots in S form a simply-laced subsystem of rank ≤ 7.
   - A rank-7 simply-laced system with more than 42 lines is E7, whose orthogonal complement in E8 is a root. Root hyperplanes hold 55 live roots (live root) or 52 (removed root). A rank-6 system has ≤ 36 lines.
   - verify.py also confirms min |T_z| = 49 over all 13,320 lattice vectors of norm ≤ 8 and over all 3,360 A2 continua.
2. **Which z can matter.** The sets with |T_z| ≤ 72 come from hyperplanes holding ≥ 32 root lines (E7, D7, E6+A1), or from z in the A2-plane orthogonal to an E6.
   - In the unimodular lattice E8, det(L) = det(L^⊥) for primitive L. So the normals of E7, D7 and E6+A1 hyperplanes have norms 2, 4 and 6, and the norm ≤ 8 enumeration covers them.
   - The E6 complements are the A2 continua.
3. **Two-plane lemma.** T_z ⊆ W(P) ∪ W(P') for two planes P, P' iff T_z = T_ρ for a removed root ρ. In that case z ∝ ρ: the 52 live roots orthogonal to ρ span ρ^⊥. This is checked for all 4,016 candidate sets against all 43,435 planes.

**Lemma 7.** No three removed roots lie in a 2-plane. Two roots of the same frame are orthogonal. For x ∈ F0 and y ∈ F1 non-orthogonal, x ± y is non-orthogonal to x and to y, while every other member of each frame is orthogonal to its own frame-mate, so x ± y is not removed. This is checked.

## 3. Theorems

**Theorem A.** If R_8 = 8, i.e. the contents landed at lines and at 0 span Q^8, then cost ≥ 56.

*Proof.* This is the master inequality with R_8 = 8. ∎

**Theorem B.** If R_8 = 7, then cost ≥ 56.

*Proof.*
1. Let δ = max{d ≤ 7 : R_d = 8}, and y a functional vanishing on V'_8 = V'_{δ+1}.
2. For p ∈ T_y, Lemma 4 gives a fresh vector b ∈ B_F with F ⊆ q_p^⊥ and y(b) ≠ 0. Then b ∉ V'_{δ+1}, so level(F) ≤ δ.
3. A level-1 frame inside q_p^⊥ is q_p^⊥ itself, so each level-1 fresh frame serves one port. A frame at level 2..δ serves ≤ m_δ ports.
4. So m_δ·F_sh + f_1 ≥ |T_y| ≥ 49.
5. The master inequality gives cost ≥ 49 + max(0, δ − 2) + ⌈49/m_δ⌉ for δ ≥ 2, and ≥ 98 for δ = 1. For δ = 7, 6, 5, 4, 3, 2 this is 56, 56, 57, 58, 63, 74. ∎

**Theorem C.** If R_8 = 6 and R_7 = 8, then cost ≥ 56.

*Proof.* The master inequality gives cost ≥ 42 + 10 + F_sh + f_1, so it suffices to show F_sh + f_1 ≥ 4. Let Z be the plane of functionals vanishing on V'_8, and project the shallow fresh vectors to Q^8/V'_8. Suppose at most three have nonzero projection.
1. Their projections cannot all be parallel: a z ∈ Z vanishing on all of them would have no active frame for T_z ≠ ∅.
2. Take z_1 ∈ Z vanishing on one projection. At most two active fresh vectors remain, and they must cover T_{z_1} (≥ 49).
3. One frame-world holds ≤ 36 ports, and a level-1 frame holds one. So exactly two remain, at distinct frames of levels 2..7, not parallel to the first.
4. Applying this to each of the three vectors: the projections are pairwise independent, the three frames are distinct, and each T_{z_j} is covered by two frames of level ≤ 7, hence by two plane-worlds.
5. By Lemma 6, each z_j ∝ a removed root. These are three pairwise independent removed roots in the 2-plane Z, contradicting Lemma 7. ∎

## 4. What is not covered

**The open profiles.** Theorems A–C leave exactly:
- at least three directions that never reach a line (R_8 ≤ 5);
- exactly two, with at least one ending above the planes (R_8 = 6 and R_7 ≤ 7).

**Why covering alone does not close them.** In those profiles the master inequality reduces the question to the **shallow covering problem**. Let Z = (V'_8)^⊥, with dim Z = s.
- Choose pairs (β_j ∈ Z, frame F_j at level ≤ 7, with β_j in the flag piece of its level).
- For every nonzero z ∈ Z, these must satisfy T_z ⊆ ⋃_{⟨z, β_j⟩ ≠ 0} W(F_j).
- This is equivalent to proj_Z(r_p) ∈ span{β_j : p ∈ W(F_j)} for every port.
- One must show the number of pairs is ≥ Σ_i (9 − max(δ_i, 2)) over the shallow directions i.
- Choosing z orthogonal to s − 1 of the β's gives only the Singleton-type bound s + 1, while 2s is needed even when every δ_i = 7.

**Recorded data.**
- A live root's T is not covered by three 3-dimensional frames.
- A removed root's T is.

**Scope.**
- Only root-code centre contents are covered. Non-root-code contents (rank-9 realisations, factors outside the root code) are out of scope.
- The bound is on descent mass, not on the r^t moment.
- No supplier, κ or selected-result change is claimed.
