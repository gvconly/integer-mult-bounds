# Paid floors for q-conjugated exits: GLOBAL120 termination and the odd-E8 negative decoder (no κ claim)

Two open lines end with the same missing piece: a decoder for a bank whose code has been conjugated by the central Gaussian Q. Each line asks whether a paid decoder fits its budget. This package proves exact lower bounds on the number of paid calls, of any rank, that such a decoder needs.

The counting follows [coded-active-query-bound](../coded-active-query-bound/README.md) and [encoded-exit-obstructions](../encoded-exit-obstructions/README.md). It is stated here in an output-coset form (Lemma A), together with a suffix-projection form (Lemma B). Two additions are new: an IDLE-commutation analysis (Lemmas D and E) and a path lemma for three fresh stages (Lemma G).

## Results

| Target | Paid calls needed | Available or allowed |
|---|---:|---:|
| GLOBAL120 exit q⁻¹J⁻¹q on one 120-field bank | **≥ 120** (240 suffice) | — |
| GLOBAL120 nominated retirement cut (4)–(5) | **≥ 240** as nominated (C_P or C_I); **≥ 232** if Y0 was decoded earlier | 112-port query bank + 8 mode calls = 120 |
| GLOBAL120 retirements of pairs 0–6 (no scheduled paid call on Y_k after contact k) | **≥ 840** new calls, or existing calls rerouted from other banks | 916 additional rank-1 jobs at C* = 0.00330007 |
| GLOBAL120 retirements of pairs 7–14 | ≥ 120 non-IDLE calls each (in-sector, code-chart assumption); their suffix IDLEs commute with the code | — |
| GLOBAL120, all 15 retirements (same conditions) | **≥ 1,800**, new or rerouted | 916 |
| Odd-E8 negative post-decoder J = Q⁻¹E1⁻¹E2E3⁻¹Q | **≥ v = 112³ = 1,404,928** per emitted bank (2v suffice) | — |
| The same, per port per local call | **≥ 1/3** | margin 0.04382 (exceeded **7.6×**) |

### GLOBAL120

- The 120 owners' centre lines in the balanced code's normalized chart are pairwise distinct, under all four chart conventions checked. So every output row of the q-conjugated exit has private non-free coset coordinates.
- **The nominated cut is not buildable from the retained calls.** The exit's coordinates belong only to the retiring bank. The continuing neighbour needs a differently coded contact in the same Q-cosets. The retained 120 calls cannot serve both, so the cut needs at least 112 more.
- **Pairs 0–6 cost at least 840.** For each contact k, the controller and the neutral literal invoice were checked for paid calls landing on Y_k afterwards. For pairs 0–6 there are none. Lemma B then forces at least 7·120 = 840 new or rerouted paid calls. Unless those calls are rerouted, that is 91.7% of the optimistic allowance, before any other open GLOBAL120 cost.
- **Scheduled IDLEs commute with the conjugated code.** A stage Gaussian commutes with a code living in another stage's chart. Only a layer on the code's own stage, or a q, flips the code to an ordinary one: D_t·barJ_t = (X⁻¹J_tX)·D_t (Lemma D).
- **Every IDLE misses the code's stage.** The verifier checks that every IDLE supplies only stages its pair does not run itself. The package assumes, as the chart does at contact 0, that the retiring code lives in the chart of the pair's last FULL stage, 2k+1. Under that assumption no scheduled IDLE contains the code's stage.
- **What that forces, in-sector.** Suppose every operation on retired-bank data commutes with the IDLE's Gaussian. Then the IDLE factors out, and each retirement needs at least 120 non-IDLE calls (Lemma E). That totals at least 1,800 additional or rerouted calls, against 916. In that sector, the optimistic allowance cannot fund the Gamma termination unless at least 884 existing calls are rerouted, or a changed native FULL or contact places a stage-(2k+1) or central layer after each retirement.
- **What remains open.** The in-sector condition excludes the lineage's proposed IDLE pre/post codes in the IDLE stages' charts. Without it, only the 840 for pairs 0–6 is proven, and those codes stay open for pairs 7–14.

### Odd-E8

- For three fresh stages, the path j → d → d → d has a binary signature that no other path from input j produces. So every owner row of J has its own non-free coordinate.
- The witness coefficient is 1 · (1541/104663) · (K+P)_{aa′}/2. Every base row of K+P has 48 off-diagonal nonzeros, so the coefficient is nonzero for every owner.
- Hence any post-decoder of the emitted negative bank needs at least one paid call per owner: one third of an H1 per port per local call. The whole conditional margin of the source screen is 0.0438.
- The N3 word has no IDLE and no paid call after the emission.
- If the same premises hold for full120 (not re-verified here), its declared margin of 0.0893 is exceeded 3.7×.

The written argument, premises and scope are in [PROOF.md](PROOF.md).

## Scope

- **What is covered.** Post-decoders that see the emitted bank and independent workspace, with paid calls in the monomial and phase-word algebra (Gaussians of any rank), in the width-uniform model of the cited packages.
- **What is not covered.**
  - correlated or pre-return histories (encoded-exit-obstructions prices one such route at 8 H8 against 112 H1);
  - changed targets or contact types;
  - IDLE pre/post codes in the IDLE stages' charts, and cross-stage re-aiming (outside Lemma E's sector);
  - changed native FULLs;
  - GLOBAL120 designs that reroute other banks' scheduled paid calls through retired banks.
- **Not tight.** The bounds are counts, and the elementary decoders use twice as many calls.
- **No κ claim.** No supplier, κ or selected-result change is claimed.

## Reproduce

Python 3.12, standard library only, from the repository root. Takes about 20 s.

```sh
python3.12 -I -B research/q-conjugated-exit-paid-bounds/verify.py
python3.12 -B scripts/research_package.py --package research/q-conjugated-exit-paid-bounds
```

The verifier:
- reads the eight pinned live inputs listed in VERIFY.json;
- recomputes every number above;
- compares the replay with RESULT.json;
- refuses to run under `-O`, because its checks are assertions.

Controls:
- a common centre line, or a diagonal code, gives rank 0;
- the uncoded contacts reproduce the retained 112 and 120 calls;
- the old 112-port code reproduces coded-active-query-bound's 112;
- removing cross-cube couplings removes every odd-E8 witness.

## Open next questions

1. **Levers outside Lemma E.** Can IDLE pre/post codes in an IDLE stage's chart, or a changed native FULL that defers part of its own-stage Gaussian past the contact, carry the decode without extra paid mass?
   - Why it matters: these are the GLOBAL120 levers that Lemma E's sector excludes, besides rerouting.
   - First step: write pair 7's retirement with a stage-17 pre/post code around its suffix IDLE, and test it against all-input targets.
   - Blocked by: the missing 120-wide RIGHT, endpoint and allocation interface.
2. **Correlated histories.** Can a correlated or pre-emission copy of the negative bank go below one paid call per owner?
   - First step: price one explicit pre-emission copy route for the odd-E8 bank, with its retained-history invoice.
3. **Exact complexity.** Where between v and 2v, and between 120 and 240 for one GLOBAL120 bank, does the true count lie? This is the same gap that coded-active-query-bound leaves at 112 and 224.

## Credit

- **Counting lemmas.** Dillon's coded-active-query-bound, with its credited lineage (Dugongue #4929; Chris R/Clear1nsight #4932; Cinna Davis #4694/#4698; Cao Yifan #4610; zxv). Joel Pulikkan's encoded-exit-obstructions (Lemmas 1–2).
- **Odd-E8.** JoyboyBrian's PR5508 (odd-e8-negative-exit-algebra) supplies the premises, target and receipts. THU-NMRC's PR5464 (odd-e8-lowrank-e-native) supplies the three-call word and margin. Ali Romman supplies the odd-E8 geometry and the fresh-stage separation lineage (PR4933/5065). The full120 package is THU-NMRC's.
- **GLOBAL120.** The chart, controller, invoice and screen come from global120-optimistic-budget. That package credits Ali Romman's controller, Dillon's coupled provider, Chris R's neutral completion and the native and compiler predecessors. The balanced-J native source is credited as in its package.
- **New here:**
  - the output-coset and suffix-projection forms of the counting lemmas;
  - the cut and retirement bounds, the IDLE-commutation analysis, and the invoice-checked controller audit;
  - the three-stage private-witness lemma;
  - the budget comparisons.

Prepared by G. V. Chudnovsky with Anthropic Claude assistance. Apache-2.0.
