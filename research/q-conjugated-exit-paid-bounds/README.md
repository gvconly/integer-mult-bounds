# Paid floors for q-conjugated exits: GLOBAL120 termination and the odd-E8 negative decoder (no κ claim)

Two open lines end with the same missing piece: a decoder for a bank whose code has been conjugated by the central Gaussian Q. Each line asks whether a paid decoder fits its budget. This package proves exact lower bounds on the number of paid calls, of any rank, that such a decoder needs.

The counting follows [coded-active-query-bound](../coded-active-query-bound/README.md) and [encoded-exit-obstructions](../encoded-exit-obstructions/README.md). It is stated here in an output-coset form (Lemma A), together with a suffix-projection form (Lemma B). A new path lemma covers three fresh stages (Lemma D).

## Results

| Target | Paid calls needed | Available or allowed |
|---|---:|---:|
| GLOBAL120 exit q⁻¹J⁻¹q on one 120-field bank | **≥ 120** (240 suffice) | — |
| GLOBAL120 nominated retirement cut (4)–(5) | **≥ 240** as nominated (C_P or C_I); **≥ 232** if Y0 was decoded earlier | 112-port query bank + 8 mode calls = 120 |
| GLOBAL120 retirements of pairs 0–6 (no scheduled paid call on Y_k after contact k) | **≥ 840** new calls, or existing calls rerouted from other banks | 916 additional rank-1 jobs at C* = 0.00330007 |
| GLOBAL120 retirements of pairs 7–14 | ≥ 120 each, met in count by each pair's suffix IDLE (112 main + 8 mode calls) | — |
| Odd-E8 negative post-decoder J = Q⁻¹E1⁻¹E2E3⁻¹Q | **≥ v = 112³ = 1,404,928** per emitted bank (2v suffice) | — |
| The same, per port per local call | **≥ 1/3** | margin 0.04382 (exceeded **7.6×**) |

### GLOBAL120

- The 120 owners' centre lines in the balanced code's normalized chart are pairwise distinct, under all four chart conventions checked. So every output row of the q-conjugated exit has private non-free coset coordinates.
- **The nominated cut is not buildable from the retained calls.** The exit's coordinates belong only to the retiring bank. The continuing neighbour needs a differently coded contact in the same Q-cosets. The retained 120 calls cannot serve both, so the cut needs at least 112 more.
- **Pairs 0–6 cost at least 840.** For each contact k, the controller and the neutral literal invoice were checked for paid calls landing on Y_k afterwards. For pairs 0–6 there are none. Lemma B then forces at least 7·120 = 840 new or rerouted paid calls. Unless those calls are rerouted, that is 91.7% of the optimistic allowance, before any other open GLOBAL120 cost.
- **Pairs 7–14 are open.** Their suffix IDLEs already apply the 120 calls Lemma B requires. Whether an IDLE can also carry the decode is the open coded-caller IDLE question. Counting does not decide it, and it does not exclude rerouting.

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
  - IDLE-borne decodes;
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

1. **IDLE-borne decodes.** Can an already-paid IDLE carry the q-conjugated decode, for GLOBAL120 pairs 7–14 and for coded callers generally?
   - Why it matters: counting allows it, and it is the main remaining lever for both lines.
   - First step: one frozen pair with its suffix IDLE, tested against all-input targets.
   - Blocked by: the 120-wide RIGHT, endpoint and allocation interface listed in the GLOBAL120 chart.
2. **Correlated histories.** Can a correlated or pre-emission copy of the negative bank go below one paid call per owner?
   - First step: price one explicit pre-emission copy route for the odd-E8 bank, with its retained-history invoice.
3. **Exact complexity.** Where between v and 2v, and between 120 and 240 for one GLOBAL120 bank, does the true count lie? This is the same gap that coded-active-query-bound leaves at 112 and 224.

## Credit

- **Counting lemmas.** Dillon's coded-active-query-bound, with its credited lineage (Dugongue #4929; Chris R/Clear1nsight #4932; Cinna Davis #4694/#4698; Cao Yifan #4610; zxv). Joel Pulikkan's encoded-exit-obstructions (Lemmas 1–2).
- **Odd-E8.** JoyboyBrian's PR5508 (odd-e8-negative-exit-algebra) supplies the premises, target and receipts. THU-NMRC's PR5464 (odd-e8-lowrank-e-native) supplies the three-call word and margin. Ali Romman supplies the odd-E8 geometry and the fresh-stage separation lineage (PR4933/5065). The full120 package is THU-NMRC's.
- **GLOBAL120.** The chart, controller, invoice and screen come from global120-optimistic-budget. That package credits Ali Romman's controller, Dillon's coupled provider, Chris R's neutral completion and the native and compiler predecessors. The balanced-J native source is credited as in its package.
- **New here:**
  - the output-coset and suffix-projection forms of the counting lemmas;
  - the cut and retirement bounds, and the invoice-checked controller audit;
  - the three-stage private-witness lemma;
  - the budget comparisons.

Prepared by G. V. Chudnovsky with Anthropic Claude assistance. Apache-2.0.
