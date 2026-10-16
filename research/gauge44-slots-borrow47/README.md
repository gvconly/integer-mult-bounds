# 47 anticipatory dirty-borrowing aliases on the seven-slot k8 word at T = 8: conditional κ = 0.001687525611032835

This package applies **Cinna Davis's anticipatory dirty borrowing (PR2226, `research/anticipatory-dirty-borrowing`)** to
the selected seven-slot k8 word (the PR1042/1094 word; its canonical form is the `CANONICAL-WORD.json.gz` pinned by
`research/gauge44-slots-arcflow/data/PLAN.json`). It lowers the result natively and hosts it at **T = 8** with a new
whole-owner map. The global code is the borrow-aware code of `research/gmx-t74-borrow43-handoff`, which runs at eight
replicas from its SPEC geometry. **Prepared by G. V. Chudnovsky with Anthropic Claude assistance.** Credits are in
NOTICE.md.

- **κ = 0.001687525611032835.** The complex lower endpoint is 0.001690378167482793552316, and the complex side binds;
  PR589's bit supplier gives 0.001733427039589844234749. The conservative κ is 0.001687525608148634.
- **Gains:** +0.6956% over research/gauge44-slots-mixed-arcflow (0.001675867534384272, #1175), +0.8030% over
  research/gauge44-slots-arcflow (0.001674083482880778, same word without borrowing), and +0.8887% over the selected
  record (#1094, 0.001672660197076918).
- **Local word.** PR2226 aliases 80 retired helper registers. This package keeps the **47 "full" aliases**: every alias
  that is not a low (cold-lender) alias and whose lender is a helper register. The rule drops the 32 low aliases and the
  source-lender alias of register 653. Each kept alias is a rank-zero change of logical reference K_G of its lender
  register. The 47 removed helpers are registers 595–703. Results:
  - R = 480 → **433** helpers, widths `{4:56,5:32,7:30,8:122,9:240}` → `{4:37,5:4,7:25,8:127,9:240}`;
  - local paid profile `{1:2412,2:418,3:432,4:142,5:82,8:8}` (mass 5586) → `{1:2424,2:418,3:427,4:140,5:42,8:8}` (mass 5375);
  - 73 private FULL-frame zero buffers in the inverse cleanup, charged as physical stock. Physical is 7684 = 7168 data
    + 8 centres + 433 helpers + 73 zero buffers + 2.
  
  None of PR2226's 64 initial corrections involves only kept aliases, so the word has none. With the same aggregate
  owner model, all 80 aliases price at about 0.0016782 and the 47 full aliases at 0.0016905 (complex lower endpoint).
- **Owner map.** 31,176 whole-helper roles (72 calls × 433 helpers) on arcflow's unchanged 72-call schedule, owner
  windows and cap-8/16/32/48 idle stock. The idle remainder is `{1:2664, 3:989, 5:515, 7:134, 48:2595}`, and the global
  profile is `{1:177192, 2:30096, 3:31733, 4:10080, 5:3539, 7:134, 8:576, 48:2595}`. Rank mass is 520704 and the deficit
  2560, exactly as before.
  - The windowed aggregate LP bound is `{1:2664, 3:988.05, 5:519.75, 7:124.5, 48:2595.95}` (complex endpoint
    0.0016904699). The shipped map is the integer optimum of the pattern-time model on the LP support. Widening the
    support to every pattern with remainder ≤ 3 gives the same optimum.
  - Unlike arcflow's word, the windows bind here. The time-free relaxation would allow about 0.0016917.
  - The strict assembly prints a `scope` string inherited from its code ("The owner map attains the fixed-inventory
    bin bound"). That string does not describe this map; this paragraph states what the map attains.
- **Evidence:** `RESULT.json` is the strict replay's `REPLAY.json`. `evidence/` holds the receipts of a complete fresh
  `verify_borrow47.py` run: borrow-word audit and controls, the five native steps, emit/audit with 19 mutation controls,
  invoice, assembly, 69 admission controls, 12 bit-supplier controls and the independent price.
- **Control:** the same code, run on research/gauge44-slots-arcflow's own data and locals with this package's T = 8
  geometry, reproduces arcflow's pinned κ 0.001674083482880778 exactly (`evidence/CONTROL-ARCFLOW-T8.json`).

## Verification

```sh
python -B research/gauge44-slots-borrow47/verify_borrow47.py --main . --work /tmp/borrow47
# PASS_BORROW47_COMPLETE_T8_VERIFICATION 0.001687525611032835
```

This takes about 5 minutes; the native phase lowering takes about 3 minutes. It needs Python 3.12+ and NumPy, as
`research/gmx-t74-borrow43-handoff` does. Every stage starts from pinned main sources (`construction/borrow/BORROW-PINS.json`,
`data/UPSTREAM-PINS.json`):

1. `construction/borrow/derive_witness.py` runs PR2226's own `replay.py` and `operators.py` on the rule-selected 47
   aliases. The only change is that the whole-witness count assertion (80 pairs, 32 low) is removed, as an exact
   one-line substitution. The output must equal `WITNESS-f47.json` except for timing fields.
2. `construction/borrow/build_borrow47.py` checks that the selected word is arcflow's pinned canonical word. It then runs
   gmx's unchanged `new/borrow/` tools:
   - `build_word.py` (word content `4cf8a00c…` pinned);
   - `audit_local.py --base` (`PASS_INDEPENDENT_EXACT_BORROWING_WORD`);
   - `borrow_controls.py` (every mutant rejected for its stated reason);
   - the PR2478 cut-centre and PR2508 helper-index guards, with their negative controls.

   Finally it binds the word's residual widths to `data/PLAN.json`.
3. gmx's `new/run_native_borrow.py` runs the five native steps (transport, ledgers, independent local replay, wrappers,
   phases with the independent warp-recipe recheck). `--native-from DIR` reuses a completed run of the same word.
4. `construction/owners/make_plan.py` regenerates `data/PLAN.json` from arcflow's PLAN, the native provider and
   `data/OWNERS.json.gz`. `construction/owners/recount_owners.py` (arcflow's) recounts the map using only the standard
   library.
5. `verify_global.py` is arcflow's strict replay driver with the gmx code. It runs restore-parent, bind, emit, audit with
   mutation controls, invoice, assembly, admission controls and bit-supplier controls. Then arcflow's
   `independent/price_profile.py` reprices `data/PLAN.json`.

Re-deriving inputs (not run by the verifier):
- `construction/owners/t8pack.py PLAN OUT` (NumPy, SciPy/HiGHS) regenerates a map at the same remainder cost.
  HiGHS is not deterministic across versions, so the frozen `data/OWNERS.json.gz` is the admitted object.
- `construction/spec/derive_spec.py` repeats discovery → `new/unpin_spec.py` → `freeze_spec.py`. The base SPEC adds
  arcflow's T = 8 parent geometry (11704/12218/16128, 504 parent helpers, payload 512) and outer refinement
  (24/18/28 digits, 12 steps) to arcflow's SPEC.
- `freeze_spec.py` is gmx's top-level freezer with one correction: with a single provider, `code/assembly_controls.py`
  drops its four partner-only substitutions, so the expected admission-control count is 69 instead of 73. `code/` is
  `research/gmx-t74-borrow43-handoff/code` byte for byte; its unused `code/freeze_spec.py` is a stale T74 variant and is
  not run.

## Open next questions

1. **Borrowing on #1175's mixed providers.**
   - *Relevance:* research/gauge44-slots-mixed-arcflow gained about +0.1% over arcflow by mixing two provider words. The
     same aliasing on those words could stack with this package's +0.8%.
   - *First step:* run PR2226's alias search, or gmx's `new/borrow/alias_enumerate.py`, on the two mixed words. Then
     price the full aliases with the aggregate owner model.
   - *Blocker:* no alias witness exists yet for those words, and the two providers must agree on the helper count R.
2. **Windows.**
   - *Relevance:* on this word the owner windows cost about 1.2·10⁻⁶ in the complex endpoint (windowed 0.0016905 vs
     time-free 0.0016917).
   - *First step:* search the 72-call order with the window LP as objective.
   - *Blocker:* a new order changes the windows and the controller chronology, which emit/audit must re-admit.
3. **The 33 dropped aliases.**
   - *Relevance:* the low and source-lender aliases are about break-even one at a time but lose when all 80 are taken.
     Their initial corrections and lender widths interact with the owner packing.
   - *First step:* price small low-alias subsets on top of the 47 with the window LP.
4. **Admission.** Like the other T = 8 packages, this is a conditional finite witness until the maintainers' pinned
   scoped execution.

This is a finite conditional witness. It inherits every assumption of the #1005/#1018/#1094 chain, PR2226 and the gmx
borrow lowering: all-size compiler, common-label/mutable-source realization, whole-call hosting, precision/tape, exact
recovery, analytic transfer and the bit H1–H5 interfaces. No optimality beyond the stated T = 8 whole-helper family on
arcflow's schedule is claimed.
