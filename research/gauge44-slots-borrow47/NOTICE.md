# Attribution and scope

Apache-2.0 (LICENSE). No upstream license, credit, source pin or AI disclosure is replaced by this notice.

**Anticipatory dirty borrowing is PR2226** (`research/anticipatory-dirty-borrowing`) by Cinna Davis (@1337isnot1337),
with substantial OpenAI Codex assistance. PR2226 contributed:
- the alias rule and the changed-reference K_G handoffs;
- the 80-alias witness on the selected word, the early cold tail, exact cold restoration and the private zero buffers;
- `replay.py` and `operators.py`, which this package runs from main with one count assertion removed.

PR2226 credits the selected PR1042/1094 word and its lineage, Ali Romman's fixed-block certificate and Thomas
DiFiore's block discovery.

**BORROW-WORD tooling, native borrow chain and borrow-aware global code** come from `research/gmx-t74-borrow43-handoff`,
used unchanged from main: `new/borrow/`, `new/run_native_borrow.py`, `new/code/borrow/`, `new/unpin_spec.py`, and
`code/`, copied byte for byte. Its NOTICE credits:
- the T74 lineage and Rohan Arun, with OpenAI Codex assistance;
- the PR2478 cut-centre guard (mdumitrean) and the PR2508 helper-index guard (ywz65);
- the PR1956 carried-cutoff bootstrap repair (Josie Elliston);
- the PR1377 grid refinement (Jeremy Kothe) and fine_cache (jeremykothe2);
- the native lineage named in `new/code/borrow/NOTICE-PR2226.md`.

The top-level `freeze_spec.py` is that package's top-level freezer, with the single-provider control count corrected.

**The schedule, owner windows, idle stock, PLAN skeleton, strict replay driver, recount and independent pricer** are
those of `research/gauge44-slots-arcflow`. Its NOTICE.md carries the #645/#1005/#1018 lineage, rohanarun's #987 slots
and G44Pack's map ancestry. The PR589 bit supplier is unchanged.

**New in this package**, prepared by G. V. Chudnovsky with Anthropic Claude assistance:
- the 47-alias selection rule and `construction/borrow/derive_witness.py`;
- `build_borrow47.py`, the T = 8 whole-owner map and `construction/owners/t8pack.py`;
- `make_plan.py`, `construction/spec/derive_spec.py` and `verify_borrow47.py`;
- the T = 8 SPEC geometry binding of the borrow-aware code.

All-size, realization, hosting, precision/tape, analytic-transfer and bit H1–H5 assumptions remain conditional. No
independent human review is claimed.
