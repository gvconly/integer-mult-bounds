# Independent check of #5440 (p19 coded-reverse chain, κ candidate 0.0033999964)

**Finding: the shipped replay does not tie the coded chain to the invoiced controller, and with the invoiced
semantics the coded chain is not exact.** No new κ. Prepared by G. V. Chudnovsky with Anthropic Claude assistance;
Apache-2.0. This branch is #5440's head (`c0f52b5b1`, on main `32890bc5d`) plus this directory; run from its root:

```sh
python3 -B research/pr5440-independent-check/validate.py .           # simulator == PR5440's own physical_code / signed_map / gaussian / replay (p=2,3)
python3 -B research/pr5440-independent-check/their_code_checks.py .  # uses only PR5440's modules: letter D, and standard word vs controller()
python3 -B research/pr5440-independent-check/formal_check.py .       # which listing semantics reproduce controller() (collapsed and fresh-block)
python3 -B research/pr5440-independent-check/coded_formal.py .       # coded listing with controller()'s semantics vs controller()
python3 -B research/pr5440-independent-check/freshtest.py .          # coded listing vs standard, collapsed vs fresh-block geometry
python3 -B research/pr5440-independent-check/enum128.py .            # p=2 fresh-block: 128 block assignments of the code gates (~45 min)
python3 -B research/pr5440-independent-check/tolerance.py .          # extra rank-one jobs the p=19 invoice can absorb
```

Each shipped `.log` is the stdout of the script of the same name; `enum128.py` and `tolerance.py` also write their
`.json`. Standard library only; arithmetic mod 65537 with i = 256, as in #5440. `fresh.py` simulates the listing in
the Walsh chart: Gaussian gates are diagonal phases there, a normalizer K = Z_z P_G moves frequency u to Gu ⊕ z, and
the mode code mixes owners at a common frequency. `validate.py` checks every primitive and the whole p=2 and p=3
replays against #5440's address-domain code, entry for entry.

## 1. The 0-residual replay compares two copies of a different operator

`p19_word.py` compares the coded listing with `standard()`. Both use the same three conventions, and each differs
from `coded-reverse-controller`'s `controller()`, the operator whose diagonal census the invoice prices and which
`p19_bank.py` evaluates:

| | listing and `standard()` | `controller()` and the cited primitive |
|---|---|---|
| reverse-stage letter | R = D·Q on both roles | D only (`X ← D(X − QCY)`, `Y ← DY`, PHYSICAL.md) |
| idles | R^k on all nine axes, including the thin middle idles | D over the selected retained stages; middle idles thinned to 1 |
| per-pair odd bridges | `role 2j+1 += ∓ role 2j+3` | `role 2j+3 += ∓ role 2j+1` |

The table's claims were checked in two ways:
- **With #5440's own modules** (p=2, collapsed nine-axis geometry), `standard()` differs from `controller()`,
  evaluated with `p19_bank.py`'s monomial rule, on **229,372 of 229,376** entries.
- **With `controller()`'s conventions put back**, the uncoded listing reproduces `controller()` exactly. The
  residual is 0 at p=2 and p=3, in both geometries.

So the gate structure itself is fine. The deviations are the letter, the idles and the bridges.

## 2. With the invoiced letter D the coded listing is not exact

The listing with D letters (#5440's own `replay(..., letter='D')`) differs from its standard word on **114,688**
entries at p=2. With all three conventions of `controller()`, the coded listing differs from `controller()` by
4,804 entries (p=2) and 11,015 (p=3) in the collapsed geometry (`coded_formal.py`). This is the commutation fact
`physical_gauge_lift.py` itself records: the code commutes with R = DQ, but not with D.

## 3. The invoiced fresh-block geometry

The invoice prices central axis + one 8-axis block per retained stage (ambient rank 153). In that geometry the
coded chain was only ever replayed in collapsed form: `p19_word.py` puts every stage on one nine-axis field, and
`p19_bank.py` checks the *uncoded* trace. In fresh-block geometry the coded listing fails for every policy tested:
- codes on each gate's own block, or all codes on block 0;
- #5440's letter and idle conventions, or `controller()`'s.

At p=2 the two reverse-stage codes were fixed to their own stage's block, since the native call at stage j carries
stage j's normalizers. The other seven code gates (two loads, three forward reads, two unloads) were run over both
blocks, 128 assignments, with a 30-second cap per run (`ENUM128.json`):

| semantics | assignments completed | exact |
|---|---|---|
| `controller()`'s (letter D, formal idles, formal bridges) | 102 of 128 (26 hit the cap) | **0** of 102; minimum residual 3,741 |
| #5440's own (letter R, nine-axis idles, its bridges) | 128 of 128 | **0** of 128; minimum residual 4,436 |

The 26 capped runs are undecided, not shown inexact. This is consistent with `research/fresh-stage-code-separation`.

## 4. The cost of the missing commutation exceeds the margin

To obtain R from the primitive's D, each reverse stage needs a central Q on both roles. That is 2 × 112 = 224
rank-one Gaussian jobs, the same count `research/physical-inverse-pair` derives for its exit Q⁻¹CQ. At p=19 it
comes to **4,256 per call**. The p=19 invoice (`p19_word.invoice`, profile P_NEW) has a total rank-mass surplus of
**2,016**. Its root stays ≥ 0.0034 only up to **1,103** extra rank-one jobs, and above the selected 0.003027288235
only up to **1,202**. With 4,256 extra jobs the root is 0 (`tolerance.py`, `TOLERANCE.json`).

Even if the Q jobs were paid, the R-letter chain computes a different operator from the one whose census is
priced (section 1).

## 5. Also

At #5440's head, `verify.py` exits with `ValueError: both gains are recorded`. `P19_CLAIM.json` records a 12.312%
gain against the newly selected #5413, while `verify.py` still requires more than 12.5%.

## What would settle it

A coded realization of `controller()` itself (letter D, fresh blocks) whose exactness is shown in the invoiced
geometry, with every Q/Q⁻¹ job and recoding charged. Alternatively, an exact gauge repair costing fewer than about
1,100 rank-one jobs per p=19 call.
