# Recorded search (not run by verify.py)

`pairmix.c` is an exhaustive C search. Its input file `mixdata.txt` is written by `export_mixdata.py`, which runs from the repository root under Python 3.12.

It asks whether one mixed rotation a e_j + b e_k can replace the two pure rotations of directions j and k in an R_8 = 0 centre. Every other direction is left unchanged. The search ranges over:
- every pair j < k;
- every common line l ⊆ h_j^⊥ ∩ h_k^⊥;
- every plane P_j ⊃ l inside h_j^⊥ and every plane P_k ⊃ l inside h_k^⊥;
- every other plane Pq of P_j + P_k.

That is 285,653,480 configurations. For each, it checks the pair-local service conditions of PROOF.md §7: every needer missing x_j sees P_k or misses both, symmetrically for k, every missing needer sees Pq, and all double missers share one ratio.

Recorded result (gcc -O2, about 50 s):
- **no configuration serves every needer**;
- exactly one configuration per pair fails at a single port, for the 10 pairs listed in `checks.MIX_NEAR`. verify.py re-checks those 10 exactly.

```sh
python3.12 -B research/p104-centre-outside-root-code/evidence/export_mixdata.py /tmp/mixdata.txt
gcc -O2 -o /tmp/pairmix research/p104-centre-outside-root-code/evidence/pairmix.c
(cd /tmp && ./pairmix | tail -1)      # configs 285653480 solutions 0 best 1
```

## Chained repairs (`repair_shapes.py`, about 90 s, arithmetic mod 2^61 - 1)

For the four matched pair-mix designs (cost 52, four ports unserved), the script asks whether **three** unit-step repairs can serve those ports, re-using earlier landings.

Recorded result, for all four matchings:
- no vector of Q^8 suits three of the four ports;
- every shape in which each port is served by a single repair is impossible: S1, S211 and S2 (definitions in the script) all give 0.

The check lets an earlier step's content range over its whole source space, so it over-approximates and these zeros are sound. Still open: a port served *jointly* by two repairs, neither suiting it alone.

```sh
python3.12 -B research/p104-centre-outside-root-code/evidence/repair_shapes.py
```
