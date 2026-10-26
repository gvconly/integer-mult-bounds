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
