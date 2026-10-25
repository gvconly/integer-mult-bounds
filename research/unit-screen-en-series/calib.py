"""Calibration: the screener must reproduce the repository's published numbers before it scores anything new.
  format-limits-corrected A.3: bit P111 (x = 63/111, h = 8) 6.0005e-3; complex P112 (64/112, h = 9) 5.1863e-3
  label-field-census: E8 Hamming frame x_lb = 8/15 -> 8.009e-3; E7 x = 7/9 -> 1.847e-3
  e8-decoder-floors-and-cold14-notes: E8 decoder 72 on 120 lines; the paper's triples at h = 24: X* ~ 7.4e-3 and 630
  orthogonal neighbours per line."""
import itertools, time
import numpy as np
from fractions import Fraction as Fr
from screen import screen, show, ceiling
from e8units import hamming_e8
t0 = time.time()
for name, h, x, want in (('bit P111 (63/111, h=8)', 8, Fr(63, 111), 6.0005e-3), ('complex P112 (64/112, h=9)', 9, Fr(64, 112), 5.1863e-3),
                         ('E8 Hamming frame x_lb 8/15', 8, Fr(8, 15), 8.009e-3), ('E7 7/9', 7, Fr(7, 9), 1.847e-3)):
    b = ceiling(h, x)
    print('%-30s X=%.4e@s%d  (published %.4e)' % (name, float(b[0]), b[1], want))
L, G = hamming_e8()
sup = [[i for i in range(120) if int(L[i] @ f) != 0] for f in L[:8]]
show(screen('E8 root lines', G, sup))
for hh in (12, 18, 24):
    T = list(itertools.combinations(range(hh), 3))
    M = np.zeros((len(T), hh), dtype=np.int64)
    for i, t in enumerate(T):
        M[i, list(t)] = 1
    show(screen('paper triples h=%d (G = |S n T| - 1)' % hh, M @ M.T - 1, None, cliques=False))
print('%.0f s' % (time.time() - t0))
