"""The polynomial-method bound for compatible units, and what it allows on the ceiling yardstick.

Claim.  Let G be a real symmetric v x v matrix of rank h with G_ii = 2 and G_ij in {-1, 0, 1} for i != j (any
signature).  Then v <= C(h+2, 3).
Proof.  Write G = R Q R^T with R real v x h and Q = diag(+-1) (spectral decomposition, keep the h nonzero
eigenvalues).  Put f(t) = t^3 - t, which vanishes on {-1, 0, 1} and has f(2) = 6.  Because r_i^T Q r_i = 2 for every
i, the entrywise matrix F_ij = G_ij^3 - G_ij G_ii G_jj / 4 equals f(G_ij), so F = 6 I has rank v.  Each entry of F is
a form of degree 3 in r_i and degree 3 in r_j, so rank F <= C(h+2, 3).
E8 attains the bound (120 = C(10, 3); its root lines are a tight spherical 7-design).  The same argument with the
odd polynomial vanishing on {0, +-1, +-2} (norm 4) gives C(h+4, 5), which the Leech lattice attains (98280 = C(28, 5)).

This script checks F = 6I on E8, prints the ceiling of a unit sitting exactly at the bound with label rank c = h
(E8, the triples and the series all have c >= h; units with c < h would score higher on the proxy), and the number
of lines a rank-h unit with c = h + 1 needs to beat the best series member (T_20, 8.17e-3)."""
import math
import numpy as np
from fractions import Fraction as Fr
from screen import ceiling
from e8units import hamming_e8

L, G = hamming_e8()
G = G.astype(np.int64)
d = np.diag(G)
F = G ** 3 - G * np.outer(d, d) // 4
assert (F == 6 * np.eye(120, dtype=np.int64)).all()
print('E8: F = G^3 - G (d d^T)/4 = 6 I checked; v = 120 = C(10,3), so E8 attains the bound')
print('\n h   v_max   x_min    cap X  (unit at the bound, c = h)')
for h in (8, 9, 10, 11, 12, 13, 14, 16, 18, 20, 22, 24, 28, 32):
    vmax = math.comb(h + 2, 3); x = Fr(h * h, vmax); b = ceiling(h, x)
    print('%2d %7d  %.4f  %.3e' % (h, vmax, float(x), float(b[0])), flush=True)
print('\nlines needed at rank h (c = h + 1) to beat 8.17e-3:')
target = Fr(817, 100000)
for h in (9, 10, 11, 12, 13, 14, 16):
    bound = math.comb(h + 2, 3)
    top = ceiling(h, Fr(h * (h + 1), bound))
    if not top or top[0] < target:
        print('h=%2d  bound=%4d  even at the bound X=%.3e' % (h, bound, float(top[0]) if top else 0), flush=True)
        continue
    lo, hi = h * (h + 1), bound
    while hi - lo > 1:
        mid = (lo + hi) // 2
        b = ceiling(h, Fr(h * (h + 1), mid))
        if b and b[0] >= target:
            hi = mid
        else:
            lo = mid
    print('h=%2d  bound=%4d  need v >= %4d (fill %.3f)' % (h, bound, hi, hi / bound), flush=True)
