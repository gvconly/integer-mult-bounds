"""The E_n level-0/1 series T_n: the A_{n-1} roots (e_i - e_j; 0) together with the paper's triples (t_S; 1) in Z^{n,1}
(form x.y - m m').  All lie in rho^perp, rho = (1^n; 3), and are the roots of levels 0 and 1 for the A_{n-1} grading
of the E_n lattice:
  n = 9:  rho is isotropic and rho^perp / rho = E8; T_9 is the whole E8 root system (A8 + the 84 triples of a 9-set),
          the 'Lorentzian E8 chart' noted in research/c5-dyadic-decoder-gate;
  n = 10: rho has norm 1, so rho^perp is even unimodular of signature (9,1), i.e. II_{9,1} (the E10 lattice);
  n > 10: Lorentzian.
The paper's unit is the level-1 part (triples) only.  For each n: v, rank h, label rank c, exact cost of the
coordinate-total decoder over the n+1 coordinates of Z^{n,1}, ceilings with window h (bit) and c (complex), the ceiling
at every odd s in 3..13, orthogonal lines, cover numbers, and (n <= 12) maximality among all roots of height <= 2."""
import itertools, math, sys, json, time
from pathlib import Path
import numpy as np
from fractions import Fraction as Fr
from screen import rank_q, ceiling, S_LIST
from e8units import f2_rank
from skel import X
from covers import cover
from lorentz import roots


def T(n):
    A = []
    for i, j in itertools.combinations(range(n), 2):
        x = [0] * (n + 1); x[i] = 1; x[j] = -1; A.append(x)
    for t in itertools.combinations(range(n), 3):
        x = [0] * (n + 1)
        for a in t:
            x[a] = 1
        x[n] = 1; A.append(x)
    return np.array(A, dtype=np.int64)


def main(ns):
    rows = []
    for n in ns:
        t0 = time.time()
        R = T(n); Q = np.diag([1] * n + [-1]).astype(np.int64); G = R @ Q @ R.T
        v = len(R); off = ~np.eye(v, dtype=bool)
        assert set(np.diag(G)) == {2} and set(G[off].tolist()) <= {-1, 0, 1}
        assert (R @ Q @ np.array([1] * n + [3]) == 0).all()          # all in rho^perp
        h = rank_q(G.tolist()) if v <= 1400 else int(np.linalg.matrix_rank(G.astype(float)))
        B = ((1 + G) % 2).astype(np.uint8); c = f2_rank(B)
        cost = sum(f2_rank(B[np.ix_(S, S)]) for S in (np.flatnonzero(R[:, a]) for a in range(n + 1)))
        x = Fr(cost, v)
        ev = np.linalg.eigvalsh(G.astype(float))
        r = dict(n=n, v=v, h=h, c=c, bound=math.comb(h + 2, 3), cost=cost, x=round(float(x), 4),
                 orth=int(((G == 0) & off).sum(1).max()), negative_eigenvalues=int((ev < -1e-9).sum()))
        r['fill'] = round(v / r['bound'], 3)
        bh, bc = ceiling(h, x), ceiling(c, x)
        r['X_window_h'] = (float(bh[0]), bh[1]); r['X_window_c'] = (float(bc[0]), bc[1])
        r['X_by_s_window_h'] = {s: float(X(s, h, 1, x)['hi']) for s in S_LIST if X(s, h, 1, x)}
        for kind in ('bit', 'complex'):
            cv = cover(R, Q, kind); r['cover_' + kind] = dict(cover=cv['cover'], two_max=cv['two_max'], readers=cv['readers'])
        if n <= 12:
            Rall = roots(n, (0, 1, 2)); ok = (np.abs(Rall @ Q @ R.T) <= 1).all(1)
            inR = {tuple(y) for y in R} | {tuple(-y) for y in R}
            r['roots_height_le2_compatible_with_all'] = int(sum(1 for k in np.flatnonzero(ok) if tuple(Rall[k]) not in inR))
        rows.append(r)
        print(json.dumps(r), '%.0fs' % (time.time() - t0), flush=True)
    return rows


if __name__ == '__main__':
    ns = [int(a) for a in sys.argv[1:]] or [9, 10, 11, 12, 14, 16, 18, 20, 22, 24]
    rows = main(ns)
    if not sys.argv[1:]:          # only the full default run refreshes the shipped receipt
        json.dump(rows, open(Path(__file__).resolve().parent / 'en_series.json', 'w'), indent=1)
