"""Unit screener: score any finite line configuration on the community's unit yardstick.

A unit is a set of v lines with an integer Gram matrix G (after an optional shift G - gamma*J, as in the upstream
triple construction, whose bit labels use the form I - J/9).  Following research/e8-decoder-floors-and-cold14-notes
(ade_screen.py) and research/label-field-census:
  h      window = rank_Q(G)
  B      label pairing (J + G) mod 2, which must have unit diagonal (labels nondegenerate)
  c      rank_F2(B)
  x_prx  c*h/v                         (ade_screen proxy)
  x_frm  sum_f rank_F2 B[S_f, S_f] / v (exact cost of the coordinate-total decoder on an orthogonal frame f, when one
                                       is supplied and the Gram identity G = sum_f (Lf)(Lf)^T/|f|^2 holds)
  x_lb   sum_f clique(S_f) / v         (census lower bound for that decoder; needs cliques, small v only)
  orth   number of lines orthogonal to a line (helper-cost driver: E8 63, triples at h=24: 630)
and the certified B_K format ceiling X(s, h, 1, x) of research/format-limits-corrected (skel.X), maximised over odd
s in 3..13.  A ceiling is an upper bound for the old B_K format, used here only as the common relative yardstick.
"""
import sys
from fractions import Fraction as Fr
from pathlib import Path
import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / 'research/format-limits-corrected/scripts'))
from skel import X                      # noqa: E402
from e8units import f2_rank, first_clique  # noqa: E402

S_LIST = (3, 5, 7, 9, 11, 13)


def ceiling(h, x, k=1):
    """best certified upper end of the root over S_LIST, and the s attaining it (None if no positive saving)."""
    best = None
    for s in S_LIST:
        r = X(s, h, k, Fr(x))
        if r is None:
            continue
        if best is None or r['hi'] > best[0]:
            best = (r['hi'], s)
    return best


def rank_q(G):
    """exact rank over Q of an integer matrix (fraction-free Bareiss on Python ints)."""
    A = [[int(a) for a in row] for row in G]
    n, m = len(A), len(A[0]) if A else 0
    r, prev = 0, 1
    for col in range(m):
        piv = next((i for i in range(r, n) if A[i][col]), None)
        if piv is None:
            continue
        A[r], A[piv] = A[piv], A[r]
        for i in range(r + 1, n):
            if A[i][col]:
                a, b = A[r][col], A[i][col]
                A[i] = [(a * A[i][j] - b * A[r][j]) // prev for j in range(m)]
            else:
                A[i] = [(A[r][col] * A[i][j]) // prev for j in range(m)]
        prev = A[r][col]
        r += 1
        if r == n:
            break
    return r


def max_clique(adj, upto=40):
    k = 1
    while k < upto and first_clique(adj, k + 1):
        k += 1
    return k


def screen(name, G, frame_supports=None, cliques=True, note=''):
    G = np.asarray(G, dtype=np.int64)
    v = len(G)
    off = ~np.eye(v, dtype=bool)
    B = ((1 + G) % 2).astype(np.uint8)
    out = dict(name=name, v=v, note=note)
    out['diag_ok'] = bool((np.diag(B) == 1).all())
    h = rank_q(G.tolist()) if v <= 3000 else int(np.linalg.matrix_rank(G.astype(float)))
    c = f2_rank(B)
    orth = (G == 0) & off
    od = orth.sum(1)
    out.update(h=h, c=c, orth_min=int(od.min()), orth_max=int(od.max()),
               inner=sorted(set(G[off].tolist())), norms=sorted(set(np.diag(G).tolist())))
    xs = {'prx': Fr(c * h, v)}
    if frame_supports is not None:
        cost = sum(f2_rank(B[np.ix_(S, S)]) for S in frame_supports)
        xs['frm'] = Fr(cost, v)
        out['frame_cost'] = cost
        if cliques:
            adj = ((G != 0) & off)
            lb = 0
            for S in frame_supports:
                sub = adj[np.ix_(S, S)].tolist()
                lb += max_clique(sub)
            xs['lb'] = Fr(lb, v)
    out['x'] = {k: str(val) for k, val in xs.items()}
    out['X'] = {}
    for k, val in xs.items():
        b = ceiling(h, val)
        out['X'][k] = None if b is None else (float(b[0]), b[1])
    return out


def show(o):
    xs = '  '.join('x_%s=%.4f X=%s' % (k, float(Fr(v)), ('%.3e@s%d' % o['X'][k]) if o['X'][k] else 'none')
                   for k, v in o['x'].items())
    print('%-34s v=%5d h=%3d c=%3d orth=%d..%d diag_ok=%s  %s %s' % (
        o['name'], o['v'], o['h'], o['c'], o['orth_min'], o['orth_max'], o['diag_ok'], xs, o['note']), flush=True)
