"""Cover numbers of a source's readers by hyperplane sections, both reader types:
  bit-type readers     = lines orthogonal to the source (format-limits D: E8 -> E7's 63 lines, cover 3)
  complex-type readers = lines NOT orthogonal to the source (format-limits D: P112 -> 52 readers, cover 2)
Candidate hyperplanes w^perp: the configuration's lines, the coordinate vectors e_i and the differences e_i - e_j.
Greedy set cover gives an upper bound on the cover number; 'two_max' is the largest union of two sections among the
300 largest candidate sections (so two_max < readers means these candidates need at least three)."""
import itertools, sys, json
import numpy as np
import screen  # noqa: F401  (puts the repository engine on sys.path)
from e8units import hamming_e8

def cands(L, n):
    E = np.eye(n, dtype=np.int64)
    D = np.array([E[i] - E[j] for i, j in itertools.combinations(range(n), 2)], dtype=np.int64)
    return np.vstack([L, E, D])

def cover(L, Q, kind, src=0):
    ip = (L @ Q @ L[src]).ravel()
    sel = (ip == 0) if kind == 'bit' else (ip != 0)
    sel[src] = False
    R = L[sel]; C = cands(L, L.shape[1])
    Z = (R @ Q @ C.T == 0)
    Z = Z[:, Z.sum(0) < len(R)]                       # proper sections only
    rem = np.ones(len(R), dtype=bool); sizes = []
    while rem.any() and len(sizes) < 30:
        g = (Z & rem[:, None]).sum(0); j = int(g.argmax())
        if g[j] == 0: sizes.append('uncoverable:%d' % rem.sum()); break
        sizes.append(int(g[j])); rem &= ~Z[:, j]
    cnt = Z.sum(0); top = np.argsort(-cnt)[:300]
    two = max(int((Z[:, a] | Z[:, b]).sum()) for a in top for b in top)
    return dict(readers=int(len(R)), largest=int(cnt.max()), greedy=sizes, cover=len(sizes), two_max=two)

if __name__ == '__main__':
    out = {}
    L, G = hamming_e8(); L = L.astype(np.int64); Q = np.eye(8, dtype=np.int64)
    for kind in ('bit', 'complex'):
        r = cover(L, Q, kind); out['E8 ' + kind] = r; print('E8      %-7s %s' % (kind, r), flush=True)
    for hh in (16, 20, 24):
        T = list(itertools.combinations(range(hh), 3)); M = np.zeros((len(T), hh + 1), dtype=np.int64)
        for i, t in enumerate(T): M[i, list(t)] = 1
        M[:, hh] = 1; Qt = np.diag([1] * hh + [-1]).astype(np.int64)
        for kind in ('bit', 'complex'):
            r = cover(M, Qt, kind); out['triples%d %s' % (hh, kind)] = r; print('triples%d %-7s %s' % (hh, kind, r), flush=True)
