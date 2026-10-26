"""Source-bound geometry of the P104 centre, recomputed from the published word.

Copied unchanged (apart from this note) from research/p104-multigen-centre-floor/common.py (same author), so that this
package replays on its own.

Input (live, unchanged): research/p104-productive-controller/WORD.json.gz (schema p104-productive-source/1).

From the word alone this module recovers
  * the eight owned centre lifetimes: copy frame F_j, descent chain, landing line [g_j], content u_j on the ports;
  * an integral E8 realisation of the 104 ports: u_j(s) = <x_j, r_s>/4 with the x_j a Z-basis of E8 (det Gram = 1),
    unique up to scale from the orthogonality pattern of the port labels;
  * all 120 E8 lines, the 16 removed lines (two orthogonal frames) and their F2 labels.
All arithmetic is exact (integers and Fractions)."""
import gzip
import json
from fractions import Fraction as Fr
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
WORD = 'research/p104-productive-controller/WORD.json.gz'


def par(x):
    return x.bit_count() & 1


def load_word(rel=WORD):
    return json.loads(gzip.decompress((ROOT / rel).read_bytes()))


def f2span(vs):
    basis = []
    for v in vs:
        for b in basis:
            v = min(v, v ^ b)
        if v:
            basis.append(v)
    return basis


def centre(W):
    """Replay the word; return the eight owned lifetimes with their port contents."""
    fr = W['frames']
    live = list(W['live'])
    idx = {r: i for i, r in enumerate(live)}
    rows = [{r: Fr(1)} for r in range(743)]
    cur = list(W['initial']) + [None] * 8
    life = {}
    for t, e in enumerate(W['program']):
        op, r = e[0], e[1]
        if op == 'MOVE':
            _, r, a, b, private = e
            if r >= 735:
                assert private and len(fr[b]) < len(fr[a])
                L = life[r]
                content = {idx[s]: c for s, c in rows[r].items() if s in idx}
                L['steps'].append((len(fr[a]), len(fr[b])))
                L['contents'].append(content)
                L['landing'] = fr[b]
            cur[r] = b
        elif op == 'INIT':
            cur[r] = e[2]
            rows[r] = {}
            life[r] = {'init_frame': fr[e[2]], 'steps': [], 'contents': []}
        elif op == 'CLEAR':
            cur[r] = None
            rows[r] = {}
        elif op == 'SCALE':
            c = Fr(e[-2])
            rows[r] = {k: c * v for k, v in rows[r].items()}
        else:
            s, c = e[2], Fr(e[3])
            if op == 'COPY':
                rows[r] = {k: c * v for k, v in rows[s].items()}
                life[r]['copy_frame'] = fr[e[-1]]
            else:
                for k, v in rows[s].items():
                    nv = rows[r].get(k, 0) + c * v
                    if nv:
                        rows[r][k] = nv
                    else:
                        rows[r].pop(k, None)
    out = []
    for r in sorted(life):
        L = life[r]
        assert all(c == L['contents'][0] for c in L['contents']), 'centre content changes during descent'
        out.append({'register': r, 'copy_frame': L['copy_frame'], 'landing': L['landing'], 'steps': L['steps'],
                    'cost': sum(a - b for a, b in L['steps']), 'u': L['contents'][0]})
    return out


def _nullspace(rows, n):
    rows = [list(r) for r in rows]
    piv, r = [], 0
    for c in range(n):
        p = next((i for i in range(r, len(rows)) if rows[i][c] != 0), None)
        if p is None:
            continue
        rows[r], rows[p] = rows[p], rows[r]
        pv = rows[r][c]
        rows[r] = [x / pv for x in rows[r]]
        for i in range(len(rows)):
            if i != r and rows[i][c] != 0:
                f = rows[i][c]
                rows[i] = [a - f * b for a, b in zip(rows[i], rows[r])]
        piv.append(c)
        r += 1
    basis = []
    for f in (c for c in range(n) if c not in piv):
        v = [Fr(0)] * n
        v[f] = Fr(1)
        for i, c in enumerate(piv):
            v[c] = -rows[i][f]
        basis.append(v)
    return basis


def _inverse(A):
    n = len(A)
    M = [[Fr(x) for x in row] + [Fr(int(i == j)) for j in range(n)] for i, row in enumerate(A)]
    for c in range(n):
        p = next(r for r in range(c, n) if M[r][c] != 0)
        M[c], M[p] = M[p], M[c]
        pv = M[c][c]
        M[c] = [x / pv for x in M[c]]
        for r in range(n):
            if r != c and M[r][c] != 0:
                f = M[r][c]
                M[r] = [a - f * b for a, b in zip(M[r], M[c])]
    return [row[n:] for row in M]


def _det(A):
    M = [[Fr(x) for x in row] for row in A]
    n, d = len(M), Fr(1)
    for c in range(n):
        p = next((i for i in range(c, n) if M[i][c] != 0), None)
        if p is None:
            return Fr(0)
        if p != c:
            M[c], M[p] = M[p], M[c]
            d = -d
        d *= M[c][c]
        for i in range(c + 1, n):
            f = M[i][c] / M[c][c]
            M[i] = [a - f * b for a, b in zip(M[i], M[c])]
    return d


class Geometry:
    """Integral E8 coordinates in the basis of the eight centre contents."""

    def __init__(self, W=None):
        W = W or load_word()
        self.W = W
        self.Q = list(W['ports'])
        self.P = len(self.Q)
        self.lives = centre(W)
        assert len(self.lives) == 8
        P, Q = self.P, self.Q
        # R[s][j] = 4 u_j(s) in {0, +-1}
        R = [[self.lives[j]['u'].get(s, Fr(0)) * 4 for j in range(8)] for s in range(P)]
        assert all(x in (-1, 0, 1) for row in R for x in row)
        idx = [(a, b) for a in range(8) for b in range(a, 8)]
        eqs = [[R[s][a] * R[t][b] + (R[s][b] * R[t][a] if a != b else 0) for a, b in idx]
               for s in range(P) for t in range(s + 1, P) if par(Q[s] & Q[t])]
        ns = _nullspace(eqs, len(idx))
        self.gram_solution_dimension = len(ns)
        assert len(ns) == 1, 'orthogonality pattern must fix the Gram up to scale'
        Y = [[Fr(0)] * 8 for _ in range(8)]
        for (a, b), v in zip(idx, ns[0]):
            Y[a][b] = Y[b][a] = v
        n0 = sum(R[0][a] * Y[a][b] * R[0][b] for a in range(8) for b in range(8))
        Y = [[y * 2 / n0 for y in row] for row in Y]
        Xf = _inverse(Y)
        assert all(x.denominator == 1 for row in Xf for x in row)
        self.X = [[int(x) for x in row] for row in Xf]
        self.det = _det(self.X)
        assert self.det == 1
        roots = []
        for s in range(P):
            c = [sum(Y[a][b] * R[s][b] for b in range(8)) for a in range(8)]
            assert all(x.denominator == 1 for x in c)
            roots.append(tuple(int(x) for x in c))
        self.roots = roots
        assert all(self.ip(r, r) == 2 for r in roots)
        # check u_j(s) = <x_j, r_s>/4, x_j = e_j
        for j in range(8):
            e = tuple(int(i == j) for i in range(8))
            assert all(self.ip(e, roots[s]) == R[s][j] for s in range(P))
        self._close()

    def ip(self, a, b):
        X = self.X
        return sum(a[i] * X[i][j] * b[j] for i in range(8) if a[i] for j in range(8) if b[j])

    @staticmethod
    def line(v):
        for x in v:
            if x:
                return tuple(v) if x > 0 else tuple(-y for y in v)
        raise ValueError('zero vector')

    def _close(self):
        allr = set()
        for r in list(self.roots) + [tuple(int(i == j) for i in range(8)) for j in range(8)]:
            allr.add(r)
            allr.add(tuple(-x for x in r))
        changed = True
        while changed:
            changed = False
            cur = list(allr)
            for a in cur:
                for b in cur:
                    if self.ip(a, b) == -1:
                        c = tuple(x + y for x, y in zip(a, b))
                        if c not in allr:
                            allr.add(c)
                            changed = True
        assert len(allr) == 240
        self.lines = sorted({self.line(v) for v in allr})
        live = {self.line(r) for r in self.roots}
        self.removed = [v for v in self.lines if v not in live]
        assert len(self.lines) == 120 and len(self.removed) == 16
        self.removed_labels = [self.label(v) for v in self.removed]
        # label dictionary: <q_p, q_s> = 1 iff r_p orthogonal to r_s, or p = s
        Q = self.Q
        assert all(par(q) == 1 for q in Q) and len(set(Q)) == self.P
        assert all(par(Q[a] & Q[b]) == int(a == b or self.ip(self.roots[a], self.roots[b]) == 0)
                   for a in range(self.P) for b in range(self.P))

    def label(self, rho):
        """The F2 label l with <l, q_p> = [rho orthogonal to r_p] for every port (unique: the labels span F2^9)."""
        piv = []
        for p in range(self.P):
            m, b = self.Q[p], int(self.ip(rho, self.roots[p]) == 0)
            for pm, pb, pc in piv:
                if m >> pc & 1:
                    m ^= pm
                    b ^= pb
            if m:
                c = m.bit_length() - 1
                piv = [(pm ^ m, pb ^ b, pc) if pm >> c & 1 else (pm, pb, pc) for pm, pb, pc in piv]
                piv.append((m, b, c))
            else:
                assert b == 0, 'no consistent label'
        assert len(piv) == 9
        lab = 0
        for pm, pb, pc in piv:
            if pb:
                lab |= 1 << pc
        return lab
