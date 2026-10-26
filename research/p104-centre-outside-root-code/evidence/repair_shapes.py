#!/usr/bin/env python3
"""RECORDED evidence (not run by verify.py; about 15-20 minutes): chained repairs of total cost 3 for the four matched
pair-mix designs of PROOF.md section 7.

Each matched design costs 52 and leaves four ports unserved.  A repair is a unit step: source frame F, landing G of
codimension 1, content in the span of contents landed at frames inside F.  Linear algebra is mod p = 2^61 - 1, on the
four unserved ports only.

1. No vector of Q^8 suits three of the four ports.  A vector suits port t if it lies in S_t + Q r_t and outside S_t.
2. Standalone unit steps never serve two ports (verify.py, mixing_census).

So with three unit steps, if every port is fixed by a single step, some step serves two ports by re-using an earlier
landing G_z inside its source.  The shapes of total cost 3 are:
   S1   z -> x serves {a,b}, z -> y serves {c,d}
   S211 z serves one port c alone, z -> x serves {a,b}, one more step for the last port
   S2   z -> x serves {a,b} -> y serves {c,d} from a source containing G_x
The search lets the earlier content range over its whole source content space, an over-approximation.  So "0" is a sound
negative for these shapes.  Not covered: a port fixed jointly by two steps, neither suiting it alone.

  python3.12 -B research/p104-centre-outside-root-code/evidence/repair_shapes.py
Recorded output: for all four matchings, "suits three: False | S1 0 S211 0 S2 0"."""
import itertools
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(HERE))
from common import Geometry, par, f2span  # noqa: E402
from checks import MIX_NEAR, span_elems  # noqa: E402

p = 2 ** 61 - 1
G = Geometry()
P, Q = G.P, G.Q
roots = [tuple(int(c) for c in r) for r in G.roots]


def red(vs):
    out = []
    for v in ([x % p for x in w] for w in vs):
        for (c, b) in out:
            if v[c]:
                f = v[c]
                v = [(a - f * bb) % p for a, bb in zip(v, b)]
        c = next((i for i in range(8) if v[i]), None)
        if c is None:
            continue
        inv = pow(v[c], p - 2, p)
        v = [x * inv % p for x in v]
        out = [(cc, [(a - bb[c] * x) % p for a, x in zip(bb, v)] if bb[c] else bb) for (cc, bb) in out]
        out.append((c, v))
    return [b for (_, b) in out]


def inside(v, B):
    return len(red(B + [v])) == len(red(B))


def intersect(A, B):
    if not A or not B:
        return []
    M = A + [[(-x) % p for x in b] for b in B]
    n = len(M)
    rows = [[M[i][c] for i in range(n)] for c in range(8)]
    piv, r = [], 0
    for c in range(n):
        q = next((i for i in range(r, 8) if rows[i][c]), None)
        if q is None:
            continue
        rows[r], rows[q] = rows[q], rows[r]
        inv = pow(rows[r][c], p - 2, p)
        rows[r] = [x * inv % p for x in rows[r]]
        for i in range(8):
            if i != r and rows[i][c]:
                f = rows[i][c]
                rows[i] = [(a - f * b) % p for a, b in zip(rows[i], rows[r])]
        piv.append(c)
        r += 1
    vecs = []
    for fc in (c for c in range(n) if c not in piv):
        co = [0] * n
        co[fc] = 1
        for i, c in enumerate(piv):
            co[c] = (-rows[i][fc]) % p
        vecs.append([sum(co[i] * A[i][d] for i in range(len(A))) % p for d in range(8)])
    return red(vecs)


def elems(b):
    return frozenset(span_elems(f2span(list(b))))


def wok(t, Yb):
    return all(par(Q[t] & v) == 0 for v in Yb)


def perp_to(key, h):
    return all(par(e & h) == 0 for e in key)


def main():
    near = {}
    for (j, k, l, a, b) in MIX_NEAR:
        need_j = {t for t in range(P) if roots[t][j]}
        Wj = {t for t in range(P) if par(Q[t] & l) == 0 and par(Q[t] & a) == 0}
        D = sorted(need_j - Wj)
        rat = {}
        for t in D:
            key = (roots[t][j] * roots[t][k] > 0)
            rat.setdefault(key, []).append(t)
        sign = 1 if len(rat.get(True, [])) == 27 else -1
        near[(j, k)] = ([l, a], [l, b], [a, b], sign)
    matchings = [c for c in itertools.combinations(sorted(near), 4) if len({x for e in c for x in e}) == 8]
    for M in matchings:
        frames = []
        for (j, k) in M:
            Pj, Pk, Pq, sign = near[(j, k)]
            frames.append((Pj, [int(i == j) for i in range(8)]))
            frames.append((Pk, [int(i == k) for i in range(8)]))
            v = [0] * 8
            v[j], v[k] = 1, sign
            frames.append((Pq, v))
        S, bad = {}, []
        for t in range(P):
            St = red([v for (fb, v) in frames if wok(t, f2span(fb))])
            if not inside(list(roots[t]), St):
                bad.append(t)
                S[t] = St

        def cons(key):
            return red([v for (fb, v) in frames if elems(fb) <= key])

        def suits(V, T):
            A = V
            for t in T:
                A = intersect(A, S[t] + [list(roots[t])])
                if not A:
                    return False
            return all(not all(inside(a, S[t]) for a in A) for t in T)
        full = red([[int(i == j) for i in range(8)] for j in range(8)])
        triple = any(suits(full, T) for T in itertools.combinations(bad, 3))
        sums = {}
        for r in range(1, len(frames) + 1):
            for sub in itertools.combinations(range(len(frames)), r):
                Xb = f2span([x for i in sub for x in frames[i][0]])
                sums.setdefault(elems(Xb), Xb)
        H = {(a, b): Q[a] ^ Q[b] for a, b in itertools.permutations(bad, 2)}

        def landings(Xb, h):
            key = elems(Xb)
            Gin = f2span([e for e in key if par(e & h) == 0])
            if len(Gin) == len(Xb) - 1:
                return [Gin]
            if len(Gin) == len(Xb):
                out, seen = [], set()
                for w in range(1, 512):
                    if any(par(w & v) for v in Xb):
                        Gb = f2span([e for e in key if par(e & w) == 0])
                        if elems(Gb) not in seen:
                            seen.add(elems(Gb))
                            out.append(Gb)
                return out
            return []

        def two_fix_after(Gz, Cz, pair):
            a, b = pair
            h = H[pair]
            hits = []
            planes_in = [i for i, (fb, v) in enumerate(frames) if all(par(x & h) == 0 for x in fb)]
            for r in range(0, len(planes_in) + 1):
                for sub in itertools.combinations(planes_in, r):
                    Fb = f2span(list(Gz) + [x for i in sub for x in frames[i][0]])
                    key = elems(Fb)
                    if not perp_to(key, h):
                        continue
                    Gx = f2span([e for e in key if par(e & Q[a]) == 0 and par(e & Q[b]) == 0])
                    if len(Fb) - len(Gx) != 1:
                        continue
                    V = red(cons(key) + Cz)
                    if suits(V, (a, b)):
                        hits.append((Fb, Gx, V))
            return hits
        res = {'S1': 0, 'S211': 0, 'S2': 0}
        pairings = [((bad[0], bad[1]), (bad[2], bad[3])), ((bad[0], bad[2]), (bad[1], bad[3])),
                    ((bad[0], bad[3]), (bad[1], bad[2]))]
        for (Pab, Pcd) in pairings + [(y, x) for (x, y) in pairings]:
            for key, Xb in sums.items():
                Cz = cons(key)
                for Gz in landings(Xb, H[Pab]):
                    gk = elems(Gz)
                    first = two_fix_after(Gz, Cz, Pab)
                    if perp_to(gk, H[Pcd]) and first and two_fix_after(Gz, Cz, Pcd):
                        res['S1'] += 1
                    for c in Pcd:
                        if perp_to(gk, Q[c]) and suits(Cz, (c,)) and first:
                            res['S211'] += 1
                    for (Fb, Gx, V) in first:
                        if perp_to(elems(Gx), H[Pcd]) and two_fix_after(Gx, V, Pcd):
                            res['S2'] += 1
        print('matching', M, 'bad', bad, '| suits three:', triple,
              '| S1 %d S211 %d S2 %d' % (res['S1'], res['S211'], res['S2']), flush=True)


if __name__ == '__main__':
    main()
