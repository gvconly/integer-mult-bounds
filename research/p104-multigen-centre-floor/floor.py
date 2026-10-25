"""Finite certificates used by PROOF.md (all exact).

capacities  m_d = max #port labels in a d-dim subspace of F2^9, d = 1..8, by exhaustive enumeration.
sources     RC first-generation content exists only at F2^9 and the 16 removed hyperplanes (Lemma 2).
good        live roots orthogonal to r_p span r_p^perp for every port (Lemma 1).
tsets       |T_z| >= 49 for all lattice z of norm <= 8 and all E6 continua (with the classification in PROOF.md
            this covers every z); the candidate sets with |T_z| <= 72.
two_plane   T_z is covered by two plane-worlds iff T_z = T_rho for a removed root rho (Lemma 6).
no_a2       no three removed roots span a 2-plane."""
from itertools import combinations
from common import par, f2span


def subspaces(k, n=9):
    """RREF bases of all k-dim subspaces of F2^n."""
    for piv in combinations(range(n), k):
        free = [[c for c in range(n) if c > piv[i] and c not in piv] for i in range(k)]
        nfree = sum(len(f) for f in free)
        for bits in range(1 << nfree):
            rows, b = [], bits
            for i in range(k):
                r = 1 << piv[i]
                for c in free[i]:
                    if b & 1:
                        r |= 1 << c
                    b >>= 1
                rows.append(r)
            yield rows


def perp_masks(G):
    """mask[w] = bitmask of ports p with <q_p, w> = 0."""
    return [sum(1 << p for p in range(G.P) if par(G.Q[p] & w) == 0) for w in range(512)]


def capacities(G):
    mask = perp_masks(G)
    full = (1 << G.P) - 1
    m, counts = {}, {}
    for k in range(1, 9):
        best, cnt = 0, 0
        for rows in subspaces(k):
            acc = full
            for r in rows:
                acc &= mask[r]
            c = acc.bit_count()
            cnt += 1
            if c > best:
                best = c
        m[9 - k] = best
        counts[9 - k] = cnt
    return m, counts


def qrank(G, vecs):
    rows = [list(v) for v in vecs]
    r = 0
    from fractions import Fraction as Fr
    rows = [[Fr(x) for x in v] for v in rows]
    for c in range(8):
        p = next((i for i in range(r, len(rows)) if rows[i][c] != 0), None)
        if p is None:
            continue
        rows[r], rows[p] = rows[p], rows[r]
        for i in range(len(rows)):
            if i != r and rows[i][c] != 0:
                f = rows[i][c] / rows[r][c]
                rows[i] = [a - f * b for a, b in zip(rows[i], rows[r])]
        r += 1
    return r


def good(G):
    return sum(qrank(G, [G.roots[s] for s in range(G.P) if s != p and G.ip(G.roots[s], G.roots[p]) == 0]) == 7
               for p in range(G.P))


def sources(G):
    """For each hyperplane w^perp: rank of the roots of ports outside it.  Removed labels: 7, others: 8."""
    removed = set(G.removed_labels)
    bad = []
    for w in range(1, 512):
        rk = qrank(G, [G.roots[p] for p in range(G.P) if par(G.Q[p] & w) == 1])
        if (w in removed and rk != 7) or (w not in removed and rk != 8):
            bad.append((w, rk))
    return bad


def lattice8(G):
    """All E8 vectors of norm <= 8, one per sign pair (13,320 = (240+2160+6720+17520)/2)."""
    R2 = []
    for v in G.lines:
        R2 += [v, tuple(-x for x in v)]
    level = {G.line(v): v for v in R2}
    frontier = dict(level)
    for _ in range(3):
        new = {}
        for v in frontier.values():
            for r in R2:
                w = tuple(a + b for a, b in zip(v, r))
                n = G.ip(w, w)
                if 0 < n <= 8:
                    k = G.line(w)
                    if k not in level and k not in new:
                        new[k] = w
        level.update(new)
        frontier = new
    return level


def tmask(G, z):
    return sum(1 << p for p in range(G.P) if G.ip(z, G.roots[p]) != 0)


def tsets(G):
    vecs = lattice8(G)
    by_norm = {}
    cands = {}
    tmin = 10 ** 9
    for k, z in vecs.items():
        t = tmask(G, z)
        n = G.ip(z, z)
        c = t.bit_count()
        by_norm.setdefault(n, set()).add(c)
        tmin = min(tmin, c)
        if c <= 72:
            cands.setdefault(t, set()).add('norm%d' % n)
    # E6 continua: generic z in span(a, b) for an A2 {a, b}: T = ports with <a,r> != 0 or <b,r> != 0
    e6min = 10 ** 9
    nA2 = 0
    for i, a in enumerate(G.lines):
        for b in G.lines[i + 1:]:
            if G.ip(a, b) in (1, -1):
                nA2 += 1
                t = tmask(G, a) | tmask(G, b)
                e6min = min(e6min, t.bit_count())
                if t.bit_count() <= 72:
                    cands.setdefault(t, set()).add('E6-continuum')
    return {'vectors_norm_le_8': len(vecs), 'T_sizes_by_norm': {str(n): sorted(s) for n, s in sorted(by_norm.items())},
            'min_T_lattice': tmin, 'A2_pairs': nA2, 'min_T_E6_continuum': e6min}, cands


def two_plane(G, cands):
    mask = perp_masks(G)
    worlds = []
    for a in range(1, 512):
        for b in range(a + 1, 512):
            if (a ^ b) > b:
                worlds.append(mask[a] & mask[b])
    assert len(worlds) == 43435
    inv = [~w for w in worlds]
    covered = []
    for t in cands:
        for iw in inv:
            rest = t & iw
            if rest.bit_count() > 36:
                continue
            labs = [G.Q[p] for p in range(G.P) if rest >> p & 1]
            if len(f2span(labs)) <= 7:
                covered.append(t)
                break
    removed_T = {tmask(G, rho) for rho in G.removed}
    return len(worlds), covered, removed_T


def no_a2(G):
    rem = set(G.removed) | {tuple(-x for x in v) for v in G.removed}
    hits = []
    for i, a in enumerate(G.removed):
        for b in G.removed[i + 1:]:
            if G.ip(a, b) in (1, -1):
                for s in (1, -1):
                    if tuple(x + s * y for x, y in zip(a, b)) in rem:
                        hits.append((a, b))
    return hits
