"""Exact finite checks for PROOF.md (contents outside the root code).

Every function recomputes its numbers from the published word (common.Geometry) with integers, Fractions
and F2 bit arithmetic.  Model C pricing uses the exact rank formula of the dense table:
    rank_F2 B[R, A] = dim U_A - dim(U_A cap U_R^perp),  U_X = span of the labels q_t, t in X,
cross-checked here against a direct F2 elimination of the matrix B[R, A] = (<q_t, q_s>)."""
import random
from itertools import combinations
from fractions import Fraction as Fr
from common import par, f2span


# ---------------------------------------------------------------- F2 helpers (9-bit vectors)

def span_elems(basis):
    el = [0]
    for b in basis:
        el += [e ^ b for e in el]
    return el


def dim(vs):
    return len(f2span(list(vs)))


def perp(vs, n=9):
    """Basis of the orthogonal complement (standard dot product) of span(vs) in F2^n."""
    vs = list(vs)
    return f2span([w for w in range(1, 1 << n) if all(par(w & v) == 0 for v in vs)])


def meet(A, B):
    """Basis of span(A) cap span(B)."""
    ea = set(span_elems(f2span(A)))
    return f2span([e for e in span_elems(f2span(B)) if e in ea and e])


def rank_formula(Q, R, A):
    UA = f2span([Q[s] for s in A])
    if not UA or not R:
        return 0
    URp = perp([Q[t] for t in R])
    return len(UA) - len(meet(UA, URp))


def rank_direct(Q, R, A):
    rows = []
    A = list(A)
    for t in R:
        m = 0
        for i, s in enumerate(A):
            if par(Q[t] & Q[s]):
                m |= 1 << i
        rows.append(m)
    return len(f2span(rows))


# ---------------------------------------------------------------- isometry with #968's unit

def isometry_968(G, units):
    """Find a linear isometry from the word's E8 (x-basis coordinates, form G.ip) onto #2338's E8 lines
    (doubled coordinates) carrying the 16 removed lines onto F0 + H4 (#968's P104 = E8 - F0 - H4)."""
    E8 = units.E8
    rem968 = [E8[i] for i in units.F0 + units.H4F1]
    cand = rem968 + [tuple(-x for x in v) for v in rem968]
    ipd = lambda a, b: Fr(sum(x * y for x, y in zip(a, b)), 4)       # doubled coordinates -> root norm 2
    X = [tuple(int(i == j) for i in range(8)) for j in range(8)]
    Gx = [[G.ip(a, b) for b in X] for a in X]
    found = []

    def ext(ch):
        k = len(ch)
        if k == 8:
            found.append(list(ch))
            return
        for v in cand:
            if ipd(v, v) == Gx[k][k] and all(ipd(v, ch[i]) == Gx[k][i] for i in range(k)):
                ext(ch + [v])
                if found:
                    return
    ext([])
    assert found, 'no isometric image of the x-basis among the removed roots of #968'
    B = found[0]
    img = lambda c: tuple(sum(c[j] * B[j][i] for j in range(8)) for i in range(8))
    line = lambda v: max(tuple(v), tuple(-x for x in v))
    L968 = {line(v) for v in E8}
    R968 = {line(v) for v in rem968}
    live_img = {line(img(r)) for r in G.roots}
    rem_img = {line(img(v)) for v in G.removed}
    assert live_img | rem_img <= L968 and rem_img == R968 and live_img == L968 - R968
    return {'live_lines_matched': len(live_img), 'removed_lines_matched': len(rem_img),
            'image_of_x_basis_is_removed': True}


# ---------------------------------------------------------------- Collapse: controls on actual designs

def unit(j):
    return tuple(int(i == j) for i in range(8))


def needers(G, j):
    """Ports whose root has a nonzero x_j-coefficient (the frame-pair users of content j)."""
    return [t for t in range(G.P) if G.roots[t][j] != 0]


def support(G, z):
    return [s for s in range(G.P) if G.ip(z, G.roots[s]) != 0]


def collapse_selected(G):
    """The selected centre: 8 pure lifetimes (fresh content x_j copied at h_j^perp, 7 unit steps to the line g_j).
    Collapsed decoder = the frame-pair decoder; Model C cost of each class = 7 = its number of steps."""
    out = []
    for j, L in enumerate(G.lives):
        # the x_j form a basis, so r_t has a unique expansion: the users of content j are exactly its needers
        R, A = needers(G, j), support(G, unit(j))
        cf = set(span_elems(f2span(L['copy_frame'])))
        assert all(G.Q[s] in cf for s in A), 'content support lies in its copy frame'
        hperp = {w for w in range(512) if par(w & G.label(unit(j))) == 0}
        assert set(span_elems(f2span(L['copy_frame']))) == hperp, 'copy frame = h_j^perp'
        rf, rd = rank_formula(G.Q, R, A), rank_direct(G.Q, R, A)
        assert rf == rd
        steps = L['cost']
        # Union lemma bound: dim F_root - dim(F_root cap landing)
        Froot = f2span(L['copy_frame'])
        bound = len(Froot) - len(meet(Froot, f2span(L['landing'])))
        assert rf <= bound <= steps
        assert all(par(G.Q[t] & L['landing'][0]) == 0 for t in R)        # every needer sees the landing
        out.append((rf, bound, steps))
    return {'classes': len(out), 'model_C_costs': [o[0] for o in out], 'union_bounds': [o[1] for o in out],
            'steps': [o[2] for o in out], 'collapsed_cost': sum(o[0] for o in out), 'centre_cost': sum(o[2] for o in out)}


def collapse_all_planes(G, ALL_PLANES):
    """The R_8 = 0 centre of p104-multigen-centre-floor: x_j copied at h_j^perp descends 6 steps to the plane P_j; the
    landed content is carried by a public register up to P_j + P'_j (3-dim) and a second owned lifetime re-descends it
    1 step to the plane P'_j.  With the carrier holding only that content (helper endpoints free, as in the predecessor)
    every step is pure, and the class collapses to one Model C term of cost 7."""
    out = []
    for j, (Pj, P2) in enumerate(ALL_PLANES):
        h = G.label(unit(j))
        Froot = f2span([w for w in range(1, 512) if par(w & h) == 0])
        e1, e2 = f2span(list(Pj)), f2span(list(P2))
        g = G.lives[j]['landing'][0]
        assert len(e1) == len(e2) == 2 and len(f2span(e1 + e2)) == 3
        assert all(par(v & h) == 0 for v in span_elems(e1))                    # P_j inside the copy frame h_j^perp
        assert set(span_elems(e1)) & set(span_elems(e2)) == {0, g}              # P_j cap P'_j = [g_j]
        R = [t for t in needers(G, j) if all(par(G.Q[t] & v) == 0 for v in e1) or all(par(G.Q[t] & v) == 0 for v in e2)]
        assert R == needers(G, j), 'every needer of x_j sees P_j or P2_j'
        A = support(G, unit(j))
        rf = rank_formula(G.Q, R, A)
        assert rf == rank_direct(G.Q, R, A)
        X = meet(meet(Froot, e1), e2)
        bound = len(Froot) - len(X)
        steps = (8 - 2) + (len(f2span(e1 + e2)) - 2)
        assert rf <= bound <= steps
        out.append((rf, bound, steps))
    return {'classes': len(out), 'model_C_costs': [o[0] for o in out], 'union_bounds': [o[1] for o in out],
            'steps': [o[2] for o in out], 'collapsed_cost': sum(o[0] for o in out), 'centre_cost': sum(o[2] for o in out)}


# ---------------------------------------------------------------- randomized tests of Lemma S and the Union lemma

def _rand_sub(rng, d, n=9):
    while True:
        b = f2span([rng.randrange(1, 1 << n) for _ in range(d)])
        if len(b) == d:
            return b


def _rand_extension(rng, base, d):
    """Random d-dim subspace containing span(base)."""
    b = f2span(base)
    while len(b) < d:
        b = f2span(b + [rng.randrange(1, 512)])
    return b


def _rand_hyperplane(rng, F):
    """Random codim-1 subspace of span(F)."""
    F = f2span(F)
    while True:
        w = rng.randrange(1, 512)
        if any(par(w & v) for v in F):
            return f2span([v for v in span_elems(F) if v and par(w & v) == 0])


def union_lemma_random(G, trials, seed):
    """Random pure classes: fresh content (random support inside a random frame) and carried unit steps from
    frames containing earlier landings.  Check rank B[R, A] <= dim F_root - dim(F_root cap all landings) <= #steps,
    with R = every port seeing some landing of the class (the worst case for the users)."""
    rng = random.Random(seed)
    Q, P = G.Q, G.P
    worst_slack, tested = None, 0
    for _ in range(trials):
        d = rng.randrange(2, 10)
        Froot = _rand_sub(rng, d)
        el = set(span_elems(Froot))
        A = [s for s in range(P) if Q[s] in el and rng.random() < 0.7]
        if not A:
            continue
        lands = [_rand_hyperplane(rng, Froot)]
        steps = 1
        for _k in range(rng.randrange(0, 8)):
            parent = rng.choice(lands)
            F = _rand_extension(rng, parent, min(9, len(parent) + rng.randrange(1, 3)))
            lands.append(_rand_hyperplane(rng, F))
            steps += 1
        R = [t for t in range(P) if any(all(par(Q[t] & v) == 0 for v in Lb) for Lb in lands)]
        X = Froot
        for Lb in lands:
            X = meet(X, Lb) if X else []
        bound = len(Froot) - len(X)
        r = rank_formula(Q, R, A)
        assert r == rank_direct(Q, R, A)
        assert r <= bound <= steps, (r, bound, steps)
        tested += 1
        s = steps - r
        worst_slack = s if worst_slack is None else min(worst_slack, s)
    return {'attempted': trials, 'tested': tested, 'min_slack_steps_minus_rank': worst_slack}


def lemma_s_random(G, trials, seed):
    """Random centres WITH mixing.  For a random final step, check
    |ancestry| >= dim(sigma + G) - dim G >= rank B[world(G), supp], with sigma the label span of the union of the
    supports of all fresh parts in the ancestry (generic coefficients: no cancellation, the worst case)."""
    rng = random.Random(seed)
    Q, P = G.Q, G.P
    min_slack, tested = None, 0
    for _ in range(trials):
        steps = []          # (F, G, parents, fresh_support)
        n = rng.randrange(2, 12)
        for i in range(n):
            avail = list(range(i))
            par_ = rng.sample(avail, min(len(avail), rng.randrange(0, 3))) if avail else []
            base = [v for k in par_ for v in steps[k][1]]
            d = min(9, max(len(f2span(base)) + 1, rng.randrange(2, 10)))
            F = _rand_extension(rng, base, d)
            el = set(span_elems(F))
            fresh = [s for s in range(P) if Q[s] in el and rng.random() < 0.5] if (not par_ or rng.random() < 0.5) else []
            if not par_ and not fresh:
                fresh = [s for s in range(P) if Q[s] in el][:1]
            steps.append((F, _rand_hyperplane(rng, F), par_, fresh))
        i = n - 1
        anc, stack = set(), [i]
        while stack:
            k = stack.pop()
            if k in anc:
                continue
            anc.add(k)
            stack += steps[k][2]
        supp = sorted({s for k in anc for s in steps[k][3]})
        if not supp:
            continue
        Gi = steps[i][1]
        sigma = f2span([Q[s] for s in supp])
        lhs = len(anc)
        top = len(f2span([v for k in anc for v in steps[k][0]])) - len(Gi)
        mid = len(f2span(sigma + Gi)) - len(Gi)
        assert lhs >= top >= mid, (lhs, top, mid)
        W = [t for t in range(P) if all(par(Q[t] & v) == 0 for v in Gi)]
        r = rank_formula(Q, W, supp)
        assert r == rank_direct(Q, W, supp)
        assert lhs >= mid >= r, (lhs, mid, r)
        tested += 1
        s = lhs - r
        min_slack = s if min_slack is None else min(min_slack, s)
    return {'attempted': trials, 'tested': tested, 'min_slack_ancestry_minus_rank': min_slack}


def mixing_gap(G):
    """Why Collapse needs mixing-freeness: x_a and x_b land on a common line g (7 + 7 unit steps) and their sum is
    carried g -> 0 (1 step).  The centre pays 15 and gives the ports of world(g) span{x_a, x_b} and every other port
    x_a + x_b.  Model C prices: the three deliveries term by term, and the two natural decoders giving every port at
    least the same span, (x_a, world(g)) + (x_a + x_b, all) and (x_a, all) + (x_b, all).  Every price is cross-checked
    by direct elimination.  Lemma S holds for the sum (ancestry 15 >= its Model C price)."""
    allp = list(range(G.P))

    def price(R, A):
        r = rank_formula(G.Q, R, A)
        assert r == rank_direct(G.Q, R, A)
        return r
    rows = []
    for a, b in combinations(range(8), 2):
        ha, hb = G.label(unit(a)), G.label(unit(b))
        for g in (g for g in range(1, 512) if par(g & ha) == 0 and par(g & hb) == 0):
            W = [t for t in allp if par(G.Q[t] & g) == 0]
            Aa, Ab = support(G, unit(a)), support(G, unit(b))
            As = support(G, tuple(int(i in (a, b)) for i in range(8)))
            ca, cb, cs = price(W, Aa), price(W, Ab), price(allp, As)
            alt1 = price(W, Aa) + price(allp, As)
            alt2 = price(allp, Aa) + price(allp, Ab)
            assert ca == cb == 7
            rows.append((min(alt1, alt2), ca + cb + cs, a, b, g, len(W), alt1, alt2, cs))
    best_alt = min(r[0] for r in rows)
    ex = min(rows)
    assert best_alt == 16 and ex[1] == 23 and 15 >= ex[8]
    return {'pairs_and_lines': len(rows), 'centre_cost': 15, 'term_by_term_model_C': ex[1],
            'best_natural_model_C_decoder': best_alt, 'example': {'contents': [ex[2], ex[3]], 'line': ex[4],
            'world_size': ex[5], 'alternatives': [ex[6], ex[7]]}, 'lemma_S_price_of_sum': ex[8]}


# ---------------------------------------------------------------- evidence 1: truncation

def _qrank(rows):
    R = [list(r) for r in rows]
    rk, n = 0, len(R[0]) if R else 0
    for c in range(n):
        p = next((i for i in range(rk, len(R)) if R[i][c] != 0), None)
        if p is None:
            continue
        R[rk], R[p] = R[p], R[rk]
        for i in range(len(R)):
            if i != rk and R[i][c] != 0:
                f = R[i][c] / R[rk][c]
                R[i] = [x - f * y for x, y in zip(R[i], R[rk])]
        rk += 1
    return rk


def truncation_census(G):
    """Truncate one frame-pair content x_j to the ports with labels in a 7-dim subframe h_j^perp cap k^perp that
    still contains g_j (saving one unit).  Every port re-solves its combination; count ports with no clean vector."""
    P, Q = G.P, G.Q
    lands = [L['landing'][0] for L in G.lives]
    h = [G.label(unit(j)) for j in range(8)]
    u = [[Fr(G.ip(unit(j), G.roots[s])) for s in range(P)] for j in range(8)]
    M = [[s for s in range(P) if s == t or G.ip(G.roots[t], G.roots[s]) == 0] for t in range(P)]
    sees = [[j for j in range(8) if par(Q[t] & lands[j]) == 0] for t in range(P)]

    def clean(t, contents):
        cols = M[t]
        rows = [[c[s] for s in cols] for c in contents]
        tgt = [Fr(int(s == t)) for s in cols]
        return _qrank(rows) == _qrank(rows + [tgt])
    assert all(clean(t, [u[j] for j in sees[t]]) for t in range(P))
    per_content, tested = [], 0
    for j in range(8):
        full = sum(1 for s in range(P) if u[j][s] != 0)
        seen, mins = set(), None
        for k in range(1, 512):
            if k == h[j] or par(k & lands[j]):
                continue
            kept = frozenset(s for s in range(P) if u[j][s] != 0 and par(Q[s] & k) == 0)
            if len(kept) == full or kept in seen:
                continue
            seen.add(kept)
            tested += 1
            uj = [u[j][s] if s in kept else Fr(0) for s in range(P)]
            broken = sum(1 for t in range(P) if j in sees[t] and not clean(t, [uj if i == j else u[i] for i in sees[t]]))
            mins = broken if mins is None else min(mins, broken)
        per_content.append(mins)
    users = [len(needers(G, j)) for j in range(8)]
    return {'truncations_tested': tested, 'min_broken_ports_per_content': per_content, 'users_per_content': users}


# ---------------------------------------------------------------- evidence 2: removed-root directions, overcomplete

def removed_direction_family(G):
    """One-sided decoders whose directions are removed roots (each with the -1 label bonus).  Port t chooses a set
    D_t of removed roots spanning r_t; direction rho is delivered to U_rho = {t : rho in D_t} at cost
    dim span(q_U) - [l_rho in span(q_U)].  Options per port: all minimal spanning sets of size <= 4 and the
    frame-pair support.  Check that the frame-pair (cost 56) has no improving single-port change."""
    P, Q = G.P, G.Q
    rem = [tuple(int(c) for c in v) for v in G.removed]
    lab = G.removed_labels
    supp = [support(G, rho) for rho in rem]
    for i in range(16):        # cost_dir below uses U_{supp rho} = l_rho^perp (8-dim)
        assert set(span_elems(f2span([Q[s] for s in supp[i]]))) == {w for w in range(512) if par(w & lab[i]) == 0}

    def qr(vs):
        return _qrank([[Fr(x) for x in v] for v in vs]) if vs else 0
    opts = []
    for t in range(P):
        rt = tuple(int(c) for c in G.roots[t])
        mins = []
        for k in range(1, 5):
            for D in combinations(range(16), k):
                if any(set(m) <= set(D) for m in mins):
                    continue
                vs = [rem[i] for i in D]
                if qr(vs) == k and qr(vs + [rt]) == k:
                    mins.append(D)
        opts.append(mins)
    xidx = []
    for j in range(8):
        e = unit(j)
        xidx.append(rem.index(e) if e in rem else rem.index(tuple(-x for x in e)))
    choice = []
    for t in range(P):
        want = tuple(sorted(xidx[j] for j in range(8) if G.roots[t][j] != 0))
        if want not in opts[t]:
            opts[t].append(want)
        choice.append(opts[t].index(want))

    def cost_dir(U, i):
        if not U:
            return 0
        b = f2span([Q[t] for t in U])
        x = lab[i]
        for v in sorted(b, reverse=True):
            x = min(x, x ^ v)
        return len(b) - (1 if x == 0 else 0)

    def total(ch):
        U = [[] for _ in range(16)]
        for t, c in enumerate(ch):
            for i in opts[t][c]:
                U[i].append(t)
        return sum(cost_dir(U[i], i) for i in range(16))
    base = total(choice)

    def users(ch):
        U = [[] for _ in range(16)]
        for t, c in enumerate(ch):
            for i in opts[t][c]:
                U[i].append(t)
        return U
    for i, U in enumerate(users(choice)):   # cross-check the direction price against both rank computations
        assert cost_dir(U, i) == rank_formula(Q, U, supp[i]) == rank_direct(Q, U, supp[i])
    best_after_move, moves = None, 0
    for t in range(P):
        for c in range(len(opts[t])):
            if c == choice[t]:
                continue
            ch = list(choice)
            ch[t] = c
            v = total(ch)
            moves += 1
            best_after_move = v if best_after_move is None else min(best_after_move, v)
    sizes = {}
    for o in opts:
        m = min(len(x) for x in o)
        sizes[m] = sizes.get(m, 0) + 1
    return {'options': 'every minimal spanning set of removed roots of size <= 4, plus the frame-pair support',
            'options_total': sum(len(o) for o in opts),
            'frame_pair_cost': base, 'single_port_moves': moves, 'min_cost_after_one_move': best_after_move,
            'ports_by_min_spanning_size': {str(k): sizes[k] for k in sorted(sizes)}}


# ---------------------------------------------------------------- evidence 3: why frame counting cannot prove 56

def live_frames(G):
    P = G.P
    orth = [[p != s and G.ip(G.roots[p], G.roots[s]) == 0 for s in range(P)] for p in range(P)]
    frames = []

    def ext(cl, cand):
        if len(cl) == 8:
            frames.append(cl)
            return
        for i, c in enumerate(cand):
            ext(cl + [c], [d for d in cand[i + 1:] if orth[c][d]])
    ext([], list(range(P)))
    return frames


FRAME_TOUCH_WITNESSES = [(7, [5, 6]), (6, [353, 10, 12]), (5, [497, 258, 132, 24])]


def frame_touch(G, frames):
    """Every live 8-frame must meet the worlds of at least 8 terms (their directions span Q^8).  A counting proof
    of 56 would need every term of cost c to meet at most 759c/7 frames.  These landing frames L (level 9 - dim L,
    contained in a removed hyperplane, so cost level - 1) meet more."""
    out = []
    for level, rows in FRAME_TOUCH_WITNESSES:
        b = f2span(rows)
        assert 9 - len(b) == level
        eps = any(all(par(l & v) == 0 for v in b) for l in G.removed_labels)
        assert eps
        W = {t for t in range(G.P) if all(par(G.Q[t] & v) == 0 for v in b)}
        tau = sum(1 for F in frames if W & set(F))
        cost = level - 1
        assert 7 * tau > 759 * cost
        out.append({'level': level, 'cost': cost, 'world_size': len(W), 'frames_met': tau,
                    'counting_allowance': round(759 * cost / 7, 2)})
    return out


# ---------------------------------------------------------------- evidence 4: mixing in the R_8 = 0 profile

# (j, k, l, a, b): P_j = <l, a>, P_k = <l, b>, mixed rotation plane Pq = <a, b> (the frame-pair lines g_j = a, g_k = b)
MIX_NEAR = [(0, 5, 56, 7, 448), (0, 6, 443, 7, 67), (1, 4, 13, 25, 491), (1, 5, 38, 25, 448), (1, 7, 367, 25, 137),
            (2, 6, 392, 52, 67), (2, 7, 322, 52, 137), (3, 4, 70, 82, 491), (3, 6, 494, 82, 67), (3, 7, 292, 82, 137)]


def _planes_in(vs):
    vs = sorted(set(v for v in vs if v))
    out = set()
    for i, a in enumerate(vs):
        for b in vs[i + 1:]:
            out.add(tuple(sorted((a, b, a ^ b))))
    return sorted(out)


def _qbasis(vs):
    out = []
    for v in vs:
        if _qrank(out + [v]) > len(out):
            out.append(list(v))
    return out


def _qintersect(A, B):
    """basis of span(A) cap span(B) in Q^8 (exact)."""
    if not A or not B:
        return []
    M = [list(a) for a in A] + [[-x for x in b] for b in B]
    n = len(M)
    rows = [[Fr(M[i][c]) for i in range(n)] for c in range(8)]
    piv, r = [], 0
    for c in range(n):
        p = next((i for i in range(r, 8) if rows[i][c] != 0), None)
        if p is None:
            continue
        rows[r], rows[p] = rows[p], rows[r]
        pv = rows[r][c]
        rows[r] = [x / pv for x in rows[r]]
        for i in range(8):
            if i != r and rows[i][c] != 0:
                f = rows[i][c]
                rows[i] = [a - f * b for a, b in zip(rows[i], rows[r])]
        piv.append(c)
        r += 1
    vecs = []
    for fcol in (c for c in range(n) if c not in piv):
        v = [Fr(0)] * n
        v[fcol] = Fr(1)
        for i, c in enumerate(piv):
            v[c] = -rows[i][fcol]
        vecs.append([sum(v[i] * A[i][d] for i in range(len(A))) for d in range(8)])
    return _qbasis(vecs)


def mixing_census(G):
    """R_8 = 0 centres: x_j lands on a plane P_j inside h_j^perp (6 steps), rotations cost 1 (a plane inside a 3-dim frame
    holding the landed planes of the contents it combines).  Checks (exact):
      * a plane inside h_j^perp sees at most 30 of the 52 needers of x_j (any plane: 32);
      * the ten near-solutions: one mixed rotation e_j +- e_k at <g_j, g_k> serves all needers of x_j and x_k except one
        sign-mismatched double misser;
      * the four perfect matchings of these pairs give centres of cost 48 + 4 = 52 serving every port but four;
      * standalone repairs (source frame = a sum of landed planes, cheapest target, any admissible content) fix at most
        two of the four ports, and two only at cost 2: so every matched design needs >= 4 more units (56)."""
    P, Q = G.P, G.Q
    pm = [sum(1 << t for t in range(P) if par(Q[t] & w) == 0) for w in range(512)]
    W = lambda pl: pm[pl[0]] & pm[pl[1]]
    roots = [tuple(int(c) for c in r) for r in G.roots]
    need = [sum(1 << t for t in range(P) if roots[t][j]) for j in range(8)]
    h = [G.label(unit(j)) for j in range(8)]
    g = [L['landing'][0] for L in G.lives]
    allpl = _planes_in(range(512))
    assert len(allpl) == 43435
    seen_in, seen_any = [], []
    for j in range(8):
        inside = _planes_in([w for w in range(512) if par(w & h[j]) == 0])
        assert len(inside) == 10795
        seen_in.append(max((W(p) & need[j]).bit_count() for p in inside))
        seen_any.append(max((W(p) & need[j]).bit_count() for p in allpl))
    assert seen_in == [30] * 8 and seen_any == [32] * 8
    pl = lambda x, y: tuple(sorted((x, y, x ^ y)))
    near = {}
    for (j, k, l, a, b) in MIX_NEAR:
        assert a == g[j] and b == g[k]
        Pj, Pk, Pq = pl(l, a), pl(l, b), pl(a, b)
        assert all(par(v & h[j]) == 0 for v in Pj) and all(par(v & h[k]) == 0 for v in Pk)
        Mj, Mk = need[j] & ~W(Pj), need[k] & ~W(Pk)
        Wq = W(Pq)
        assert Mj == Mk and Mj.bit_count() == 28 and Mj & ~Wq == 0
        D = [t for t in range(P) if Mj >> t & 1]
        ratios = {}
        for t in D:
            ratios.setdefault(Fr(roots[t][j], roots[t][k]), []).append(t)
        assert sorted(len(v) for v in ratios.values()) == [1, 27]
        main = next(r for r, v in ratios.items() if len(v) == 27)
        badp = next(v[0] for v in ratios.values() if len(v) == 1)
        near[(j, k)] = (Pj, Pk, Pq, main, badp)
    matchings = [c for c in combinations(sorted(near), 4) if len({x for e in c for x in e}) == 8]
    assert len(matchings) == 4
    out = []
    for M in matchings:
        landed = []
        for (j, k) in M:
            Pj, Pk, Pq, main, _ = near[(j, k)]
            landed += [(Pj, [Fr(int(i == j)) for i in range(8)]), (Pk, [Fr(int(i == k)) for i in range(8)])]
            v = [Fr(0)] * 8
            v[j], v[k] = Fr(main.numerator), Fr(main.denominator)
            landed.append((Pq, v))
        S = [_qbasis([v for (p, v) in landed if W(p) >> t & 1]) for t in range(P)]
        bad = [t for t in range(P) if _qrank(S[t] + [list(roots[t])]) > len(S[t])]
        assert sorted(bad) == sorted(near[e][4] for e in M)
        best, frames = {}, set()
        for r in range(1, len(landed) + 1):
            for sub in combinations(range(len(landed)), r):
                Xb = f2span([x for i in sub for x in landed[i][0]])
                key = frozenset(span_elems(Xb))
                if key in frames:
                    continue
                frames.add(key)
                C = _qbasis([v for (p, v) in landed if set(p) <= key])
                for nb in range(1, 5):
                    for B in combinations(bad, nb):
                        Yb = f2span([e for e in key if all(par(e & Q[t]) == 0 for t in B)])
                        cost = len(Xb) - len(Yb)
                        if nb in best and best[nb] <= cost:
                            continue
                        A = C
                        for t in B:
                            A = _qintersect(A, S[t] + [list(roots[t])])
                            if not A:
                                break
                        if A and all(_qrank(S[t] + A) > len(S[t]) for t in B):
                            best[nb] = cost
        assert 3 not in best and 4 not in best and best.get(1) == 1 and best.get(2, 2) == 2
        out.append({'pairs': [list(e) for e in M], 'cost_before_repairs': 52, 'unserved_ports': sorted(bad),
                    'source_frames': len(frames), 'cheapest_repair_by_ports_fixed': {str(k): v for k, v in sorted(best.items())}})
    return {'needers_seen_by_one_plane': {'inside_copy_frame': 30, 'any_plane': 32, 'needers': 52},
            'near_solutions': [[j, k] for (j, k, *_r) in MIX_NEAR], 'double_missers_per_pair': 28,
            'sign_mismatched_ports_per_pair': 1, 'matched_designs': out}
