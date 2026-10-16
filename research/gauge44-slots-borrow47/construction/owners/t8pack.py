"""T8 whole-owner fill for a single-provider plan: window LP -> MILP on LP support -> exact time-incidence MILP -> Konig colouring -> roles.
usage: t8pack.py PLAN.json OUT-OWNERS.json.gz"""
import collections, gzip, json, math, sys, time
import numpy as np, scipy.sparse as sp
from scipy.optimize import linprog, milp, LinearConstraint, Bounds
from scipy.sparse import coo_matrix
from scipy.sparse.csgraph import maximum_flow
plan = json.load(open(sys.argv[1])); NT = len(plan['order'])
order = [tuple(x) for x in plan['order']]
pcw = {(z['stage'], z['replica']): z['widths'] for z in plan['per_call_widths']}
widths_list = pcw[order[0]]; assert all(pcw[c] == widths_list for c in order)
INV = collections.Counter(widths_list); W = sorted(INV)
L = lambda r: math.log(73 / r)
def patterns(cap):
    out = []
    def rec(i, rem, cur):
        if i == len(W): out.append(tuple(cur)); return
        for k in range(rem // W[i] + 1): rec(i + 1, rem - k * W[i], cur + [k])
    rec(0, cap, []); return out
PATS = {c: patterns(c) for c in (8, 16, 32, 48)}
rem = lambda cap, p: cap - sum(a * b for a, b in zip(p, W))
members = collections.defaultdict(list)
for key, (cap, times) in plan['windows'].items():
    rep, ch = map(int, key.split(',')); members[(cap, tuple(sorted(times)))].append((rep, ch))
cl = sorted(members); NB = {k: 112 * len(members[k]) for k in cl}

def model(keep):
    xi = {}; yi = {}; cost = []
    for ci, (cap, S) in enumerate(cl):
        for p in keep[ci]:
            r = rem(cap, p); xi[ci, p] = len(cost); cost.append(0.0 if r == 0 else r * (L(r) - L(48)))
    nx = len(cost)
    for ci, (cap, S) in enumerate(cl):
        for w in W:
            for t in S: yi[ci, w, t] = nx + len(yi)
    I = []; J = []; V = []; lo = []; hi = []
    def row(es, a, b):
        r = len(lo)
        for j, v in es: I.append(r); J.append(j); V.append(v)
        lo.append(a); hi.append(b)
    bytime = collections.defaultdict(list)
    for ci, (cap, S) in enumerate(cl):
        row([(xi[ci, p], 1.0) for p in keep[ci]], NB[cl[ci]], NB[cl[ci]])
        for wi, w in enumerate(W): row([(xi[ci, p], float(p[wi])) for p in keep[ci] if p[wi]] + [(yi[ci, w, t], -1.0) for t in S], 0, 0)
        for t in S: bytime[t].append(ci); row([(yi[ci, w, t], 1.0) for w in W], 0, NB[cl[ci]])
    for t in range(NT):
        for w in W: row([(yi[ci, w, t], 1.0) for ci in bytime[t]], INV[w], INV[w])
    A = sp.csr_matrix((V, (I, J)), shape=(len(lo), nx + len(yi)))
    return xi, yi, np.array(cost + [0.0] * len(yi)), A, np.array(lo, float), np.array(hi, float), nx

t0 = time.time()
full = {ci: PATS[k[0]] for ci, k in enumerate(cl)}
xi, yi, c, A, lo, hi, nx = model(full)
eq = lo == hi
lpres = linprog(c, A_ub=sp.vstack([A[~eq], -A[~eq]]).tocsr(), b_ub=np.concatenate([hi[~eq], -lo[~eq]]), A_eq=A[eq], b_eq=lo[eq], bounds=(0, None), method='highs')
print('LP', lpres.status, 'obj/call %.4f' % (lpres.fun / NT), 'secs %.0f' % (time.time() - t0), flush=True)
lpidle = collections.Counter()
for (ci, p), j in xi.items():
    r = rem(cl[ci][0], p)
    if r and lpres.x[j] > 1e-9: lpidle[r] += lpres.x[j]
print('LP idle', {k: round(v, 2) for k, v in sorted(lpidle.items())}, flush=True)
sup = {ci: [p for p in full[ci] if lpres.x[xi[ci, p]] > 1e-7] for ci in range(len(cl))}
del A
res = None
for slack in tuple(int(x) for x in __import__("os").environ.get("SLACKS", "0,1,2,4").split(",")):
    keep = {ci: sup[ci] + [p for p in full[ci] if p not in sup[ci] and rem(cl[ci][0], p) <= slack and rem(cl[ci][0], p) > 0] for ci in range(len(cl))}
    xi, yi, c, A, lo, hi, nx = model(keep)
    integ = np.zeros(len(c)); integ[:nx] = 1
    res = milp(c, constraints=LinearConstraint(A, lo, hi), integrality=integ, bounds=Bounds(0, np.inf), options=dict(time_limit=900, mip_rel_gap=1e-6))
    print('MILP slack', slack, 'nx', nx, 'status', res.status, 'obj/call', None if res.x is None else '%.4f' % (res.fun / NT), 'secs %.0f' % (time.time() - t0), flush=True)
    if res.x is not None: break
bins = {}
for ci, k in enumerate(cl):
    lst = []
    for p in keep[ci]: lst += [p] * int(round(res.x[xi[ci, p]]))
    assert len(lst) == NB[k], (k, len(lst)); bins[ci] = lst
idle = collections.Counter(rem(k[0], p) for ci, k in enumerate(cl) for p in bins[ci] if rem(k[0], p))
print('idle', dict(sorted(idle.items())), flush=True)

# exactly liftable per-pattern time incidences x(c,p,t,w): sum_t x = n*mult_p(w), sum_w x(c,p,t,.) <= n, per-call inventories
from scipy.sparse.csgraph import maximum_bipartite_matching
NP = {}
for ci in range(len(cl)):
    for p in bins[ci]: NP[ci, p] = NP.get((ci, p), 0) + 1
xv = {}; I = []; J = []; V = []; lo = []; hi = []
def row(es, a, b):
    r = len(lo)
    for j_, v_ in es: I.append(r); J.append(j_); V.append(v_)
    lo.append(a); hi.append(b)
bytw = collections.defaultdict(list)
for (ci, p), n in NP.items():
    cap, S = cl[ci]
    for wi, w in enumerate(W):
        if p[wi]:
            for t in S: xv[ci, p, t, w] = len(xv); bytw[t, w].append(xv[ci, p, t, w])
    for wi, w in enumerate(W):
        if p[wi]: row([(xv[ci, p, t, w], 1.0) for t in S], n * p[wi], n * p[wi])
    if any(p):
        for t in S: row([(xv[ci, p, t, w], 1.0) for wi, w in enumerate(W) if p[wi]], 0, n)
for t in range(NT):
    for w in W: row([(j_, 1.0) for j_ in bytw[t, w]], INV[w], INV[w])
A = sp.csr_matrix((V, (I, J)), shape=(len(lo), len(xv)))
r2 = milp(np.zeros(len(xv)), constraints=LinearConstraint(A, lo, hi), integrality=np.ones(len(xv)), bounds=Bounds(0, np.inf), options=dict(time_limit=1800))
print('incidence MILP', len(xv), 'vars status', r2.status, 'secs %.0f' % (time.time() - t0), flush=True)
assert r2.x is not None
X = {k: int(round(r2.x[j_])) for k, j_ in xv.items()}
assigned = []   # (t, w, ci, bin index within class)
binbase = {}; q = collections.Counter()
for ci in range(len(cl)):
    for p in sorted(set(bins[ci])): binbase[ci, p] = q[ci]; q[ci] += NP[ci, p]
for (ci, p), n in NP.items():
    cap, S = cl[ci]
    if not any(p): continue
    pos = [w for wi, w in enumerate(W) for _ in range(p[wi])]; m = len(pos); St = list(S)
    assert len(St) >= m
    # multigraph times x positions, positions of width w share x(c,p,t,w) evenly-greedily so each has degree n
    M = np.zeros((len(St), len(St)), dtype=np.int64)   # right side: m positions + (|S|-m) dummies
    for w in set(pos):
        js = [j_ for j_ in range(m) if pos[j_] == w]; need = {j_: n for j_ in js}
        k = 0
        for ti, t in enumerate(St):
            a_ = X[ci, p, t, w]
            while a_:
                while need[js[k]] == 0: k += 1
                d = min(a_, need[js[k]]); M[ti, js[k]] += d; need[js[k]] -= d; a_ -= d
    deg = M[:, :m].sum(1); assert deg.max() <= n
    k = m
    for ti in range(len(St)):
        d = n - deg[ti]
        while d:
            room = n - M[:, k].sum()
            if room == 0: k += 1; continue
            e = min(d, room); M[ti, k] += e; d -= e
    assert (M.sum(0) == n).all() and (M.sum(1) == n).all()
    for b in range(n):
        mt = maximum_bipartite_matching(sp.csr_matrix((M > 0).astype(np.int8)), perm_type='column')
        assert (mt >= 0).all(), 'no perfect matching'
        for ti, j_ in enumerate(mt):
            M[ti, j_] -= 1
            if j_ < m: assigned.append((St[ti], pos[j_], ci, binbase[ci, p] + b, j_))
print('coloured', len(assigned), flush=True)
pairs, sides = plan['pairs'], plan['sides']
hw = collections.defaultdict(list)
for h, w in enumerate(widths_list): hw[w].append(h)
binid = {}
for ci, k in enumerate(cl):
    q = 0
    for (rep, ch) in sorted(members[k]):
        for port in range(112): binid[ci, q] = (rep, ch, port); q += 1
cursor = collections.Counter(); offset = collections.Counter(); roles = []
for t, w, ci, p, j_ in sorted(assigned, key=lambda z: (z[2], z[3], z[4])):
    s, rep = order[t]; h = hw[w][cursor[t, w]]; cursor[t, w] += 1; orep, ch, port = binid[ci, p]
    roles.append(dict(child=ch, copy=112 * orep + port, family=orep * 896 + (2 * pairs[str(ch)] + sides[str(ch)]) * 112 + port, helper=h, offset=offset[orep, ch, port],
                      owner_replica=orep, port=port, replica=rep, stage=s, time=t, width=w))
    offset[orep, ch, port] += w
assert all(cursor[t, w] == INV[w] for t in range(NT) for w in W)
roles.sort(key=lambda z: (z['time'], z['helper']))
fam = collections.Counter((z['time'], z['family']) for z in roles); assert max(fam.values()) == 1
idle2 = collections.Counter()
for key, (cap, times) in plan['windows'].items():
    rep, ch = map(int, key.split(','))
    for port in range(112):
        r = cap - offset[rep, ch, port]
        if r: idle2[r] += 1
print('roles', len(roles), 'idle', dict(sorted(idle2.items())), flush=True)
open(sys.argv[2], 'wb').write(gzip.compress(json.dumps(roles, separators=(',', ':')).encode(), mtime=0))
