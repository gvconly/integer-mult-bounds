"""realize.py COUNTS.json OUTDIR [width order]: exact FILL.json.gz / ORDER.json from integer class pattern counts.
Phase 1: per-width max flows from calls to window classes (scipy maximum_flow), with per-(class, call) capacities;
phase 2: per class, units to bins by per-width bipartite flows (each bin at most one helper per call), exact ILP fallback,
and capacity learning on failures.  Helpers are then named per call and prefixes packed.  reconstruct.py checks the
result literally.  Construction-side.  Prepared by G. V. Chudnovsky with Anthropic Claude assistance. Apache-2.0."""
import collections, gzip, json, sys, random
from pathlib import Path
import numpy as np
from scipy.sparse import coo_matrix
from scipy.sparse.csgraph import maximum_flow
import orderlp as O
d = json.load(open(sys.argv[1])); out = Path(sys.argv[2]); out.mkdir(parents=True, exist_ok=True)
order = d['order']; NT = len(order); W = O.widths; INV = O.INV
windows, _ = O.mod.actual_windows(O.parent, order)
# classes -> member (replica, child) groups
members = collections.defaultdict(list)
for k, (cap, ts) in windows.items():
    r, c = map(int, k.split(',')); members[(cap, tuple(sorted(ts)))].append((r, c))
cls = [(c[0], tuple(c[1])) for c in d['classes']]
bins = {}   # class index -> list of pattern tuples (len = 112 * #members)
pc = collections.defaultdict(list)
for cap, S, p, x in d['solution']: pc[(cap, tuple(S))] += [tuple(p)] * x
for ci, key in enumerate(cls):
    b = pc[key]; assert len(b) == 112 * len(members[key]), (key, len(b), len(members[key])); bins[ci] = b
def flow(N, edges, src, sink):
    r, c, v = zip(*edges); m = coo_matrix((np.array(v, dtype=np.int64), (r, c)), shape=(N, N)).tocsr(); return maximum_flow(m, src, sink)
used = {ci: [set() for _ in bins[ci]] for ci in bins}; assigned = []
# phase 1: per-width global flows, with residual per-(class, call) capacity = nonempty bins - units already routed at that call
nc = len(cls); WORDER = [int(x) for x in (sys.argv[3] if len(sys.argv) > 3 else "4,5,7,8,9").split(",")]
nonempty = {ci: sum(1 for b in bins[ci] if any(b)) for ci in range(nc)}
withw = {(ci, w): sum(1 for b in bins[ci] if b[W.index(w)]) for ci in range(nc) for w in W}
captw = {}; captot = {}
for round_ in range(40):
    assigned = []
    routed = collections.Counter(); units = collections.defaultdict(dict)   # units[ci][(t,w)] = n
    for w in WORDER:
        j = W.index(w); sink = 1 + NT + nc; E = []
        for t in range(NT): E.append((0, 1 + t, INV[w]))
        for ci in range(nc):
            need = sum(b[j] for b in bins[ci])
            if not need: continue
            E.append((1 + NT + ci, sink, need))
            for t in cls[ci][1]:
                cap = min(withw[ci, w], min(nonempty[ci], captot.get((ci, t), 10**9)) - routed[ci, t], captw.get((ci, t, w), 10**9))
                if cap > 0: E.append((1 + t, 1 + NT + ci, cap))
        F = flow(sink + 1, E, 0, sink); print('FLOW', w, F.flow_value, NT * INV[w], flush=True)
        assert F.flow_value == NT * INV[w]
        co = F.flow.tocoo()
        for a, b, v in zip(co.row, co.col, co.data):
            if v > 0 and 1 <= a <= NT and 1 + NT <= b < sink:
                t = int(a - 1); ci = int(b - 1 - NT); units[ci][(t, w)] = int(v); routed[ci, t] += int(v)
    # phase 2: per class, (t,w) -> (p,w) [cap 1] -> sink [slots]; then repair cross-width same-call conflicts by swaps
    fails = []
    for ci in range(nc):
        U = units[ci]
        if not U: continue
        nb = len(bins[ci]); need = sum(U.values()); ok = False
        import itertools
        for worder in ([9, 7, 8, 5, 4], [4, 5, 9, 7, 8], [7, 9, 8, 4, 5], [8, 9, 7, 4, 5], [5, 4, 7, 9, 8]):
            slot = collections.defaultdict(list); usedt = [set() for _ in range(nb)]; good = True
            for w in worder:
                wi = W.index(w); ks = [k for k in U if k[1] == w]
                if not ks: continue
                tt = [k[0] for k in ks]; T0 = 1 + nb; sk = T0 + len(tt); E = []
                for p, bb in enumerate(bins[ci]):
                    if bb[wi]:
                        E.append((0, 1 + p, bb[wi]))
                        for k_, t in enumerate(tt):
                            if t not in usedt[p]: E.append((1 + p, T0 + k_, 1))
                for k_, t in enumerate(tt): E.append((T0 + k_, sk, U[(t, w)]))
                ff = flow(sk + 1, E, 0, sk)
                if ff.flow_value != sum(U[k] for k in ks): good = False; break
                G = ff.flow.tocoo()
                for a_, b_, v_ in zip(G.row, G.col, G.data):
                    if v_ > 0 and 1 <= a_ < T0 and T0 <= b_ < sk:
                        p = a_ - 1; t = tt[b_ - T0]; usedt[p].add(t); slot[p].append((t, w))
            if good: ok = True; break
        if not ok:
            # exact ILP: binary z[t,w,p] for slots; each (t,w) count met, each bin's slots met, each (bin, t) at most 1
            import numpy as _np, scipy.sparse as _sp
            from scipy.optimize import milp as _milp, LinearConstraint as _LC, Bounds as _B
            var = [(t, w, p) for (t, w) in U for p, bb in enumerate(bins[ci]) if bb[W.index(w)]]
            vid = {v: i for i, v in enumerate(var)}; I_ = []; J_ = []; lo_ = []; hi_ = []; r_ = 0
            grp = collections.defaultdict(list)
            for i, (t, w, p) in enumerate(var): grp[('u', t, w)].append(i); grp[('b', p, w)].append(i); grp[('pt', p, t)].append(i)
            for key, ids in grp.items():
                for i in ids: I_.append(r_); J_.append(i)
                if key[0] == 'u': lo_.append(U[(key[1], key[2])]); hi_.append(U[(key[1], key[2])])
                elif key[0] == 'b': v_ = bins[ci][key[1]][W.index(key[2])]; lo_.append(v_); hi_.append(v_)
                else: lo_.append(0); hi_.append(1)
                r_ += 1
            Am = _sp.csr_matrix((_np.ones(len(I_)), (I_, J_)), shape=(r_, len(var)))
            rr = _milp(_np.zeros(len(var)), constraints=_LC(Am, lo_, hi_), integrality=_np.ones(len(var)), bounds=_B(0, 1), options=dict(time_limit=120))
            if rr.x is not None:
                slot = collections.defaultdict(list)
                for i, (t, w, p) in enumerate(var):
                    if rr.x[i] > 0.5: slot[p].append((t, w))
                ok = True; print('    ILP rescued class', ci, flush=True)
            else:
                tot = collections.Counter()
                for (t, w), n in U.items(): tot[t] += n
                for t, n in tot.most_common(3): captot[(ci, t)] = n - 1
                fails.append((ci, -1, need)); continue
        for p, lst in slot.items():
            for t, w in lst: assigned.append((t, w, ci, p))
    print('round', round_, 'phase-2 failures', len(fails), fails[:4], flush=True)
    if not fails: break
if fails: raise SystemExit(3)
# name helpers and bins
prov = O.rd(O.F / 'new/BORROW-PROVIDERS.json'); widths = prov['gauge44']['widths']; assert all(v['widths'] == widths for v in prov.values())
hw = collections.defaultdict(list)
for h, w in enumerate(widths): hw[w].append(h)
binid = {}
for ci, key in enumerate(cls):
    k = 0
    for (r, c) in sorted(members[key]):
        for port in range(112): binid[ci, k] = (r, c, port); k += 1
cursor = collections.Counter(); offset = collections.Counter(); rows = []
for t, w, ci, p in sorted(assigned):
    h = hw[w][cursor[t, w]]; cursor[t, w] += 1; r, c, port = binid[ci, p]; rows.append([t, h, r, c, port, offset[r, c, port]]); offset[r, c, port] += w
assert all(cursor[t, w] == INV[w] for t in range(NT) for w in W)
shipped = O.rd(O.F / 'new/borrow-fill/FILL.json.gz')
fill = dict(columns=shipped['columns'], dispatch_overrides=shipped['dispatch_overrides'], rows=rows)
(out / 'FILL.json.gz').write_bytes(gzip.compress(json.dumps(fill, separators=(',', ':')).encode(), mtime=0))
(out / 'ORDER.json').write_text(json.dumps(order))
idle = collections.Counter()
for k, (cap, ts) in windows.items():
    r, c = map(int, k.split(','))
    for port in range(112):
        rem = cap - offset[r, c, port]
        if rem: idle[rem] += 1
print('DONE rows', len(rows), 'idle', dict(sorted(idle.items())), flush=True)
