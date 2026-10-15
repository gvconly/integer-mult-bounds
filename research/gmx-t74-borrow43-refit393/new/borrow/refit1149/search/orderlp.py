"""orderlp.py ORDER-SPEC [SOL.json]
Window-aware aggregate packing LP for a call order (the LP of window_lp.py, parametrised by the order): ORDER-SPEC is
'shipped', 'D=<d>' (key replica + 4*stage + d*[stage>=1]; d = 4 is PR2571's shifted diagonal), 'K=a,d,e' (key replica +
a*stage + d*[stage>=1] + e*[stage=8]) or a JSON order.  Prints the idle excess per call; with SOL.json also writes the
LP pattern counts per window class.  Results (per call): shipped 275.010, d=3 and d=4 272.801, d=5 272.886, d=6 273.719,
d=8 275.386, a=5 276.433, a=3 273.249; time-free 271.593.  Construction-side (numpy, scipy).
Prepared by G. V. Chudnovsky with Anthropic Claude assistance. Apache-2.0."""
import collections, gzip, importlib.util, json, math, sys, time
import numpy as np, scipy.sparse as sp
from scipy.optimize import linprog
from pathlib import Path
F = Path(__file__).resolve().parents[4]   # package root
rd = lambda p: json.loads(gzip.decompress(p.read_bytes()) if p.suffix == '.gz' else p.read_bytes())
parent = rd(F / 'inputs/GLOBAL-TEMPLATE.json.gz')
spw = importlib.util.spec_from_file_location('windows', F / 'vendor/build_parent.py'); mod = importlib.util.module_from_spec(spw); spw.loader.exec_module(mod)
INV = {4: 26, 5: 9, 7: 62, 8: 78, 9: 240}; widths = sorted(INV)
L = lambda r: math.log(73 / r)
def patterns(cap):
    out = []
    def rec(i, rem, cur):
        if i == len(widths): out.append(tuple(cur)); return
        for k in range(rem // widths[i] + 1): rec(i + 1, rem - k * widths[i], cur + [k])
    rec(0, cap, []); return out
PATS = {c: patterns(c) for c in (8, 16, 32, 48)}
def make_order(key):
    calls = [(s, r) for s in range(9) for r in range(74)]
    return [list(c) for c in sorted(calls, key=key)]
def lp(order, want_x=False):
    windows, _ = mod.actual_windows(parent, order); NT = len(order)
    classes = collections.Counter()
    for k, (cap, times) in windows.items(): classes[(cap, tuple(sorted(times)))] += 112
    cl = list(classes.items()); xi = {}; yi = {}; cost = []
    for ci, ((cap, S), n) in enumerate(cl):
        for p in PATS[cap]:
            r = cap - sum(k * w for k, w in zip(p, widths)); xi[ci, p] = len(cost); cost.append(0.0 if r == 0 else r * (L(r) - L(48)))
    nx = len(cost)
    for ci, ((cap, S), n) in enumerate(cl):
        for w in widths:
            for t in S: yi[ci, w, t] = nx + len(yi)
    cost = np.array(cost + [0.0] * len(yi)); I = []; J = []; V = []; rhs = []
    def row(es, b):
        r = len(rhs)
        for j, v in es: I.append(r); J.append(j); V.append(v)
        rhs.append(b)
    bytime = collections.defaultdict(list)
    for ci, ((cap, S), n) in enumerate(cl):
        row([(xi[ci, p], 1.0) for p in PATS[cap]], n)
        for wi, w in enumerate(widths): row([(xi[ci, p], float(p[wi])) for p in PATS[cap] if p[wi]] + [(yi[ci, w, t], -1.0) for t in S], 0.0)
        for t in S: bytime[t].append(ci)
    for t in range(NT):
        for w in widths: row([(yi[ci, w, t], 1.0) for ci in bytime[t]], float(INV[w]))
    A = sp.csr_matrix((V, (I, J)), shape=(len(rhs), len(cost)))
    res = linprog(cost, A_eq=A, b_eq=np.array(rhs), bounds=(0, None), method='highs')
    idle = collections.Counter()
    for (ci, p), j in xi.items():
        if res.x[j] > 1e-6:
            r = cl[ci][0][0] - sum(k * w for k, w in zip(p, widths))
            if r: idle[r] += float(res.x[j])
    if want_x:
        sol = [[cl[ci][0][0], list(cl[ci][0][1]), list(p), float(res.x[j])] for (ci, p), j in xi.items() if res.x[j] > 1e-7]
        return res.fun / NT, {r: round(v / NT, 3) for r, v in sorted(idle.items())}, sol
    return res.fun / NT, {r: round(v / NT, 3) for r, v in sorted(idle.items())}
if __name__ == '__main__':
    arg = sys.argv[1]; t = time.time()
    if arg.startswith('D='):
        D = int(arg[2:]); order = make_order(lambda c: (c[1] + 4 * c[0] + D * (c[0] >= 1), c[0]))
    elif arg.startswith('K='):
        a, D, E = map(int, arg[2:].split(',')); order = make_order(lambda c: (c[1] + a * c[0] + D * (c[0] >= 1) + E * (c[0] == 8), c[0]))
    elif arg == 'shipped':
        order = rd(F / 'new/borrow-fill/ORDER.json')
    else:
        order = json.load(open(arg))
    if arg == 'D=0': print('D=0 equals shipped order:', order == rd(F / 'new/borrow-fill/ORDER.json'))
    if len(sys.argv) > 2:
        v, idle, sol = lp(order, True); json.dump(dict(order=order, solution=sol), open(sys.argv[2], 'w'))
    else:
        v, idle = lp(order)
    print(arg, 'excess/call %.3f' % v, idle, 'secs %.0f' % (time.time() - t), flush=True)
