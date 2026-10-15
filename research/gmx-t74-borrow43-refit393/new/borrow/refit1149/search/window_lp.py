#!/usr/bin/env python3
"""window_lp.py [--out RESULT.json]
Window-aware aggregate packing bound for the shipped owner fill (new/borrow-fill/FILL.json.gz with ORDER.json / DISPATCH.json):
bins (replica, child, port) are grouped by (cap, window); every bin takes one whole-helper pattern; a slot of width w in a bin of
window S must be served by a helper of width w from a call whose time lies in S.  Per-call distinct-family rules are ignored, so the
LP value is a lower bound on the idle excess  sum_r idle_r r (ln(73/r) - ln(73/48))  of every fill that respects the windows of this
order, and it is compared with the shipped fill and with the time-free bound (PR2679's ceil74.py class ceiling).  On PR2679's order
and fill the three numbers were 275.010 (window LP), 275.026 (fill) and 271.593 (time-free) per call (results/WINDOW-LP-PR2679-ORDER.json:
its 170 rank-7 remainders are forced by the windows of that order); on this package's shifted order 272.801 and 272.832 (results/WINDOW-LP.json).  Construction-side evidence (numpy, scipy); not part of the verification.
Prepared by G. V. Chudnovsky with Anthropic Claude assistance. Apache-2.0."""
import argparse, collections, gzip, importlib.util, json, math, sys, time
import numpy as np, scipy.sparse as sp
from scipy.optimize import linprog
from pathlib import Path
F = Path(__file__).resolve().parents[4]          # package root (new/borrow/refit1149/search -> package)
def read(p):
    b = p.read_bytes(); return json.loads(gzip.decompress(b) if p.suffix == '.gz' else b)
def patterns(cap, widths):
    out = []
    def rec(i, rem, cur):
        if i == len(widths): out.append(tuple(cur)); return
        for k in range(rem // widths[i] + 1): rec(i + 1, rem - k * widths[i], cur + [k])
    rec(0, cap, []); return out
def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter); ap.add_argument('--out'); a = ap.parse_args()
    plan = read(F / 'inputs/PLAN.json'); parent = read(F / 'inputs/GLOBAL-TEMPLATE.json.gz'); order = read(F / 'new/borrow-fill/ORDER.json')
    spw = importlib.util.spec_from_file_location('windows', F / 'vendor/build_parent.py'); mod = importlib.util.module_from_spec(spw); spw.loader.exec_module(mod)
    windows, schedule = mod.actual_windows(parent, order)
    prov = read(F / 'new/BORROW-PROVIDERS.json'); inv = collections.Counter(next(iter(prov.values()))['widths'])
    assert all(collections.Counter(d['widths']) == inv for d in prov.values())
    L = lambda r: math.log(73 / r); widths = sorted(inv); NT = len(order)
    classes = collections.Counter()
    for key, (cap, times) in windows.items(): classes[(cap, tuple(sorted(times)))] += 112
    cl = list(classes.items()); pats = {cap: patterns(cap, widths) for cap in sorted({c for c, S in classes})}
    xi = {}; yi = {}; cost = []
    for ci, ((cap, S), n) in enumerate(cl):
        for p in pats[cap]:
            r = cap - sum(k * w for k, w in zip(p, widths)); xi[ci, p] = len(cost); cost.append(0.0 if r == 0 else r * (L(r) - L(48)))
    nx = len(cost)
    for ci, ((cap, S), n) in enumerate(cl):
        for w in widths:
            for t in S: yi[ci, w, t] = nx + len(yi)
    nv = nx + len(yi); cost = np.array(cost + [0.0] * len(yi))
    I, J, V, rhs = [], [], [], []
    def row(entries, b):
        r = len(rhs)
        for j, v in entries: I.append(r); J.append(j); V.append(v)
        rhs.append(b)
    for ci, ((cap, S), n) in enumerate(cl):
        row([(xi[ci, p], 1.0) for p in pats[cap]], n)
        for wi, w in enumerate(widths): row([(xi[ci, p], float(p[wi])) for p in pats[cap] if p[wi]] + [(yi[ci, w, t], -1.0) for t in S], 0.0)
    for t in range(NT):
        for w in widths: row([(yi[ci, w, t], 1.0) for ci, ((cap, S), n) in enumerate(cl) if t in S], float(inv[w]))
    A = sp.csr_matrix((V, (I, J)), shape=(len(rhs), nv)); t0 = time.time()
    res = linprog(cost, A_eq=A, b_eq=np.array(rhs), bounds=(0, None), method='highs'); assert res.status == 0, res.message
    idle = collections.Counter()
    for (ci, p), j in xi.items():
        if res.x[j] > 1e-6:
            r = cl[ci][0][0] - sum(k * w for k, w in zip(p, widths))
            if r: idle[r] += float(res.x[j])
    # shipped fill
    fill = read(F / 'new/borrow-fill/FILL.json.gz'); dispatch = {(s, r): pid for s, r, pid in read(F / 'new/borrow-fill/DISPATCH.json')}
    for s_, r_, p_ in fill['dispatch_overrides']: dispatch[s_, r_] = p_
    used = collections.Counter(); caps = {}
    for key, (cap, times) in windows.items(): r_, ch = map(int, key.split(',')); caps[r_, ch] = cap
    for time_, helper, rep, child, port, offset in fill['rows']:
        stage, replica = order[time_]; used[rep, child, port] += prov[dispatch[stage, replica]]['widths'][helper]
    shipped = collections.Counter()
    for (r_, ch), cap in caps.items():
        for port in range(112):
            rem = cap - used.get((r_, ch, port), 0)
            if rem: shipped[rem] += 1
    ex = lambda idl: sum(v * r * (L(r) - L(48)) for r, v in idl.items()) / NT
    out = dict(classes=len(cl), bins=sum(classes.values()), calls=NT, inventory_per_call={str(k): v for k, v in sorted(inv.items())},
               window_lp_excess_per_call=round(res.fun / NT, 4), window_lp_idle_per_call={str(r): round(v / NT, 4) for r, v in sorted(idle.items())},
               shipped_excess_per_call=round(ex(shipped), 4), shipped_idle={str(r): v for r, v in sorted(shipped.items())},
               time_free_note='time-free pattern optimum of this inventory: 271.593 per call (PR2679 new/borrow-fill/search/ceil74.py class ceiling 0.00171514043461)',
               lp_size=dict(rows=A.shape[0], columns=A.shape[1], nonzeros=int(A.nnz)), seconds=round(time.time() - t0, 1),
               reading='shipped_excess - window_lp_excess is all a fill of this order can still gain; window_lp_excess - 271.593 is what only a different order can gain')
    print(json.dumps(out, indent=1))
    if a.out: Path(a.out).write_text(json.dumps(out, indent=1) + '\n')
if __name__ == '__main__':
    main()
