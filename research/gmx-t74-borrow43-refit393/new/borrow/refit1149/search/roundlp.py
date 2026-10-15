"""roundlp.py SOL.json OUT.json: integer pattern counts per window class from the (near-integral: 11 of 666 classes
fractional) LP solution of orderlp.py; integral classes kept, fractional classes re-solved as a small MILP that keeps every
width total.  Construction-side.  Prepared by G. V. Chudnovsky with Anthropic Claude assistance. Apache-2.0."""
import json, collections, math, sys
import numpy as np
from scipy.optimize import milp, LinearConstraint, Bounds
import orderlp as O
d = json.load(open(sys.argv[1])); W = O.widths; L = O.L
windows, _ = O.mod.actual_windows(O.parent, d['order'])
classes = collections.Counter()
for k, (cap, ts) in windows.items(): classes[(cap, tuple(sorted(ts)))] += 112
fixed = collections.Counter(); free = []
byclass = collections.defaultdict(list)
for cap, S, p, x in d['solution']: byclass[(cap, tuple(S))].append((tuple(p), x))
for key, lst in byclass.items():
    if all(abs(x - round(x)) < 1e-6 for p, x in lst):
        for p, x in lst:
            if round(x): fixed[(key, p)] += round(x)
    else:
        for p, x in lst: free.append((key, p, x))
fk = sorted({k for k, p, x in free})
have = {(k, p) for k, p, x in free}
for k in fk:
    for p in O.PATS[k[0]]:
        r = k[0] - sum(a * b for a, b in zip(p, W))
        if r in (0, 1, 3, 5, k[0]) and (k, p) not in have: free.append((k, p, 0.0))
print('free classes', len(fk), 'vars', len(free))
# variables: free (key,p) counts; constraints: per class total = n - fixed in class; per width totals preserved
c = []; A = []; lo = []; hi = []
nv = len(free)
for key in fk:
    row = [1.0 if k == key else 0.0 for k, p, x in free]; n = classes[key] - sum(v for (kk, p), v in fixed.items() if kk == key)
    A.append(row); lo.append(n); hi.append(n)
for wi, w in enumerate(W):
    row = [float(p[wi]) for k, p, x in free]; tot = sum(p[wi] * x for k, p, x in free)
    A.append(row); lo.append(round(tot)); hi.append(round(tot)); print('width', w, 'fractional total', tot)
cost = [ (k[0] - sum(a * b for a, b in zip(p, W))) for k, p, x in free]
cost = [r * (L(r) - L(48)) if r else 0.0 for r in cost]
# stay near LP: bounds floor/ceil
res = milp(np.array(cost), constraints=LinearConstraint(np.array(A), lo, hi), integrality=np.ones(nv),
           bounds=Bounds(0, [classes[k] for k, p, x in free]))
print('round status', res.status, res.message)
sol = [[key[0], list(key[1]), list(p), v] for (key, p), v in fixed.items()] + [[k[0], list(k[1]), list(p), int(round(res.x[i]))] for i, (k, p, x) in enumerate(free) if round(res.x[i])]
idle = collections.Counter()
for cap, S, p, x in sol:
    r = cap - sum(a * b for a, b in zip(p, W))
    if r: idle[r] += x
NT = len(d['order']); obj = sum(x * (r * (L(r) - L(48))) for r, x in idle.items()) / NT
print('rounded idle', dict(sorted(idle.items())), 'excess/call %.4f' % obj)
json.dump(dict(order=d['order'], classes=[[k[0], list(k[1]), n] for k, n in classes.items()], solution=sol, obj=obj), open(sys.argv[2], 'w'))
