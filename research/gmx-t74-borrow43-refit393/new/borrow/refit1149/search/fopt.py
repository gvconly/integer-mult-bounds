"""Exact frame optimisation (HiGHS MILP over the local polytope of each component) with optional scalar node prices."""
import numpy as np, scipy.sparse as sp
from scipy.optimize import milp, LinearConstraint, Bounds
import framemodel as x4model

def weights(b, m=73):
    return [0.0] + [d * (m / d) ** b for d in range(1, 10)]

def solve_component(M, labels, unary, pairs, wt, price):
    rows, rhs, cols = x4model.lp(labels, unary, pairs)
    nr, nc = len(rows), len(cols)
    c = np.zeros(nc)
    for j, col in enumerate(cols):
        c[j] = sum(x * wt[d] for d, x in enumerate(col[1]))
        if col[0][0] == 'x' and price is not None: c[j] += price(col[0][1], col[0][2])
    I, J, V = [], [], []
    for j, col in enumerate(cols):
        for i, val in col[2]: I.append(i); J.append(j); V.append(val)
    A = sp.csr_matrix((V, (I, J)), shape=(nr, nc))
    res = milp(c, constraints=LinearConstraint(A, rhs, rhs), integrality=np.ones(nc), bounds=Bounds(0, 1), options=dict(time_limit=600))
    assert res.status == 0, res.message
    asg = {}
    for j, col in enumerate(cols):
        if col[0][0] == 'x' and res.x[j] > 0.5: asg[col[0][1]] = col[0][2]
    return asg

def optimise(M, b, price=None):
    wt = weights(b); asg = {}
    for comp in M.components():
        labels, unary, pairs = M.component_problem(comp)
        if len(comp) == 1:
            node = comp[0]
            best = min(labels[node], key=lambda F: sum(x * wt[d] for d, x in enumerate(unary[node][F])) + (price(node, F) if price else 0.0))
            asg[node] = best
        else:
            asg.update(solve_component(M, labels, unary, pairs, wt, price))
    return asg
