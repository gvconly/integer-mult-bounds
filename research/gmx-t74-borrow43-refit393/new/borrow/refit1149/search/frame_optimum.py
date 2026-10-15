#!/usr/bin/env python3
"""frame_optimum.py WORD.json.gz [--out RESULT.json]
Exact frame optimum of a BORROW-WORD under its own frame contract (new/borrow/audit_local.py run_borrow section 4): the
scalar program, warps, entrances, ends, centre pins and target levels are fixed; every gate frame is free in its
propagated interval [lo, hi]; the objective is the first-order excess moment  sum_d c_d * (d (73/d)^b - d (73/48)^b)
(b = 0.001715), i.e. the paid children priced against hosting the same mass as rank-48 idle.  The problem decomposes into
connected components of gates sharing registers and each component is solved exactly (HiGHS MILP over the local polytope).
Reports the optimum profile, the given profile, and the minimal set of gate-frame changes that realises the optimum
(a tie-break of 1e-4 per changed gate keeps every cost-neutral relabelling at its given frame).
Also run with --free-ends to let the end frames of restored helpers move (lowering only, inside their shipped ends) and
with --free-warps to let every B-phase warp's frames move (lender frame and target frame of equal dimension).
Construction-side evidence (needs numpy and scipy); not part of the verification.  Apache-2.0."""
import argparse, gzip, json, math, sys, time
from collections import Counter
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
import bwmodel, fopt
B = 0.001715
L = lambda r: math.log(73 / r)
fopt.weights = lambda b, m=73: [0.0] + [d * (73 / d) ** b - d * (73 / 48) ** b for d in range(1, 10)]
def excess(P): return sum(n * r * (L(r) - L(48)) for r, n in P.items())
def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('word'); ap.add_argument('--out'); ap.add_argument('--free-ends', action='store_true'); ap.add_argument('--free-warps', action='store_true')
    a = ap.parse_args()
    w = json.load(gzip.open(a.word))
    low = {int(q): bwmodel.basis(U) for q, U in w['helper_final'].items()}
    t = time.time()
    M = bwmodel.BWModel(w, free_ends=low.keys() if a.free_ends else (), end_domain=low if a.free_ends else None, free_warps=a.free_warps)
    base = M.baseline(); P0 = M.profile(base)
    assert P0 == {int(k): v for k, v in w['report']['local_profile'].items()}, 'model does not reproduce the reported profile'
    price = lambda node, F: 0.0 if F == base.get(node) else 1e-4
    asg = fopt.optimise(M, B, price); P1 = M.profile(asg)
    changes = []
    for nd in sorted(asg, key=str):
        if asg[nd] != base.get(nd):
            if nd[0] == 'G':
                sec, e = M.gates[nd[1]]; idx = nd[1] - (0 if sec == 'A' else M.nA if sec == 'B' else M.nA + M.nB)
                if sec == 'A': idx = sum(1 for i in range(nd[1]) if M.gates[i][0] == 'A')   # index among framed A gates (prelude warps excluded)
                changes.append(dict(section=sec, index=idx, gate=e[:-1], given=list(base[nd]), optimum=list(asg[nd])))
            else:
                changes.append(dict(node=str(nd), given=str(base.get(nd)), optimum=str(asg[nd])))
    comps = M.components()
    res = dict(word=Path(a.word).name, given_profile={str(k): v for k, v in sorted(P0.items())}, optimum_profile={str(k): v for k, v in sorted(P1.items())},
               given_excess=round(excess(P0), 4), optimum_excess=round(excess(P1), 4), gain_excess=round(excess(P0) - excess(P1), 4),
               free_gate_nodes=len([n for n in M.domain if n[0] == 'G']), components=len(comps), largest_component=max(len(c) for c in comps),
               free_ends=a.free_ends, free_warps=a.free_warps, inventory={str(k): v for k, v in sorted(Counter(M.widths(asg).values()).items())},
               minimal_changes=changes, seconds=round(time.time() - t, 1),
               objective='first-order excess moment sum_d c_d d((73/d)^b - (73/48)^b), b = 0.001715; exact per component (HiGHS MILP, local polytope)')
    print(json.dumps(res, indent=1))
    if a.out: Path(a.out).write_text(json.dumps(res, indent=1) + '\n')
if __name__ == '__main__':
    main()
