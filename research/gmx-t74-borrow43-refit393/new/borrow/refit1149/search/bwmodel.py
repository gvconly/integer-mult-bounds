"""Frame-assignment model of a gmx BORROW-WORD (T74 lineage), following new/borrow/audit_local.py run_borrow section 4:
time order = A (prelude warps/target writes first) -> A/B boundary pins (centres exact, targets at their levels) -> B with the
K_G warps inserted before B[p] -> restore -> final climbs to the end frames.  A warp pins the lender's frame to U just before
it (the last framed gate on the lender before the warp is fixed at U) and restarts the lender's chain at E.
End frames of restored helpers (helper_final) may be left free (dims lo..8) for the end-lowering study."""
import gzip, json
from collections import defaultdict, Counter
import framemodel as xm
from framemodel import basis, contained, join, meet, interval, FULL, H

def registers(e):
    return tuple(e[1:3]) if e[0] == 'add' else tuple(e[1:2])

class BWModel(xm.Model):
    def __init__(self, w, free_ends=(), end_domain=None, free_warps=False):
        rep = w['report']; self.v, self.R = v, R = rep['v'], rep['R']; assert rep['h'] == H
        self.n = n = 2 * v + R; self.w = w
        removed = set(w['removed']); live = [r for r in range(n) if r not in removed]; self.live = live
        end = [FULL] * v + [basis([x for x in range(1, 1 << H) if (x & p).bit_count() % 2 == 0]) for p in w['ports']] + [FULL] * R
        low = {int(q): basis(U) for q, U in w.get('helper_final', {}).items()}
        for q, U in low.items(): end[q] = U
        start = [(p,) for p in w['ports']] + [()] * (v + R)
        for q, U in w.get('start_extra', {}).items(): start[int(q)] = basis(U)
        self.start, self.end, self.low = start, end, low
        lev = {int(t): basis(U) for t, U in w.get('target_levels', {}).items()}
        pins = {q: basis(U) for q, U in zip(w['centers'], w['center_frames'])}
        for t in range(v): pins[v + t] = lev.get(t, ())
        self.ceil = {}
        self.centre_children = []
        for q, U in zip(w['centers'], w['center_frames']):
            U = basis(U); served = {lev.get(e[1] - v, ()) for e in w['cut'] if e[2] == q}
            for T in sorted(served, key=lambda x: -len(x)):
                assert contained(T, U)
                if T != U: self.centre_children.append(len(U) - len(T))
                U = T
        # ---- time-ordered events: ('gate', sec, e) framed, ('warp', r, U, E), ('pins',)
        ev = []
        npre = w.get('prelude_events', 0)
        for e in w['A'][:npre]:
            if e[0] == 'warp': ev.append(('warp', e[1], basis(e[3]), basis(e[4])))
            else: ev.append(('gate', 'A', e))
        for e in w['A'][npre:]: ev.append(('gate', 'A', e))
        ev.append(('pins',))
        byp = defaultdict(list)
        for op, rec in enumerate(w['warps'], 1): byp[rec[0]].append(op)
        B = w['B']
        for p in range(len(B) + 1):
            for op in byp.get(p, []):
                _, lender, cols, U, E = w['warps'][op - 1]; ev.append(('warp', lender, basis(U), basis(E)))
            if p < len(B): ev.append(('gate', 'B', B[p]))
        for e in w['restore']: ev.append(('gate', 'T', e))
        self.events = ev
        gates = []; self.given = []; self.regs = []
        fixed = {}; chains = defaultdict(list); order = []
        for q in live: fixed[('S', q)] = start[q]; chains[q].append(('S', q)); order.append(('S', q))
        last_gate_node = {}   # register -> its last gate node so far (to pin at a warp)
        warpcount = 0; breaks = set()   # (pred, W) pairs: no edge, the lender's frame is re-set by the warp
        couplings = []; wgiven = {}; sec_of = {}; wfree = set(); wsucc = []
        for t, item in enumerate(ev):
            if item[0] == 'pins':
                for q, U in pins.items():
                    fixed[('P', q)] = U; chains[q].append(('P', q)); order.append(('P', q)); last_gate_node[q] = None
                continue
            if item[0] == 'warp':
                _, r, U, E = item
                prev = chains[r][-1]
                node = ('W', warpcount); warpcount += 1
                if prev[0] == 'G' and free_warps and sec_of[prev] == 'B':
                    # free warp: lender's last gate frame U' and the warp target E' are variables with dim U' = dim E';
                    # E' is any dim-equal subspace of the borrower's first gate frame (pseudo label = its dimension)
                    couplings.append((prev, node)); wgiven[node] = E; wgiven[prev] = U; wfree.add(node)
                else:
                    if prev[0] == 'G':
                        if prev in fixed: assert fixed[prev] == U, ('warp pin conflict', prev, fixed[prev], U)
                        fixed[prev] = U
                    else:
                        assert fixed[prev] == U, ('warp from a non-gate node with wrong frame', prev, fixed[prev], U)
                    fixed[node] = E
                breaks.add((prev, node)); chains[r].append(node); order.append(node)
                continue
            _, sec, e = item
            i = len(gates); gates.append((sec, e)); self.given.append(basis(e[-1])); rs = registers(e); self.regs.append(rs)
            assert len(set(rs)) == len(rs)
            node = ('G', i); order.append(node); sec_of[node] = sec
            for q in rs: chains[q].append(node)
        self.gates = gates
        self.nA = sum(1 for s, e in gates if s == 'A'); self.nB = sum(1 for s, e in gates if s == 'B')
        self.free_ends = set(free_ends)
        for q in live:
            chains[q].append(('E', q)); order.append(('E', q))
            if q not in self.free_ends: fixed[('E', q)] = end[q]
        # chains must be ordered consistently with `order`
        for q, ch in chains.items():
            for a, b in zip(ch, ch[1:]):
                if a in wfree: breaks.add((a, b)); wsucc.append((a, b))
        self.wsucc = wsucc; self.wfree = wfree
        self.chains = chains; self.fixed_frames = fixed
        position = {node: k for k, node in enumerate(order)}
        pred = defaultdict(list); succ = defaultdict(list)
        for q, ch in chains.items():
            for a, b in zip(ch, ch[1:]):
                assert position[a] < position[b], (a, b)
                if (a, b) in breaks: continue
                succ[a].append((q, b)); pred[b].append((q, a))
        lo, hi = {}, {}
        for node in order:
            if node in fixed:
                lo[node] = fixed[node]
                for q, a in pred[node]:
                    assert contained(lo[a], fixed[node]), ('fixed frame below forward join', node, a, lo[a], fixed[node])
            else:
                lo[node] = join(lo[a] for q, a in pred[node]) if pred[node] else ()
        for node in reversed(order):
            if node in fixed: hi[node] = fixed[node]
            else:
                U = FULL
                for q, b in succ[node]: U = meet(U, hi[b])
                if node[0] == 'E' and end_domain is not None:
                    U = meet(U, end_domain.get(node[1], FULL))
                hi[node] = U
            assert contained(lo[node], hi[node]), ('empty domain', node)
        self.lo, self.hi = lo, hi
        self.couplings = couplings; self.wgiven = wgiven
        self.domain = {}
        gsucc = dict(wsucc)
        for (prev, node) in couplings:
            g = gsucc[node]
            dims = [d for d in range(len(lo[prev]), len(hi[prev]) + 1) if d <= len(hi[g])]
            self.domain[node] = [('d',) * d for d in dims]
            assert self.domain[node], ('empty W domain', node)
        for node in order:
            if node in fixed or node in self.domain: continue
            if lo[node] != hi[node]:
                if False:
                    # enumerate only dims of the full quotient would explode; restrict to subspaces built from lo + up to 4 extra
                    self.domain[node] = interval(lo[node], hi[node]) if len(hi[node]) - len(lo[node]) <= 4 else None
                else:
                    self.domain[node] = interval(lo[node], hi[node])
        mult = defaultdict(int)
        for q, ch in chains.items():
            for a, b in zip(ch, ch[1:]):
                if (a, b) in breaks: continue
                mult[(a, b)] += 1
        self.edges = dict(mult)

    def components(self):
        parent = {node: node for node in self.domain}
        def find(x):
            while parent[x] != x:
                parent[x] = parent[parent[x]]; x = parent[x]
            return x
        for (a, b) in list(self.edges) + list(self.couplings) + list(self.wsucc):
            if a in self.domain and b in self.domain: parent[find(a)] = find(b)
        groups = defaultdict(list)
        for node in self.domain: groups[find(node)].append(node)
        comps = [sorted(g, key=lambda x: (x[0], x[1])) for g in groups.values()]
        return sorted(comps, key=lambda c: c[0])

    def component_problem(self, comp):
        labels, unary, pairs = super().component_problem(comp)
        members = set(comp); zero = [0] * (H + 1)
        for (a, b) in self.couplings:
            if a in members and b in members:
                pairs[(a, b)] = {(Fa, Fb): list(zero) for Fa in labels[a] for Fb in labels[b] if len(Fa) == len(Fb)}
                assert pairs[(a, b)], ('no dim-compatible warp pair', a, b)
        for (a, b) in self.wsucc:
            if a in members:
                Fbs = labels[b] if b in members else [self.frame_of(b, None)]
                tab = {}
                for Fa in labels[a]:
                    for Fb in Fbs:
                        if len(Fa) <= len(Fb):
                            c = list(zero); d = len(Fb) - len(Fa)
                            if d: c[d] += 1
                            tab[(Fa, Fb)] = c
                if b in members: pairs[(a, b)] = tab
                else:
                    for Fa in labels[a]:
                        if (Fa, Fbs[0]) in tab:
                            for d, x in enumerate(tab[(Fa, Fbs[0])]): unary[a][Fa][d] += x
                        else:
                            unary[a][Fa][0] += 10 ** 6   # infeasible label (dimension above the forced successor)
        return labels, unary, pairs

    def profile(self, assignment):
        count = defaultdict(int, super().profile(assignment))
        for (a, b) in self.wsucc:
            d = len(self.frame_of(b, assignment)) - len(self.frame_of(a, assignment))
            assert d >= 0, ('warp target above the first gate', a, b)
            if d: count[d] += 1
        return dict(sorted(count.items()))

    def baseline(self):
        asg = {('G', i): self.given[i] for i in range(len(self.gates)) if ('G', i) in self.domain}
        for q in self.free_ends:
            if ('E', q) in self.domain: asg[('E', q)] = self.end[q]
        for (prev, node) in self.couplings:
            if node in self.domain: asg[node] = ('d',) * len(self.wgiven[node])
        return asg

    def widths(self, assignment):
        return {q: len(self.frame_of(('E', q), assignment)) - len(self.start[q]) for q in range(2 * self.v, self.n) if q in set(self.live)}

if __name__ == '__main__':
    import sys
    path = sys.argv[1]
    w = json.load(gzip.open(path))
    M = BWModel(w)
    base = M.baseline(); P = M.profile(base)
    print('paid profile (baseline frames)', P, 'reported', w['report']['local_profile'])
    print('gates', len(M.gates), 'free gate nodes', len(M.domain), 'fixed', len(M.fixed_frames))
    big = [nd for nd, d in M.domain.items() if d is None]
    print('nodes with quotient > 4:', len(big))
    hist = Counter(len(M.hi[nd]) - len(M.lo[nd]) for nd in M.domain)
    print('slack histogram', dict(sorted(hist.items())))
