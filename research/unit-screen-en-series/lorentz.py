"""Search for dense {0,+-1}-Gram units among the roots of the Lorentzian lattice I_{n,1} = Z^{n,1}.

A unit compatible with the (J+G) mod 2 label rule has norm 2 and pairwise inner products in {0,+-1} (lines).  Such
sets satisfy v <= C(h+2,3) + h at rank h (Hadamard-cube rank argument).  Both E8 (inside I_{9,1} = E8 + I_{1,1}) and the
upstream triples (vectors (t_S; 1)) live here.  Roots (x; m) with |x|^2 - m^2 = 2, height m in {0,1,2}; the
compatibility graph joins two roots when |<u,w>| <= 1; large cliques are found by randomized greedy plus a
1-for-2 swap local search on bitsets.  usage: lorentz.py n [heights] [seconds]"""
import itertools, math, random, sys, time, json
import numpy as np


def roots(n, heights):
    out = []
    for m in heights:
        target = 2 + m * m
        # integer vectors x in Z^n with |x|^2 = target, entries in [-2, 2]
        for sup_size in range(1, target + 1):
            for sup in itertools.combinations(range(n), sup_size):
                for mags in itertools.product((1, 2), repeat=sup_size):
                    if sum(a * a for a in mags) != target:
                        continue
                    for sg in itertools.product((1, -1), repeat=sup_size):
                        x = [0] * n
                        for k, a, s in zip(sup, mags, sg):
                            x[k] = a * s
                        if m == 0 and next(a for a in x if a) < 0:
                            continue                      # one representative of each +- pair at height 0
                        out.append(tuple(x) + (m,))
    return np.array(out, dtype=np.int64)


def clique_search(adj_bits, nv, seconds, seed=1):
    rng = random.Random(seed)
    best = []
    t0 = time.time()
    order = list(range(nv))
    while time.time() - t0 < seconds:
        rng.shuffle(order)
        cl, cand = [], (1 << nv) - 1
        for i in order:
            if cand >> i & 1:
                cl.append(i)
                cand &= adj_bits[i]
        # local search: try to drop one vertex and add two
        improved = True
        while improved and time.time() - t0 < seconds:
            improved = False
            S = set(cl)
            for out_v in rng.sample(cl, len(cl)):
                rest = [u for u in cl if u != out_v]
                cand = (1 << nv) - 1
                for u in rest:
                    cand &= adj_bits[u]
                cand &= ~(1 << out_v)
                adds = []
                c2 = cand
                while c2:
                    i = (c2 & -c2).bit_length() - 1
                    c2 &= c2 - 1
                    adds.append(i)
                rng.shuffle(adds)
                for a in adds:
                    rem = cand & adj_bits[a] & ~(1 << a)
                    if rem:
                        b = (rem & -rem).bit_length() - 1
                        cl = rest + [a, b]
                        cand2 = cand & adj_bits[a] & adj_bits[b]
                        while cand2:                      # greedy extend further
                            i = (cand2 & -cand2).bit_length() - 1
                            cl.append(i)
                            cand2 &= adj_bits[i] & ~(1 << i)
                        improved = True
                        break
                if improved:
                    break
        if len(cl) > len(best):
            best = list(cl)
            print('  clique %d at %.0fs' % (len(best), time.time() - t0), flush=True)
    return best


if __name__ == '__main__':
    n = int(sys.argv[1])
    heights = tuple(int(c) for c in sys.argv[2]) if len(sys.argv) > 2 else (0, 1)
    secs = float(sys.argv[3]) if len(sys.argv) > 3 else 120
    R = roots(n, heights)
    Q = np.diag([1] * n + [-1]).astype(np.int64)
    nv = len(R)
    RQ = (R @ Q).astype(np.int16)
    Rt = R.T.astype(np.int16)
    adj = []
    for s0 in range(0, nv, 2048):
        Gm = RQ[s0:s0 + 2048] @ Rt
        ok = np.abs(Gm) <= 1
        for k in range(len(ok)):
            ok[k, s0 + k] = False
        P = np.packbits(ok, axis=1, bitorder='little')
        adj += [int.from_bytes(r.tobytes(), 'little') for r in P]
    print('n=%d heights=%s roots(lines)=%d' % (n, heights, nv), flush=True)
    best = clique_search(adj, nv, secs)
    C = R[best]
    json.dump(dict(n=n, heights=heights, vectors=C.tolist()), open('lorentz_n%d_h%s.json' % (n, ''.join(map(str, heights))), 'w'))
    print('best', len(best))
