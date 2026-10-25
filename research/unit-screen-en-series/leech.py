"""Leech lattice minimal vectors (scaled by sqrt 8: norm 32, inner products 0, +-8, +-16) -> 98280 lines, rank 24.
Reports the polynomial-bound fill, orthogonal lines per line, the label-rule check, a greedy clique of pairwise
non-orthogonal lines (complex window >= it) and the bit-type cover number of a line's orthogonal readers."""
import itertools, math, random, numpy as np, time
t0 = time.time()
gg = ['100000000000110111000101','010000000000011011100011','001000000000101101110001','000100000000010110111001',
      '000010000000001011011101','000001000000000101101111','000000100000100010110111','000000010000110001011011',
      '000000001000111000101101','000000000100011100010111','000000000010101110001011','000000000001111111111110']
Gg = np.array([[int(ch) for ch in r] for r in gg])
code = np.array([np.array(c) @ Gg % 2 for c in itertools.product((0, 1), repeat=12)])
V = []
for i, j in itertools.combinations(range(24), 2):                       # (+-4, +-4, 0^22)
    for s in (4, -4):
        x = np.zeros(24, int); x[i] = 4; x[j] = s; V.append(x)           # one of each +-pair
for w in code[code.sum(1) == 8]:                                         # (+-2)^8 on octads, even # of minus
    sup = np.flatnonzero(w)
    for sg in itertools.product((1, -1), repeat=7):
        if sg.count(-1) % 2 == 0:
            x = np.zeros(24, int); x[sup] = 2 * np.array((1,) + sg); V.append(x)
for w in code:                                                           # (-3, 1^23) + 4 * codeword pattern
    if w[0] == 0:                                                        # one of each +-pair: fix sign via w[0]
        for i in range(24):
            x = np.where(w == 1, -1, 1).astype(int); x = x.copy()
            # vector: (1^24) - 2*w (signs) then position i gets -3 (sign flipped relative)
            y = x.copy(); y[i] = -3 * x[i]; V.append(y)
L = np.array(V); print('lines', len(L), 'norms', set((L * L).sum(1).tolist()))
ip = L @ L[0]; vals = sorted(set(ip.tolist()) - {32})
print('inner products with line 0:', {int(a): int((ip == a).sum()) for a in sorted(set(ip.tolist()))})
orth = int((ip == 0).sum()); print('orthogonal lines per line:', orth)
h = 24; s = len(set(abs(a) for a in vals if a != 0)) * 2 + 1
bound = math.comb(h + 4, 5)   # f = t(t^2-64)(t^2-256) is odd, so it homogenizes to degree 5 (constant norm)
print('distinct off-diagonal values', s, ' polynomial bound', bound, ' fill %.3f' % (len(L) / bound))
print('label rule (J+G) applies?', all((a // 8) % 2 == 1 or a == 0 for a in vals), '(inner products /8 must be 0 or odd)')
# greedy cliques of pairwise non-orthogonal lines, seeded inside the 4600 lines at inner product 16 with line 0
rng = random.Random(1); best = 0
cap = np.flatnonzero(np.abs(ip) == 16); cap = np.concatenate([[0], cap])
Lc = L[cap]; NZ = (Lc @ Lc.T) != 0
for _ in range(300):
    order = list(range(len(cap))); rng.shuffle(order); cl = []; cand = np.ones(len(cap), bool)
    for i in order:
        if cand[i]: cl.append(i); cand &= NZ[i]
    best = max(best, len(cl))
print('greedy clique of pairwise non-orthogonal lines >=', best)
# bit-type cover of line 0's orthogonal readers by sections w^perp, w over a sample of lines
R = L[ip == 0]; samp = L[rng.sample(range(len(L)), 3000)]
Z = (R @ samp.T == 0); rem = np.ones(len(R), bool); sizes = []
while rem.any() and len(sizes) < 40:
    g = (Z & rem[:, None]).sum(0); j = int(g.argmax())
    if g[j] == 0: sizes.append('uncovered %d' % rem.sum()); break
    sizes.append(int(g[j])); rem &= ~Z[:, j]
print('greedy cover of the %d orthogonal readers (sampled hyperplanes): %d sections, sizes %s' % (len(R), len(sizes), sizes[:12]))
print('%.0f s' % (time.time() - t0))
