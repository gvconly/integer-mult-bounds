"""Screen of lattice shells and set systems against the unit yardstick (see screen.py).
For each configuration: decoder matrix N (integer, after shift/scale; diag constant), rank h, the distinct off-diagonal
values s, the polynomial-method bound on v at that rank, lines orthogonal to a line (helper driver), a greedy lower
bound on the largest set of pairwise non-orthogonal lines (complex window >= it), and, when the (J+G) mod 2 label rule
applies (diag even, off-diagonal 0 or odd), c, x = c*h/v and the ceiling."""
import itertools, math, random, sys, time, json
import numpy as np
from fractions import Fraction as Fr
from screen import rank_q, ceiling  # noqa: E402  (also puts the engine on sys.path)
from e8units import f2_rank, hamming_e8

def poly_bound(h, vals, diag):
    """Polynomial-method bound on v at rank h.  f(t) = prod_{a in vals}(t - a) applied entrywise gives f(diag) I (rank v).
    The diagonal is constant, so a monomial t^j with j of the same parity as s = deg f lifts to degree s (multiply by
    (G_ii G_jj / diag^2)^((s-j)/2)); monomials of the other parity lift to degree s - 1.  Hence
    v <= C(h+s-1, s) + [f has a term of the other parity] * C(h+s-2, s-1)."""
    coef = np.poly1d([1.0])
    for a in vals: coef = coef * np.poly1d([1.0, -float(a)])
    cs = coef.coeffs[::-1]          # cs[i] = coefficient of t^i
    s = len(vals)
    other = any(abs(ci) > 1e-9 and (i - s) % 2 for i, ci in enumerate(cs))
    return int(math.comb(h + s - 1, s) + (math.comb(h + s - 2, s - 1) if other else 0))

def greedy_clique(N, tries=200, seed=1):
    rng = random.Random(seed); v = len(N); nz = (N != 0); np.fill_diagonal(nz, False)
    best = 0
    for _ in range(tries):
        order = list(range(v)); rng.shuffle(order)
        cl = []; cand = np.ones(v, dtype=bool)
        for i in order:
            if cand[i]:
                cl.append(i); cand &= nz[i]
        best = max(best, len(cl))
    return best

rows = []
def report(name, N, note='', tries=200):
    N = np.asarray(N); v = len(N); off = ~np.eye(v, dtype=bool)
    d = set(np.diag(N).tolist()); assert len(d) == 1, d
    diag = d.pop(); vals = sorted(set(N[off].tolist()))
    h = rank_q(N.tolist()) if v <= 1500 else int(np.linalg.matrix_rank(N.astype(float)))
    orth = int(((N == 0) & off).sum(1).max())
    om = greedy_clique(N, tries)
    pb = poly_bound(h, vals, diag)
    r = dict(name=name, v=v, h=h, s=len(vals), vals=[str(x) for x in vals], bound=pb, fill=round(v / pb, 3), orth=orth, omega_lb=om, note=note)
    isint = all(float(x).is_integer() for x in vals) and float(diag).is_integer()
    if isint and int(diag) % 2 == 0 and all(int(x) == 0 or int(x) % 2 for x in vals):
        B = ((1 + N.astype(np.int64)) % 2).astype(np.uint8); c = f2_rank(B)
        x = Fr(c * h, v); b = ceiling(h, x)
        r.update(rule='J+G', c=c, x=float(x), X=(float(b[0]) if b else 0.0), s_best=(b[1] if b else None))
    else:
        r.update(rule='dyadic only')
    rows.append(r)
    print('%-38s v=%6d h=%3d s=%d vals=%-22s bound=%7d fill=%.3f orth=%6d omega>=%4d %s' % (
        name, v, h, len(vals), ','.join(r['vals'])[:22], pb, v / pb, orth, om,
        ('c=%d x=%.3f X=%.3e' % (r['c'], r['x'], r['X'])) if 'c' in r else r['rule']), note, flush=True)

def lines_of(vecs):
    seen, out = set(), []
    for x in map(tuple, vecs):
        if x in seen or tuple(-a for a in x) in seen: continue
        seen.add(x); out.append(x)
    return np.array(out, dtype=np.int64)

def roots_D(n):
    return lines_of([[ (s if k == i else (t if k == j else 0)) for k in range(n)] for i, j in itertools.combinations(range(n), 2) for s in (1, -1) for t in (1, -1)])

t0 = time.time()
L, G = hamming_e8(); report('E8 root lines', G)
E7 = [i for i in range(120) if G[0, i] == 0]; report('E7 (lines orthogonal to a root)', G[np.ix_(E7, E7)])
for n in (8, 10, 12, 16):
    R = roots_D(n); report('D%d root lines' % n, R @ R.T)
E = np.zeros((240, 16), dtype=np.int64); E[:120, :8] = L; E[120:, 8:] = L
report('E8 + E8 (decomposable)', (E @ E.T) // 2)
# upstream triples (Lorentzian norm-2 vectors (t_S; 1) of I_{h,1})
for hh in (9, 10, 12, 16, 20, 24):
    T = list(itertools.combinations(range(hh), 3)); M = np.zeros((len(T), hh), dtype=np.int64)
    for i, t in enumerate(T): M[i, list(t)] = 1
    report('triples h=%d (I_{h,1} roots)' % hh, M @ M.T - 1, tries=40)
# Barnes-Wall BW16: x in Z^16, x mod 2 in RM(1,4), sum x = 0 mod 4; minimal norm 8 (4320 vectors)
gen = [[1]*16] + [[(i >> b) & 1 for i in range(16)] for b in range(4)]
RM = {tuple(np.array(c) @ np.array(gen) % 2) for c in itertools.product((0, 1), repeat=5)}
vecs = []
for i, j in itertools.combinations(range(16), 2):
    for s, t in itertools.product((2, -2), repeat=2):
        x = [0]*16; x[i] = s; x[j] = t; vecs.append(x)
for cw in RM:
    if sum(cw) != 8: continue
    sup = [k for k in range(16) if cw[k]]
    for sg in itertools.product((1, -1), repeat=8):
        if sg.count(-1) % 2 == 0:
            x = [0]*16
            for k, s in zip(sup, sg): x[k] = s
            vecs.append(x)
BW = lines_of(vecs); assert len(BW) == 2160, len(BW)
Gb = BW @ BW.T; assert set(np.diag(Gb)) == {8}
report('Barnes-Wall BW16 minimal lines', Gb // 2, note='(Gram/2: norm 4)', tries=60)
# Golay code octads (759 blocks of S(5,8,24)) as 0/1 vectors; decoder (G - 4J)/4 has diagonal 1
gg = ['100000000000110111000101','010000000000011011100011','001000000000101101110001','000100000000010110111001',
      '000010000000001011011101','000001000000000101101111','000000100000100010110111','000000010000110001011011',
      '000000001000111000101101','000000000100011100010111','000000000010101110001011','000000000001111111111110']
Gg = np.array([[int(ch) for ch in r] for r in gg])
octads = []
for c in itertools.product((0, 1), repeat=12):
    w = np.array(c) @ Gg % 2
    if w.sum() == 8: octads.append(w)
O = np.array(octads, dtype=np.int64); assert len(O) == 759, len(O)
Go = O @ O.T; report('Golay octads S(5,8,24)', Go - 4, note='(4N = G - 4J)', tries=60)
# Steiner triple system: lines of PG(3,2) (35 lines, 15 points), and of PG(4,2)
for n in (4, 5):
    pts = list(range(1, 2**n)); lines = sorted({tuple(sorted((a, b, a ^ b))) for a in pts for b in pts if a < b})
    M = np.zeros((len(lines), len(pts)), dtype=np.int64)
    for i, l in enumerate(lines):
        for p in l: M[i, p - 1] = 1
    report('lines of PG(%d,2) (STS)' % (n - 1), M @ M.T - 1)
json.dump(rows, open(__import__('pathlib').Path(__file__).resolve().parent / 'families.json', 'w'), indent=1)
print('%.0f s' % (time.time() - t0))
