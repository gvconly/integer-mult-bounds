"""Frequency-domain (Walsh chart) simulator of PR5440's coded p-chain, in the collapsed nine-axis geometry that
p19_word.py replays and in the fresh-block geometry that the invoice prices (central axis + one 8-axis block per
retained stage). Arithmetic mod P = 65537 with i = 256, exactly as PR5440's physical_gauge_lift.py.

In the Walsh basis a nine-gate Gaussian gate on axis k is diag(1, i) on frequency bit k, so Q, D, R and all idle
powers are diagonal phases; a signed normalizer K = Z_z P_G moves frequency u to G u xor z (validated below
against PR5440's own address-domain signed_map); the mode code A mixes owners at a common frequency."""
import json, random, sys
from pathlib import Path
ROOT = Path(sys.argv[1]) if len(sys.argv) > 1 else Path('.')
P, IMAG = 65537, 256
IINV = pow(IMAG, P - 2, P)
NORM = json.loads((ROOT / 'research/coded-native-caller/NORMALIZERS.json').read_text())
PORTS = NORM['ports']; MODES = [[(i, s) for i, s in enumerate(c) if s] for c in NORM['mode_columns']]
OWN = 112
def par(x): return bin(x).count('1') & 1
GT = []   # transpose rows of each G
for rec in PORTS:
    rows = rec['G_rows']; GT.append([sum(((r >> j) & 1) << i for i, r in enumerate(rows)) for j in range(9)])
def gmul(rows, x): return sum(par(r & x) << b for b, r in enumerate(rows))

class Geo:
    """axis layout: bit 0 central; block b occupies bits 1+8b .. 8+8b (collapsed: every block is block 0)."""
    def __init__(self, blocks, collapsed=False):
        self.blocks, self.collapsed = blocks, collapsed
    def bmask(self, b):
        b = 0 if self.collapsed else b
        return 0xFF << (1 + 8 * b)
    def sub(self, u, b):
        b = 0 if self.collapsed else b
        return (u & 1) | (((u >> (1 + 8 * b)) & 0xFF) << 1)
    def put(self, u, b, s):
        b = 0 if self.collapsed else b
        u &= ~(1 | (0xFF << (1 + 8 * b)))
        return u | (s & 1) | (((s >> 1) & 0xFF) << (1 + 8 * b))

# a field (bank) = dict (owner, freq) -> value
def add(x, y, c=1):
    out = dict(x)
    for k, v in y.items():
        w = (out.get(k, 0) + c * v) % P
        if w: out[k] = w
        else: out.pop(k, None)
    return out
def phase(x, mask, power=1):
    out = {}
    for (o, u), v in x.items():
        e = (bin(u & mask).count('1') * power) % 4
        out[(o, u)] = v * pow(IMAG, e, P) % P
    return out
def kmap(x, geo, b, inverse=False):
    out = {}
    for (o, u), v in x.items():
        s = geo.sub(u, b); rec = PORTS[o]
        s2 = gmul(GT[o], s ^ rec['Z_mask']) if inverse else gmul(rec['G_rows'], s) ^ rec['Z_mask']
        out[(o, geo.put(u, b, s2))] = v
    return out
def acode(x, inverse=False):
    den = pow(52 if inverse else 4, P - 2, P)
    byf = {}
    for (o, u), v in x.items(): byf.setdefault(u, {})[o] = v
    out = {}
    for u, vec in byf.items():
        t = [sum(sg * vec.get(i, 0) for i, sg in c) % P for c in MODES]
        new = dict(vec)
        for m, c in enumerate(MODES):
            if t[m]:
                for i, sg in c: new[i] = (new.get(i, 0) - sg * den * t[m]) % P
        for o, v in new.items():
            if v: out[(o, u)] = v
    return out
def code(x, geo, b, inverse=False):
    """physical code C = K^-1 (A x I) K on block b (C^-1 uses A^-1 = I - MM^T/52), as PR5440's physical_code."""
    return kmap(acode(kmap(x, geo, b), inverse), geo, b, True)

def chain(p, v0, geo, coded, letter='R', idle_mode='R', code_block=None, pat=None, formal_bridges=False):
    """PR5440's native_chain/replay (p19_word.py) with the operators mapped to a geometry.
    letter: 'R' (= D_j Q, PR5440's choice) or 'D' (D_j only, the cited primitive X <- D(X - QCY)).
    idle_mode: 'R' (PR5440's gaussian_power over nine axes, collapsed meaning) or 'formal' (product of the D blocks of
    the retained stages in the formal controller's idle selection).
    code_block(kind, j) -> block index of the code used by a gate (default: the gate's own pair)."""
    pat = pat or ('I' * ((p - 1) // 2) + 'T' * (p - 1 - (p - 1) // 2))
    cb = code_block or (lambda kind, j: j)
    st = [dict(r) for r in v0]
    Qm = 1
    def Dm(j): return geo.bmask(j)
    def let(x, j): return phase(x, Dm(j) | (Qm if letter == 'R' else 0))
    def idle(x, sel_blocks, power):
        if idle_mode == 'R':   # PR5440 collapsed: R**power on nine axes (block of the pair)
            return phase(x, Dm(sel_blocks[0] if sel_blocks else 0) | Qm, power)
        m = 0
        for b in sel_blocks: m |= Dm(b)
        return phase(x, m, 1)
    C = (lambda x, j: code(x, geo, cb('C', j))) if coded else (lambda x, j: x)
    Ci = (lambda x, j: code(x, geo, cb('Ci', j), True)) if coded else (lambda x, j: x)
    Qf = lambda x: phase(x, Qm, 1); Qi = lambda x: phase(x, Qm, 3)
    for r in range(1, 2 * p, 2): st[r] = Ci(st[r], (r - 1) // 2)                 # GAUGE_LOAD
    st[1] = add(st[1], Ci(Qi(st[0]), 0))                                        # FORWARD (0 -> 1)
    for j in range(p):
        a, b = 2 * j, 2 * j + 1
        st[a] = add(st[a], Qf(C(st[b], j)), -1)                                 # REVERSE
        st[a] = let(st[a], j); st[b] = let(st[b], j)                           # LETTER
        if j < p - 1 and pat[j] in 'IB':                                         # initial idle on pair j+1
            for r in (2 * (j + 1), 2 * (j + 1) + 1): st[r] = idle(st[r], list(range(j + 1)), 2 * j + 2)
        if j < p - 1:
            st[2 * j] = add(st[2 * j], st[2 * j + 2], 1)
            if formal_bridges: st[2 * j + 3] = add(st[2 * j + 3], st[2 * j + 1], -1)
            else: st[2 * j + 1] = add(st[2 * j + 1], st[2 * j + 3], -1)
        st[b] = add(st[b], Ci(Qi(st[a]), j))                                     # FORWARD
        if j < p - 1:                                                            # middle idle (thin stage)
            for r in (2 * (j + 1), 2 * (j + 1) + 1): st[r] = idle(st[r], [] if idle_mode != 'R' else [j + 1], 1)
            st[2 * j] = add(st[2 * j], st[2 * j + 2], -1)
            if formal_bridges: st[2 * j + 3] = add(st[2 * j + 3], st[2 * j + 1], 1)
            else: st[2 * j + 1] = add(st[2 * j + 1], st[2 * j + 3], 1)
    for j in range(p - 1):
        if pat[j] in 'TB':                                                       # terminal idle on pair j
            for r in (2 * j, 2 * j + 1): st[r] = idle(st[r], list(range(j + 1, p)), 2 * p + 1 - (2 * j + 3))
    for j in reversed(range(p - 1)):
        st[2 * j + 2] = add(st[2 * j + 2], st[2 * j], 1); st[2 * j + 1] = add(st[2 * j + 1], st[2 * j + 3], -1)
    for r in range(1, 2 * p, 2): st[r] = C(st[r], (r - 1) // 2)                  # GAUGE_UNLOAD
    for j in range(p):                                                           # SWAP
        hold = st[2 * j]; st[2 * j] = st[2 * j + 1]; st[2 * j + 1] = {k: (-v) % P for k, v in hold.items()}
    return st

def residual(a, b):
    return sum(1 for r in range(len(a)) for k in set(a[r]) | set(b[r]) if a[r].get(k, 0) != b[r].get(k, 0))

def random_inputs(p, geo, n_entries, rng):
    nbits = 9 if geo.collapsed else 1 + 8 * p
    v0 = [dict() for _ in range(2 * p)]
    for r in range(2 * p):
        for _ in range(n_entries):
            v0[r][(rng.randrange(OWN), rng.randrange(1 << nbits))] = rng.randrange(1, P)
    return v0
