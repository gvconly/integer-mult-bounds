import sys, random, importlib.util, json
from pathlib import Path
ROOT = Path(sys.argv[1]); sys.argv = sys.argv[:2]
sys.path.insert(0, str(__import__('pathlib').Path(__file__).resolve().parent))
import fresh as F
spec = importlib.util.spec_from_file_location('L', ROOT/'research/coded-reverse-word-p19/physical_gauge_lift.py'); L = importlib.util.module_from_spec(spec); spec.loader.exec_module(L)
spec = importlib.util.spec_from_file_location('W', ROOT/'research/coded-reverse-word-p19/p19_word.py'); W = importlib.util.module_from_spec(spec); spec.loader.exec_module(W)
P = F.P
def wht(row):
    a = list(row); h = 1
    while h < 512:
        for i in range(0, 512, 2*h):
            for j in range(i, i+h):
                x, y = a[j], a[j+h]; a[j], a[j+h] = (x+y) % P, (x-y) % P
        h *= 2
    return a
def to_freq(bank):
    out = {}
    for o, row in enumerate(bank):
        Fr = wht(row)
        for u, v in enumerate(Fr):
            if v: out[(o, u)] = v
    return out
rng = random.Random(1)
geo = F.Geo(1, collapsed=True)
# 1. signed map
bank = [[rng.randrange(P) if rng.random() < 0.02 else 0 for _ in range(512)] for _ in range(112)]
for inv in (False, True):
    theirs = to_freq([L.signed_row(F.PORTS[o], bank[o], inv) for o in range(112)])
    mine = F.kmap(to_freq(bank), geo, 0, inv)
    print('signed map inverse=%s match:' % inv, theirs == mine)
# 2. gaussian
th = to_freq(L.gaussian(bank, range(9))); mi = F.phase(to_freq(bank), 0x1FF, 1); print('R match', th == mi)
th = to_freq(L.gaussian(bank, [0], True)); mi = F.phase(to_freq(bank), 1, 3); print('Q^-1 match', th == mi)
# 3. physical code
th = to_freq(L.physical_code(bank, F.PORTS, F.MODES, False)); mi = F.code(to_freq(bank), geo, 0, False); print('C match', th == mi)
th = to_freq(L.physical_code(bank, F.PORTS, F.MODES, True)); mi = F.code(to_freq(bank), geo, 0, True); print('C^-1 match', th == mi)
# 4. whole chain p=2 and p=3 vs PR5440's address-domain replay
for p in (2, 3):
    v0 = [[[rng.randrange(P) if rng.random() < 0.01 else 0 for _ in range(512)] for _ in range(112)] for _ in range(2*p)]
    C = lambda x: L.physical_code(x, F.PORTS, F.MODES, False); Ci = lambda x: L.physical_code(x, F.PORTS, F.MODES, True)
    out, _ = W.replay(p, W.listing(p), v0, L, C, Ci)
    ref = W.standard(p, v0, 'R', L)
    theirs = [to_freq(b) for b in out]; theirs_std = [to_freq(b) for b in ref]
    f0 = [to_freq(b) for b in v0]
    mine = F.chain(p, f0, geo, True); mine_std = F.chain(p, f0, geo, False)
    print('p=%d: my coded == their coded: %s | my standard == their standard: %s | their coded vs standard residual %d | mine %d'
          % (p, mine == theirs, mine_std == theirs_std, F.residual(theirs, theirs_std), F.residual(mine, mine_std)))
