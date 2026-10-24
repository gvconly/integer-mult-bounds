import sys, random, itertools
sys.path.insert(0, str(__import__('pathlib').Path(__file__).resolve().parent))
import fresh as F
rng = random.Random(7)
print('collapsed controls:')
for p in (2, 3):
    geo = F.Geo(p, collapsed=True)
    v0 = F.random_inputs(p, geo, 60, rng)
    for letter, idle in (('R', 'R'), ('D', 'R'), ('R', 'formal'), ('D', 'formal')):
        r = F.residual(F.chain(p, v0, geo, True, letter, idle), F.chain(p, v0, geo, False, letter, idle))
        print('  p=%d collapsed letter=%s idle=%s residual=%d' % (p, letter, idle, r))
print('fresh-block geometry (central axis + one 8-axis block per retained stage):')
pol = {'own pair block': None, 'all codes on block 0': (lambda kind, j: 0)}
for p in (2, 3):
    geo = F.Geo(p)
    v0 = F.random_inputs(p, geo, 60, rng)
    for (letter, idle), (pname, cb) in itertools.product((('R', 'R'), ('R', 'formal'), ('D', 'formal')), pol.items()):
        a = F.chain(p, v0, geo, True, letter, idle, cb); b = F.chain(p, v0, geo, False, letter, idle, cb)
        tot = sum(len(x) for x in b)
        print('  p=%d fresh letter=%s idle=%-6s codes=%-20s residual=%d (uncoded output: %d nonzero entries)' % (p, letter, idle, pname, F.residual(a, b), tot))
