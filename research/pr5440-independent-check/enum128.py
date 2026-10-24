"""p=2 fresh-block geometry, invoiced semantics (letter D, formal idles, formal bridges) and PR5440's semantics.
The two REVERSE codes are fixed to their own stage's block (the native call at stage j carries stage j's normalizers);
the seven other code gates (2 loads, 3 forward reads, 2 unloads) range over both blocks: 128 assignments.
Runs longer than 30 s are recorded as not completed."""
import sys, random, itertools, time, json, signal
sys.path.insert(0, str(__import__('pathlib').Path(__file__).resolve().parent))
import fresh as F
class TO(Exception): pass
def alarm(*a): raise TO()
signal.signal(signal.SIGALRM, alarm)
p = 2; geo = F.Geo(p); rng = random.Random(5); v0 = F.random_inputs(p, geo, 1, rng)
# gate order of code calls in fresh.chain: load r1, load r3, FORWARD(0), REVERSE j0, FORWARD j0, REVERSE j1, FORWARD j1, unload r1, unload r3
res = {}
for name, letter, idle, fb in (('invoiced semantics (D, formal idles, formal bridges)', 'D', 'formal', True),
                               ('PR5440 listing semantics (R, nine-axis idles, its bridges)', 'R', 'R', False)):
    std = F.chain(p, v0, geo, False, letter, idle, formal_bridges=fb)
    t0 = time.time(); done = zeros = timeouts = 0; best = None
    for free in itertools.product(range(p), repeat=7):
        bits = list(free[:3]) + [0, free[3], 1, free[4], free[5], free[6]]
        it = iter(bits); signal.alarm(30)
        try:
            out = F.chain(p, v0, geo, True, letter, idle, (lambda kind, j: next(it)), formal_bridges=fb); signal.alarm(0)
        except TO:
            timeouts += 1; continue
        done += 1; r = F.residual(out, std); zeros += (r == 0)
        if best is None or r < best[0]: best = (r, bits)
    res[name] = dict(assignments=128, completed=done, timeouts=timeouts, zero_residual=zeros, min_residual=best[0] if best else None, argmin=best[1] if best else None)
    print(name, res[name], '(%.0f s)' % (time.time() - t0), flush=True)
json.dump(res, open(__import__('pathlib').Path(__file__).resolve().parent / 'ENUM128.json', 'w'), indent=1)
