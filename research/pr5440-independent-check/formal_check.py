"""Is my fresh-geometry uncoded chain (letter D, formal idles) the invoiced formal controller (PR4940 controller(),
thinned), evaluated with PR5440's own p19_bank semantics (centre exponent -> Q power, stage bit -> D on its block)?"""
import sys, random, importlib.util
from pathlib import Path
from fractions import Fraction
sys.path.insert(0, str(__import__('pathlib').Path(__file__).resolve().parent))
import fresh as F
ROOT = Path(sys.argv[1])
sys.argv = sys.argv[:1]
spec = importlib.util.spec_from_file_location('c4940', ROOT/'research/coded-reverse-controller/verify.py'); mod = importlib.util.module_from_spec(spec)
import contextlib, io
with contextlib.redirect_stdout(io.StringIO()): spec.loader.exec_module(mod)
P = F.P
def formal_apply(p, pat, v0, geo):
    keep = mod.keep_of(p, pat); co, _ = mod.controller(p, keep, {s: 1 for s in range(2*p+1)})
    thick = set(range(1, 2*p, 2)); M = mod.thin(co, thick); blocks = {s: i for i, s in enumerate(sorted(thick))}
    n = 2*p; out = [dict() for _ in range(n)]
    for term, m in M.items():
        q, sel = term[0], term[1:]
        cnt = {}
        for s, bit in enumerate(sel):
            if bit and s in blocks: cnt[geo.bmask(blocks[s])] = cnt.get(geo.bmask(blocks[s]), 0) + bit
        for i in range(n):
            for j in range(n):
                c = m[i][j]
                if not c: continue
                cc = (c.numerator % P) * pow(c.denominator % P, P-2, P) % P
                img = F.phase(v0[j], 1, q % 4) if q else dict(v0[j])
                for mk, e in cnt.items(): img = F.phase(img, mk, e)
                out[i] = F.add(out[i], img, cc)
    return out
rng = random.Random(3)
for p in (2, 3):
    pat = 'I'*((p-1)//2) + 'T'*(p-1-(p-1)//2)
    for collapsed in (True, False):
        geo = F.Geo(p, collapsed); v0 = F.random_inputs(p, geo, 40, rng)
        want = formal_apply(p, pat, v0, geo)
        for letter, idle, fb in (('D', 'formal', False), ('D', 'formal', True), ('R', 'R', False), ('R', 'R', True)):
            got = F.chain(p, v0, geo, False, letter, idle, formal_bridges=fb)
            print('p=%d %-9s: uncoded listing letter=%s idle=%-6s formal_bridges=%-5s vs invoiced formal controller: residual %d' % (p, 'collapsed' if collapsed else 'fresh', letter, idle, fb, F.residual(got, want)))
