"""Confirmation with PR5440's own modules (p19_word.py, physical_gauge_lift.py), collapsed nine-axis geometry, p=2:
(a) the coded listing vs its standard word with the cited primitive's letter D (PR5440 supports letter='D');
(b) PR5440's standard word (letter R) vs the invoiced formal controller (PR4940 controller(), thinned), evaluated on the
    same 512 addresses with p19_bank.py's monomial semantics (centre exponent -> Q power, retained stage -> D)."""
import sys, random, importlib.util, contextlib, io
from pathlib import Path
ROOT = Path(sys.argv[1]); sys.argv = sys.argv[:1]
def load(n, p):
    s = importlib.util.spec_from_file_location(n, p); m = importlib.util.module_from_spec(s)
    with contextlib.redirect_stdout(io.StringIO()): s.loader.exec_module(m)
    return m
PKG = ROOT/'research/coded-reverse-word-p19'
L = load('L', PKG/'physical_gauge_lift.py'); W = load('W', PKG/'p19_word.py'); B = load('B', PKG/'p19_bank.py')
mod = load('c4940', ROOT/'research/coded-reverse-controller/verify.py')
import json
rec = json.loads((ROOT/'research/coded-native-caller/NORMALIZERS.json').read_text())
modes = [[(i, s) for i, s in enumerate(c) if s] for c in rec['mode_columns']]
C = lambda x: L.physical_code(x, rec['ports'], modes, False); Ci = lambda x: L.physical_code(x, rec['ports'], modes, True)
P = W.P; rng = random.Random(9); p = 2
v0 = [[[rng.randrange(P) if rng.random() < 0.02 else 0 for _ in range(512)] for _ in range(112)] for _ in range(2*p)]
for letter in ('R', 'D'):
    out, _ = W.replay(p, W.native_chain(p, W.pattern_of(p), letter), v0, L, C, Ci)
    ref = W.standard(p, v0, letter, L)
    print('(a) letter %s: coded listing vs its standard word, differing entries: %d' % (letter, sum(L.diff(out[r], ref[r]) for r in range(2*p))))
# (b) formal controller on 512 addresses: centre = axis 0, every retained stage's 8 axes = axes 1..8 (collapsed)
pat = W.pattern_of(p); keep = mod.keep_of(p, pat)
co, _ = mod.controller(p, keep, {s: 1 for s in range(2*p+1)}); M = mod.thin(co, set(range(1, 2*p, 2)))
def mono(row, term):
    x = row[:]; q = term[0]
    if q % 4: B.butterfly(x, 0, q % 4)
    k = sum(b for s, b in enumerate(term[1:]) if s % 2)
    for _ in range(k):
        for ax in range(1, 9): B.butterfly(x, ax, 1)
    return x
want = [[[0]*512 for _ in range(112)] for _ in range(2*p)]
for term, m in M.items():
    for i in range(2*p):
        for j in range(2*p):
            c = m[i][j]
            if not c: continue
            cc = (c.numerator % P) * pow(c.denominator % P, P-2, P) % P
            for o in range(112):
                img = mono(v0[j][o], term)
                want[i][o] = [(a + cc*b) % P for a, b in zip(want[i][o], img)]
ref = W.standard(p, v0, 'R', L)
print('(b) PR5440 standard word (letter R) vs invoiced formal controller, differing entries: %d of %d'
      % (sum(L.diff(ref[r], want[r]) for r in range(2*p)), 2*p*112*512))
