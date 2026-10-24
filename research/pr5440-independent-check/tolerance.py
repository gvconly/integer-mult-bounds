"""How many extra rank-one jobs the p=19 invoice of PR5440 can absorb. Uses PR5440's own p19_word.invoice and root
(pattern I^9 T^9, odd roles thick, profile P_NEW); extra jobs are added to the rank-one demand.
usage: python3 -B tolerance.py ROOT   (ROOT = checkout of PR5440's head)"""
import sys, json, importlib.util
from pathlib import Path
ROOT = Path(sys.argv[1]); sys.argv = sys.argv[:1]
spec = importlib.util.spec_from_file_location('W', ROOT / 'research/coded-reverse-word-p19/p19_word.py')
W = importlib.util.module_from_spec(spec); spec.loader.exec_module(W)
p = 19; pat = 'I' * 9 + 'T' * 9; thick = frozenset(range(1, 2 * p, 2))
T, D = W.invoice(p, pat, thick, W.P_NEW)
mass = lambda H: sum(r * n for r, n in H.items())
def rootx(e):
    D2 = dict(D); D2[1] = D2.get(1, 0) + e; return W.root(T, D2)
out = dict(root=W.root(T, D), surplus=mass(T) - mass(D), extra_224x19=rootx(224 * 19),
           extra_1103=rootx(1103), extra_1104=rootx(1104), extra_1202=rootx(1202), extra_1203=rootx(1203))
for k, v in out.items(): print(k, v)
json.dump(out, open(Path(__file__).resolve().parent / 'TOLERANCE.json', 'w'), indent=1)
