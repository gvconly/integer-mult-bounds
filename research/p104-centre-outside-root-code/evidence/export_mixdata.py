"""Write the input of pairmix.c: ports, labels, the copy-frame labels h_j and the x-basis coordinates of every port."""
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from common import Geometry  # noqa: E402

G = Geometry()
h = [G.label(tuple(int(i == j) for i in range(8))) for j in range(8)]
with open(sys.argv[1], 'w') as f:
    f.write('%d\n' % G.P + ' '.join(map(str, G.Q)) + '\n' + ' '.join(map(str, h)) + '\n')
    for r in G.roots:
        f.write(' '.join(str(int(c)) for c in r) + '\n')
