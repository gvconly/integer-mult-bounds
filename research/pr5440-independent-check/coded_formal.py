"""Coded listing (PR5440's gate structure: gauge loads/unloads on odd roles, coded reverse reads Q.C, forward reads
C^-1 after Q^-1) with controller()'s semantics (letter D_j, formal idles, formal bridges), against the uncoded listing
with the same semantics, which formal_check.py shows equals controller() exactly. Codes on each gate's own pair block.
usage: python3 -B coded_formal.py ROOT   (ROOT = checkout of PR5440's head)"""
import sys, random
sys.path.insert(0, str(__import__('pathlib').Path(__file__).resolve().parent))
import fresh as F
rng = random.Random(21)
for p, collapsed in ((2, True), (2, False), (3, True)):
    geo = F.Geo(p, collapsed); v0 = F.random_inputs(p, geo, 10, rng)
    std = F.chain(p, v0, geo, False, 'D', 'formal', formal_bridges=True)
    cod = F.chain(p, v0, geo, True, 'D', 'formal', formal_bridges=True)
    print('p=%d %-9s coded(controller semantics, own-pair codes) vs controller: residual %d'
          % (p, 'collapsed' if collapsed else 'fresh', F.residual(cod, std)), flush=True)
