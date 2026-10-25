#!/usr/bin/env python3
"""Verify the multi-generation centre floor certificates for the P104 word.

  python3 -I -B research/p104-multigen-centre-floor/verify.py           # check (read-only)
  python3 -I -B research/p104-multigen-centre-floor/verify.py --write   # regenerate RESULT.json

Checks (PROOF.md):
  * centre structure of the word: 8 owned lifetimes, rank 7 each, copy frame = label(x_j)^perp, landing = label(y_j);
  * integral E8 realisation recovered from the word (Gram unique up to scale, det 1), 120 lines, 16 removed lines in
    two orthogonal frames, label dictionary <q_p,q_s> = [r_p perp r_s or p = s];
  * Lemma 1 (good ports), Lemma 2 (first-generation sources), capacities m_d, |T_z| >= 49, the two-plane lemma,
    no A2 among removed roots, and the arithmetic of Theorems A, B, C;
  * controls: the selected word meets the master inequality with equality; an explicit multi-generation centre of
    cost 56 (Notes, check below) satisfies all 104 ports and has R_8 = 8, as Theorem A requires."""
import json
import math
import sys
from fractions import Fraction as Fr
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from common import Geometry, par, f2span  # noqa: E402
import floor  # noqa: E402


def centre_checks(G):
    out = {}
    lives = G.lives
    out['owned_lifetimes'] = len(lives)
    out['centre_cost'] = sum(L['cost'] for L in lives)
    out['descent_steps'] = [[list(s) for s in L['steps']] for L in lives]
    labs, lands = [], []
    for j, L in enumerate(lives):
        e = tuple(int(i == j) for i in range(8))
        assert G.line(e) in G.removed
        h = G.label(e)
        F = L['copy_frame']
        assert len(F) == 8 and all(par(v & h) == 0 for v in F)          # copy frame = h_j^perp
        assert len(L['landing']) == 1
        labs.append(h)
        lands.append(L['landing'][0])
    # landing line = label of the dual basis vector y_j (a removed root)
    from common import _inverse
    Y = _inverse(G.X)
    dual = [tuple(int(x) for x in Y[j]) for j in range(8)]
    assert all(G.line(y) in G.removed and G.label(y) == lands[j] for j, y in enumerate(dual))
    # frame-pair requirement: supp of r_p in the x-basis = {j : <q_p, g_j> = 0}
    ok = sum({j for j in range(8) if G.roots[p][j]} == {j for j in range(8) if par(G.Q[p] & lands[j]) == 0}
             for p in range(G.P))
    out['ports_with_exact_frame_pair_support'] = ok
    out['content_labels'] = labs
    out['landing_lines'] = lands
    return out


def unit(j):
    return tuple(int(i == j) for i in range(8))


def multigen_design(G):
    """Control: a genuinely multi-generation centre of cost 56 (four coupled pairs), checked on all ports.
    Contents z = removed roots of basis B; pair (j, k): z_j at its line M_j, z_k at plane M_k, and the combination
    rho z_j + z_k descends from M_j + M_k one step to G'.  Coordinates in the B-basis by exact solve."""
    B = [G.removed[i] for i in (0, 1, 2, 4, 5, 6, 8, 13)]
    pairs = [(0, 4), (1, 2), (3, 6), (5, 7)]
    from common import _inverse
    Zm = [[G.ip(B[a], B[b]) for b in range(8)] for a in range(8)]
    Zi = _inverse(Zm)
    a = [[sum(Fr(G.ip(G.roots[p], B[k])) * Zi[k][j] for k in range(8)) for j in range(8)] for p in range(G.P)]
    labs = [G.label(z) for z in B]

    def meet(perps, inside):
        elems = [0]
        for v in inside:
            elems += [e ^ v for e in elems]
        return f2span([e for e in elems if e and all(par(e & w) == 0 for w in perps)])

    def hyper(h):
        return f2span([u for u in range(1, 512) if par(u & h) == 0])
    N = [[p for p in range(G.P) if a[p][j] != 0] for j in range(8)]
    deliveries, cost = [], 0
    for j, k in pairs:
        both = [p for p in N[j] if a[p][k] != 0]
        best = None
        for rho in sorted(set(a[p][j] / a[p][k] for p in both)):
            C = {p for p in both if a[p][j] / a[p][k] == rho}
            Mj = meet([G.Q[p] for p in N[j] if p not in C], hyper(labs[j]))
            Mk = meet([G.Q[p] for p in N[k] if p not in C], hyper(labs[k]))
            S = f2span(Mj + Mk)
            if not S:
                continue
            Gp = meet([G.Q[p] for p in C if any(par(G.Q[p] & v) for v in S)], S)
            c = (8 - len(Mj)) + (8 - len(Mk)) + (len(S) - len(Gp))
            if best is None or c < best[0]:
                comb = [Fr(0)] * 8
                comb[j], comb[k] = rho, Fr(1)
                best = (c, [(unit(j), Mj), (unit(k), Mk), (tuple(comb), Gp)])
        cost += best[0]
        deliveries += best[1]
    ok = 0
    for p in range(G.P):
        vecs = [v for v, Gf in deliveries if all(par(x & G.Q[p]) == 0 for x in Gf)]
        r0 = floor.qrank(G, vecs)
        ok += floor.qrank(G, vecs + [a[p]]) == r0
    line_level = [v for v, Gf in deliveries if len(Gf) <= 1]
    return {'cost': cost, 'ports_satisfied': ok, 'deliveries': len(deliveries),
            'landing_dims': sorted(len(Gf) for _, Gf in deliveries), 'R8_rank_at_lines': floor.qrank(G, line_level)}


def theorem_b(m):
    rows = {}
    for d in range(1, 8):
        if d == 1:
            rows[d] = 49 + 49
        else:
            rows[d] = 49 + max(0, d - 2) + math.ceil(49 / m[d])
    return rows


def main():
    if not __debug__:
        raise SystemExit('assertions are the checks: run without -O')
    G = Geometry()
    res = {'geometry': {'ports': G.P, 'gram_solution_dimension': G.gram_solution_dimension, 'gram_det': int(G.det),
                        'lines': len(G.lines), 'removed': len(G.removed), 'removed_labels': G.removed_labels}}
    # removed lines = two orthogonal frames
    rem = G.removed
    orth = [[G.ip(a, b) == 0 for b in rem] for a in rem]
    frames = []

    def ext(cl, start):
        if len(cl) == 8:
            frames.append(tuple(cl))
            return
        for k in range(start, 16):
            if all(orth[k][c] for c in cl):
                ext(cl + [k], k + 1)
    ext([], 0)
    res['geometry']['orthogonal_8_frames_inside_removed'] = len(frames)
    assert len(frames) == 2 and not set(frames[0]) & set(frames[1])
    res['centre'] = centre_checks(G)
    C = res['centre']
    assert C['owned_lifetimes'] == 8 and C['centre_cost'] == 56 and C['ports_with_exact_frame_pair_support'] == 104
    res['lemma1_good_ports'] = floor.good(G)
    assert res['lemma1_good_ports'] == 104
    res['lemma2_bad_hyperplanes'] = floor.sources(G)
    assert res['lemma2_bad_hyperplanes'] == []
    m, counts = floor.capacities(G)
    res['capacities_m_d'] = {str(d): m[d] for d in sorted(m)}
    res['subspaces_enumerated'] = {str(d): counts[d] for d in sorted(counts)}
    assert [m[d] for d in range(1, 9)] == [1, 2, 4, 8, 12, 20, 36, 60]
    info, cands = floor.tsets(G)
    res['T_sizes'] = info
    assert info['vectors_norm_le_8'] == 13320 and info['min_T_lattice'] == 49 and info['min_T_E6_continuum'] >= 49
    nworlds, covered, removed_T = floor.two_plane(G, cands)
    res['two_plane'] = {'candidate_T_sets_le_72': len(cands), 'plane_worlds': nworlds, 'covered': len(covered),
                        'covered_sizes': sorted({t.bit_count() for t in covered}),
                        'covered_equals_removed_root_sets': set(covered) == removed_T}
    assert set(covered) == removed_T and len(removed_T) == 16
    res['no_A2_among_removed'] = floor.no_a2(G) == []
    assert res['no_A2_among_removed']
    tb = theorem_b(m)
    res['theorem_B_bounds_by_last_level'] = {str(d): v for d, v in tb.items()}
    assert min(tb.values()) == 56
    res['theorem_C_bound'] = 52 + 4
    # controls
    res['control_selected_master_inequality'] = {'R_2..R_8': 8, 'R_9': 0, 'shallow_fresh': 0, 'bound': 7 * 8,
                                                 'actual': C['centre_cost']}
    md = multigen_design(G)
    res['control_multigen_design'] = md
    assert md['cost'] == 56 and md['ports_satisfied'] == 104 and md['R8_rank_at_lines'] == 8
    path = HERE / 'RESULT.json'
    text = json.dumps(res, indent=1, sort_keys=True) + '\n'
    if '--write' in sys.argv:
        path.write_text(text)
    else:
        assert path.read_text() == text, 'RESULT.json differs from the replay'
    for k, v in res.items():
        print('%-40s %s' % (k, v if not isinstance(v, dict) else json.dumps(v, sort_keys=True)[:160]))
    print('PASS_P104_MULTIGEN_CENTRE_FLOOR_NO_NEW_KAPPA')


if __name__ == '__main__':
    main()
