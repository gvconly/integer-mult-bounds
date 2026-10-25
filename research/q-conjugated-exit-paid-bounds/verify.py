#!/usr/bin/env python3
"""Replay both applications, run the controls and compare with RESULT.json.

  python3 -B research/q-conjugated-exit-paid-bounds/verify.py           # check (read-only)
  python3 -B research/q-conjugated-exit-paid-bounds/verify.py --write   # regenerate RESULT.json

Controls: (1) with one common centre line q commutes with the code and the decode is free (rank 0), while
identifying only two owners' lines leaves rank 120 (rows keep other private coordinates);
(2) a diagonal code needs no paid call; (3) the uncoded contact and the old 112-port code reproduce the retained
query count and coded-active-query-bound's 112; (4) zeroing the cross-cube couplings removes every odd-E8 witness."""
import json
import sys
from fractions import Fraction as Fr
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import global120 as g  # noqa: E402
import e8  # noqa: E402
from common import rank, inverse  # noqa: E402


def controls():
    out = {}
    rows, r, c, M = g.owners()
    J = g.code(M)
    Ji = inverse(J)
    r2 = list(r)
    r2[1] = r2[0]
    out['decode_rank_with_two_owners_identified'] = rank(g.conj_rows(g.Classes(r2), Ji, 'Y', range(g.N)))
    out['decode_rank_with_one_common_centre_line'] = rank(g.conj_rows(g.Classes([r[0]] * g.N), Ji, 'Y', range(g.N)))
    I = [[Fr(int(a == b)) for b in range(g.N)] for a in range(g.N)]
    out['decode_rank_diagonal_code'] = rank(g.conj_rows(g.Classes(r), I, 'Y', range(g.N)))
    cl = g.Classes(r)
    out['uncoded_CI_contact_rank'] = rank([{(cl.q[a], bk, a): Fr(1) for bk in ('X0', 'X1')} for a in range(g.N)] * 2)
    nz, kt, M8, P, K, n = e8.base()
    KP = [[K[a][b] + P[a][b] for b in range(n)] for a in range(n)]
    cube = {p: ci for ci, cc in enumerate(kt['literal_cube_partition']) for p in cc['ports']}
    cut = [[x if cube[a] == cube[b] else 0 for b, x in enumerate(row)] for a, row in enumerate(KP)]
    out['odd_e8_witness_rows_without_cross_cube_couplings'] = sum(1 for a in range(n) if any(cut[a][b] for b in range(n) if b != a))
    return out


def main():
    if not __debug__:
        raise SystemExit('assertions are the checks: run without -O')
    res = {'global120': g.main(verbose=False), 'odd_e8': e8.main(verbose=False), 'controls': controls()}
    G, E, C = res['global120'], res['odd_e8'], res['controls']
    assert G['distinct_centre_functionals'] == G['distinct_G_e0'] == G['distinct_row0_G'] == G['distinct_G_inverse_e0'] == 120
    assert G['decode_barJinv_120'] == 120 and G['decode_old112_calibration'] == 112
    assert G['cut5_CP_nominated'] == 240 and G['cut5_CI_nominated'] == 240
    assert G['cut5_CP_Y0_already_decoded'] == 232 and G['cut5_CI_Y0_already_decoded'] == 240
    assert G['calibration_uncoded_CP_contact'] == 112 and G['retirements_without_scheduled_paid_call_on_Yk'] == 7
    assert G['retirements_whose_IDLE_matches_the_count'] == 8
    assert G['forced_new_or_rerouted_calls'] == 840 and G['allowance_rank1_jobs'] == 916
    assert G['idles_disjoint_from_own_FULL_stages'] == 15 and G['in_sector_required_calls_15_retirements'] == 1800
    assert E['paid_postdecoder_lower_bound_calls'] == 112 ** 3 and E['min_offdiagonal_support_of_K_plus_P_row'] > 0
    assert C['decode_rank_with_one_common_centre_line'] == 0 and C['decode_rank_diagonal_code'] == 0
    assert C['uncoded_CI_contact_rank'] == 120 and C['odd_e8_witness_rows_without_cross_cube_couplings'] == 0
    path = HERE / 'RESULT.json'
    text = json.dumps(res, indent=1, sort_keys=True) + '\n'
    if '--write' in sys.argv:
        path.write_text(text)
    else:
        assert path.read_text() == text, 'RESULT.json differs from the replay'
    for sec in ('global120', 'odd_e8', 'controls'):
        for k, v in res[sec].items():
            if k != 'retirements':
                print('%-9s %-52s %s' % (sec, k, v))
    print('PASS_Q_CONJUGATED_EXIT_PAID_BOUNDS_NO_NEW_KAPPA')


if __name__ == '__main__':
    main()
