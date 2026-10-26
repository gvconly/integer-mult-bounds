#!/usr/bin/env python3
"""Verify the certificates of PROOF.md: P104 centre contents outside the root code.

  python3 -I -B research/p104-centre-outside-root-code/verify.py           # check (read-only)
  python3 -I -B research/p104-centre-outside-root-code/verify.py --write   # regenerate RESULT.json

Checks (PROOF.md):
  * the word's centre: 8 owned lifetimes of rank 7, integral E8 realisation, 16 removed lines;
  * the word's unit is isometric to #968's P104 = E8 - F0 - H4 of decoder-floors-frontier-format (so the Model C
    floors for #968's unit apply verbatim);
  * Collapse controls: the selected centre and the R_8 = 0 all-planes centre are mixing-free and collapse to Model C
    decoders of cost 56, with the Union-lemma bound tight;
  * randomized exact tests of the Union lemma (pure classes) and of Lemma S (with mixing);
  * the mixing gap: a two-content mixing centre that Model C prices higher than the centre pays;
  * evidence: every single truncation of a frame-pair content breaks all of its users; in the overcomplete family
    with removed-root directions the frame-pair has no improving single-port change; landing frames whose worlds meet
    more frames than a frame-counting proof of 56 allows."""
import json
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path.insert(0, str(HERE))
sys.path.insert(1, str(ROOT / 'research' / 'decoder-floors-frontier-format' / 'lib'))
from common import Geometry, par  # noqa: E402
import checks  # noqa: E402
import units  # noqa: E402  (#2338: E8 lines, F0, H4)

# the R_8 = 0 centre of research/p104-multigen-centre-floor (verify.py ALL_PLANES): bases of P_j and P'_j
ALL_PLANES = [((1, 6), (7, 505)), ((1, 24), (25, 487)), ((1, 52), (52, 458)), ((2, 80), (82, 431)),
              ((1, 490), (21, 491)), ((1, 448), (62, 448)), ((2, 65), (67, 446)), ((3, 137), (137, 373))]


def main():
    if not __debug__:
        raise SystemExit('assertions are the checks: run without -O')
    G = Geometry()
    res = {}
    lives = G.lives
    assert len(lives) == 8 and all(L['cost'] == 7 and len(L['landing']) == 1 for L in lives)
    res['centre'] = {'owned_lifetimes': len(lives), 'cost': sum(L['cost'] for L in lives),
                     'ports': G.P, 'removed_lines': len(G.removed), 'gram_det': int(G.det)}
    res['isometry_with_968'] = checks.isometry_968(G, units)
    sel = checks.collapse_selected(G)
    assert sel['collapsed_cost'] == sel['centre_cost'] == 56 and sel['model_C_costs'] == [7] * 8
    res['collapse_selected'] = sel
    ap = checks.collapse_all_planes(G, ALL_PLANES)
    assert ap['collapsed_cost'] == ap['centre_cost'] == 56
    res['collapse_all_planes'] = ap
    res['union_lemma_random'] = checks.union_lemma_random(G, 3000, 20261026)
    res['lemma_s_random'] = checks.lemma_s_random(G, 3000, 26102026)
    gap = checks.mixing_gap(G)
    assert gap['best_natural_model_C_decoder'] > gap['centre_cost']
    res['mixing_gap'] = gap
    tr = checks.truncation_census(G)
    assert tr['min_broken_ports_per_content'] == tr['users_per_content'] == [52] * 8
    res['truncation_census'] = tr
    rd = checks.removed_direction_family(G)
    assert rd['frame_pair_cost'] == 56 and rd['min_cost_after_one_move'] > 56
    res['removed_direction_family'] = rd
    frames = checks.live_frames(G)
    assert len(frames) == 759
    res['live_frames'] = len(frames)
    res['frame_touch_witnesses'] = checks.frame_touch(G, frames)
    path = HERE / 'RESULT.json'
    text = json.dumps(res, indent=1, sort_keys=True) + '\n'
    if '--write' in sys.argv:
        path.write_text(text)
    else:
        assert path.read_text() == text, 'RESULT.json differs from the replay'
    for k, v in res.items():
        print('%-28s %s' % (k, json.dumps(v, sort_keys=True)[:170]))
    print('PASS_P104_CENTRE_OUTSIDE_ROOT_CODE_NO_NEW_KAPPA')


if __name__ == '__main__':
    main()
